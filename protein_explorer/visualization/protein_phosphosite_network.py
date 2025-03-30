# This should be placed in: protein_explorer/visualization/phosphosite_network.py

"""
Functions for generating phosphosite network visualizations.
This module provides functions to create interactive network visualizations
for phosphosites and their structural matches.
"""

def create_phosphosite_network_visualization(protein_uniprot_id, phosphosites=None, structural_matches=None):
    """
    Create a network visualization of phosphosites and their structural matches.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries (optional)
        structural_matches: Dictionary mapping site names to lists of match dictionaries (optional)
        
    Returns:
        HTML string with network visualization
    """
    # Create network visualization HTML
    html = """
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Phosphosite Structural Network</h5>
        </div>
        <div class="card-body">
            <div id="rmsd-filter-container" class="mb-3">
                <label for="rmsd-filter" class="form-label">
                    RMSD Threshold: <span id="rmsd-value">2.0 Å</span>
                </label>
                <input 
                    type="range" 
                    class="form-range" 
                    id="rmsd-filter" 
                    min="0.5" 
                    max="5.0" 
                    step="0.1" 
                    value="2.0"
                    oninput="document.getElementById('rmsd-value').textContent = this.value + ' Å'; updateNetworkFilter();"
                >
                <div class="d-flex justify-content-between">
                    <small>0.5 Å (Very similar)</small>
                    <small>5.0 Å (Less similar)</small>
                </div>
            </div>
            
            <div class="mb-3">
                <div class="d-flex align-items-center flex-wrap">
                    <div class="d-flex align-items-center me-4 mb-2">
                        <div style="width: 16px; height: 16px; background-color: #4CAF50; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Known protein sites</span>
                    </div>
                    <div class="d-flex align-items-center me-4 mb-2">
                        <div style="width: 16px; height: 16px; background-color: #FF9800; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Unknown protein sites</span>
                    </div>
                    <div class="d-flex align-items-center mb-2">
                        <div style="width: 16px; height: 16px; background-color: #9E9E9E; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Structurally similar sites</span>
                    </div>
                </div>
            </div>
            
            <div id="network-container" style="height: 500px; width: 100%; position: relative; border: 1px solid #ddd; border-radius: 5px;"></div>
            
            <p class="text-muted mt-3 mb-0">
                <small>
                    This network shows the structural relationships between phosphosites in this protein and
                    similar sites in other proteins. Edges represent structural similarity with RMSD below
                    the threshold. Hover over nodes for details and click to view the site page.
                </small>
            </p>
        </div>
    </div>

    <!-- D3.js Library (if not already included) -->
    <script src="https://d3js.org/d3.v7.min.js"></script>
    
    <!-- Network Visualization Script -->
    <script>
    document.addEventListener('DOMContentLoaded', function() {
        // Setup network visualization when DOM is loaded
        setupPhosphositeNetwork();
    });
    
    let networkSimulation = null;
    
    function setupPhosphositeNetwork() {
        const networkContainer = document.getElementById('network-container');
        if (!networkContainer) {
            console.error("Network container not found");
            return;
        }
        
        console.log("Initializing phosphosite network visualization");
        
        // Reset the network container
        networkContainer.innerHTML = '';
        
        // Create information panel
        const infoPanel = document.createElement('div');
        infoPanel.className = 'node-info-panel';
        infoPanel.style.position = 'absolute';
        infoPanel.style.top = '10px';
        infoPanel.style.right = '10px';
        infoPanel.style.width = '250px';
        infoPanel.style.backgroundColor = 'white';
        infoPanel.style.border = '1px solid #ddd';
        infoPanel.style.borderRadius = '5px';
        infoPanel.style.padding = '10px';
        infoPanel.style.boxShadow = '0 0 10px rgba(0,0,0,0.1)';
        infoPanel.style.zIndex = '100';
        infoPanel.style.fontSize = '0.9rem';
        infoPanel.style.maxHeight = '380px';
        infoPanel.style.overflowY = 'auto';
        infoPanel.innerHTML = '<p class="text-center"><em>Hover over a node to see details</em></p>';
        networkContainer.appendChild(infoPanel);
        
        // Extract protein UniProt ID
        const proteinUniprotId = '""" + protein_uniprot_id + """';
        
        // Extract network data from the page DOM
        const networkData = extractNetworkDataFromDOM(proteinUniprotId);
        
        if (!networkData || networkData.nodes.length === 0) {
            networkContainer.innerHTML = '<div class="alert alert-info m-3 mt-5">No phosphosite network data available. This could be because no structural matches were found with RMSD < 2.0.</div>';
            return;
        }
        
        console.log(`Creating network with ${networkData.nodes.length} nodes and ${networkData.links.length} links`);
        
        // Create SVG element
        const svg = d3.select(networkContainer)
            .append('svg')
            .attr('width', '100%')
            .attr('height', '100%')
            .style('position', 'absolute')
            .style('top', '0')
            .style('left', '0')
            .style('background-color', 'white');
        
        // Get dimensions
        const width = networkContainer.clientWidth;
        const height = networkContainer.clientHeight;
        
        // Create group for zoom behavior
        const g = svg.append('g');
        
        // Add zoom behavior
        const zoom = d3.zoom()
            .scaleExtent([0.2, 5])
            .on('zoom', (event) => {
                g.attr('transform', event.transform);
            });
        
        svg.call(zoom);
        
        // Helper functions for styling
        function getNodeColor(d) {
            if (d.type === 'protein') {
                return d.isKnown ? '#4CAF50' : '#FF9800'; // Green for known protein sites, Orange for unknown
            }
            return '#9E9E9E'; // Gray for matched sites
        }
        
        function getLinkColor(d) {
            if (d.rmsd < 1.0) return '#4CAF50'; // Green for very low RMSD
            if (d.rmsd < 1.5) return '#8BC34A'; // Light green for low RMSD
            if (d.rmsd < 2.0) return '#CDDC39'; // Lime for medium RMSD
            return '#FFC107'; // Yellow/amber for higher RMSD
        }
        
        function getLinkWidth(d) {
            // Scale link width based on RMSD - thicker for lower RMSD
            return Math.max(1, 4 / Math.max(0.5, d.rmsd));
        }
        
        // Setup force simulation
        networkSimulation = d3.forceSimulation(networkData.nodes)
            .force('link', d3.forceLink(networkData.links)
                .id(d => d.id)
                .distance(d => d.rmsd ? d.rmsd * 30 : 60))
            .force('charge', d3.forceManyBody()
                .strength(d => d.type === 'protein' ? -200 : -100))
            .force('center', d3.forceCenter(width / 2, height / 2))
            .force('collision', d3.forceCollide().radius(d => d.size + 5))
            .force('x', d3.forceX(width / 2).strength(0.05))
            .force('y', d3.forceY(height / 2).strength(0.05));
        
        // Draw links
        const link = g.append('g')
            .selectAll('line')
            .data(networkData.links)
            .enter()
            .append('line')
            .attr('stroke', getLinkColor)
            .attr('stroke-width', getLinkWidth)
            .attr('stroke-opacity', 0.6)
            .attr('class', 'network-link');
        
        // Draw nodes
        const node = g.append('g')
            .selectAll('circle')
            .data(networkData.nodes)
            .enter()
            .append('circle')
            .attr('r', d => d.size)
            .attr('fill', getNodeColor)
            .attr('stroke', '#fff')
            .attr('stroke-width', 1.5)
            .attr('cursor', 'pointer')
            .attr('class', 'network-node')
            .attr('data-rmsd', d => d.rmsd || 0)
            .call(d3.drag()
                .on('start', dragstarted)
                .on('drag', dragged)
                .on('end', dragended));
        
        // Add labels for protein sites
        const label = g.append('g')
            .selectAll('text')
            .data(networkData.nodes.filter(d => d.type === 'protein'))
            .enter()
            .append('text')
            .attr('text-anchor', 'middle')
            .attr('dy', d => -d.size - 5)
            .style('font-size', '10px')
            .style('font-weight', 'bold')
            .style('pointer-events', 'none')
            .style('fill', '#333')
            .attr('class', 'network-node-label')
            .text(d => d.name);
        
        // Function to update info panel
        function updateInfoPanel(d) {
            let content = '';
            if (d.type === 'protein') {
                content = `
                    <h6 class="border-bottom pb-2 mb-2">${d.name} - ${d.uniprot}</h6>
                    <p><strong>Site Type:</strong> ${d.siteType}</p>
                    <p><strong>Known Site:</strong> ${d.isKnown ? 'Yes' : 'No'}</p>
                    ${d.meanPlddt ? `<p><strong>Mean pLDDT:</strong> ${d.meanPlddt}</p>` : ''}
                    ${d.nearbyCount ? `<p><strong>Nearby Residues:</strong> ${d.nearbyCount}</p>` : ''}
                    ${d.motif ? `<p><strong>Motif:</strong> <code>${d.motif}</code></p>` : ''}
                    <div class="d-grid gap-2 mt-3">
                        <a href="/site/${d.uniprot}/${d.name}" class="btn btn-sm btn-primary">View Site Details</a>
                    </div>
                `;
            } else {
                content = `
                    <h6 class="border-bottom pb-2 mb-2">${d.name} - ${d.uniprot}</h6>
                    <p><strong>Site Type:</strong> ${d.siteType}</p>
                    <p><strong>RMSD:</strong> ${d.rmsd ? d.rmsd.toFixed(2) + ' Å' : 'N/A'}</p>
                    <div class="d-grid gap-2 mt-3">
                        <a href="https://www.uniprot.org/uniprotkb/${d.uniprot}" class="btn btn-sm btn-outline-primary" target="_blank">View on UniProt</a>
                        <a href="/site/${d.uniprot}/${d.name}" class="btn btn-sm btn-primary">View Site Details</a>
                    </div>
                `;
            }
            
            infoPanel.innerHTML = content;
        }
        
        // Node interactions
        node
            .on('mouseover', function(event, d) {
                // Highlight node on hover
                d3.select(this)
                    .transition()
                    .duration(200)
                    .attr('r', d.size * 1.4);
                
                // Update info panel
                updateInfoPanel(d);
            })
            .on('mouseout', function(event, d) {
                // Restore node size
                d3.select(this)
                    .transition()
                    .duration(200)
                    .attr('r', d.size);
            })
            .on('click', function(event, d) {
                // Navigate to site details
                window.location.href = `/site/${d.uniprot}/${d.name}`;
            });
        
        // Tick function for the simulation
        networkSimulation.on('tick', () => {
            link
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y);
            
            node
                .attr('cx', d => d.x)
                .attr('cy', d => d.y);
            
            label
                .attr('x', d => d.x)
                .attr('y', d => d.y);
        });
        
        // Store nodes and links in global scope for filter function
        window.networkNodes = node;
        window.networkLinks = link;
        window.networkLabels = label;
        
        // Drag functions
        function dragstarted(event, d) {
            if (!event.active) networkSimulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
        }
        
        function dragged(event, d) {
            d.fx = event.x;
            d.fy = event.y;
        }
        
        function dragended(event, d) {
            if (!event.active) networkSimulation.alphaTarget(0);
            d.fx = null;
            d.fy = null;
        }
    }
    
    function extractNetworkDataFromDOM(proteinUniprotId) {
        try {
            // Create arrays for nodes and links
            const nodes = [];
            const links = [];
            const nodeMap = new Map(); // Track unique nodes
            
            // Get table with phosphosite data
            const phosphositeTable = document.querySelector('table.phosphosite-table');
            if (!phosphositeTable) {
                console.error("Phosphosite table not found");
                return null;
            }
            
            // Process rows to get phosphosites
            const rows = phosphositeTable.querySelectorAll('tbody tr');
            console.log(`Found ${rows.length} phosphosite rows in the table`);
            
            rows.forEach((row, index) => {
                // Get site info
                const siteCell = row.querySelector('td:first-child');
                if (!siteCell) return;
                
                const siteLink = siteCell.querySelector('a');
                const siteName = siteLink ? siteLink.textContent.trim() : siteCell.textContent.trim();
                
                if (!siteName) return;
                
                // Get site attributes
                const knownAttr = row.getAttribute('data-known');
                console.log(`Row ${index}, Site: ${siteName}, data-known attribute: "${knownAttr}"`);
                
                // Case-insensitive check for true/false strings, and also handle "1" or "yes" values
                const isKnown = knownAttr === 'true' || knownAttr === 'True' || knownAttr === 'TRUE' || 
                            knownAttr === '1' || knownAttr === 'yes' || knownAttr === 'Yes';
                
                const siteType = row.getAttribute('data-type') || siteName[0];
                const resno = row.getAttribute('data-resno') || 
                            (siteName.match(/\d+/) ? siteName.match(/\d+/)[0] : '0');
                
                // Get other metrics if available
                const meanPlddt = row.getAttribute('data-plddt') || 
                                row.querySelector('td:nth-child(3)') ? 
                                row.querySelector('td:nth-child(3)').textContent.trim() : '';
                const nearbyCount = row.getAttribute('data-nearby') || 
                                row.querySelector('td:nth-child(5)') ? 
                                row.querySelector('td:nth-child(5)').textContent.trim() : '';
                
                // Check if the row has a "Yes" in the Known column
                const knownCell = row.querySelector('td:nth-child(7)'); // 7th column is "Known"
                const knownText = knownCell ? knownCell.textContent.trim() : '';
                const isKnownByText = knownText === 'Yes';
                
                // Use either attribute or text determination for isKnown
                const finalIsKnown = isKnown || isKnownByText;
                
                console.log(`Site ${siteName} isKnown: ${finalIsKnown} (attribute: ${isKnown}, text: ${isKnownByText})`);
                
                // Get motif
                const motifCell = row.querySelector('td:nth-child(2)');
                const motif = motifCell ? motifCell.textContent.trim() : '';
                
                // Create node ID
                const nodeId = `${proteinUniprotId}_${siteName}`;
                
                // Add node if not already in the map
                if (!nodeMap.has(nodeId)) {
                    const node = {
                        id: nodeId,
                        name: siteName,
                        uniprot: proteinUniprotId,
                        type: 'protein',
                        isKnown: finalIsKnown,
                        siteType: siteType,
                        resno: resno,
                        meanPlddt: meanPlddt,
                        nearbyCount: nearbyCount,
                        motif: motif,
                        size: 10
                    };
                    
                    nodes.push(node);
                    nodeMap.set(nodeId, node);
                }
            });
            
            console.log(`Created ${nodes.length} protein nodes, known sites: ${nodes.filter(n => n.isKnown).length}`);
            
            // Find match tables
            document.querySelectorAll('.match-card').forEach(matchCard => {
                // Get the site name from the header
                const header = matchCard.querySelector('.card-header h5');
                if (!header) return;
                
                const headerText = header.textContent.trim();
                const siteMatch = headerText.match(/Site: ([^ ]+) Matches/);
                if (!siteMatch) return;
                
                const siteName = siteMatch[1];
                const sourceNodeId = `${proteinUniprotId}_${siteName}`;
                
                // Make sure source node exists
                if (!nodeMap.has(sourceNodeId)) return;
                
                // Get match table
                const matchTable = matchCard.querySelector('table');
                if (!matchTable) return;
                
                // Process match rows
                const matchRows = matchTable.querySelectorAll('tbody tr');
                matchRows.forEach(matchRow => {
                    const cells = matchRow.querySelectorAll('td');
                    if (cells.length < 3) return;
                    
                    // Get target info
                    const targetUniprotCell = cells[0];
                    const targetUniprotLink = targetUniprotCell.querySelector('a');
                    const targetUniprot = targetUniprotLink ? 
                                        targetUniprotLink.textContent.trim() : 
                                        targetUniprotCell.textContent.trim();
                    
                    const targetSite = cells[1].textContent.trim();
                    const rmsd = parseFloat(cells[2].textContent.trim());
                    
                    // Skip if RMSD is above threshold (initial filter, can be changed later by UI)
                    const currentThreshold = parseFloat(document.getElementById('rmsd-filter').value);
                    if (rmsd > currentThreshold) return;
                    
                    // Create target node ID
                    const targetNodeId = `${targetUniprot}_${targetSite}`;
                    
                    // Skip self-references
                    if (sourceNodeId === targetNodeId) return;
                    
                    // Add target node if doesn't exist
                    if (!nodeMap.has(targetNodeId)) {
                        const targetNode = {
                            id: targetNodeId,
                            name: targetSite,
                            uniprot: targetUniprot,
                            type: 'match',
                            isKnown: false,
                            siteType: targetSite[0],
                            rmsd: rmsd,
                            size: 8
                        };
                        
                        nodes.push(targetNode);
                        nodeMap.set(targetNodeId, targetNode);
                    }
                    
                    // Add link
                    links.push({
                        source: sourceNodeId,
                        target: targetNodeId,
                        rmsd: rmsd
                    });
                });
            });
            
            // Return network data if we have nodes
            if (nodes.length > 0) {
                console.log(`Final network: ${nodes.length} nodes, ${links.length} links`);
                return { nodes, links };
            }
            
            return null;
        } catch (error) {
            console.error("Error extracting network data:", error);
            return null;
        }
    }
    
    // Function to filter network by RMSD threshold
    function updateNetworkFilter() {
        if (!window.networkNodes || !window.networkLinks) return;
        
        const threshold = parseFloat(document.getElementById('rmsd-filter').value);
        
        // Filter links by RMSD
        window.networkLinks.style('display', function(d) {
            return d.rmsd <= threshold ? null : 'none';
        });
        
        // Show nodes only if they have visible connections
        window.networkNodes.style('display', function(d) {
            // Always show protein sites
            if (d.type === 'protein') return null;
            
            // For other nodes, check if they have any visible connections
            const hasVisibleConnection = window.networkLinks.data().some(link => 
                (link.source.id === d.id || link.target.id === d.id) && link.rmsd <= threshold
            );
            
            return hasVisibleConnection ? null : 'none';
        });
        
        // Update labels visibility to match nodes
        if (window.networkLabels) {
            window.networkLabels.style('display', function(d) {
                // Get the corresponding node
                const node = window.networkNodes.filter(n => n.__data__.id === d.id).node();
                if (node) {
                    // Check if the node is visible
                    return window.getComputedStyle(node).display !== 'none' ? null : 'none';
                }
                return null;
            });
        }
    }
    </script>
    """
    
    return html