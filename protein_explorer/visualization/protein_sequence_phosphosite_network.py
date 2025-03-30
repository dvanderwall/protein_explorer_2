# In protein_explorer/visualization/protein_sequence_phosphosite_network.py

def create_sequence_network_visualization(protein_uniprot_id, phosphosites=None, sequence_matches=None):
    """
    Create a network visualization of phosphosites and their sequence similarity matches.
    Similar to the structural network, but based on sequence similarity.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries (optional)
        sequence_matches: Dictionary mapping site names to lists of match dictionaries (optional)
        
    Returns:
        HTML string with network visualization
    """
    # Check if we have any sequence matches to display
    has_matches = False
    if sequence_matches:
        for matches in sequence_matches.values():
            if matches:
                has_matches = True
                break
    
    if not has_matches:
        return """
        <div class="card mb-4">
            <div class="card-header">
                <h5 class="mb-0">Phosphosite Sequence Similarity Network</h5>
            </div>
            <div class="card-body">
                <div class="alert alert-info">
                    No sequence similarity matches found for this protein. Try adjusting similarity thresholds or analyzing a different protein.
                </div>
            </div>
        </div>
        """
    
    # Create network visualization HTML
    html = """
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Phosphosite Sequence Similarity Network</h5>
        </div>
        <div class="card-body">
            <div id="similarity-filter-container" class="mb-3">
                <label for="similarity-filter" class="form-label">
                    Similarity Threshold: <span id="similarity-value">0.6</span>
                </label>
                <input 
                    type="range" 
                    class="form-range" 
                    id="similarity-filter" 
                    min="0.4" 
                    max="0.9" 
                    step="0.05" 
                    value="0.6"
                    oninput="document.getElementById('similarity-value').textContent = this.value; updateSequenceNetworkFilter();"
                >
                <div class="d-flex justify-content-between">
                    <small>0.4 (Less similar)</small>
                    <small>0.9 (Very similar)</small>
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
                        <div style="width: 16px; height: 16px; background-color: #9C27B0; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Sequence-similar sites</span>
                    </div>
                </div>
            </div>
            
            <div id="sequence-network-container" style="height: 500px; width: 100%; position: relative; border: 1px solid #ddd; border-radius: 5px;"></div>
            
            <p class="text-muted mt-3 mb-0">
                <small>
                    This network shows the sequence relationships between phosphosites in this protein and
                    similar sites in other proteins. Edges represent sequence similarity above
                    the threshold. Hover over nodes for details and click to view the site page.
                </small>
            </p>
        </div>
    </div>

    <!-- Network Visualization Script -->
    <script>
    // Use an IIFE to avoid global namespace pollution
    (function() {
        // Add a load event listener to make sure D3 is available
        document.addEventListener('DOMContentLoaded', function() {
            // Check if D3 is loaded
            if (typeof d3 === 'undefined') {
                console.error('D3.js is not loaded! Loading it now...');
                var script = document.createElement('script');
                script.src = 'https://d3js.org/d3.v7.min.js';
                script.onload = function() {
                    console.log('D3.js loaded successfully, initializing sequence network');
                    setupSequenceNetwork();
                };
                document.head.appendChild(script);
            } else {
                console.log('D3.js is already loaded, initializing sequence network');
                // Wait a bit to ensure everything is loaded
                setTimeout(setupSequenceNetwork, 500);
            }
        });
        
        function setupSequenceNetwork() {
            console.log('Setting up sequence network visualization');
            const networkContainer = document.getElementById('sequence-network-container');
            if (!networkContainer) {
                console.error("Sequence network container not found");
                return;
            }
            
            // Reset the network container
            networkContainer.innerHTML = '';
            
            // Create information panel
            const infoPanel = document.createElement('div');
            infoPanel.className = 'sequence-node-info-panel';
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
            
            // Extract network data from the page DOM for sequence similarity
            const networkData = extractSequenceNetworkDataFromDOM(proteinUniprotId);
            
            if (!networkData || networkData.nodes.length === 0) {
                networkContainer.innerHTML = '<div class="alert alert-info m-3 mt-5">No phosphosite sequence similarity network data available. This could be because no sequence similarity matches were found with similarity > 0.6.</div>';
                return;
            }
            
            console.log(`Creating sequence network with ${networkData.nodes.length} nodes and ${networkData.links.length} links`);
            
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
                return '#9C27B0'; // Purple for sequence-similar sites
            }
            
            function getLinkColor(d) {
                if (d.similarity > 0.8) return '#4CAF50'; // Green for very high similarity
                if (d.similarity > 0.7) return '#8BC34A'; // Light green for high similarity
                if (d.similarity > 0.6) return '#CDDC39'; // Lime for medium similarity
                return '#FFC107'; // Yellow/amber for lower similarity
            }
            
            function getLinkWidth(d) {
                // Scale link width based on similarity - thicker for higher similarity
                return Math.max(1, 4 * d.similarity);
            }
            
            // Setup force simulation
            const sequenceNetworkSimulation = d3.forceSimulation(networkData.nodes)
                .force('link', d3.forceLink(networkData.links)
                    .id(d => d.id)
                    .distance(d => 150 * (1 - (d.similarity || 0.5))))
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
                .attr('class', 'sequence-network-link');
            
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
                .attr('class', 'sequence-network-node')
                .attr('data-similarity', d => d.similarity || 0)
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
                .attr('class', 'sequence-network-node-label')
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
                        <p><strong>Similarity:</strong> ${d.similarity ? (d.similarity * 100).toFixed(1) + '%' : 'N/A'}</p>
                        ${d.motif ? `<p><strong>Motif:</strong> <code>${d.motif}</code></p>` : ''}
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
            sequenceNetworkSimulation.on('tick', () => {
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
            
            // Store nodes and links for filter function
            window.sequenceNetworkNodes = node;
            window.sequenceNetworkLinks = link;
            window.sequenceNetworkLabels = label;
            
            // Drag functions
            function dragstarted(event, d) {
                if (!event.active) sequenceNetworkSimulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            }
            
            function dragged(event, d) {
                d.fx = event.x;
                d.fy = event.y;
            }
            
            function dragended(event, d) {
                if (!event.active) sequenceNetworkSimulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            }
        }
        
        function extractSequenceNetworkDataFromDOM(proteinUniprotId) {
            try {
                // Create arrays for nodes and links
                const nodes = [];
                const links = [];
                const nodeMap = new Map(); // Track unique nodes
                
                // Get sequence match cards to extract data
                const sequenceMatchCards = document.querySelectorAll('.sequence-match-card');
                if (!sequenceMatchCards || sequenceMatchCards.length === 0) {
                    console.error("No sequence match cards found");
                    // Try to get data from the hidden data element
                    const dataElement = document.getElementById('sequence-match-data');
                    if (dataElement) {
                        try {
                            const matchesJson = dataElement.getAttribute('data-matches');
                            if (matchesJson) {
                                const matchesData = JSON.parse(matchesJson);
                                // Process the matches data if available
                                if (matchesData && matchesData.site_matches) {
                                    return buildNetworkFromData(matchesData.site_matches, proteinUniprotId);
                                }
                            }
                        } catch (e) {
                            console.error("Error parsing sequence match data:", e);
                        }
                    }
                    return null;
                }
                
                // Get phosphosite table to extract protein nodes
                const phosphositeTable = document.querySelector('table.phosphosite-table');
                if (!phosphositeTable) {
                    console.error("Phosphosite table not found");
                    return null;
                }
                
                // Process rows to get phosphosites
                const rows = phosphositeTable.querySelectorAll('tbody tr');
                rows.forEach(row => {
                    // Get site info
                    const siteCell = row.querySelector('td:first-child');
                    if (!siteCell) return;
                    
                    const siteLink = siteCell.querySelector('a');
                    const siteName = siteLink ? siteLink.textContent.trim() : siteCell.textContent.trim();
                    
                    if (!siteName) return;
                    
                    // Get site attributes
                    const isKnown = row.getAttribute('data-known') === 'true';
                    const siteType = row.getAttribute('data-type') || siteName[0];
                    const resno = row.getAttribute('data-resno') || 
                                siteName.match(/\\d+/)[0];
                    
                    // Get other metrics if available
                    const meanPlddt = row.getAttribute('data-plddt') || 
                                    row.querySelector('td:nth-child(3)') ? 
                                    row.querySelector('td:nth-child(3)').textContent.trim() : '';
                    const nearbyCount = row.getAttribute('data-nearby') || 
                                    row.querySelector('td:nth-child(5)') ? 
                                    row.querySelector('td:nth-child(5)').textContent.trim() : '';
                    
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
                            isKnown: isKnown,
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
                
                // Process sequence match cards
                sequenceMatchCards.forEach(matchCard => {
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
                        const similarity = parseFloat(cells[2].textContent.trim()) / 100; // Convert from percentage
                        
                        // Skip if similarity is below threshold (initial filter, can be changed later by UI)
                        const currentThreshold = parseFloat(document.getElementById('similarity-filter').value);
                        if (similarity < currentThreshold) return;
                        
                        // Create target node ID
                        const targetNodeId = `${targetUniprot}_${targetSite}`;
                        
                        // Skip self-references
                        if (sourceNodeId === targetNodeId) return;
                        
                        // Get motif if available
                        let motif = '';
                        if (cells.length > 3) {
                            const motifCell = cells[3];
                            motif = motifCell ? motifCell.textContent.trim() : '';
                        }
                        
                        // Add target node if doesn't exist
                        if (!nodeMap.has(targetNodeId)) {
                            const targetNode = {
                                id: targetNodeId,
                                name: targetSite,
                                uniprot: targetUniprot,
                                type: 'match',
                                isKnown: false,
                                siteType: targetSite[0],
                                similarity: similarity,
                                motif: motif,
                                size: 8
                            };
                            
                            nodes.push(targetNode);
                            nodeMap.set(targetNodeId, targetNode);
                        }
                        
                        // Add link from source to target
                        links.push({
                            source: sourceNodeId,
                            target: targetNodeId,
                            similarity: similarity
                        });
                    });
                });
                
                // Return network data if we have nodes
                if (nodes.length > 0) {
                    return { nodes, links };
                }
                
                return null;
            } catch (error) {
                console.error("Error extracting sequence network data:", error);
                return null;
            }
        }
        
        // Helper function to build network from raw data
        function buildNetworkFromData(siteMatches, proteinUniprotId) {
            const nodes = [];
            const links = [];
            const nodeMap = new Map();
            
            // First add nodes for each site in the protein
            for (const [site_id, matches] of Object.entries(siteMatches)) {
                // Parse site_id to get site name
                const parts = site_id.split('_');
                if (parts.length < 2) continue;
                
                const siteName = parts[1];
                
                // Create protein node
                const nodeId = `${proteinUniprotId}_${siteName}`;
                
                if (!nodeMap.has(nodeId)) {
                    const node = {
                        id: nodeId,
                        name: siteName,
                        uniprot: proteinUniprotId,
                        type: 'protein',
                        isKnown: false, // Default
                        siteType: siteName[0] || 'S',
                        size: 10
                    };
                    
                    nodes.push(node);
                    nodeMap.set(nodeId, node);
                }
                
                // Add match nodes and links
                for (const match of matches) {
                    const targetId = match.target_id;
                    const targetUniprot = match.target_uniprot;
                    const targetSite = match.target_site;
                    const similarity = match.similarity;
                    
                    // Skip if already added
                    if (nodeMap.has(targetId)) continue;
                    
                    // Create match node
                    const matchNode = {
                        id: targetId,
                        name: targetSite,
                        uniprot: targetUniprot,
                        type: 'match',
                        isKnown: false,
                        siteType: targetSite[0] || 'S',
                        similarity: similarity,
                        motif: match.motif,
                        size: 8
                    };
                    
                    nodes.push(matchNode);
                    nodeMap.set(targetId, matchNode);
                    
                    // Add link
                    links.push({
                        source: nodeId,
                        target: targetId,
                        similarity: similarity
                    });
                }
            }
            
            return { nodes, links };
        }
        
        // Function to filter network by similarity threshold
        window.updateSequenceNetworkFilter = function() {
            if (!window.sequenceNetworkNodes || !window.sequenceNetworkLinks) return;
            
            const threshold = parseFloat(document.getElementById('similarity-filter').value);
            
            // Filter links by similarity
            window.sequenceNetworkLinks.style('display', function(d) {
                return d.similarity >= threshold ? null : 'none';
            });
            
            // Show nodes only if they have visible connections
            window.sequenceNetworkNodes.style('display', function(d) {
                // Always show protein sites
                if (d.type === 'protein') return null;
                
                // For other nodes, check if they have any visible connections
                const hasVisibleConnection = window.sequenceNetworkLinks.data().some(link => 
                    (link.source.id === d.id || link.target.id === d.id) && link.similarity >= threshold
                );
                
                return hasVisibleConnection ? null : 'none';
            });
            
            // Update labels visibility to match nodes
            if (window.sequenceNetworkLabels) {
                window.sequenceNetworkLabels.style('display', function(d) {
                    // Get the corresponding node
                    const node = window.sequenceNetworkNodes.filter(n => n.__data__.id === d.id).node();
                    if (node) {
                        // Check if the node is visible
                        return window.getComputedStyle(node).display !== 'none' ? null : 'none';
                    }
                    return null;
                });
            }
        };
    })();
    </script>
    """
    
    return html