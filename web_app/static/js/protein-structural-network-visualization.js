// Save this to: web_app/static/js/protein-structural-network-visualization.js

/**
 * Phosphosite Structural Network Visualization
 * Creates an interactive visualization of phosphosite structural similarity networks
 * - Filters out phosphosites with Tyrosine (Y)
 * - Displays nodes with known kinases in purple
 * - Shows structural similarity relationships via RMSD
 */

function setupPhosphositeStructuralNetwork(proteinUniprotId) {
    console.log('Initializing structural network visualization for', proteinUniprotId);
    
    const networkContainer = document.getElementById('network-container');
    if (!networkContainer) {
        console.error("Structural network container not found");
        return;
    }
    
    // Reset the network container
    networkContainer.innerHTML = '';
    
    // Create information panel
    const infoPanel = document.createElement('div');
    infoPanel.className = 'node-info-panel';
    infoPanel.style.position = 'absolute';
    infoPanel.style.top = '10px';
    infoPanel.style.right = '10px';
    infoPanel.style.width = '280px';
    infoPanel.style.backgroundColor = 'white';
    infoPanel.style.border = '1px solid #ddd';
    infoPanel.style.borderRadius = '5px';
    infoPanel.style.padding = '10px';
    infoPanel.style.boxShadow = '0 0 10px rgba(0,0,0,0.1)';
    infoPanel.style.zIndex = '100';
    infoPanel.style.fontSize = '0.9rem';
    infoPanel.style.maxHeight = '400px';
    infoPanel.style.overflowY = 'auto';
    infoPanel.innerHTML = '<p class="text-center"><em>Hover over a node to see details</em></p>';
    networkContainer.appendChild(infoPanel);
    
    // Extract network data from the page DOM
    const networkData = extractStructuralNetworkData(proteinUniprotId);
    
    if (!networkData || networkData.nodes.length === 0) {
        networkContainer.innerHTML = '<div class="alert alert-info m-3 mt-5">No phosphosite structural network data available. This could be because no structural matches were found with RMSD < 2.0.</div>';
        return;
    }
    
    // Create a map of node IDs for quick lookup
    const nodeIdMap = new Map();
    networkData.nodes.forEach(node => {
        nodeIdMap.set(node.id, true);
    });
    
    // Filter links to include only those where both source and target nodes exist
    const validLinks = networkData.links.filter(link => {
        // Check if both source and target nodes exist
        const sourceExists = nodeIdMap.has(link.source) || 
                            (typeof link.source === 'object' && link.source.id && nodeIdMap.has(link.source.id));
        const targetExists = nodeIdMap.has(link.target) || 
                            (typeof link.target === 'object' && link.target.id && nodeIdMap.has(link.target.id));
        
        return sourceExists && targetExists;
    });
    
    console.log(`Creating structural network with ${networkData.nodes.length} nodes and ${validLinks.length} links (${networkData.links.length - validLinks.length} invalid links removed)`);
    
    // Update the links in the network data
    networkData.links = validLinks;
    
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
        // First check if node has a known kinase (regardless of type)
        if (d.known_kinase) {
            return '#9C27B0'; // Purple for any site with a known kinase
        }
        
        // Then check node type
        if (d.type === 'protein') {
            // Check is_known_phosphosite (from database)
            if (d.is_known_phosphosite === true || 
                d.is_known_phosphosite === 1 || 
                d.is_known_phosphosite === "1" || 
                d.is_known_phosphosite === 1.0 ||
                d.isKnown === true) {
                return '#4CAF50'; // Green for known protein sites
            }
            
            return '#FF9800'; // Orange for unknown protein sites
        }
        
        // For match nodes without a known kinase
        return '#9E9E9E'; // Gray for structurally similar sites without known kinases
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
    
    // Setup force simulation with the validated network data
    networkSimulation = d3.forceSimulation(networkData.nodes)
        .force('link', d3.forceLink(networkData.links)
            .id(d => d.id)
            .distance(d => d.rmsd ? d.rmsd * 30 : 60))
        .force('charge', d3.forceManyBody()
            .strength(d => d.type === 'protein' ? -50 : -25))
        .force('center', d3.forceCenter(width / 2, height / 2))
        .force('collision', d3.forceCollide().radius(d => d.size + 3))
        .force('x', d3.forceX(width / 2).strength(0.1))
        .force('y', d3.forceY(height / 2).strength(0.1));
    
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
    // Function to update info panel
    function updateInfoPanel(d) {
        let content = '';
        if (d.type === 'protein') {
            content = `
                <h6 class="border-bottom pb-2 mb-2">${d.name} - ${d.uniprot}</h6>
                <p><strong>Site Type:</strong> ${d.siteType}</p>
                <p><strong>Known Site:</strong> ${d.isKnown ? 'Yes' : 'No'}</p>
                ${d.known_kinase ? `<p><strong>Known Kinase${d.known_kinase.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(d.known_kinase)}</p>` : ''}
                ${d.meanPlddt ? `<p><strong>Mean pLDDT:</strong> ${d.meanPlddt}</p>` : ''}
                ${d.nearbyCount ? `<p><strong>Nearby Residues:</strong> ${d.nearbyCount}</p>` : ''}
                ${d.surface_accessibility ? `<p><strong>Surface Access:</strong> ${typeof d.surface_accessibility === 'number' ? d.surface_accessibility.toFixed(1) : d.surface_accessibility}%</p>` : ''}
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
                ${d.known_kinase ? `<p><strong>Known Kinase${d.known_kinase.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(d.known_kinase)}</p>` : ''}
                ${d.mean_plddt || d.site_plddt ? `<p><strong>Mean pLDDT:</strong> ${d.mean_plddt || d.site_plddt}</p>` : ''}
                ${d.nearby_count ? `<p><strong>Nearby Residues:</strong> ${d.nearby_count}</p>` : ''}
                ${d.surface_accessibility ? `<p><strong>Surface Access:</strong> ${typeof d.surface_accessibility === 'number' ? d.surface_accessibility.toFixed(1) : d.surface_accessibility}%</p>` : ''}
                ${d.motif ? `<p><strong>Motif:</strong> <code>${d.motif}</code></p>` : ''}
                <div class="d-grid gap-2 mt-3">
                    <a href="https://www.uniprot.org/uniprotkb/${d.uniprot}" class="btn btn-sm btn-outline-primary" target="_blank">View on UniProt</a>
                    <a href="/site/${d.uniprot}/${d.name}" class="btn btn-sm btn-primary">View Site Details</a>
                </div>
            `;
        }
        
        infoPanel.innerHTML = content;
    }

    // Helper function to format kinase list with badges
    function formatKinaseList(kinaseString) {
        if (!kinaseString) return '';
        
        const kinases = kinaseString.split(',').map(k => k.trim());
        return kinases.map(kinase => 
            `<span class="badge bg-primary me-1">${kinase}</span>`
        ).join(' ');
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

function extractStructuralNetworkData(proteinUniprotId) {
    try {
        console.log("Extracting structural network data for", proteinUniprotId);
        
        // Initialize node and link collections
        const nodes = [];
        const links = [];
        const nodeMap = new Map(); // Track unique nodes
        
        // Step 1: Try to get phosphosites data from the hidden data element
        const phosphositesDataElement = document.getElementById('phosphosites-data');
        let phosphositesData = [];
        
        if (phosphositesDataElement) {
            try {
                const phosphositesJson = phosphositesDataElement.getAttribute('data-sites');
                if (phosphositesJson) {
                    phosphositesData = JSON.parse(cleanJson(phosphositesJson));
                    console.log(`Found ${phosphositesData.length} phosphosites from hidden data element`);
                }
            } catch (e) {
                console.error("Error parsing phosphosites data:", e);
            }
        }
        
        // Step 2: If no hidden data, try to extract from the table
        if (phosphositesData.length === 0) {
            const phosphositeTable = document.querySelector('table.phosphosite-table');
            if (phosphositeTable) {
                const rows = phosphositeTable.querySelectorAll('tbody tr');
                
                if (rows.length > 0) {
                    console.log(`Found ${rows.length} rows in phosphosite table`);
                    
                    rows.forEach((row, index) => {
                        const siteCell = row.querySelector('td:first-child');
                        if (!siteCell) return;
                        
                        // Get site name either from the link or directly from cell
                        const siteLink = siteCell.querySelector('a');
                        const siteName = siteLink ? siteLink.textContent.trim() : siteCell.textContent.trim();
                        if (!siteName) return;
                        
                        // Skip Tyrosine (Y) sites
                        if (siteName[0] === 'Y') {
                            console.log(`Skipping Tyrosine site: ${siteName}`);
                            return;
                        }
                        
                        // Extract site data
                        const knownAttr = row.getAttribute('data-known');
                        const isKnown = knownAttr === 'true' || knownAttr === 'True' || knownAttr === 'TRUE' || 
                                        knownAttr === '1' || knownAttr === 'yes' || knownAttr === 'Yes' ||
                                        row.querySelector('td:nth-child(7)').textContent.trim() === 'Yes';
                                        
                        const siteType = row.getAttribute('data-type') || siteName[0];
                        const resno = row.getAttribute('data-resno') || 
                                    (siteName.match(/\d+/) ? siteName.match(/\d+/)[0] : '0');
                        
                        // Get motif if available
                        let motif = '';
                        const motifCell = row.querySelector('td:nth-child(2) code');
                        if (motifCell) motif = motifCell.textContent.trim();
                        
                        // Get mean pLDDT if available
                        let meanPlddt = '';
                        const pLDDTCell = row.querySelector('td:nth-child(3)');
                        if (pLDDTCell) meanPlddt = pLDDTCell.textContent.trim();
                        
                        // Get nearby count if available
                        let nearbyCount = '';
                        const nearbyCell = row.querySelector('td:nth-child(5)');
                        if (nearbyCell) nearbyCount = nearbyCell.textContent.trim();
                        
                        // Check for kinase information
                        let knownKinase = null;
                        const kinaseCell = row.querySelector('td.kinase-col, td:nth-child(8)');
                        if (kinaseCell && kinaseCell.textContent.trim() && 
                            kinaseCell.textContent.trim() !== '—' && 
                            kinaseCell.textContent.trim() !== '-') {
                            knownKinase = kinaseCell.textContent.trim();
                        }
                        
                        // Also check data attributes for kinases
                        const dataKinase = row.getAttribute('data-kinase');
                        if (!knownKinase && dataKinase && dataKinase !== 'null' && dataKinase !== '') {
                            knownKinase = dataKinase;
                        }
                        
                        // Create site data object
                        phosphositesData.push({
                            site: siteName,
                            resno: parseInt(resno),
                            siteType: siteType,
                            is_known: isKnown,
                            motif: motif,
                            meanPlddt: meanPlddt,
                            nearbyCount: nearbyCount,
                            known_kinase: knownKinase
                        });
                    });
                    
                    console.log(`Extracted ${phosphositesData.length} phosphosites from table`);
                }
            }
        }
        
        // Step 3: Add all protein site nodes (excluding Tyrosine sites)
        for (const site of phosphositesData) {
            // Skip Tyrosine (Y) sites
            if (site.site && site.site[0] === 'Y') {
                console.log(`Skipping Tyrosine site: ${site.site}`);
                continue;
            }
            
            // Create node ID - IMPORTANT: Using resno without letter
            const resno = site.resno || parseInt(site.site.substring(1));
            const nodeId = `${proteinUniprotId}_${resno}`;
            
            // Skip if already added
            if (nodeMap.has(nodeId)) continue;
            
            // Extract known kinase information if available
            const knownKinase = site.known_kinase || site.knownKinase;
            
            // Create node with all available data
            const node = {
                id: nodeId,
                name: site.site,
                uniprot: proteinUniprotId,
                type: 'protein',
                isKnown: site.is_known === true || site.is_known === 1,
                is_known_phosphosite: site.is_known_phosphosite || (site.is_known ? 1.0 : 0.0),
                siteType: site.siteType || site.site[0],
                motif: site.motif || '',
                known_kinase: knownKinase,
                meanPlddt: site.meanPlddt || site.mean_plddt || '',
                nearbyCount: site.nearbyCount || site.nearby_count || '',
                size: 10
            };
            
            // Add node to collections
            nodes.push(node);
            nodeMap.set(nodeId, node);
        }
        
        console.log(`Created ${nodes.length} protein nodes`);
        
        // Step 4: Get structural matches and create links
        
        // First try to get from hidden data element
        const matchesDataElement = document.getElementById('structural-match-data');
        if (matchesDataElement) {
            try {
                const matchesJson = matchesDataElement.getAttribute('data-matches');
                if (matchesJson) {
                    const matchesData = JSON.parse(cleanJson(matchesJson));
                    
                    if (matchesData) {
                        console.log(`Processing matches from data element for ${Object.keys(matchesData).length} sites`);
                        
                        // Process each site's matches
                        for (const [siteName, matches] of Object.entries(matchesData)) {
                            // Skip if site name starts with 'Y' (Tyrosine)
                            if (siteName[0] === 'Y') continue;
                            
                            // Skip sites without matches
                            if (!Array.isArray(matches) || matches.length === 0) continue;
                            
                            // Extract site number without letter for node ID
                            const siteResno = parseInt(siteName.replace(/[A-Z]/g, ''));
                            // Create source node ID with just the number
                            const sourceNodeId = `${proteinUniprotId}_${siteResno}`;
                            
                            for (const match of matches) {
                                // Skip matches with Y sites
                                if (match.target_site && match.target_site[0] === 'Y') continue;
                                
                                // Get target info
                                const targetUniprot = match.target_uniprot;
                                const targetSite = match.target_site;
                                const rmsd = match.rmsd;
                                
                                // Skip if missing essential info
                                if (!targetUniprot || !targetSite || !rmsd) continue;
                                
                                // Extract target residue number without letter
                                const targetResno = parseInt(targetSite.replace(/[A-Z]/g, ''));
                                // Create target node ID with just the number
                                const targetNodeId = `${targetUniprot}_${targetResno}`;
                                
                                // Skip self-references
                                if (sourceNodeId === targetNodeId) continue;
                                
                                // Extract known kinase information
                                const knownKinase = match.known_kinase;
                                
                                // Add target node if not already present
                                if (!nodeMap.has(targetNodeId)) {
                                    const targetNode = {
                                        id: targetNodeId,
                                        name: targetSite,
                                        uniprot: targetUniprot,
                                        type: 'match',
                                        isKnown: false,  // Default for matches
                                        siteType: targetSite[0],
                                        rmsd: rmsd,
                                        known_kinase: knownKinase,
                                        motif: match.motif || '',
                                        size: 8  // Slightly smaller for match sites
                                    };
                                    
                                    nodes.push(targetNode);
                                    nodeMap.set(targetNodeId, targetNode);
                                }
                                
                                // Create link
                                links.push({
                                    source: sourceNodeId,
                                    target: targetNodeId,
                                    rmsd: rmsd
                                });
                            }
                        }
                    }
                }
            } catch (e) {
                console.error("Error processing structural match data:", e);
            }
        }
        
        // If we didn't find matches from the data element, try from DOM elements
        if (links.length === 0) {
            document.querySelectorAll('.match-card').forEach(matchCard => {
                const header = matchCard.querySelector('.card-header h5');
                if (!header) return;
                
                const headerText = header.textContent.trim();
                const siteMatch = headerText.match(/Site: ([^ ]+) Matches/);
                if (!siteMatch) return;
                
                const siteName = siteMatch[1];
                
                // Skip Y sites
                if (siteName[0] === 'Y') return;
                
                // Extract site number without letter
                const siteResno = parseInt(siteName.replace(/[A-Z]/g, ''));
                // Create source node ID with just the number
                const sourceNodeId = `${proteinUniprotId}_${siteResno}`;
                
                // Make sure source node exists
                if (!nodeMap.has(sourceNodeId)) {
                    // If the source node doesn't exist, try to create it
                    console.log(`Source node ${sourceNodeId} not found, creating it`);
                    
                    nodes.push({
                        id: sourceNodeId,
                        name: siteName,
                        uniprot: proteinUniprotId,
                        type: 'protein',
                        isKnown: false,
                        siteType: siteName[0],
                        size: 10
                    });
                    
                    nodeMap.set(sourceNodeId, { id: sourceNodeId });
                }
                
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
                    
                    // Skip Y sites
                    if (targetSite[0] === 'Y') return;
                    
                    let rmsd;
                    try {
                        rmsd = parseFloat(cells[2].textContent.trim());
                    } catch (e) {
                        console.warn(`Failed to parse RMSD: ${cells[2].textContent.trim()}`);
                        return;
                    }
                    
                    if (isNaN(rmsd)) {
                        console.warn(`Invalid RMSD value: ${cells[2].textContent.trim()}`);
                        return;
                    }
                    
                    // Skip if RMSD is above threshold (initial filter, can be changed later by UI)
                    const rmsdFilter = document.getElementById('rmsd-filter');
                    const currentThreshold = rmsdFilter ? parseFloat(rmsdFilter.value) : 2.0;
                    if (rmsd > currentThreshold) return;
                    
                    // Extract target residue number without letter
                    const targetResno = parseInt(targetSite.replace(/[A-Z]/g, ''));
                    // Create target node ID with just the number
                    const targetNodeId = `${targetUniprot}_${targetResno}`;
                    
                    // Skip self-references
                    if (sourceNodeId === targetNodeId) return;
                    
                    // Check for known kinase information
                    let knownKinase = null;
                    
                    // Look for kinase information in additional cells
                    if (cells.length >= 4) {
                        const kinaseCell = cells[3];
                        if (kinaseCell && kinaseCell.textContent.trim() && 
                            kinaseCell.textContent.trim() !== 'None' && 
                            kinaseCell.textContent.trim() !== 'Unknown') {
                            knownKinase = kinaseCell.textContent.trim();
                        }
                    }
                    
                    // Add target node if not already present
                    if (!nodeMap.has(targetNodeId)) {
                        const targetNode = {
                            id: targetNodeId,
                            name: targetSite,
                            uniprot: targetUniprot,
                            type: 'match',
                            isKnown: false,  // Assume false for matches
                            siteType: targetSite[0],
                            rmsd: rmsd,
                            known_kinase: knownKinase,
                            size: 8
                        };
                        
                        nodes.push(targetNode);
                        nodeMap.set(targetNodeId, { id: targetNodeId });
                    }
                    
                    // Add link
                    links.push({
                        source: sourceNodeId,
                        target: targetNodeId,
                        rmsd: rmsd
                    });
                });
            });
        }
        
        // Return network data if we have nodes
        if (nodes.length > 0) {
            console.log(`Final structural network: ${nodes.length} nodes, ${links.length} links`);
            return { nodes, links };
        }
        
        return null;
    } catch (error) {
        console.error("Error extracting structural network data:", error);
        return null;
    }
}

// Function to filter network by RMSD threshold
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
        const hasVisibleConnection = window.networkLinks.data().some(link => {
            const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
            const targetId = typeof link.target === 'object' ? link.target.id : link.target;
            
            return ((sourceId === d.id || targetId === d.id) && link.rmsd <= threshold);
        });
        
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

// Helper function to clean JSON before parsing
function cleanJson(jsonString) {
    if (!jsonString) return "{}";
    
    return jsonString
        .replace(/NaN/g, 'null')
        .replace(/Infinity/g, 'null')
        .replace(/undefined/g, 'null')
        .replace(/\bnan\b/g, 'null')
        .replace(/\binfinity\b/g, 'null');
}

// Extract kinase information from site data
function extractKinaseInfo(data) {
    // Create an array to collect all kinases
    const kinases = [];
    
    // Check all possible locations for kinase information
    
    // Direct known_kinase field
    if (data.known_kinase && typeof data.known_kinase === 'string' && data.known_kinase !== 'unlabeled') {
        kinases.push(data.known_kinase);
    }
    
    // PhosphositePlus format KINASE_1 through KINASE_5
    for (let i = 1; i <= 5; i++) {
        const kinaseField = `KINASE_${i}`;
        if (data[kinaseField] && typeof data[kinaseField] === 'string' && data[kinaseField] !== 'unlabeled') {
            kinases.push(data[kinaseField]);
        }
    }
    
    // Target known kinase field
    if (data.target_known_kinase && typeof data.target_known_kinase === 'string' && data.target_known_kinase !== 'unlabeled') {
        kinases.push(data.target_known_kinase);
    }
    
    // Check for kinases array
    if (data.kinases && Array.isArray(data.kinases) && data.kinases.length > 0) {
        // Add all valid kinases from the array
        for (const kinase of data.kinases) {
            if (kinase && typeof kinase === 'string' && kinase !== 'unlabeled') {
                kinases.push(kinase);
            }
        }
    }
    
    // Check if top_kinases is available with scores
    if (data.top_kinases && Array.isArray(data.top_kinases) && data.top_kinases.length > 0) {
        // Add top kinase
        const topKinase = data.top_kinases[0];
        const kinaseName = topKinase.kinase || topKinase.name;
        if (kinaseName) {
            kinases.push(kinaseName);
        }
    }
    
    // Remove duplicates and return as a comma-separated string
    return [...new Set(kinases)].join(', ') || null;
}