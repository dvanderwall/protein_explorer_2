/**
 * Sequence Network Visualization
 * Creates an interactive visualization of phosphosite sequence similarity networks
 */

function sequenceNetworkVisualization(proteinUniprotId) {
    console.log('Initializing sequence network visualization for', proteinUniprotId);
    
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
    
    // Extract network data from the page DOM
    const networkData = extractSequenceNetworkData(proteinUniprotId);
    
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
    
    console.log('Sequence network visualization setup complete');
}

// Extract network data from the page
function extractSequenceNetworkData(proteinUniprotId) {
    try {
        console.log("Extracting sequence network data for", proteinUniprotId);
        
        // Create arrays for nodes and links
        const nodes = [];
        const links = [];
        const nodeMap = new Map(); // Track unique nodes
        
        // First attempt: Try to get data from the hidden data element
        const dataElement = document.getElementById('sequence-match-data');
        if (dataElement) {
            console.log("Found sequence-match-data element, attempting to parse JSON data");
            try {
                // Get the JSON data from the data-matches attribute
                const matchesJson = dataElement.getAttribute('data-matches');
                if (matchesJson) {
                    console.log("Found matches JSON data");
                    const matchesData = JSON.parse(matchesJson);
                    
                    // Process the matches data
                    if (matchesData) {
                        console.log("Processing matches data:", Object.keys(matchesData));
                        
                        // First create nodes for each site in the protein
                        for (const [siteName, matches] of Object.entries(matchesData)) {
                            // Create protein node
                            const nodeId = `${proteinUniprotId}_${siteName}`;
                            
                            // Get site data from the phosphosite table if available
                            let isKnown = false;
                            let siteType = siteName[0] || 'S';
                            let motif = '';
                            
                            // Try to get more info from the DOM
                            const siteRow = document.querySelector(`tr[data-site="${siteName}"]`);
                            if (siteRow) {
                                isKnown = siteRow.getAttribute('data-known') === 'true';
                                siteType = siteRow.getAttribute('data-type') || siteName[0];
                                const motifCell = siteRow.querySelector('td:nth-child(2) code');
                                if (motifCell) motif = motifCell.textContent.trim();
                                
                                // Try to get additional metrics
                                const meanPlddt = siteRow.getAttribute('data-plddt') || siteRow.querySelector('td:nth-child(3)')?.textContent.trim();
                                const nearbyCount = siteRow.getAttribute('data-nearby') || siteRow.querySelector('td:nth-child(5)')?.textContent.trim();
                                
                                if (!nodeMap.has(nodeId)) {
                                    const node = {
                                        id: nodeId,
                                        name: siteName,
                                        uniprot: proteinUniprotId,
                                        type: 'protein',
                                        isKnown: isKnown,
                                        siteType: siteType,
                                        motif: motif,
                                        meanPlddt: meanPlddt,
                                        nearbyCount: nearbyCount,
                                        size: 10
                                    };
                                    
                                    nodes.push(node);
                                    nodeMap.set(nodeId, node);
                                    console.log(`Added protein node: ${nodeId}, isKnown: ${isKnown}`);
                                }
                            } else {
                                // If no row exists, just create a basic node
                                if (!nodeMap.has(nodeId)) {
                                    const node = {
                                        id: nodeId,
                                        name: siteName,
                                        uniprot: proteinUniprotId,
                                        type: 'protein',
                                        isKnown: false, // Default to false if we can't determine
                                        siteType: siteType,
                                        motif: '',
                                        size: 10
                                    };
                                    
                                    nodes.push(node);
                                    nodeMap.set(nodeId, node);
                                    console.log(`Added basic protein node: ${nodeId}`);
                                }
                            }
                            
                            // Add match nodes and links
                            for (const match of matches) {
                                // Skip if similarity is below threshold (use current slider value)
                                const currentThreshold = parseFloat(document.getElementById('similarity-filter').value || 0.6);
                                if (match.similarity < currentThreshold) continue;
                                
                                const targetId = match.target_id;
                                const targetUniprot = match.target_uniprot;
                                const targetSite = match.target_site;
                                const similarity = match.similarity;
                                
                                // Skip self-references
                                if (nodeId === targetId) continue;
                                
                                // Add target node if doesn't exist
                                if (!nodeMap.has(targetId)) {
                                    const targetNode = {
                                        id: targetId,
                                        name: targetSite,
                                        uniprot: targetUniprot,
                                        type: 'match', // Sequence match
                                        isKnown: false, // Default for matches
                                        siteType: targetSite[0] || 'S',
                                        similarity: similarity,
                                        motif: match.motif || '',
                                        size: 8
                                    };
                                    
                                    nodes.push(targetNode);
                                    nodeMap.set(targetId, targetNode);
                                    console.log(`Added match node: ${targetId}, sim: ${similarity}`);
                                }
                                
                                // Add link from protein to match
                                links.push({
                                    source: nodeId,
                                    target: targetId,
                                    similarity: similarity
                                });
                                console.log(`Added link: ${nodeId} -> ${targetId}, sim: ${similarity}`);
                            }
                        }
                        
                        console.log(`Created network with ${nodes.length} nodes and ${links.length} links from JSON data`);
                        return { nodes, links };
                    }
                }
            } catch (e) {
                console.error("Error parsing sequence match data from element:", e);
            }
        }
        
        // Fallback: Use DOM-based extraction
        // ... (rest of the extraction code similar to what we had before)
        
        // If we get here with no data, return null
        console.error("No sequence match data found");
        return null;
    } catch (error) {
        console.error("Error extracting sequence network data:", error);
        return null;
    }
}

// Function to filter network by similarity threshold
function updateSequenceNetworkFilter() {
    if (!window.sequenceNetworkNodes || !window.sequenceNetworkLinks) {
        console.error("Network elements not available for filtering");
        return;
    }
    
    console.log("Updating sequence network filter");
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
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', function() {
    console.log("DOM loaded, checking for sequence network containers");
    if (document.getElementById('sequence-network-container')) {
        // Try to determine protein ID
        let proteinId = '';
        
        // Try to get from URL path
        const pathMatch = window.location.pathname.match(/\/protein\/([^\/]+)/);
        if (pathMatch) {
            proteinId = pathMatch[1];
        }
        
        // Or try to get from a known element on the page
        if (!proteinId) {
            const proteinHeader = document.querySelector('.card-header h5');
            if (proteinHeader) {
                const headerText = proteinHeader.textContent;
                const match = headerText.match(/\(([A-Z0-9]+)\)/);
                if (match) {
                    proteinId = match[1];
                }
            }
        }
        
        if (proteinId) {
            console.log("Calling sequence network visualization for", proteinId);
            sequenceNetworkVisualization(proteinId);
        } else {
            console.error("Could not determine protein ID for visualization");
        }
    }
});