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
        networkContainer.innerHTML = '<div class="alert alert-info m-3 mt-5">No phosphosite sequence similarity network data available. This could be because no sequence similarity matches were found with similarity > 0.4.</div>';
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
            // Debug output to see what properties are available
            console.log(`Node color for ${d.id}: is_known=${d.is_known}, is_known_phosphosite=${d.is_known_phosphosite}`);
            
            // First check is_known_phosphosite (from database)
            if (d.is_known_phosphosite === true || d.is_known_phosphosite === 1 || d.is_known_phosphosite === "1" || d.is_known_phosphosite === 1.0) {
                return '#4CAF50'; // Green for known protein sites
            }
            
            // Then check is_known (fallback property)
            if (d.is_known === true || d.is_known === "true" || d.is_known === "True" || d.is_known === 1) {
                return '#4CAF50'; // Green for known protein sites
            }
            
            // Default to unknown (orange)
            return '#FF9800'; // Orange for unknown protein sites
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
                <p><strong>Known Site:</strong> ${(d.is_known_phosphosite === 1 || d.is_known_phosphosite === 1.0 || d.is_known === true) ? 'Yes' : 'No'}</p>
                ${d.mean_plddt ? `<p><strong>Mean pLDDT:</strong> ${d.mean_plddt}</p>` : ''}
                ${d.nearby_count ? `<p><strong>Nearby Residues:</strong> ${d.nearby_count}</p>` : ''}
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
                <p><strong>Similarity:</strong> ${d.similarity ? (d.similarity * 100).toFixed(1) + '%' : 'N/A'}</p>
                ${d.motif ? `<p><strong>Motif:</strong> <code>${d.motif}</code></p>` : ''}
                ${d.mean_plddt ? `<p><strong>Mean pLDDT:</strong> ${d.mean_plddt}</p>` : ''}
                ${d.nearby_count ? `<p><strong>Nearby Residues:</strong> ${d.nearby_count}</p>` : ''}
                ${d.surface_accessibility ? `<p><strong>Surface Access:</strong> ${typeof d.surface_accessibility === 'number' ? d.surface_accessibility.toFixed(1) : d.surface_accessibility}%</p>` : ''}
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
    
    // Apply initial filter threshold
    updateSequenceNetworkFilter();
    
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
        
        // First: Try direct DOM approach to get phosphosite data
        // This approach doesn't rely on hidden JSON data
        const allPhosphositeRows = document.querySelectorAll('.phosphosite-table tbody tr, table.table-striped tbody tr');
        
        if (allPhosphositeRows.length > 0) {
            console.log(`Found ${allPhosphositeRows.length} rows in phosphosite table`);
            
            // First pass: Create all protein site nodes directly from the table
            allPhosphositeRows.forEach((row, index) => {
                const siteCell = row.querySelector('td:first-child');
                if (!siteCell) return;
                
                // Get site name either from the link or directly from cell
                const siteLink = siteCell.querySelector('a');
                const siteName = siteLink ? siteLink.textContent.trim() : siteCell.textContent.trim();
                if (!siteName) return;
                
                // Create the node ID
                const nodeId = `${proteinUniprotId}_${siteName}`;
                
                // Check all cells to find "Known" column which might be in different positions
                const cells = row.querySelectorAll('td');
                let knownCell = null;
                
                // Try to identify the "Known" column by checking column headers
                const headers = document.querySelectorAll('th');
                let knownColumnIndex = -1;
                headers.forEach((header, i) => {
                    const headerText = header.textContent.trim();
                    if (headerText === 'Known in PhosphositePlus' || headerText === 'Known') {
                        knownColumnIndex = i;
                        console.log(`Found "Known" column at index ${i}`);
                    }
                });
                
                // If we found the Known column index, use it
                if (knownColumnIndex >= 0 && cells.length > knownColumnIndex) {
                    knownCell = cells[knownColumnIndex];
                } else {
                    // Fallback: try common positions for the Known column
                    knownCell = cells[4] || cells[5] || cells[6];
                }
                
                const isKnown = knownCell ? knownCell.textContent.trim() === 'Yes' : false;
                
                // Get site type and other data
                let siteType = siteName[0] || 'S';
                
                // Get motif if available
                let motif = '';
                const motifCell = row.querySelector('td:nth-child(2) code');
                if (motifCell) motif = motifCell.textContent.trim();
                
                // Skip if we already have this node
                if (nodeMap.has(nodeId)) return;
                
                // Create node with is_known property correctly set
                const node = {
                    id: nodeId,
                    name: siteName,
                    uniprot: proteinUniprotId,
                    type: 'protein',
                    is_known: isKnown,  // Use the value we found
                    is_known_phosphosite: isKnown ? 1.0 : 0.0,  // Derive from isKnown
                    siteType: siteType,
                    motif: motif,
                    size: 10
                };
                
                // Add node to our collection
                nodes.push(node);
                nodeMap.set(nodeId, node);
                
                if (isKnown) {
                    console.log(`Added known protein site: ${nodeId} from table`);
                }
            });
            
            console.log(`Created ${nodes.length} protein nodes directly from table, ${nodes.filter(n => n.is_known).length} known sites`);
        }
        
        // Try to get sequence matches from the hidden data element
        const dataElement = document.getElementById('sequence-match-data');
        if (dataElement) {
            console.log("Found sequence-match-data element, attempting to parse JSON data");
            try {
                // Get the JSON data from the data-matches attribute
                const matchesJson = dataElement.getAttribute('data-matches');
                if (matchesJson) {
                    console.log("Found matches JSON string, length:", matchesJson.length);
                    
                    // First, sanitize the JSON data to handle NaN, Infinity, undefined, etc.
                    const sanitizedJson = matchesJson
                        .replace(/NaN/g, 'null')
                        .replace(/Infinity/g, 'null')
                        .replace(/undefined/g, 'null')
                        .replace(/\bnan\b/g, 'null')
                        .replace(/\binfinity\b/g, 'null');
                    
                    console.log("Sanitized JSON, first 100 chars:", sanitizedJson.substring(0, 100));
                    
                    // Parse the sanitized JSON
                    const matchesData = JSON.parse(sanitizedJson);
                    console.log("Successfully parsed matches data");
                    
                    // Process the matches data
                    if (matchesData) {
                        console.log("Processing matches data:", Object.keys(matchesData));
                        
                        // Process each site's matches
                        for (const [siteName, matches] of Object.entries(matchesData)) {
                            // Check if we have matches
                            if (!Array.isArray(matches) || matches.length === 0) {
                                console.log(`No matches for site ${siteName}`);
                                continue;
                            }
                            
                            // Create protein node if not already created
                            const nodeId = `${proteinUniprotId}_${siteName}`;
                            
                            if (!nodeMap.has(nodeId)) {
                                // Get site data from direct row again as fallback
                                const siteRow = document.querySelector(`tr[data-site="${siteName}"]`);
                                let isKnown = false;
                                let siteType = siteName[0] || 'S';
                                let motif = '';
                                
                                if (siteRow) {
                                    // Get is_known from row
                                    if (siteRow.hasAttribute('data-known')) {
                                        const knownValue = siteRow.getAttribute('data-known');
                                        isKnown = knownValue === 'true' || knownValue === '1';
                                    }
                                    
                                    // Get site type
                                    siteType = siteRow.getAttribute('data-type') || siteName[0];
                                    
                                    // Get motif
                                    const motifCell = siteRow.querySelector('td:nth-child(2) code');
                                    if (motifCell) motif = motifCell.textContent.trim();
                                }
                                
                                // Create node
                                const node = {
                                    id: nodeId,
                                    name: siteName,
                                    uniprot: proteinUniprotId,
                                    type: 'protein',
                                    is_known: isKnown,
                                    is_known_phosphosite: isKnown ? 1.0 : 0.0,
                                    siteType: siteType,
                                    motif: motif,
                                    size: 10
                                };
                                
                                nodes.push(node);
                                nodeMap.set(nodeId, node);
                                console.log(`Added protein node from matches: ${nodeId}`);
                            }
                            
                            // Add match nodes and links for this site
                            console.log(`Processing ${matches.length} matches for site ${siteName}`);
                            
                            for (const match of matches) {
                                // Make sure we have a target_id
                                if (!match.target_id) {
                                    console.warn("Match missing target_id:", match);
                                    continue;
                                }
                                
                                // MODIFIED: Use constant lower threshold of 0.4 instead of slider value
                                // Only completely exclude matches with similarity below 0.4
                                const similarity = parseFloat(match.similarity || 0);
                                if (similarity < 0.4) continue;
                                
                                const targetId = match.target_id;
                                const targetUniprot = match.target_uniprot;
                                const targetSite = match.target_site;
                                
                                // Skip self-references
                                if (nodeId === targetId) continue;
                                
                                // Add target node if doesn't exist
                                if (!nodeMap.has(targetId)) {
                                    // Extract site type from target_site
                                    const targetSiteType = targetSite && targetSite[0].match(/[A-Z]/) ? 
                                                        targetSite[0] : 'S';  // Default to 'S'
                                    
                                    const targetNode = {
                                        id: targetId,
                                        name: targetSite || 'Unknown',
                                        uniprot: targetUniprot || 'Unknown',
                                        type: 'match', // Sequence match
                                        is_known: false,  // Default for matches
                                        is_known_phosphosite: 0.0,
                                        siteType: targetSiteType,
                                        similarity: similarity,
                                        motif: match.motif || '',
                                        size: 8
                                    };
                                    
                                    nodes.push(targetNode);
                                    nodeMap.set(targetId, targetNode);
                                    console.log(`Added match node: ${targetId}, similarity: ${similarity}`);
                                }
                                
                                // Add link from protein to match
                                links.push({
                                    source: nodeId,
                                    target: targetId,
                                    similarity: similarity
                                });
                                console.log(`Added link: ${nodeId} -> ${targetId}, similarity: ${similarity}`);
                            }
                        }
                        
                        console.log(`Created network with ${nodes.length} nodes and ${links.length} links from matches`);
                    }
                }
            } catch (e) {
                console.error("Error parsing sequence match data from element:", e);
                console.error("Raw matches JSON (first 100 chars):", dataElement.getAttribute('data-matches').substring(0, 100));
            }
        }
        
        // Final network stats
        const knownProteinNodes = nodes.filter(n => n.type === 'protein' && (n.is_known || n.is_known_phosphosite === 1 || n.is_known_phosphosite === 1.0)).length;
        const unknownProteinNodes = nodes.filter(n => n.type === 'protein' && !n.is_known && n.is_known_phosphosite !== 1 && n.is_known_phosphosite !== 1.0).length;
        const matchNodes = nodes.filter(n => n.type === 'match').length;
        
        console.log(`Final network summary:`);
        console.log(`- Known protein sites (green): ${knownProteinNodes}`);
        console.log(`- Unknown protein sites (orange): ${unknownProteinNodes}`);
        console.log(`- Match sites (purple): ${matchNodes}`);
        console.log(`- Total nodes: ${nodes.length}`);
        console.log(`- Total links: ${links.length}`);
        
        return { nodes, links };
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
    console.log(`Similarity threshold set to: ${threshold}`);
    
    // Count visible links before and after
    const initialVisibleLinks = window.sequenceNetworkLinks.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    // Filter links by similarity
    window.sequenceNetworkLinks.style('display', function(d) {
        const visible = d.similarity >= threshold;
        return visible ? null : 'none';
    });
    
    // Count links after filtering
    const filteredVisibleLinks = window.sequenceNetworkLinks.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    console.log(`Links: ${initialVisibleLinks} → ${filteredVisibleLinks} after filtering`);
    
    // Show nodes only if they have visible connections
    const initialVisibleNodes = window.sequenceNetworkNodes.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    // Update node visibility
    window.sequenceNetworkNodes.style('display', function(d) {
        // Always show protein sites
        if (d.type === 'protein') {
            return null; // Always visible
        }
        
        // For match nodes, check if they have any visible connections
        const linkData = window.sequenceNetworkLinks.data();
        const hasVisibleConnection = linkData.some(link => {
            const isConnected = (
                (link.source.id === d.id || link.target.id === d.id) && 
                link.similarity >= threshold
            );
            return isConnected;
        });
        
        return hasVisibleConnection ? null : 'none';
    });
    
    // Count nodes after filtering
    const filteredVisibleNodes = window.sequenceNetworkNodes.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    console.log(`Nodes: ${initialVisibleNodes} → ${filteredVisibleNodes} after filtering`);
    
    // Update labels visibility to match nodes
    if (window.sequenceNetworkLabels) {
        window.sequenceNetworkLabels.style('display', function(d) {
            // Find the corresponding node
            const matchingNode = window.sequenceNetworkNodes.filter(function(n) {
                return n.__data__.id === d.id;
            });
            
            if (matchingNode.size() > 0) {
                const node = matchingNode.node();
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
                const match = headerText.match(/\(([A-Z0-9_]+)\)/);
                if (match) {
                    proteinId = match[1];
                }
            }
        }
        
        // Try another approach if still not found
        if (!proteinId) {
            // Look for protein ID in page content or URL
            const contentMatch = document.body.textContent.match(/UniProt ID: ([A-Z0-9_]+)/);
            if (contentMatch) {
                proteinId = contentMatch[1];
            } else {
                // Check if it's in the URL hash
                const hashMatch = window.location.hash.match(/#protein=([A-Z0-9_]+)/);
                if (hashMatch) {
                    proteinId = hashMatch[1];
                }
            }
        }
        
        if (proteinId) {
            console.log("Calling sequence network visualization for", proteinId);
            
            // IMPORTANT: Set initial threshold to 0.4 before visualization
            const thresholdSlider = document.getElementById('similarity-filter');
            if (thresholdSlider) {
                thresholdSlider.value = 0.4;
                // Update the displayed value
                const valueDisplay = document.getElementById('similarity-value');
                if (valueDisplay) valueDisplay.textContent = '0.4';
            }
            
            sequenceNetworkVisualization(proteinId);
        } else {
            console.error("Could not determine protein ID for visualization");
        }
    }
});