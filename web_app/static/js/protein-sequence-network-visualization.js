
/**
 * Enhanced Sequence Network Visualization
 * Creates an interactive visualization of phosphosite sequence similarity networks
 * - Filters out phosphosites with Tyrosine (Y)
 * - Ensures motifs are displayed in tooltips
 * - Shows kinase information when available
 * - Properly handles coloring nodes with known kinases
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
    const networkData = extractSequenceNetworkData(proteinUniprotId);
    
    if (!networkData || networkData.nodes.length === 0) {
        networkContainer.innerHTML = '<div class="alert alert-info m-3 mt-5">No phosphosite sequence similarity network data available. This could be because no sequence similarity matches were found with similarity > 0.4.</div>';
        return;
    }
    
    console.log(`Creating sequence network with ${networkData.nodes.length} nodes and ${networkData.links.length} links before filtering`);
    
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
            return '#9C27B0'; // Magenta for any site with a known kinase
        }
        
        // Then check node type
        if (d.type === 'protein') {
            // Check is_known_phosphosite (from database)
            if (d.is_known_phosphosite === true || 
                d.is_known_phosphosite === 1 || 
                d.is_known_phosphosite === "1" || 
                d.is_known_phosphosite === 1.0) {
                return '#4CAF50'; // Green for known protein sites
            }
            
            // Check is_known (fallback property)
            if (d.is_known === true || 
                d.is_known === "true" || 
                d.is_known === "True" || 
                d.is_known === 1) {
                return '#4CAF50'; // Green for known protein sites
            }
            
            // Default to unknown (orange)
            return '#FF9800'; // Orange for unknown protein sites
        }
        
        // For match nodes without a known kinase
        return '#E91E63'; // Pink for sequence-similar sites without known kinases
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
    
    // Setup force simulation - store it in window for access by filter function
    window.sequenceNetworkSimulation = d3.forceSimulation(networkData.nodes)
        .force('link', d3.forceLink(networkData.links)
            .id(d => d.id)
            .distance(d => 100 * (1 - (d.similarity || 0.5))))
        .force('charge', d3.forceManyBody()
            .strength(d => d.type === 'protein' ? -50 : -25))
        .force('center', d3.forceCenter(width / 2, height / 2))
        .force('collision', d3.forceCollide().radius(d => d.size + 3))
        .force('x', d3.forceX(width / 2).strength(0.1))
        .force('y', d3.forceY(height / 2).strength(0.1));
    
    // Draw links - store in window for access by filter function
    window.sequenceNetworkLinks = g.append('g')
        .selectAll('line')
        .data(networkData.links)
        .enter()
        .append('line')
        .attr('stroke', getLinkColor)
        .attr('stroke-width', getLinkWidth)
        .attr('stroke-opacity', 0.6)
        .attr('class', 'sequence-network-link');
    
    // Draw nodes - store in window for access by filter function
    window.sequenceNetworkNodes = g.append('g')
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
    
    // Add labels for protein sites - store in window for access by filter function
    window.sequenceNetworkLabels = g.append('g')
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
                ${d.known_kinase ? `<p><strong>Known Kinase${d.known_kinase.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(d.known_kinase)}</p>` : ''}
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
                ${d.known_kinase ? `<p><strong>Known Kinase${d.known_kinase.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(d.known_kinase)}</p>` : ''}
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
    
    // Helper function to format kinase list with badges
    function formatKinaseList(kinaseString) {
        if (!kinaseString) return '';
        
        const kinases = kinaseString.split(',').map(k => k.trim());
        return kinases.map(kinase => 
            `<span class="badge bg-primary me-1">${kinase}</span>`
        ).join(' ');
    }
    
    // Node interactions
    window.sequenceNetworkNodes
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
    window.sequenceNetworkSimulation.on('tick', () => {
        window.sequenceNetworkLinks
            .attr('x1', d => d.source.x)
            .attr('y1', d => d.source.y)
            .attr('x2', d => d.target.x)
            .attr('y2', d => d.target.y);
        
        window.sequenceNetworkNodes
            .attr('cx', d => d.x)
            .attr('cy', d => d.y);
        
        window.sequenceNetworkLabels
            .attr('x', d => d.x)
            .attr('y', d => d.y);
    });
    
    // Apply initial filter threshold
    updateSequenceNetworkFilter();
    
    // Add a legend for the node colors
    addNetworkLegend(svg, width);
    
    // Update the legend in the container's HTML to match updated color scheme
    updateNetworkLegendHTML();
    
    // Drag functions
    function dragstarted(event, d) {
        if (!event.active) window.sequenceNetworkSimulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }
    
    function dragged(event, d) {
        d.fx = event.x;
        d.fy = event.y;
    }
    
    function dragended(event, d) {
        if (!event.active) window.sequenceNetworkSimulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }
    
    console.log('Sequence network visualization setup complete');
}

// 1. Update the addNetworkLegend function to include all node types with correct colors
function addNetworkLegend(svg, width) {
    const legend = svg.append("g")
        .attr("class", "legend")
        .attr("transform", `translate(20, 20)`);
    
    const legendData = [
        { color: "#4CAF50", label: "Known protein sites" },
        { color: "#FF9800", label: "Unknown protein sites" },
        { color: "#9C27B0", label: "Sites with known kinase" },
        { color: "#E91E63", label: "Sequence-similar sites" }
    ];
    
    // Add background for better visibility
    legend.append("rect")
        .attr("x", -5)
        .attr("y", -5)
        .attr("width", 190)
        .attr("height", legendData.length * 20 + 10)
        .attr("fill", "white")
        .attr("opacity", 0.8)
        .attr("rx", 5)
        .attr("ry", 5);
    
    const legendItems = legend.selectAll(".legend-item")
        .data(legendData)
        .enter()
        .append("g")
        .attr("class", "legend-item")
        .attr("transform", (d, i) => `translate(0, ${i * 20})`);
    
    // Add colored circles (not squares - fix this issue)
    legendItems.append("circle")
        .attr("r", 6)
        .attr("cx", 6)
        .attr("cy", 4)
        .attr("fill", d => d.color);
    
    // Add labels
    legendItems.append("text")
        .attr("x", 15)
        .attr("y", 8)
        .text(d => d.label)
        .attr("font-size", "12px");
}

// Extract network data from the page with enhanced motif and kinase handling
function extractSequenceNetworkData(proteinUniprotId) {
    try {
        console.log("Extracting sequence network data for", proteinUniprotId);
        
        // Create arrays for nodes and links
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
            const allPhosphositeRows = document.querySelectorAll('.phosphosite-table tbody tr, table.table-striped tbody tr');
            
            if (allPhosphositeRows.length > 0) {
                console.log(`Found ${allPhosphositeRows.length} rows in phosphosite table`);
                
                allPhosphositeRows.forEach((row, index) => {
                    const siteCell = row.querySelector('td:first-child');
                    if (!siteCell) return;
                    
                    // Get site name either from the link or directly from cell
                    const siteLink = siteCell.querySelector('a');
                    const siteName = siteLink ? siteLink.textContent.trim() : siteCell.textContent.trim();
                    if (!siteName) return;
                    
                    // Skip Tyrosine (Y) sites as requested
                    if (siteName[0] === 'Y') {
                        console.log(`Skipping Tyrosine site: ${siteName}`);
                        return;
                    }
                    
                    // Extract site data
                    const isKnown = row.getAttribute('data-known') === 'true';
                    let siteType = siteName[0];
                    
                    // Get motif if available
                    let motif = '';
                    const motifCell = row.querySelector('td:nth-child(2) code');
                    if (motifCell) motif = motifCell.textContent.trim();
                    
                    // Create site data object
                    phosphositesData.push({
                        site: siteName,
                        resno: parseInt(siteName.substring(1)),
                        siteType: siteType,
                        is_known: isKnown,
                        motif: motif
                    });
                });
                
                console.log(`Extracted ${phosphositesData.length} phosphosites from table`);
            }
        }
        
        // Step 3: Add all protein site nodes (excluding Tyrosine sites)
        for (const site of phosphositesData) {
            // Skip Tyrosine (Y) sites as requested
            if (site.site && site.site[0] === 'Y') {
                console.log(`Skipping Tyrosine site: ${site.site}`);
                continue;
            }
            
            const nodeId = `${proteinUniprotId}_${site.resno}`;
            
            // Skip if already added
            if (nodeMap.has(nodeId)) continue;
            
            // Extract known kinase information if available
            const knownKinase = extractKinaseInfo(site);
            
            // Create node with all available data
            const node = {
                id: nodeId,
                name: site.site || `${site.siteType || 'S'}${site.resno}`,
                uniprot: proteinUniprotId,
                type: 'protein',
                is_known: site.is_known === true || site.is_known === 1,
                is_known_phosphosite: site.is_known_phosphosite || (site.is_known ? 1.0 : 0.0),
                siteType: site.siteType || site.site[0] || 'S',
                motif: site.motif || '',
                known_kinase: knownKinase,
                mean_plddt: site.mean_plddt || site.meanPLDDT || '',
                nearby_count: site.nearby_count || site.nearbyCount || '',
                surface_accessibility: site.surface_accessibility || site.surfaceAccessibility || '',
                size: 10
            };
            
            // Add node to our collection
            nodes.push(node);
            nodeMap.set(nodeId, node);
        }
        
        console.log(`Created ${nodes.length} protein nodes`);
        
        // Step 4: Get sequence matches from the hidden data element
        const matchesDataElement = document.getElementById('sequence-match-data');
        if (matchesDataElement) {
            console.log("Found sequence-match-data element");
            try {
                // Get the JSON data from the data-matches attribute
                const matchesJson = matchesDataElement.getAttribute('data-matches');
                if (matchesJson) {
                    // Clean and parse the JSON
                    const sanitizedJson = cleanJson(matchesJson);
                    const matchesData = JSON.parse(sanitizedJson);
                    
                    // Process the matches data
                    if (matchesData) {
                        console.log("Processing matches data:", Object.keys(matchesData));
                        
                        // Cache of supplementary data to avoid duplicated requests
                        const suppDataCache = new Map();
                        
                        // Process each site's matches
                        for (const [siteName, matches] of Object.entries(matchesData)) {
                            // Skip if the site name starts with 'Y' (Tyrosine)
                            if (siteName && siteName[0] === 'Y') {
                                console.log(`Skipping matches for Tyrosine site: ${siteName}`);
                                continue;
                            }
                            
                            // Skip any sites without matches
                            if (!Array.isArray(matches) || matches.length === 0) {
                                console.log(`No matches for site ${siteName}`);
                                continue;
                            }
                            
                            // Get site info
                            const siteType = siteName[0] || 'S';
                            const siteNumber = parseInt(siteName.substring(1));
                            
                            // Skip invalid sites
                            if (isNaN(siteNumber)) continue;
                            
                            // Create node ID
                            const nodeId = `${proteinUniprotId}_${siteNumber}`;
                            
                            // Find the protein node for this site
                            let proteinNode = nodes.find(n => n.id === nodeId);
                            
                            // If we don't have this protein node yet (possibly because it wasn't in phosphositesData)
                            if (!proteinNode && siteType !== 'Y') {
                                // Find any existing motif info in the matches
                                let motifInfo = '';
                                for (const match of matches) {
                                    if (match.query_motif) {
                                        motifInfo = match.query_motif;
                                        break;
                                    }
                                }
                                
                                // Create the protein node
                                proteinNode = {
                                    id: nodeId,
                                    name: siteName,
                                    uniprot: proteinUniprotId,
                                    type: 'protein',
                                    is_known: false,  // Default to unknown
                                    is_known_phosphosite: 0,
                                    siteType: siteType,
                                    motif: motifInfo,
                                    size: 10
                                };
                                
                                // Add to nodes collection
                                nodes.push(proteinNode);
                                nodeMap.set(nodeId, proteinNode);
                                console.log(`Added protein node from matches: ${nodeId}`);
                            }
                            
                            // Skip if we still don't have a valid protein node (e.g., for Tyrosine sites)
                            if (!proteinNode) continue;
                            
                            // Process the matches for this site
                            console.log(`Processing ${matches.length} matches for site ${siteName}`);
                            
                            for (const match of matches) {
                                // Skip if missing target_id
                                if (!match.target_id) continue;
                                
                                // Skip matches with Tyrosine site type
                                if (match.site_type === 'Y' || 
                                    (match.target_site && match.target_site[0] === 'Y')) {
                                    console.log(`Skipping Tyrosine match: ${match.target_site}`);
                                    continue;
                                }
                                
                                // Get similarity and validate
                                const similarity = parseFloat(match.similarity || 0);
                                if (similarity < 0.4) continue;
                                
                                const targetId = match.target_id;
                                const targetUniprot = match.target_uniprot;
                                const targetSite = match.target_site;
                                
                                // Skip self-references
                                if (nodeId === targetId) continue;
                                
                                // Extract or retrieve motif information
                                let motif = match.motif || '';
                                
                                // If no motif is available, try to extract from supplementary data
                                if (!motif && !suppDataCache.has(targetId)) {
                                    // Attempt to find motif elsewhere in the page
                                    motif = findMotifInPage(targetId, targetSite, targetUniprot);
                                    suppDataCache.set(targetId, { motif });
                                } else if (suppDataCache.has(targetId)) {
                                    // Use cached data
                                    const cachedData = suppDataCache.get(targetId);
                                    motif = motif || cachedData.motif;
                                }
                                
                                // Extract known kinase information - enhanced to check more sources
                                const knownKinase = extractKinaseInfo(match);
                                
                                // Check if this is a match within the same protein
                                const isSameProtein = targetUniprot === proteinUniprotId;
                                // For matches within the same protein, we still want to keep the match type
                                // This ensures we don't incorrectly classify nodes
                                
                                // Add target node if not already present
                                if (!nodeMap.has(targetId)) {
                                    // Extract site type from target_site
                                    const targetSiteType = match.site_type || 
                                                         (targetSite && targetSite[0].match(/[A-Z]/) ? 
                                                          targetSite[0] : 'S');
                                    
                                    // Skip if it's a Tyrosine site
                                    if (targetSiteType === 'Y') continue;
                                    
                                    const targetNode = {
                                        id: targetId,
                                        name: targetSite || 'Unknown',
                                        uniprot: targetUniprot || 'Unknown',
                                        // Keep 'match' type even if it's the same protein
                                        type: 'match',
                                        is_known: match.is_known === true || match.is_known === 1,
                                        is_known_phosphosite: match.is_known_phosphosite || 0,
                                        siteType: targetSiteType,
                                        similarity: similarity,
                                        motif: motif,
                                        known_kinase: knownKinase,
                                        mean_plddt: match.mean_plddt || match.site_plddt || '',
                                        nearby_count: match.nearby_count || '',
                                        surface_accessibility: match.surface_accessibility || '',
                                        size: 8
                                    };
                                    
                                    nodes.push(targetNode);
                                    nodeMap.set(targetId, targetNode);
                                }
                                
                                // Add link between protein site and match
                                links.push({
                                    source: nodeId,
                                    target: targetId,
                                    similarity: similarity
                                });
                            }
                        }
                    }
                }
            } catch (e) {
                console.error("Error processing sequence match data:", e);
            }
        }
        
        // Final filtering step - remove any remaining Tyrosine nodes/links
        const filteredNodes = nodes.filter(node => 
            node.siteType !== 'Y' && 
            (node.name && node.name[0] !== 'Y')
        );
        
        const filteredLinks = links.filter(link => {
            const sourceNode = nodeMap.get(link.source);
            const targetNode = nodeMap.get(link.target);
            if (!sourceNode || !targetNode) return false;
            
            return sourceNode.siteType !== 'Y' && 
                   targetNode.siteType !== 'Y' && 
                   sourceNode.name[0] !== 'Y' && 
                   targetNode.name[0] !== 'Y';
        });
        
        console.log(`Network after filtering Tyrosine sites:
                    Nodes: ${nodes.length} → ${filteredNodes.length}
                    Links: ${links.length} → ${filteredLinks.length}`);
        
        return { nodes: filteredNodes, links: filteredLinks };
        
    } catch (error) {
        console.error("Error extracting sequence network data:", error);
        return { nodes: [], links: [] };
    }
}

// Helper function to clean JSON before parsing
function cleanJson(jsonString) {
    return jsonString
        .replace(/NaN/g, 'null')
        .replace(/Infinity/g, 'null')
        .replace(/undefined/g, 'null')
        .replace(/\bnan\b/g, 'null')
        .replace(/\binfinity\b/g, 'null');
}

// Enhanced function to extract kinase information from site data
// Looks for kinase data in multiple possible locations
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

// Function to search for motif information elsewhere on the page
function findMotifInPage(siteId, siteName, uniprotId) {
    // Try different patterns to find motif information
    
    // First, check for dedicated motif/sequence cells in tables
    const tables = document.querySelectorAll('table');
    for (const table of tables) {
        // Check if this table contains site information
        const siteRows = table.querySelectorAll('tr');
        for (const row of siteRows) {
            // Check if this row is for our site
            let matchesOurSite = false;
            const cells = row.querySelectorAll('td');
            
            for (const cell of cells) {
                // Check if cell contains site identifier
                if (cell.textContent.includes(siteName) || 
                    cell.textContent.includes(siteId) ||
                    (cell.innerHTML.includes(uniprotId) && cell.innerHTML.includes(siteName))) {
                    matchesOurSite = true;
                    break;
                }
            }
            
            // If this row is for our site, look for motif in a code element
            if (matchesOurSite) {
                const motifCell = row.querySelector('td code.motif-sequence, td:nth-child(2) code, td:nth-child(4) code');
                if (motifCell && motifCell.textContent) {
                    return motifCell.textContent.trim();
                }
            }
        }
    }
    
    // If we can't find the motif in tables, try the sequence match data
    const allSequenceMatchCards = document.querySelectorAll('.sequence-match-card');
    for (const card of allSequenceMatchCards) {
        const rows = card.querySelectorAll('tr');
        for (const row of rows) {
            const idCell = row.querySelector('td:first-child a');
            const siteCell = row.querySelector('td:nth-child(2)');
            
            // Check if this row matches our target
            if ((idCell && idCell.textContent.trim() === uniprotId) &&
                (siteCell && siteCell.textContent.trim() === siteName)) {
                const motifCell = row.querySelector('td code.motif-sequence');
                if (motifCell && motifCell.textContent) {
                    return motifCell.textContent.trim();
                }
            }
        }
    }
    
    // If still not found, look for similarity match data in the page content
    const sequenceMatchesElement = document.getElementById('sequence-match-data');
    if (sequenceMatchesElement) {
        try {
            const matchesJson = sequenceMatchesElement.getAttribute('data-matches');
            if (matchesJson) {
                const matchesData = JSON.parse(cleanJson(matchesJson));
                
                // Search through all matches
                for (const [siteName, matches] of Object.entries(matchesData)) {
                    if (Array.isArray(matches)) {
                        for (const match of matches) {
                            if (match.target_id === siteId || 
                                (match.target_uniprot === uniprotId && match.target_site === siteName)) {
                                if (match.motif) {
                                    return match.motif;
                                }
                            }
                        }
                    }
                }
            }
        } catch (e) {
            console.log("Error parsing match data for motif search:", e);
        }
    }
    
    return ''; // Return empty string if no motif found
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
    
    // First show all nodes and links to reset any previous filtering
    window.sequenceNetworkNodes.style('display', null);
    window.sequenceNetworkLinks.style('display', null);
    
    // Now count all nodes and links before applying the new filter
    const initialVisibleLinks = window.sequenceNetworkLinks.size();
    const initialVisibleNodes = window.sequenceNetworkNodes.size();
    
    // Filter links by similarity
    window.sequenceNetworkLinks.style('display', function(d) {
        return d.similarity >= threshold ? null : 'none';
    });
    
    // Create a set to track nodes that should be visible
    const visibleNodeIds = new Set();
    
    // Add all protein nodes as they should always be visible
    window.sequenceNetworkNodes.each(function(d) {
        if (d.type === 'protein') {
            visibleNodeIds.add(d.id);
        }
    });
    
    // Add match nodes that have connections above the threshold
    window.sequenceNetworkLinks.each(function(d) {
        if (d.similarity >= threshold) {
            // Handle both object and string source/target
            const sourceId = typeof d.source === 'object' ? d.source.id : d.source;
            const targetId = typeof d.target === 'object' ? d.target.id : d.target;
            visibleNodeIds.add(sourceId);
            visibleNodeIds.add(targetId);
        }
    });
    
    // Update node visibility based on our visibility set
    window.sequenceNetworkNodes.style('display', function(d) {
        return visibleNodeIds.has(d.id) ? null : 'none';
    });
    
    // Update label visibility
    if (window.sequenceNetworkLabels) {
        window.sequenceNetworkLabels.style('display', function(d) {
            return visibleNodeIds.has(d.id) ? null : 'none';
        });
    }
    
    // Count nodes and links after filtering
    const filteredVisibleLinks = window.sequenceNetworkLinks.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    const filteredVisibleNodes = window.sequenceNetworkNodes.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    console.log(`Links: ${initialVisibleLinks} → ${filteredVisibleLinks} after filtering`);
    console.log(`Nodes: ${initialVisibleNodes} → ${filteredVisibleNodes} after filtering`);
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', function() {
    console.log("DOM loaded, checking for sequence network containers");
    const networkContainer = document.getElementById('sequence-network-container');
    if (networkContainer) {
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
            sequenceNetworkVisualization(proteinId);
        } else {
            console.error("Could not determine protein ID for visualization");
        }
    }
});

// 2. Update the updateNetworkLegendHTML function to ensure all legend items are present
function updateNetworkLegendHTML() {
    const legendContainer = document.querySelector('.d-flex.align-items-center.flex-wrap');
    if (legendContainer) {
        legendContainer.innerHTML = `
            <div class="d-flex align-items-center me-4 mb-2">
                <div style="width: 16px; height: 16px; background-color: #4CAF50; border-radius: 50%; margin-right: 6px;"></div>
                <span class="small">Known protein sites</span>
            </div>
            <div class="d-flex align-items-center me-4 mb-2">
                <div style="width: 16px; height: 16px; background-color: #FF9800; border-radius: 50%; margin-right: 6px;"></div>
                <span class="small">Unknown protein sites</span>
            </div>
            <div class="d-flex align-items-center me-4 mb-2">
                <div style="width: 16px; height: 16px; background-color: #9C27B0; border-radius: 50%; margin-right: 6px;"></div>
                <span class="small">Sites with known kinase</span>
            </div>
            <div class="d-flex align-items-center mb-2">
                <div style="width: 16px; height: 16px; background-color: #E91E63; border-radius: 50%; margin-right: 6px;"></div>
                <span class="small">Sequence-similar sites</span>
            </div>
        `;
    }
}


