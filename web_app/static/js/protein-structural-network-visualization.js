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
    
    // Update RMSD slider range
    const rmsdFilter = document.getElementById('rmsd-filter');
    if (rmsdFilter) {
        rmsdFilter.min = "0.1";
        rmsdFilter.max = "10.0";
        if (parseFloat(rmsdFilter.value) > 10.0) {
            rmsdFilter.value = "10.0";
            document.getElementById('rmsd-value').textContent = "10.0 Å";
        }
    }
    
    // Extract network data from the page DOM
    const networkData = extractStructuralNetworkData(proteinUniprotId);
    
    if (!networkData || networkData.nodes.length === 0) {
        networkContainer.innerHTML = '<div class="alert alert-info m-3 mt-5">No phosphosite structural network data available. This could be because no structural matches were found with RMSD < 10.0.</div>';
        return;
    }
    
    console.log(`Creating structural network with ${networkData.nodes.length} nodes and ${networkData.links.length} links`);
    
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
        return '#E91E63'; // Gray for structurally similar sites without known kinases
    }
    
    function getLinkColor(d) {
        if (d.rmsd < 1.0) return '#4CAF50'; // Green for very low RMSD
        if (d.rmsd < 2.0) return '#8BC34A'; // Light green for low RMSD
        if (d.rmsd < 4.0) return '#CDDC39'; // Lime for medium RMSD
        if (d.rmsd < 7.0) return '#FFC107'; // Yellow/amber for higher RMSD
        return '#FF9800'; // Orange for highest RMSD values
    }
    
    function getLinkWidth(d) {
        // Scale link width based on RMSD - thicker for lower RMSD
        return Math.max(0.5, 5 / Math.max(1, d.rmsd));
    }
    



    // Setup force simulation with reduced physics
    window.networkSimulation = d3.forceSimulation(networkData.nodes)
        .force('link', d3.forceLink(networkData.links)
            .id(d => d.id)
            .distance(d => d.rmsd ? d.rmsd * 40 : 150))
        .force('charge', d3.forceManyBody()
            .strength(d => d.type === 'protein' ? -300 : -50))
        .force('center', d3.forceCenter(width / 2, height / 2))
        .force('collision', d3.forceCollide().radius(d => d.size + 5))
        .force('x', d3.forceX(width / 2).strength(0.01))
        .force('y', d3.forceY(height / 2).strength(0.01))
        .alphaDecay(0.05)
        .velocityDecay(0.6);
    
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
        if (!d) return;
        
        console.log("Showing info for node:", d); // Debug log
        
        let content = '';
        if (d.type === 'protein') {
            // For source protein sites
            content = `
                <h6 class="border-bottom pb-2 mb-2">${d.name || ''} - ${d.uniprot || proteinUniprotId}</h6>
                <p><strong>Site Type:</strong> ${d.siteType || d.name[0] || 'S'}</p>
                <p><strong>Known Site:</strong> ${d.isKnown || d.is_known_phosphosite === 1 ? 'Yes' : 'No'}</p>`;
                    
            // Only add kinase info if available
            if (d.known_kinase) {
                content += `<p><strong>Known Kinase${d.known_kinase.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(d.known_kinase)}</p>`;
            }
            
            // Only add fields that exist
            if (d.meanPLDDT || d.mean_plddt) {
                content += `<p><strong>Mean pLDDT:</strong> ${d.meanPLDDT || d.mean_plddt}</p>`;
            }
            
            if (d.nearbyCount || d.nearby_count) {
                content += `<p><strong>Nearby Residues:</strong> ${d.nearbyCount || d.nearby_count}</p>`;
            }
            
            if (d.surface_accessibility || d.surfaceAccessibility) {
                const accessValue = d.surface_accessibility || d.surfaceAccessibility;
                content += `<p><strong>Surface Access:</strong> ${typeof accessValue === 'number' ? accessValue.toFixed(1) : accessValue}%</p>`;
            }
            
            if (d.motif) {
                content += `<p><strong>Motif:</strong> <code>${d.motif}</code></p>`;
            }
            
            content += `
                <div class="d-grid gap-2 mt-3">
                    <a href="/site/${d.uniprot || proteinUniprotId}/${d.name}" class="btn btn-sm btn-primary">View Site Details</a>
                </div>
            `;
        } else {
            // For match nodes - major improvements to data extraction
            
            // Access the raw match data if available for maximum data availability
            const matchData = d.matchData || d;
            
            // Extract display values with robust fallbacks
            const displayName = d.name || '';
            const displayUniprot = d.uniprot || '';
            const displaySiteType = d.siteType || (displayName && displayName[0]) || 'S';
            
            content = `
                <h6 class="border-bottom pb-2 mb-2">${displayName} - ${displayUniprot}</h6>
                <p><strong>Site Type:</strong> ${displaySiteType}</p>`;
                    
            // Add RMSD with safety check - critical for match nodes
            if (d.rmsd !== undefined && d.rmsd !== null && !isNaN(d.rmsd)) {
                content += `<p><strong>RMSD:</strong> ${parseFloat(d.rmsd).toFixed(2)} Å</p>`;
            }
            
            // Check for kinase info in all possible locations
            const kinaseInfo = d.known_kinase || 
                              matchData.known_kinase || 
                              matchData.target_known_kinase ||
                              null;
                              
            if (kinaseInfo) {
                content += `<p><strong>Known Kinase${kinaseInfo.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(kinaseInfo)}</p>`;
            }
            
            // Comprehensive check for pLDDT values across all possible field names
            const plddt = d.mean_plddt || d.site_plddt || d.meanPLDDT || 
                         matchData.mean_plddt || matchData.site_plddt || matchData.pLDDT ||
                         null;
            if (plddt !== null && plddt !== undefined && plddt !== '') {
                content += `<p><strong>Mean pLDDT:</strong> ${plddt}</p>`;
            }
            
            // Check all possible fields for nearby residue count
            const nearbyCount = d.nearby_count || d.nearbyCount || d.NeighborCount || 
                              matchData.nearby_count || matchData.nearbyCount || matchData.NeighborCount ||
                              null;
            if (nearbyCount !== null && nearbyCount !== undefined && nearbyCount !== '') {
                content += `<p><strong>Nearby Residues:</strong> ${nearbyCount}</p>`;
            }
            
            // Comprehensive check for surface accessibility values
            const accessValue = d.surface_accessibility || d.surfaceAccessibility || 
                              matchData.surface_accessibility || matchData.surfaceAccessibility ||
                              (matchData.HydroxylExposure ? matchData.HydroxylExposure * 100 : null);
            
            if (accessValue !== null && accessValue !== undefined && accessValue !== '') {
                // Format percentage properly
                content += `<p><strong>Surface Access:</strong> ${typeof accessValue === 'number' ? accessValue.toFixed(1) : accessValue}%</p>`;
            }
            
            // Check for motif in both the node and raw match data
            const motif = d.motif || matchData.motif || null;
            if (motif) {
                content += `<p><strong>Motif:</strong> <code>${motif}</code></p>`;
            }
            
            content += `
                <div class="d-grid gap-2 mt-3">
                    <a href="https://www.uniprot.org/uniprotkb/${displayUniprot}" class="btn btn-sm btn-outline-primary" target="_blank">View on UniProt</a>
                    <a href="/site/${displayUniprot}/${displayName}" class="btn btn-sm btn-primary">View Site Details</a>
                </div>
            `;
        }
        
        infoPanel.innerHTML = content;
    }

    // Helper function to format kinase list with badges
    function formatKinaseList(kinaseString) {
        if (!kinaseString) return '';
        
        const kinases = kinaseString.split(',').map(k => k.trim()).filter(k => k);
        return kinases.map(kinase => 
            `<span class="badge bg-primary me-1">${kinase}</span>`
        ).join(' ');
    }
    
    // Node interactions
    node
        .on('mouseover', function(event, d) {
            // Highlight node on hover
            // Log data for debugging
            console.log("Hovering node:", d);

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
            window.location.href = `/site/${d.uniprot || proteinUniprotId}/${d.name}`;
        });
    
    // Tick function for the simulation
    window.networkSimulation.on('tick', () => {
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
    
    // Set up legend
    addNetworkLegend(svg, width);
    
    // Drag functions
    function dragstarted(event, d) {
        if (!event.active) window.networkSimulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }
    
    function dragged(event, d) {
        d.fx = event.x;
        d.fy = event.y;
    }
    
    function dragended(event, d) {
        if (!event.active) window.networkSimulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }
}

// Add a legend to the network visualization
function addNetworkLegend(svg, width) {
    const legend = svg.append("g")
        .attr("class", "legend")
        .attr("transform", `translate(20, 20)`);
    
    const legendData = [
        { color: "#4CAF50", label: "Known protein sites" },
        { color: "#FF9800", label: "Unknown protein sites" },
        { color: "#9C27B0", label: "Sites with known kinase" },
        { color: "#E91E63", label: "Structurally similar sites" }
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
    
    // Add colored circles
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

function extractStructuralNetworkData(proteinUniprotId) {
    try {
        console.log("Extracting structural network data for", proteinUniprotId);
        
        // Initialize node and link collections
        const nodes = [];
        const links = [];
        const nodeMap = new Map(); // Track unique nodes by ID
        
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
                        let meanPLDDT = '';
                        const pLDDTCell = row.querySelector('td:nth-child(3)');
                        if (pLDDTCell) meanPLDDT = pLDDTCell.textContent.trim();
                        
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
                        
                        // Create site data object matching Python naming conventions
                        phosphositesData.push({
                            site: siteName,
                            resno: parseInt(resno),
                            siteType: siteType,
                            is_known: isKnown,
                            is_known_phosphosite: isKnown ? 1.0 : 0.0,
                            motif: motif,
                            meanPLDDT: meanPLDDT,
                            mean_plddt: meanPLDDT,
                            nearbyCount: nearbyCount,
                            nearby_count: nearbyCount,
                            known_kinase: knownKinase,
                            surface_accessibility: row.getAttribute('data-surface') || ''
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
            
            // Create node ID - for source protein sites
            const resno = site.resno || parseInt(site.site.substring(1));
            if (isNaN(resno)) {
                console.warn(`Invalid residue number for site: ${site.site}`);
                continue;
            }
            
            const nodeId = `${proteinUniprotId}_${resno}`;
            
            // Skip if already added
            if (nodeMap.has(nodeId)) continue;
            
            // Extract known kinase information if available
            const knownKinase = site.known_kinase || site.knownKinase;
            
            // Create node with all available data using Python naming conventions
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
                meanPLDDT: site.meanPLDDT || site.mean_plddt || '',
                mean_plddt: site.meanPLDDT || site.mean_plddt || '',
                nearbyCount: site.nearbyCount || site.nearby_count || '',
                nearby_count: site.nearbyCount || site.nearby_count || '',
                surface_accessibility: site.surface_accessibility || '',
                surfaceAccessibility: site.surface_accessibility || '',
                size: 10,
                hasConnections: false  // Track whether this node has connections
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
                    console.log("Found structural match data");
                    const matchesData = JSON.parse(cleanJson(matchesJson));
                    
                    if (matchesData) {
                        console.log(`Processing matches from data element for ${Object.keys(matchesData).length} sites`);
                        
                        // Process each site's matches
                        for (const [siteName, matches] of Object.entries(matchesData)) {
                            // Skip if site name starts with 'Y' (Tyrosine)
                            if (siteName[0] === 'Y') continue;
                            
                            // Skip sites without matches
                            if (!Array.isArray(matches) || matches.length === 0) continue;
                            
                            // Extract site number for consistent node IDs
                            let sourceNodeId;
                            if (siteName.includes('_')) {
                                // Format: UniProtID_ResidueNumber
                                sourceNodeId = siteName;
                            } else {
                                // Format might be just the site name (e.g., "S102")
                                const siteMatch = siteName.match(/^([STY])(\d+)$/);
                                if (siteMatch) {
                                    const siteType = siteMatch[1];
                                    const siteNumber = parseInt(siteMatch[2]);
                                    sourceNodeId = `${proteinUniprotId}_${siteNumber}`;
                                } else {
                                    // Try simple extraction of digits if pattern doesn't match
                                    const siteResno = parseInt(siteName.replace(/[A-Z]/g, ''));
                                    if (isNaN(siteResno)) {
                                        console.warn(`Invalid residue number for site: ${siteName}`);
                                        continue;
                                    }
                                    sourceNodeId = `${proteinUniprotId}_${siteResno}`;
                                }
                            }
                            
                            console.log(`Processing ${matches.length} matches for site ${siteName} (${sourceNodeId})`);
                            
                            // Ensure source node exists
                            if (!nodeMap.has(sourceNodeId)) {
                                // If the source node doesn't exist, create it
                                console.log(`Source node ${sourceNodeId} not found, creating it`);
                                
                                // Determine site type
                                let siteType;
                                if (siteName && siteName.match(/^[STY]/)) {
                                    siteType = siteName[0];
                                } else {
                                    siteType = 'S'; // Default to Serine if unknown
                                }
                                
                                // Create new source node
                                const sourceNode = {
                                    id: sourceNodeId,
                                    name: siteName,
                                    uniprot: proteinUniprotId,
                                    type: 'protein',
                                    isKnown: false,
                                    siteType: siteType,
                                    size: 10,
                                    hasConnections: false
                                };
                                
                                nodes.push(sourceNode);
                                nodeMap.set(sourceNodeId, sourceNode);
                            }
                            
                            // Mark this source node as having connections
                            const sourceNode = nodeMap.get(sourceNodeId);
                            sourceNode.hasConnections = true;
                            
                            // Process each match
                            for (const match of matches) {
                                // Skip matches with Y sites
                                if ((match.target_site && match.target_site[0] === 'Y') || 
                                    match.site_type === 'Y' || 
                                    match.ResidueType === 'Y') continue;
                                
                                // Get target info
                                const targetUniprot = match.target_uniprot;
                                let targetSite = match.target_site;
                                const rmsd = parseFloat(match.rmsd);
                                
                                // Skip if missing essential info or invalid RMSD
                                if (!targetUniprot || !targetSite || isNaN(rmsd)) {
                                    console.warn(`Missing or invalid match data for ${siteName}:`, match);
                                    continue;
                                }
                                
                                // Determine target site type (multiple possible sources)
                                let targetSiteType = match.site_type || match.ResidueType || 
                                                  (targetSite && targetSite[0].match(/[STY]/) ? targetSite[0] : 'S');
                                
                                // Skip Tyrosine sites
                                if (targetSiteType === 'Y') continue;
                                
                                // If target site doesn't have the site type prefix, add it
                                if (targetSite && !targetSite.match(/^[STY]/)) {
                                    targetSite = `${targetSiteType}${targetSite}`;
                                }
                                
                                // Extract target residue number 
                                let targetResno;
                                if (targetSite.match(/\d+/)) {
                                    targetResno = parseInt(targetSite.match(/\d+/)[0]);
                                } else if (match.ResidueNumber) {
                                    targetResno = parseInt(match.ResidueNumber);
                                } else {
                                    // Try to extract from target_id if available
                                    const targetIdParts = match.target_id ? match.target_id.split('_') : [];
                                    if (targetIdParts.length >= 2) {
                                        targetResno = parseInt(targetIdParts[1]);
                                    }
                                    
                                    if (isNaN(targetResno)) {
                                        console.warn(`Cannot determine target residue number: ${targetSite}`);
                                        continue;
                                    }
                                }
                                
                                // Create target node ID
                                const targetNodeId = `${targetUniprot}_${targetResno}`;
                                
                                // Skip self-references
                                if (sourceNodeId === targetNodeId) continue;
                                
                                // Extract kinase information 
                                const knownKinase = extractKinaseInfo(match);
                                
                                // Add target node if not already present
                                if (!nodeMap.has(targetNodeId)) {
                                    // Use the raw match data to ensure all fields are available
                                    const targetNode = {
                                        id: targetNodeId,
                                        name: targetSite,
                                        uniprot: targetUniprot,
                                        type: 'match',
                                        isKnown: match.is_known_phosphosite === 1 || match.is_known === true,
                                        is_known_phosphosite: match.is_known_phosphosite || 0,
                                        siteType: targetSiteType,
                                        rmsd: rmsd,
                                        // Comprehensive extraction of all possible fields
                                        known_kinase: knownKinase,
                                        motif: match.motif || '',
                                        mean_plddt: match.mean_plddt || match.site_plddt || '',
                                        meanPLDDT: match.mean_plddt || match.site_plddt || '',
                                        nearby_count: match.nearby_count || match.NeighborCount || '',
                                        nearbyCount: match.nearby_count || match.NeighborCount || '',
                                        surface_accessibility: match.surface_accessibility || 
                                                             (match.HydroxylExposure ? match.HydroxylExposure * 100 : ''),
                                        // Store the full original match data for complete access
                                        matchData: match,
                                        size: 8,  // Slightly smaller for match sites
                                        hasConnections: true
                                    };
                                    
                                    nodes.push(targetNode);
                                    nodeMap.set(targetNodeId, targetNode);
                                }
                                
                                // Create link with explicit source and target
                                links.push({
                                    source: sourceNodeId,
                                    target: targetNodeId,
                                    rmsd: rmsd
                                });
                                
                                console.log(`Added link: ${sourceNodeId} -> ${targetNodeId}, RMSD: ${rmsd}`);
                            }
                        }
                    }
                }
            } catch (e) {
                console.error("Error processing structural match data:", e);
                console.error(e.stack); // Add stack trace for better debugging
            }
        } else {
            console.warn("Structural match data element not found");
        }
        
        // If we didn't find matches from the data element, try from DOM elements
        if (links.length === 0) {
            console.log("No links from data element, trying to extract from DOM");
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
                if (isNaN(siteResno)) {
                    console.warn(`Invalid residue number for site: ${siteName}`);
                    return;
                }
                
                // Create source node ID with just the number
                const sourceNodeId = `${proteinUniprotId}_${siteResno}`;
                
                // Make sure source node exists
                if (!nodeMap.has(sourceNodeId)) {
                    // If the source node doesn't exist, create it
                    console.log(`Source node ${sourceNodeId} not found, creating it`);
                    
                    const sourceNode = {
                        id: sourceNodeId,
                        name: siteName,
                        uniprot: proteinUniprotId,
                        type: 'protein',
                        isKnown: false,
                        siteType: siteName[0],
                        size: 10,
                        hasConnections: false
                    };
                    
                    nodes.push(sourceNode);
                    nodeMap.set(sourceNodeId, sourceNode);
                }
                
                // Mark this source node as having connections
                const sourceNode = nodeMap.get(sourceNodeId);
                sourceNode.hasConnections = true;
                
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
                    if (rmsdFilter) {
                        rmsdFilter.min = "0.1";
                        rmsdFilter.max = "10.0";
                        if (!rmsdFilter.value) {
                            rmsdFilter.value = "4.0"; // Default value
                            const rmsdValueElement = document.getElementById('rmsd-value');
                            if (rmsdValueElement) {
                                rmsdValueElement.textContent = "4.0 Å";
                            }
                        }
                    }
                    const currentThreshold = rmsdFilter ? parseFloat(rmsdFilter.value) : 10.0;
                    if (rmsd > currentThreshold) return;
                    
                    // Extract target residue number without letter
                    const targetResno = parseInt(targetSite.replace(/[A-Z]/g, ''));
                    if (isNaN(targetResno)) {
                        console.warn(`Invalid target residue number: ${targetSite}`);
                        return;
                    }
                    
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
                            size: 8,
                            hasConnections: true
                        };
                        
                        nodes.push(targetNode);
                        nodeMap.set(targetNodeId, targetNode);
                    }
                    
                    // Add link with explicit source and target as strings
                    links.push({
                        source: sourceNodeId,
                        target: targetNodeId,
                        rmsd: rmsd
                    });
                    
                    console.log(`Added link from DOM: ${sourceNodeId} -> ${targetNodeId}, RMSD: ${rmsd}`);
                });
            });
        }
        
        // Only include nodes that have connections
        const connectedNodes = nodes.filter(node => node.hasConnections || node.type === 'match');
        
        // Log final counts and verify data
        console.log(`Final structural network: ${connectedNodes.length} nodes, ${links.length} links`);
        
        // Verify that all link sources and targets exist as nodes
        const nodeIds = new Set(connectedNodes.map(node => node.id));
        const invalidLinks = links.filter(link => {
            const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
            const targetId = typeof link.target === 'object' ? link.target.id : link.target;
            return !nodeIds.has(sourceId) || !nodeIds.has(targetId);
        });
        
        if (invalidLinks.length > 0) {
            console.warn(`Found ${invalidLinks.length} invalid links with missing nodes`);
            // Filter out invalid links
            const validLinks = links.filter(link => {
                const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
                const targetId = typeof link.target === 'object' ? link.target.id : link.target;
                return nodeIds.has(sourceId) && nodeIds.has(targetId);
            });
            console.log(`Filtered to ${validLinks.length} valid links`);
            
            return { nodes: connectedNodes, links: validLinks };
        }
        
        if (connectedNodes.length > 0 && links.length > 0) {
            return { nodes: connectedNodes, links };
        }
        
        return null;
    } catch (error) {
        console.error("Error extracting structural network data:", error);
        console.error(error.stack); // Add stack trace for better debugging
        return null;
    }
}

// Function to filter network by RMSD threshold
// Function to filter network by RMSD threshold
function updateNetworkFilter() {
    if (!window.networkNodes || !window.networkLinks) {
        console.error("Network elements not available for filtering");
        return;
    }
    
    const filterElement = document.getElementById('rmsd-filter');
    
    // If filter element doesn't exist, use default threshold
    const threshold = filterElement ? parseFloat(filterElement.value) : 4.0;
    console.log(`RMSD threshold set to: ${threshold} Å`);
    
    // First show all nodes and links to reset any previous filtering
    window.networkNodes.style('display', null);
    window.networkLinks.style('display', null);
    
    // Filter links by RMSD
    window.networkLinks.style('display', function(d) {
        return d.rmsd <= threshold ? null : 'none';
    });
    
    // Create a set to track nodes that should be visible
    const visibleNodeIds = new Set();
    
    // Add all protein nodes as they should always be visible
    window.networkNodes.each(function(d) {
        if (d.type === 'protein') {
            visibleNodeIds.add(d.id);
        }
    });
    
    // Add match nodes that have connections below the threshold
    window.networkLinks.each(function(d) {
        if (d.rmsd <= threshold) {
            // Handle both object and string source/target
            const sourceId = typeof d.source === 'object' ? d.source.id : d.source;
            const targetId = typeof d.target === 'object' ? d.target.id : d.target;
            visibleNodeIds.add(sourceId);
            visibleNodeIds.add(targetId);
        }
    });
    
    // Update node visibility based on our visibility set
    window.networkNodes.style('display', function(d) {
        return visibleNodeIds.has(d.id) ? null : 'none';
    });
    
    // Update label visibility
    if (window.networkLabels) {
        window.networkLabels.style('display', function(d) {
            return visibleNodeIds.has(d.id) ? null : 'none';
        });
    }
    
    // Count nodes and links after filtering
    const filteredVisibleLinks = window.networkLinks.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    const filteredVisibleNodes = window.networkNodes.filter(function() {
        return window.getComputedStyle(this).display !== 'none';
    }).size();
    
    console.log(`Links: ${window.networkLinks.size()} → ${filteredVisibleLinks} after filtering`);
    console.log(`Nodes: ${window.networkNodes.size()} → ${filteredVisibleNodes} after filtering`);
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