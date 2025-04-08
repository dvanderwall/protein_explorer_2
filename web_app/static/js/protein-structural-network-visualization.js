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
        return '#9E9E9E'; // Gray for structurally similar sites without known kinases
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
    networkSimulation = d3.forceSimulation(networkData.nodes)
        .force('link', d3.forceLink(networkData.links)
            .id(d => d.id)
            .distance(d => d.rmsd ? d.rmsd * 50 : 100))
        .force('charge', d3.forceManyBody()
            .strength(d => d.type === 'protein' ? -50 : -25))
        .force('center', d3.forceCenter(width / 2, height / 2))
        .force('collision', d3.forceCollide().radius(d => d.size + 2))
        .force('x', d3.forceX(width / 2).strength(0.2))
        .force('y', d3.forceY(height / 2).strength(0.2))
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
        let content = '';
        if (d.type === 'protein') {
            content = `
                <h6 class="border-bottom pb-2 mb-2">${d.name} - ${d.uniprot}</h6>
                <p><strong>Site Type:</strong> ${d.siteType}</p>
                <p><strong>Known Site:</strong> ${d.isKnown ? 'Yes' : 'No'}</p>
                ${d.known_kinase ? `<p><strong>Known Kinase${d.known_kinase.includes(',') ? 's' : ''}:</strong> ${formatKinaseList(d.known_kinase)}</p>` : ''}
                ${d.meanPLDDT || d.mean_plddt ? `<p><strong>Mean pLDDT:</strong> ${d.meanPLDDT || d.mean_plddt}</p>` : ''}
                ${d.nearbyCount || d.nearby_count ? `<p><strong>Nearby Residues:</strong> ${d.nearbyCount || d.nearby_count}</p>` : ''}
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
            
            // Create node ID - IMPORTANT: Using resno without letter for consistency
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
                            
                            // Extract site number without letter consistently
                            const siteResno = parseInt(siteName.replace(/[A-Z]/g, ''));
                            if (isNaN(siteResno)) {
                                console.warn(`Invalid residue number for site: ${siteName}`);
                                continue;
                            }
                            
                            // Create source node ID with just the number
                            const sourceNodeId = `${siteName}`;

                            console.log(`Processing ${matches.length} matches for site ${siteName}`);
                            
                            // Ensure source node exists
                            if (!nodeMap.has(sourceNodeId)) {
                                // If the source node doesn't exist, create it
                                console.log(`Source node ${sourceNodeId} not found, creating it`);
                                
                                // Create minimal source node
                                const sourceNode = {
                                    id: sourceNodeId,
                                    name: siteName,
                                    uniprot: proteinUniprotId,
                                    type: 'protein',
                                    isKnown: false,
                                    siteType: siteName[0],
                                    size: 10
                                };
                                
                                nodes.push(sourceNode);
                                nodeMap.set(sourceNodeId, sourceNode);
                            }
                            
                            // Process each match
                            for (const match of matches) {
                                // Skip matches with Y sites
                                if (match.target_site && match.target_site[0] === 'Y') continue;
                                
                                // Get target info
                                const targetUniprot = match.target_uniprot;
                                const targetSite = match.target_site;
                                const rmsd = parseFloat(match.rmsd);
                                
                                // Skip if missing essential info or invalid RMSD
                                if (!targetUniprot || !targetSite || isNaN(rmsd)) {
                                    console.warn(`Missing or invalid match data for ${siteName}:`, match);
                                    continue;
                                }
                                
                                // Extract target residue number without letter
                                const targetResno = parseInt(targetSite.replace(/[A-Z]/g, ''));
                                if (isNaN(targetResno)) {
                                    console.warn(`Invalid target residue number: ${targetSite}`);
                                    continue;
                                }
                                
                                // Create target node ID with just the number
                                const targetNodeId = `${targetUniprot}_${targetResno}`;
                                
                                // Skip self-references
                                if (sourceNodeId === targetNodeId) continue;
                                
                                // Extract kinase information using Python's naming convention
                                const knownKinase = match.known_kinase;
                                
                                // Add target node if not already present
                                if (!nodeMap.has(targetNodeId)) {
                                    const targetNode = {
                                        id: targetNodeId,
                                        name: targetSite,
                                        uniprot: targetUniprot,
                                        type: 'match',
                                        isKnown: match.is_known_phosphosite === 1 || match.is_known === true,
                                        siteType: targetSite[0],
                                        rmsd: rmsd,
                                        known_kinase: knownKinase,
                                        motif: match.motif || '',
                                        mean_plddt: match.mean_plddt || match.site_plddt || '',
                                        meanPLDDT: match.mean_plddt || match.site_plddt || '',
                                        nearby_count: match.nearby_count || '',
                                        nearbyCount: match.nearby_count || '',
                                        surface_accessibility: match.surface_accessibility || '',
                                        surfaceAccessibility: match.surface_accessibility || '',
                                        size: 8  // Slightly smaller for match sites
                                    };
                                    
                                    nodes.push(targetNode);
                                    nodeMap.set(targetNodeId, targetNode);
                                }
                                
                                // Create link with explicit source and target as strings
                                // This is critical for D3 to properly connect nodes
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
                        size: 10
                    };
                    
                    nodes.push(sourceNode);
                    nodeMap.set(sourceNodeId, sourceNode);
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
                    if (rmsdFilter) {
                        rmsdFilter.min = "0.1";
                        rmsdFilter.max = "10.0";
                        rmsdFilter.value = "4.0"; // Default value
                        document.getElementById('rmsd-value').textContent = "4.0 Å";
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
                            size: 8
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
        
        // Log final counts and verify data
        console.log(`Final structural network: ${nodes.length} nodes, ${links.length} links`);
        
        // Verify that all link sources and targets exist as nodes
        const invalidLinks = links.filter(link => {
            return !nodeMap.has(link.source) || !nodeMap.has(link.target);
        });
        
        if (invalidLinks.length > 0) {
            console.warn(`Found ${invalidLinks.length} invalid links with missing nodes`);
            // Filter out invalid links
            const validLinks = links.filter(link => nodeMap.has(link.source) && nodeMap.has(link.target));
            console.log(`Filtered to ${validLinks.length} valid links`);
            links = validLinks;
        }
        
        if (nodes.length > 0 && links.length > 0) {
            // Create a set of all connected node IDs
            const connectedNodeIds = new Set();
            
            // Add all nodes that appear in links
            links.forEach(link => {
                const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
                const targetId = typeof link.target === 'object' ? link.target.id : link.target;
                connectedNodeIds.add(sourceId);
                connectedNodeIds.add(targetId);
            });
            
            // Filter nodes to only include those with connections
            // Exception: if there are no links at all, show all protein nodes
            const filteredNodes = links.length > 0 
                ? nodes.filter(node => connectedNodeIds.has(node.id)) 
                : nodes.filter(node => node.type === 'protein');
            
            console.log(`Filtered from ${nodes.length} nodes to ${filteredNodes.length} connected nodes`);
            
            return { nodes: filteredNodes, links };
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
    if (!window.networkNodes || !window.networkLinks) return;
    
    const threshold = parseFloat(document.getElementById('rmsd-filter').value);
    console.log(`Filtering network with RMSD threshold: ${threshold}`);
    
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