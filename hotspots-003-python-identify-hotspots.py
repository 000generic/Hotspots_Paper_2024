#! python
input_positions = open('output/1-gtfs-parsed-for-genome-positions-of-genes', 'r')
input_keepers = open( 'output/2-species16_ALL_X_blast-hits-better-than-1e-60-KEEP', 'r')
output_hotspots = open('output/3-species16-hotspots-1e-60-window200', 'w')
output_hotspots_short = open('output/3-species16-hotspots-1e-60-window200-SHORT', 'w')
output_gene_errors = open('output/3-errors-gene-hotspots-1e-60-window200', 'w')
output_position_errors = open('output/3-position-errors-hotspots-1e-60-window200', 'w')
output_human_gtf_issues =  open( 'output/3-errors-gtf-vs-proteome-human_gtf_issues', 'w' )

genes_positions = {}
positions_genes = {}
chromosomes_positions = {}

GSpp_gspp = {}
GSpp_gspp[ 'HSap' ] = 'Homo_sapiens'
GSpp_gspp[ 'DMel' ] = 'Drosophila_melanogaster'
GSpp_gspp[ 'CEle' ] = 'Caenorhabditis_elegans'
GSpp_gspp[ 'OBim' ] = 'Octopus_bimaculoides'
GSpp_gspp[ 'PPel' ] = 'Patella_pellucida'
GSpp_gspp[ 'PCae' ] = 'Patella_caerulea'
GSpp_gspp[ 'PVul' ] = 'Patella_vulgata'
GSpp_gspp[ 'CSqu' ] = 'Chrysomallon_squamiferum'
GSpp_gspp[ 'GAeg' ] = 'Gigantopelta_aegis'
GSpp_gspp[ 'ACal' ] = 'Aplysia_californica'
GSpp_gspp[ 'MAre' ] = 'Mya_arenaria'
GSpp_gspp[ 'PMax' ] = 'Pecten_maximus'
GSpp_gspp[ 'MTro' ] = 'Mytilus_trossulus'
GSpp_gspp[ 'OEdu' ] = 'Ostrea_edulis'
GSpp_gspp[ 'CVir' ] = 'Crassostrea_virginica'
GSpp_gspp[ 'CGig' ] = 'Crassostrea_gigas'

# ../../../species16/gffs/output/CVir_simplified.gff
# CVir_NC_035780.1        LOC111110729    43111   66897   -
# Patella_caerulea__SNE40_000001

chopped_to_full = {}  # human gtf issue !!!
for genus_positions in input_positions:  # 
    input_genus_positions = open( genus_positions[ :-1 ], 'r' )
    GSpp = genus_positions.split( '/' )[ -1 ].split( '.' )[ 0 ] # GSpp not genus
    gspp = GSpp_gspp[ GSpp ]
    for position in input_genus_positions:  

        if position[ 0 ] != '#':
            info = position.strip().split( '\t' )
            if len( info ) > 1:
                chromosome = '_'.join( info[ 0 ].split( '_' )[ 1: ] )
                start = info[ 2 ]
                stop = info[ 3 ]
                strand = info[ 4 ]
                gene = info[ 1 ]

                
                if gspp == 'Homo_sapiens': # human gtf issue!!!!
                    
                    gene_chopped = '-'.join( gene.split( '-')[ :-1 ] )
                    chopped_to_full[ gspp + '__' + gene_chopped ] = gspp + '__' +gene
                    

                gspp_gene = gspp + '__' + gene
                gspp_chromosome = gspp + '__' + chromosome
                position_stats = gspp + '\t' + chromosome + '\t' + start + '\t' + stop

                if gspp_gene not in genes_positions:
                    genes_positions[ gspp_gene ] = position_stats
                else:
                    error = 'gene ' + gspp_gene + ' already in genes_positions'
                    print( error )
                    print( info )
                if position_stats not in positions_genes:
                    positions_genes[ position_stats ] = gspp_gene
                    if gspp_chromosome not in chromosomes_positions:
                        chromosomes_positions[ gspp_chromosome ] = []
                    chromosomes_positions[ gspp_chromosome ].append( position_stats )
                else:
                    error = 'position ' + position_stats + ' already in positions_genes'
                    print( error  )
                    print( info )
genus_hot_spots = {}
hot_spot_id = 0
#Patella_caerulea__SNE40_000001
for keeper in input_keepers:  # Aplysia_californica     Aplysia_californica__LOC118478065__XP_035826908.1       Aplysia_californica__LOC118478065__XP_035826908.1       2.13e-106
    
    keeper = keeper.strip()
    info = keeper.split('\t')
    genus = info[ 0 ]
    query_gene_pep = info[1]
    query = genus + '__' + query_gene_pep.split( '__' )[ 1 ]
    hit_gene_pep = info[2] 
    hit = genus + '__' + hit_gene_pep.split( '__' )[ 1 ]
    evalue = info[3]
    
    if genus not in genus_hot_spots:
        genus_hot_spots[ genus ] = {}

    if query in genes_positions.keys():
        query_position_stats = genes_positions[ query ]
    else:
        if genus == 'Homo_sapiens':  # this is dealing with what I think are lncRNAs in the GTF...need to understand better in the future
            # query = genus + '__' + '-'.join( query_gene_pep.split( '__' )[ 1 ].split( '-' )[ :-1 ] )
            query = chopped_to_full[ query ]
            query_position_stats = genes_positions[ query ]
            output_human_gtf_issues.write( keeper )
        else:
            print( 'Error in gene id of ' + genus ) # this should not be happening for most specices - only human seems to have dash issues so far
            
    query_chromosome = query_position_stats.split( '\t' )[ 1 ]
    query_genus_chromosome = genus + '__' + query_chromosome
    query_index = chromosomes_positions[ query_genus_chromosome ].index( query_position_stats )

    if hit in genes_positions.keys():
        hit_position_stats = genes_positions[ hit ]
    else:
        if genus == 'Homo_sapiens':  # this is dealing with what I think are lncRNAs in the GTF...need to understand better in the future
            #print( hit + '\n' )
            #hit = genus + '__' + '-'.join( hit_gene_pep.split( '__' )[ 1 ].split( '-' )[ :-1 ] )
            hit = chopped_to_full[ hit ]
            hit_position_stats = genes_positions[ hit ]
            output_human_gtf_issues.write( keeper )
        else:
            print( 'Error in gene id of ' + genus ) # this should not be happening for most specices - only human seems to have dash issues so far
    upstream_gene = ''
    down_stream_gene = ''

    # Check for homolog within 10 or 100 or other number of genes to either side of query gene on scaffold
    first_gene_position = 0
    final_gene_on_chromosome = chromosomes_positions[ query_genus_chromosome ][ -1 ]  # Achatina	chr1_arrow_pilon	118799724	118817081
    final_gene_position = chromosomes_positions[ query_genus_chromosome ].index( final_gene_on_chromosome )

    if query_index - 100 < first_gene_position:
        hotspot_range_start = first_gene_position
    else:
        hotspot_range_start = query_index - 100
    if query_index + 100 > final_gene_position:
        hotspot_range_stop = final_gene_position
    else:
        hotspot_range_stop = query_index + 100
    hotspot_range = range( hotspot_range_start, hotspot_range_stop )

    # Evaluate genes neighboring query on scaffold for query blast hits at evalue threshold or better
    if genus == 'Acanthopleura':  # TODO need to sort out errors in processing Acanthopleura
        continue

    if len( chromosomes_positions[ query_genus_chromosome ] ) == 1:
        continue

    for gene_index in hotspot_range:
        if gene_index == query_index:
            continue
        else:
            gene_position_stats = chromosomes_positions[ query_genus_chromosome ][ gene_index ]
            if hit_position_stats == gene_position_stats:
                add_to_hot_spot = True
                for hot_spot in genus_hot_spots[ genus ].keys():
                        if query in genus_hot_spots[ genus ][ hot_spot ]:
                            genus_hot_spots[ genus ][ hot_spot ].append( hit )
                            genus_hot_spots[ genus ][ hot_spot ] = list( set( genus_hot_spots[ genus ][ hot_spot ] ) )
                            add_to_hot_spot = False

                        elif hit in genus_hot_spots[ genus ][ hot_spot ]:
                            genus_hot_spots[ genus ][ hot_spot ].append( query )
                            genus_hot_spots[ genus ][ hot_spot] = list( set( genus_hot_spots[ genus ][ hot_spot ] ) )
                            add_to_hot_spot = False
                        else:
                            continue
                if add_to_hot_spot:
                    genus_hot_spots[ genus ][ hot_spot_id ] = [ query, hit ]
                    hot_spot_id += 1

genus_hotspots_unique = {}
for genus in genus_hot_spots.keys():
    hot_spot_count = 0
    remove_hotspots = []

    for hot_spot_probe in genus_hot_spots[ genus ].keys():
        for gene in genus_hot_spots[ genus ][ hot_spot_probe ]:
            for hot_spot_test in genus_hot_spots[ genus ].keys():
                if hot_spot_test in remove_hotspots:
                    continue
                else:
                    if hot_spot_test != hot_spot_probe:
                        if gene in genus_hot_spots[ genus ][ hot_spot_test ]:
                            for gene_test in genus_hot_spots[ genus ][ hot_spot_test ]:
                                if gene_test not in genus_hot_spots[ genus ][ hot_spot_probe ]:
                                    genus_hot_spots[ genus ][ hot_spot_probe ].append( gene_test )
                            remove_hotspots.append( hot_spot_test )

    for remove_hotspot in list( set( remove_hotspots ) ):
        del( genus_hot_spots[ genus ][ remove_hotspot  ] )

    seeds_hotspots = {}
    for hot_spot in sorted( genus_hot_spots[genus].keys() ):
        random_order = []
        for genus_gene in genus_hot_spots[genus][ hot_spot ]:
            gene = genus_gene.split( '__' )[ 1 ]
            position = '__'.join( genes_positions[ genus_gene ].split( '\t' ) )
            position_gene = position + '__' + gene
            random_order.append( position_gene )
        sorted_order = sorted( random_order )

        seed = sorted_order[ 0 ]  # Aplysia__NW_004797388.1__000516733__000521123__LOC101848961
        seeds_hotspots[ seed ] = sorted_order

    for seed in sorted( seeds_hotspots.keys() ):
        ordered_hot_spot = sorted( seeds_hotspots[ seed ] )

        chromosome = ordered_hot_spot[ 0 ].split( '__' )[ 1 ]
        hotspot_start = int( ordered_hot_spot[ 0 ].split( '__' )[ 2 ] )
        hotspot_stop = int( ordered_hot_spot[ -1 ].split( '__' )[ 3 ] )
        hot_spot_location =  'Location_' + chromosome + '_' + str( hotspot_start ) + '_' + str( hotspot_stop )

        hot_spot_count += 1
        hot_spot_id = 'hotspot_e60_w200_' + genus + '_' + str( hot_spot_count )
        hot_spot_size = 'Paralogs_' + str(len( ordered_hot_spot ))

        output = genus + '\t' + hot_spot_id + '\t' + hot_spot_size + '\t' + hot_spot_location + '\t'
        for position_gene in ordered_hot_spot:
            output = output + position_gene + ', '
        output = output[:-2] + '\n'
        output_hotspots.write(output)

        output = genus + '\t' + hot_spot_id + '\t' + hot_spot_size + '\t' + hot_spot_location + '\t'
        for position_gene in ordered_hot_spot:
            gene = position_gene.split( '__' )[ -1 ]
            output = output + gene + ', '
        output = output[:-2] + '\n'
        output_hotspots_short.write(output)

input_positions.close()
input_keepers.close()
output_hotspots.close()
output_gene_errors.close()
output_position_errors.close()
output_hotspots_short.close()
