#! python

input_reports = open( 'output/4-list-blast-reports-Cyclina', 'r' )
output_clusters = open( 'output/5-all-species-all-clusters-all-homologs-Cyclina', 'w' )
output_map = open( 'output/5-map-final-clusterids-to-initial-cluster-ids-Cyclina', 'w' )

# output/Cyclina-granulata_X_Cyclina-granulata

print( 'Input reports' )

gspp_gene_hits = {}
for next_report in input_reports:

    report_path = next_report[ :-1 ]

    print( report_path )
    
    info_report = report_path.split( '_X_' )
    query_species = info_report[ 0 ].split( '/' )[ -1 ]
    subject_species = info_report[ -1 ]

    if query_species not in gspp_gene_hits.keys():
        gspp_gene_hits[ query_species ] = {}
                        
    # Metazoa-Annelida-Polychaeta-Order_unclassified14-Capitellidae-Capitella-teleta-gigantic20220627164247_seq0000244112aa   Metazoa-Annelida-Polychaeta-Order_unclassified14-Capitellidae-Capitella-teleta-gigantic20220627164247_seq0000255669aa   1.04e-21
                        
    input_report = open( report_path, 'r' )
    for next_line in input_report:
        
        info = next_line[ :-1 ].split( '\t' )
        query = info[ 0 ]
        hit = info[ 1 ]
        evalue = info[ 2 ] # can be 1) 0.0 or 2) 0.001 or 3) 2.10e-136
        info_evalue_e = evalue.split( 'e-' )
        
        keeper = False

        if evalue == '0.0': # keep evalues less than e-180 I think it is
            keeper = True
        
        elif len( info_evalue_e ) > 1:
            evalue_e_integer = int( info_evalue_e[ -1 ] )
            if evalue_e_integer > 9:
                keeper = True
            
        #elif len( info_evalue_e ) == 1:
        #    info_evalue_dot = evalue.split( '.')
        #    evalue_dot_string = info_evalue_dot[ -1 ]
        #    if len( evalue_dot_string ) > 2:
        #        keeper = True

        else:
            print( evalue )

        if keeper == True:
            if query not in gspp_gene_hits[ query_species ].keys():
                gspp_gene_hits[ query_species ][ query ] = []
                gspp_gene_hits[ query_species ][ query ].append( hit )
            else:
                gspp_gene_hits[ query_species ][ query ].append( hit )

print( 'Process hits' ) 

counter = 0
genes_clusters = {}
clusters_genes = {}
for next_gspp in gspp_gene_hits.keys(): # read through species in dictionary

    print( next_gspp )

    for next_gene in gspp_gene_hits[ next_gspp ]: # read through genes for species in dictionary

        print( next_gene )
        
        previous_clusterids = []
        counter = counter + 1
        new_cluster = 'initial_Cyclina_homolog_cluster_' + str( counter )
        clusters_genes[ new_cluster ] = []
        
        for next_hit in gspp_gene_hits[ next_gspp ][ next_gene ]: # read through blast report hits per gene per species

            if len( next_hit.split( 'Cyclina' ) ) > 1:
            
                if next_hit in genes_clusters.keys(): # grab previous cluster id if there is one for the hit
                    previous_cluster = genes_clusters[ next_hit ]
                    previous_clusterids.append( previous_cluster )
                genes_clusters[ next_hit ] = new_cluster # update cluster id for the hit to new cluster id
            clusters_genes[ new_cluster ].append( next_hit ) # add the hit gene id to the list of gene ids for the new cluster id
                
        all_previous_clusters = list( set( previous_clusterids ) ) # make the set of all previous cluster ids for all hits to the gene in the species
        for next_cluster in all_previous_clusters:
            if next_cluster == new_cluster: # protect new_cluster just in case
                pass
            else:
                for next_gene in clusters_genes[ next_cluster ]:
                    if len( next_gene.split( 'Cyclina' ) ) > 1:
                        genes_clusters[ next_gene ] = new_cluster # update all genes in all previous clusters containing a hit to the gene in the species to the new cluster id
                    clusters_genes[ new_cluster ].append( next_gene ) # add the hit gene ideas of prevous clusters to the new cluster id's list of genes
        clusters_genes[ new_cluster ] = list( set( clusters_genes[ new_cluster ] ) ) # make the list of genes for the new cluster id a non-redundant list

        for next_cluster in all_previous_clusters:
            if next_cluster == new_cluster: # protect new_cluster just in case
                pass
            else:
                clusters_genes.pop( next_cluster ) #remove previous cluster ids now that their genes have been updated to the new cluster id

print( 'Process clusters' )

counter = 0
for next_cluster in clusters_genes:

    counter = counter + 1

    final_cluster = 'AAAB_cluster_Cyclina_' + str( counter )
    output = final_cluster + '\t' + next_cluster + '\n'
    output_map.write( output )
    
    homolog_count = len( clusters_genes[ next_cluster ] )
    output = final_cluster + '\t' + str( homolog_count ) + '\t'
    for next_homolog in  clusters_genes[ next_cluster ]:
        output = output + next_homolog + ', '
    output = output[ :-2 ] + '\n'
    output_clusters.write( output )

input_reports.close()
output_clusters.close()
output_map.close()
