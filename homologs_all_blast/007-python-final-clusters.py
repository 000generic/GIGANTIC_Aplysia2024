#! python

input_clusters_list = open( 'output/6-list-aaabs-per-species-evalue-minus10', 'r' )
output_sequences = open( 'output/7-output-evalue-minus10/7-blast-clusters-species19-sequences', 'w' )
output_clusters = open( 'output/7-output-evalue-minus10/7-blast-clusters-species19-clusters', 'w' )

# output/5-all-species-all-clusters-all-homologs-Homo
clusters_sequences = {}
sequences_clusters = {}
counter = 0
for next_clusters in input_clusters_list:

    print( next_clusters )
    
    input_clusters = open( next_clusters[ :-1 ], 'r' )
    for next_cluster in input_clusters:

        counter = counter + 1
        print( str( counter ) )
        clusterid = 'bc' + str( counter ).zfill( 10 )
        
        info = next_cluster[ :-1 ].split( '\t' )
        seqids = info[ -1 ].split( ', ' )
        outside_clusters = []
        clusters_sequences[ clusterid ] = []
        for next_seq in seqids:
            clusters_sequences[ clusterid ].append( next_seq )
    
            if next_seq in sequences_clusters.keys():
                clusterid_outside = sequences_clusters[ next_seq ]
                if clusterid_outside != clusterid: 
                    outside_clusters.append( clusterid_outside )
                    
        outside_clusters = list( set( outside_clusters ) )

        for next_cluster_outside in outside_clusters:
            clusters_sequences[ clusterid ] = list( set( clusters_sequences[ clusterid ] + clusters_sequences[ next_cluster_outside ] ) )
            clusters_sequences.pop( next_cluster_outside )

        for next_seq in clusters_sequences[ clusterid ]:
            sequences_clusters[ next_seq ] = clusterid

species = [ 'Homo', 'Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Nautilus', 'Argonauta', 'Octopus', 'Gigantopelta', 'Lottia', 'Pomacea', 'Achatina', 'Elysia', 'Aplysia' ]
counter = 0
species_clusters_sequences = {}
for next_species in species:
    species_clusters_sequences[ next_species ] = {}

for next_cluster in clusters_sequences.keys():

    counter = counter + 1
    clusterid = 'bc' + str( counter ).zfill( 10 )
    
    seqids = clusters_sequences[ next_cluster ]
    total_sequences = str( len( seqids ) )
    cluster_sequences = ', '.join( seqids )
    
    output = clusterid + '\t' + total_sequences + '\t' + cluster_sequences + '\n'
    output_clusters.write( output )
    
    for next_seq in seqids:
        info = next_seq.split( '-' )
        genus = info[ 5 ]
        species_clusters_sequences[ genus ][ clusterid ] = seqids
        output = next_seq + '\t' + clusterid + '\t' + total_sequences +  '\n'
        output_sequences.write( output )

#species_output = [ 'Aplysia', 'Octopus' ]
for next_species in sorted( species_clusters_sequences.keys() ):
    new_file = 'output/7-blast-clusters-species19-sequences-' + next_species
    output_sequences = open( new_file, 'w' )
    for next_cluster in species_clusters_sequences[ next_species ]:
        seqids = species_clusters_sequences[ next_species ][ next_cluster ]
        cluster_sequences = ', '.join( seqids )
        for next_seq in seqids:
            if len( next_seq.split( next_species ) ) > 1:
                output = next_seq + '\t' + next_cluster +  '\n'
                output_sequences.write( output )
    output_sequences.close()
    
input_clusters_list.close()
output_clusters.close()
output_sequences.close()
