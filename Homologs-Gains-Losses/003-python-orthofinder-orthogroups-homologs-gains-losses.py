#! python

input_map_nodes_to_clades = open( 'list-nodes-clades', 'r' )
input_map_clades_to_ancestors = open( 'list-clades-ancestors', 'r' )
input_map_clades_to_names = open( 'list-clades-names', 'r' )
input_nodes = open( 'output/2-list-Ns', 'r' )
output_counts = open( 'output/3-orthofinder-orthogroups-gains-losses-per-species-tree-node', 'w' )

nodes_clades = {}
clades_names = {}
clades_ancestors = {}
clades_orthogroups_conserved = {}
clades_orthogroups_gains = {}
clades_orthogroups_losses = {}
clades_homologs = {}

# read in node to clade map
for next_line in input_map_nodes_to_clades:
    info = next_line[ :-1 ].split( 't')
    node = info[ 0 ]
    clade = info[ 1 ]
    nodes_clades[ node ] = clade

# read in clade to ancestor map
for next_line in input_map_clades_to_ancestors:
    info = next_line[ :-1 ].split( 't')
    clade = info[ 0 ]
    ancestor = info[ 1 ]
    clades_ancestors[ clade ] = ancestor

# read in clade to name map
for next_line in input_map_clades_to_names:
    info = next_line[:-1].split('t')
    clade = info[ 0 ]
    name = info[ 1 ]
    clades_names[ clade ] = name

# set up lists for each clade for conserved, gained, and lost orthogroups per OrthoFinder species tree node clustering
for next_clade in clades_names.keys():
    clades_orthogroups_total[ next_clade ] = []
    clades_orthogroups_gains[ next_clade ] = []
    clades_orthogroups_losses[ next_clade ] = []

# read in orthogroups per OrthoFinder species tree node clustering into dictionary
nodes_orthogroups = {}
for next_node in input_nodes:
    next_node = next_node[ :-1 ]
    node_id = next_node.split( '.' )[ 0 ]
    clade_id = nodes_clades[ node_id ]
    clade_name = clade_name[ clade_id ]
    ancestor_id = clades_ancestors[ clade_id ]

    input_node = open( next_node, 'r' )
    for next_line in input_node:
        next_line = next_line[ :-1 ]
        info = next_line.split( '\t' )
        orthogroup_id = info[ 1 ]
        seqids = ', '.join( info[ 3: ] )
        info_seqids = seqids.split( ', ' ) # dictionary of seqids that form the orthogroup

        # ancestral Bilateria N0 orthogroups C37

        if node_id == 'N0':
            clade_1 = [ 'Homo' ] # C19
            clade_2 = [ 'Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                        'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta',
                        'Pomacea', 'Achatina', 'Elysia', 'Aplysia'] # C36
            test_clade_1 = False
            test_clade_2 = False

            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:
                clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Protostomia N1 orthogroups C36
        if node_id == 'N1':
            clade_1 = ['Drosophila', 'Caenorhabditis']  # C35
            clade_2 = ['Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus', 'Argonauta', 'Octopus',
                       'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea', 'Achatina',
                       'Elysia', 'Aplysia']  # C34
            test_clade_1 = False
            test_clade_2 = False

            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True

            if test_clade_1 == True and test_clade_2 == True:  # orthogroup conserved: test if orthogroup in both clade 1 and 2
                clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                if orthogroup_id not in clades_orthogroups_total[ ancestor_id ]: # orthogroup gain: test if orthogroup in clade 1 and 2 but not ancestor
                    clades_orthogroups_gains[ clade_id ].append( orthogroup_id )

            if test_clade_1 == True and test_clade_2 = False:  # orthogroup conserved: test if orthogroup is present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )

            if test_clade_1 == True and test_clade_2 = False:  # orthogroup conserved: test if orthogroup is present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )

            for next_orthogroup in clades_orthogroups_total[ ancestor_id ]:  # orthogroup loss: test if orthogroup present in ancestor but not clades 1 and 2
                if next_orthogroup not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_losses[ clade_id ].append( orthogroup_id )


output = 'Clade ID' + '\t' + 'Clade Name' + '\t' + 'Orthogroups Total' + '\t' + 'Orthogroups Conserved' + '\t' + 'Orthogroups Gained' + '\t' + 'Orthogroups Lost' + '\n'
output_counts.write( output )
for next_clade in clades_orthogroups_total.keys():
    ancestral_orthogroups = str( len( list( set( clades_orthogroups_total[ next_clade ] ) ) ) )
    gained_orthogroups =  str( len( list( set( clades_orthogroups_gains[ next_clade ] ) ) ) )
    conserved_orthogroups = str( len( list( set( clades_orthogroups_total[ next_clade ] ) ) )  -  len( list( set( clades_orthogroups_gains[ next_clade ] ) ) ) )
    lost_orthogroups = str( len( list( set( clades_orthogroups_losses[ next_clade ] ) ) ) )
    clade_name = clades_names[ next_clade ]
    output = next_clade + '\t' + clade_name  + '\t' +  ancestral_orthogroups + '\t' conserved_orthogroups + '\t' + gained_orthogroups + '\t' + lost_orthogroups + '\n'
    output_counts.write( output )

input_map_nodes_to_clades.close()
input_map_clades_to_names.close()
input_nodes.close()
output_counts.close()
