#! python
# Code Eric Edsinger
# Identify and count homologs, gains, and losses per species in OrthoFinder Species19 genomes clustering

input_map_nodes_to_clades = open( 'input/list-nodes-clades', 'r' )
input_map_clades_to_ancestors = open( 'input/list-clades-ancestors', 'r' )
input_map_clades_to_names = open( 'input/list-clades-names', 'r' )
input_nodes = open( 'output/2-list-Ns', 'r' )
output_counts = open( 'output/quick3-orthofinder-orthogroups-gains-losses-per-species-tree-node', 'w' )

######### Process input files

nodes_clades = {}
clades_names = {}
clades_ancestors = {}
clades_orthogroups_total = {}
clades_orthogroups_gains = {}
clades_orthogroups_losses = {}

print( 'Reading non-node files into dictionaries!')

# read in node to clade map
for next_line in input_map_nodes_to_clades:
    info = next_line[ :-1 ].split( '\t')
    node = info[ 0 ]
    clade = info[ 1 ]
    nodes_clades[ node ] = clade

# read in clade to ancestor map
for next_line in input_map_clades_to_ancestors:
    info = next_line[ :-1 ].split( '\t')
    clade = info[ 0 ]
    ancestor = info[ 1 ]
    clades_ancestors[ clade ] = ancestor

# read in clade to name map
for next_line in input_map_clades_to_names:
    info = next_line[:-1].split('\t')
    clade = info[ 0 ]
    name = info[ 1 ]
    clades_names[ clade ] = name

# set up lists for each clade for conserved, gained, and lost orthogroups per OrthoFinder species tree node clustering
for next_clade in clades_names.keys():
    clades_orthogroups_total[ next_clade ] = []
    clades_orthogroups_gains[ next_clade ] = []
    clades_orthogroups_losses[ next_clade ] = []

# read in orthogroups per OrthoFinder species tree node clustering into dictionary
for next_node in input_nodes:
    next_node = next_node[ :-1 ]
    node_id = next_node.split( '.' )[ 0 ]
    clade_id: str = nodes_clades[ node_id ]
    clade_name = clades_names[ clade_id ]
    if clade_id == 'C37':
        ancestor_id = 'NA'
    else:
        ancestor_id = clades_ancestors[ clade_id ]

########## Begin processing OrthoFinder clustering of species tree nodes N0-N17

    input_node = open( next_node, 'r' )
    print( 'Processing: ' + next_node )

    for next_line in input_node:
        next_line = next_line[ :-1 ]
        info = next_line.split( '\t' )
        orthogroup_id: str = info[ 1 ]
        seqids = ', '.join( info[ 3: ] )
        info_seqids = seqids.split( ', ' )  # dictionary of seqids that form the orthogroup

        # ancestral Bilateria N0 orthogroups C37

        if node_id == 'N0':

            # C37 Bilateria
            clade_1 = [ 'Homo' ]  # C19
            clade_2 = [ 'Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                        'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta',
                        'Pomacea', 'Achatina', 'Elysia', 'Aplysia']  # C36
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[clade_id].append(orthogroup_id)
                clades_orthogroups_total[ 'C19' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2
                clades_orthogroups_total[ 'C19' ].append( orthogroup_id )

########## Process Single-species lineage gains

            # C01 Aplysia
            clade_1 = ['Aplysia']  # C01
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta',
                       'Pomacea', 'Achatina', 'Elysia', 'Homo']  # All but Aplysia
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C01 Euopisthobranchia-Aplysia
                clades_orthogroups_total[ 'C01' ].append( orthogroup_id )

            # C02 Elysia
            clade_1 = ['Elysia']  # C02
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta',
                       'Pomacea', 'Achatina', 'Aplysia', 'Homo']  # All but Elysia
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C02 Panpulmonata-Sacoglossa-Elysia
                clades_orthogroups_total[ 'C02' ].append( orthogroup_id )

            # C03 Lissachatina
            clade_1 = ['Achatina']  # C03
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta',
                       'Pomacea', 'Elysia', 'Aplysia', 'Homo']  # All but Lissachatina
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C03 Panpulmonata-Stylommatomorpha-Lissachatina
                clades_orthogroups_total[ 'C03' ].append( orthogroup_id )

            # C04 Pomacea
            clade_1 = ['Pomacea']  # C04
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Pomacea
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C04 Caenogastropoda-Pomacea
                clades_orthogroups_total[ 'C04' ].append( orthogroup_id )

            # C05 Gigantopelta
            clade_1 = ['Gigantopelta']  # C05
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Gigantopelta
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C05 Vetigastropoda-Gigantopelta
                clades_orthogroups_total[ 'C05' ].append( orthogroup_id )

            # C06 Lottia
            clade_1 = ['Lottia']  # C06
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Lottia
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C06 Patellogastropoda-Lottia
                clades_orthogroups_total[ 'C06' ].append( orthogroup_id )

            # C07 Mizuhopecten
            clade_1 = ['Mizuhopecten']  # C07
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Mizuhopecten
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C07 Pectinidae-Mizuhopecten
                clades_orthogroups_total[ 'C07' ].append( orthogroup_id )

            # C08 Crassostrea
            clade_1 = ['Crassostrea']  # C08
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Cyclina', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Crassostrea
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C08 Ostreidae-Crassostrea
                clades_orthogroups_total[ 'C08' ].append( orthogroup_id )

            # C09 Cyclina
            clade_1 = ['Cyclina']  # C09
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Octopus', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Cyclina
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C09 Heteroconchia-Cyclina
                clades_orthogroups_total[ 'C09' ].append( orthogroup_id )

            # C10 Octopus
            clade_1 = ['Octopus']  # C10
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Argonauta', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Octopus
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C10 Octopodoidea-Octopus
                clades_orthogroups_total[ 'C10' ].append( orthogroup_id )

            # C11 Argonauta
            clade_1 = ['Argonauta']  # C11
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Argonauta
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C11 Argonautoidea-Argonauta
                clades_orthogroups_total[ 'C11' ].append( orthogroup_id )

            # C12 Nautilus
            clade_1 = ['Nautilus']  # C12
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Acanthopleura', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Nautilus
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C12 Nautiloidea-Nautilus
                clades_orthogroups_total[ 'C12' ].append( orthogroup_id )

            # C13 Acanthopleura
            clade_1 = ['Acanthopleura']  # C13
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Lingula', 'Nautilus', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Acanthopleura
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C13 Polyplacophora-Acanthopleura
                clades_orthogroups_total[ 'C13' ].append( orthogroup_id )

            # C14 Lingula
            clade_1 = ['Lingula']  # C14
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Phoronis', 'Acanthopleura', 'Nautilus', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Lingula
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C14 Brachiopoda-Lingula
                clades_orthogroups_total[ 'C14' ].append( orthogroup_id )

            # C15 Phoronis
            clade_1 = ['Phoronis']  # C15
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Capitella', 'Lingula', 'Acanthopleura', 'Nautilus', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Phoronis
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C15 Phoronida-Phoronis
                clades_orthogroups_total[ 'C15' ].append( orthogroup_id )

            # C16 Capitella
            clade_1 = ['Capitella']  # C16
            clade_2 = ['Drosophila', 'Caenorhabditis', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Phoronis
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C16 Annelida-Capitella
                clades_orthogroups_total[ 'C16' ].append( orthogroup_id )

            # C17 Drosophila
            clade_1 = ['Drosophila']  # C17
            clade_2 = ['Capitella', 'Caenorhabditis', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Drosophila
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C17 Arthropoda-Drosophila
                clades_orthogroups_total[ 'C17' ].append( orthogroup_id )

            # C18 Caenorhabditis
            clade_1 = ['Caenorhabditis']  # C18
            clade_2 = ['Capitella', 'Drosophila', 'Phoronis', 'Lingula', 'Acanthopleura', 'Nautilus', 'Argonauta',
                       'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea',
                       'Achatina', 'Elysia', 'Aplysia', 'Homo']  # All but Caenorhabditis
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len(next_seqid.split(next_clade_2_species)) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == False:   # process single-species lineage C18 Nematoda-Caenorhabditis
                clades_orthogroups_total[ 'C18' ].append( orthogroup_id )

########## Continue processing OrthoFinder clustering of species tree nodes N0-N17

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
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[clade_id]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Ecdysozoa N2 orthogroups C35
        if node_id == 'N2':
            clade_1 = [ 'Drosophila' ]  # C17
            clade_2 = [ 'Caenorhabditis' ]  # C18
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C17' ].append( orthogroup_id )
                clades_orthogroups_total[ 'C18' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C17' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[clade_id]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C18' ].append( orthogroup_id )

        # ancestral Spiralia N3 orthogroups C34
        if node_id == 'N3':
            clade_1 = [ 'Capitella', 'Phoronis', 'Lingula', ]  # C33
            clade_2 = [ 'Acanthopleura', 'Nautilus', 'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea', 'Achatina', 'Elysia', 'Aplysia' ]  # C31
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Spiralian_Clade_4 N4 orthogroups C33
        if node_id == 'N4':
            clade_1 = [ 'Capitella' ]  # C16
            clade_2 = [ 'Phoronis', 'Lingula' ]  # C32
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C16' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C16' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Tetraneuralia-Mollusca N5 orthogroups C31
        if node_id == 'N5':
            clade_1 = [ 'Acanthopleura' ]  # C13
            clade_2 = [ 'Nautilus', 'Argonauta', 'Octopus', 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea', 'Achatina', 'Elysia', 'Aplysia'  ]  # C30
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C13' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C13' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Spiralian_Clade_2-Lophophorata N6 orthogroups C32
        if node_id == 'N6':
            clade_1 = [ 'Lingula' ]  # C14
            clade_2 = [ 'Phoronis'  ]  # C15
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C14' ].append( orthogroup_id )
                clades_orthogroups_total[ 'C15' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C14' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C15' ].append( orthogroup_id )

        # ancestral Conchifera-Molluscan_Clade_1 N7 orthogroups C30
        if node_id == 'N7':
            clade_1 = [ 'Nautilus', 'Argonauta', 'Octopus' ]  # C29
            clade_2 = [ 'Cyclina', 'Crassostrea', 'Mizuhopecten', 'Lottia', 'Gigantopelta', 'Pomacea', 'Achatina', 'Elysia', 'Aplysia' ]  # C27
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Cephalopoda N8 orthogroups C29
        if node_id == 'N8':
            clade_1 = [ 'Nautilus' ]  # C12
            clade_2 = [ 'Argonauta', 'Octopus' ]  # C28
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C12' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C12' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Molluscan_Clade_2 N9 orthogroups C27
        if node_id == 'N9':
            clade_1 = [ 'Cyclina', 'Crassostrea', 'Mizuhopecten' ]  # C26
            clade_2 = [ 'Lottia', 'Gigantopelta', 'Pomacea', 'Achatina', 'Elysia', 'Aplysia'  ]  # C24
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Coleoidea/Octopoda N10 orthogroups C28
        if node_id == 'N10':
            clade_1 = [ 'Argonauta' ]  # C11
            clade_2 = [ 'Octopus'  ]  # C10
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C11' ].append( orthogroup_id )
                clades_orthogroups_total[ 'C10' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C11' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C10' ].append( orthogroup_id )

        # ancestral Bivalvia N11 orthogroups C26
        if node_id == 'N11':
            clade_1 = [ 'Cyclina' ]  # C09
            clade_2 = [  'Crassostrea', 'Mizuhopecten'  ]  # C25
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C09' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C09' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Gastropoda N12 orthogroups C24
        if node_id == 'N12':
            clade_1 = [  'Lottia', 'Gigantopelta' ]  # C23
            clade_2 = [  'Pomacea', 'Achatina', 'Elysia', 'Aplysia'  ]  # C22
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Pteriomorphia N13 orthogroups C25
        if node_id == 'N13':
            clade_1 = [ 'Crassostrea'  ]  # C08
            clade_2 = [  'Mizuhopecten'  ]  # C07
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C08' ].append( orthogroup_id )
                clades_orthogroups_total[ 'C07' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C08' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C07' ].append( orthogroup_id )

        # ancestral Psilogastropoda N14 orthogroups C23
        if node_id == 'N14':
            clade_1 = [ 'Lottia' ]  # C06
            clade_2 = [  'Gigantopelta'   ]  # C05
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C06' ].append( orthogroup_id )
                clades_orthogroups_total[ 'C05' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C06' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C05' ].append( orthogroup_id )

        # ancestral Angiogastropoda N15 orthogroups C22
        if node_id == 'N15':
            clade_1 = [ 'Pomacea' ]  # C04
            clade_2 = [ 'Achatina', 'Elysia', 'Aplysia' ]  # C21
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C04' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C04' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Heterobranchia-Euthyneura N16 orthogroups C21
        if node_id == 'N16':
            clade_1 = [ 'Aplysia']  # C01
            clade_2 = [ 'Achatina', 'Elysia',  ]  # C20
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C01' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C01' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )

        # ancestral Panpulmonata N17 orthogroups C20
        if node_id == 'N17':
            clade_1 = [ 'Achatina' ]  # C03
            clade_2 = [ 'Elysia' ]  # C02
            test_clade_1 = False
            test_clade_2 = False
            for next_seqid in info_seqids:
                for next_clade_1_species in clade_1:
                    if len( next_seqid.split( next_clade_1_species ) ) > 1:
                        test_clade_1 = True
                for next_clade_2_species in clade_2:
                    if len( next_seqid.split( next_clade_2_species ) ) > 1:
                        test_clade_2 = True
            if test_clade_1 == True and test_clade_2 == True:   # test if orthogroup conserved: orthogroup in both clade 1 and 2
                if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                    clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                clades_orthogroups_total[ 'C03' ].append( orthogroup_id )
                clades_orthogroups_total[ 'C02' ].append( orthogroup_id )
            if test_clade_1 == True and test_clade_2 == False:   # test if orthogroup conserved: orthogroup present in clades 1 not 2 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C03' ].append( orthogroup_id )
            if test_clade_1 == False and test_clade_2 == True:   # test if orthogroup conserved: orthogroup present in clades 2 not 1 and ancestor A
                if orthogroup_id in clades_orthogroups_total[ ancestor_id ]:
                    if orthogroup_id not in clades_orthogroups_total[ clade_id ]:
                        clades_orthogroups_total[ clade_id ].append( orthogroup_id )
                    clades_orthogroups_total[ 'C02' ].append( orthogroup_id )


########## Read counts to output
print( 'Calculating gains and losses...')

output = 'Clade ID' + '\t' + 'Clade Name' + '\t' + 'Orthogroups Total' + '\t' + 'Orthogroups Conserved' + '\t' + 'Orthogroups Gained' + '\t' + 'Orthogroups Lost' + '\n'
output_counts.write( output )

for next_clade in sorted( clades_orthogroups_total.keys() ):

    clade_name = clades_names[ next_clade ]

    total_orthogroups = str( len( list( set( clades_orthogroups_total[ next_clade ] ) ) ) )
    clade_orthogroups = list( set( clades_orthogroups_total[ next_clade ] ) )

    if next_clade != 'C37':
        ancestor_id = clades_ancestors[ next_clade ]
        ancestor_orthogroups = list( set(  clades_orthogroups_total[ ancestor_id ] ) )

        conserved_orthogroups = 0
        gained_orthogroups = 0
        for next_orthogroup in clade_orthogroups:
            if next_orthogroup in ancestor_orthogroups:
                conserved_orthogroups = conserved_orthogroups + 1
            else:
                gained_orthogroups = gained_orthogroups + 1

        lost_orthogroups = 0
        for next_orthogroup in ancestor_orthogroups:
            if next_orthogroup not in clade_orthogroups:
                lost_orthogroups = lost_orthogroups + 1

    else:
        conserved_orthogroups = 'NA'
        gained_orthogroups = 'NA'
        lost_orthogroups = 'NA'

    output = next_clade + '\t' + clade_name  + '\t' + total_orthogroups + '\t' + str( conserved_orthogroups ) + '\t' + str( gained_orthogroups ) + '\t' + str( lost_orthogroups ) + '\n'
    output_counts.write( output )

input_map_nodes_to_clades.close()
input_map_clades_to_ancestors.close()
input_map_clades_to_names.close()
input_nodes.close()
output_counts.close()
