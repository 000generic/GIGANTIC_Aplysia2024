#! python
# DeepLoc script

# deeploc2 --fasta --output output --model Fast --plot

input_fastas = open( 'output/1-list-fastas', 'r' )
output_commands = open('003-deeploc2-species19', 'w')

for next_fasta in input_fastas:
        fasta = next_fasta[ :-1 ]
        output = 'deeploc2 --fasta ' + fasta + ' --output output --model Fast\n'
        output_commands.write( output )

input_fastas.close()
output_commands.close()
