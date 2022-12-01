#! python
# GIGANTIC species tree script

import os
import re
from os import path

class Origin:
    def __init__(self, filename):
        self.id = filename.split('-')[-1]
        self.info = {}
        seq = []
        with open('./output/4-clipkit/' + filename,'r') as f:
            file=f.read()
            for species in file.split('>'):
                if species.strip():
                    i=species.strip().split('\n')
                    header=i[0]
                    j=header.strip().split('-')
                    name='-'.join(j[5:7])
                    seq=species.strip().split('\n')[1:]
                    self.info[name]=seq
          
        self.fill=[]
        block=[]
        for i in seq:
            self.fill.append('X' * len(i))
        for i in self.fill:
            block.append(len(i))
            self.len=sum(block)

    def write_data(self, all_species):
        for i in all_species: 
            if i in self.info:
                
                # Create files for appending 
                with open('./output/5-catsequences/'+i, 'a') as f:
                    for output in self.info[i]:
                        f.write(output)
                       # f.write('\n')

            # Fill space with 'X's if a gene is absent in certain species.
            else:
                with open('./output/5-catsequences/'+ i,'a') as f:
                    for output in self.fill:
                        f.write(output)
                        #f.write('\n')


# -redo 
for i in os.listdir('./output/5-catsequences/'):
  os.remove('./output/5-catsequences/'+i)

#Pass files to Class Origin
origin_file=[]
with open('output/1-list/busco-list','r') as f:
    id=f.read()
    for i in id.split('\n'):
        if os.path.exists('output/4-clipkit/clipkit-'+ i):
            origin_file.append(Origin('clipkit-'+i))

# Create a list of species names  
all_species=[]
for i in origin_file:
    with open('output/1-list/species-list','r') as f:
        all_species=f.read().strip().split('\n')
   # Remove empty strings from the list 
      #  all_species = '-'.join(all_species).split()       
all_species =sorted (list(all_species))

for i in all_species:
    with open('./output/5-catsequences/'+i, 'a') as f:
        f.write('>' + i +'\n')
for i in origin_file:
    i.write_data(all_species)


