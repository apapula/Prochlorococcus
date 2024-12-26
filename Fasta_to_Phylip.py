from Bio import AlignIO
from pathlib import Path
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter
import re

if 1:
    groupname='C1'#HLII'#BB_NoClose'
    genenum=1627

    file = 'Gene%i_%s.fasta' %(genenum,groupname)
    #file = '../input_data/Alignments/Gene%i_%s.fasta' %(genenum,groupname)
    pathway = '../input_data/Alignments/'
if 0:
    file='Kashtan_Clique1_Genomes.fasta'
    pathway='../../New_Kashtan/input_data/'
alignment = AlignIO.read(Path(pathway, file), 'fasta')



# write out in relaxed phylip
fileoutphy = re.sub('fasta', 'phy', file)
with open(Path(pathway, fileoutphy), 'w') as output_handle:        
    SequentialPhylipWriter(output_handle).write_alignment(alignment, id_width=30)
