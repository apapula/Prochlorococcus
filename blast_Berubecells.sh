testedit
#!/bin/bash                                                                                                                                 
#SBATCH --job-name=Berube_cells_BLAST_Against_4refs                                                                                                          
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                  
#SBATCH --mail-user=apaps@stanford.edu # Where to send mail                                                                                  
#SBATCH --time=23:59:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)                                                       
#SBATCH --ntasks=1                                                                                                                           
#SBATCH --cpus-per-task=1                                                                                                                    
#SBATCH --mem-per-cpu=2G                                                       

#SBATCH --array=1-789 # Array range                                                                                                        
ml perl
ml biology
ml ncbi-blast+


# Data dir paths                                                                                                                              
results_dir="/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/Sept22_AllCells_New/Blast_Output/"

data_path="/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/Genomes"


readarray query_list <  NT_Genes_Filelist.txt
myquery=$(echo "${query_list[$(($SLURM_ARRAY_TASK_ID - 1))]}" | sed 's/[[:space:]]//g') # Get sample ID and remove spaces at end if necessary

echo $myquery

subject_file="/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/Genomes/MIT9301_nt.fna"

thisquery=$data_path/$myquery


queryshortname=$(echo ${myquery:0:10})

echo $queryshortname

makeblastdb -dbtype nucl -parse_seqids -in ${thisquery}

outfile=${results_dir}/${queryshortname}_MIT9301nt.out

blastn -db ${thisquery}  -query ${subject_file} -qcov_hsp_perc 55 -word_size 8 -out ${outfile}  -outfmt "6 qseqid sseqid pident qcovhsp length" -task blastn

sed -i '1s/^/qseqid\tsseqid\tpident\tqcovhsp\tlength\n/' $outfile

