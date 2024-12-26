#!/bin/bash                                                                                                                                 
#SBATCH --job-name=Orthologous_gene_grouper                                                                                                            
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                  
#SBATCH --mail-user=apaps@stanford.edu # Where to send mail                                                                                  
#SBATCH --time=23:59:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)                                                       
#SBATCH --ntasks=1                                                                                                                           
#SBATCH --cpus-per-task=1                                                                                                                    
#SBATCH --mem-per-cpu=2G                                                       

#SBATCH --array=1-999# Array range                                                                                                     
ml perl
ml biology
ml ncbi-blast+

ml python/3.6.1

# Data dir paths                                                                                                                              
results_dir="/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/Flexible_Orthogroups_Unaligned"
aligned_results_dir="/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/Flexible_Orthogroups_Aligned"

filename_path="/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/Flexible_Gene_Name_Files"



readarray query_list <  ${filename_path}/Flexible_Gene_File_Name_List.txt 
mygenes=$(echo "${query_list[$(($SLURM_ARRAY_TASK_ID - 1))]}" | sed 's/[[:space:]]//g') # Get sample ID and remove spaces at end if necessary


echo $mygenes


geneindex=$(($SLURM_ARRAY_TASK_ID+2000))

#define genenamesfile; let that be the list iterated through above

geneseqfile=FlexibleGene${geneindex}.fasta

#here already have gene names so don't need to do?

#grep "$mygenes" ${data_path}/* | cut -f2| uniq > $genenamesfile

genenamesfile=${filename_path}/$mygenes

head "/scratch/groups/dsfisher/Prochlorococcus/Data_Berube/All_Berube_Genes.fasta"

blastdbcmd -db "/scratch/groups/dsfisher/Prochlorococcus/Data_Berube/All_Berube_Genes.fasta" -dbtype nucl -entry_batch $genenamesfile -outfmt '%f' -out $results_dir/$geneseqfile


#rm $genenamesfile


if [[ -s $results_dir/$geneseqfile ]]; then
    
   
    
    python3 codon_aware_align_alana.py -i $results_dir/${geneseqfile} -m $GROUP_HOME/bin/muscle -o ${aligned_results_dir}/${geneseqfile}
fi;

