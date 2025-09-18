#!/bin/bash
#SBATCH --job-name=Orthologous_gene_grouper                                                                                                            
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                  
#SBATCH --mail-user=apaps@stanford.edu # Where to send mail                                                                                  
#SBATCH --time=23:59:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)                                                       
#SBATCH --ntasks=1                                                                                                                           
#SBATCH --cpus-per-task=1                                                                                                                    
#SBATCH --mem-per-cpu=1G                                                       

#SBATCH --array=1-999# Array range                                                                                                     
ml perl
ml biology
ml ncbi-blast+

ml python/3.6.1

# Data dir paths

yourpath='***'

results_dir=${yourpath}"/Unaligned_Gene_Groups"
aligned_results_dir=${yourpath}"/Aligned_Gene_Groups"
data_path=${yourpath}"/New_New_Filters"



readarray query_list <  MIT9301_nt_genelist.txt
mygenes=$(echo "${query_list[$(($SLURM_ARRAY_TASK_ID - 1))]}" | sed 's/[[:space:]]//g') # Get sample ID and remove spaces at end if necessary

echo $mygenes


geneindex=$SLURM_ARRAY_TASK_ID

genenamesfile=Gene${geneindex}_Names_kash.txt
geneseqfile=Gene${geneindex}_Kashtan_and_HLII.fasta


grep "$mygenes" ${data_path}/* | cut -f2| uniq > $genenamesfile



blastdbcmd -db /scratch/groups/dsfisher/Prochlorococcus/Databases/Berube_Kashtan_HLII_Database.fna -dbtype nucl -entry_batch $genenamesfile -outfmt '%f' -out $results_dir/$geneseqfile


rm $genenamesfile



if [[ -s $results_dir/$geneseqfile ]]; then
    
      
    python3 codon_aware_align_alana.py -i $results_dir/${geneseqfile} -m $GROUP_HOME/bin/muscle -o ${aligned_results_dir}/${geneseqfile}
    
fi;

