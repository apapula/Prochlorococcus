'''
for each c1 cell: rehead the fasta to have >orfnum_cellname
blsat again composite reference core to get non C1 core ORFs for that cell
then take all cells non core and concatenate
blast all against all?
make orthologous gene groups, using edges and connected components

write the groups to a file.
then, for close cells, for ORFs that don't hit each other (non core), color by whether they are in a group with at least 1 other cell
'''

import C1_AllNT_Close_Pair_BLASTing as blastrepo
import numpy as np
import matplotlib.pyplot as plt
import statistics
import ORFs_Per_SAG as orfrepo
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import csv
import math
from os.path import exists
import subprocess
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Cell_Lists import C1_without_ribotype2,clique1_without_ribotype,hliilist2,BB_list,AllKashBATS_without_ribotype, BB_closegrp_1,C3,hliilist2,Atlantic,N_Pacific,S_Pacific,cl2_without_ribotype,cl4_without_ribotype,cl5_without_ribotype,cl6_without_ribotype, C2_without,C3_without
import Load_Sequences_Repo as loadrepo
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import poisson
import ast
from scipy.stats import linregress, pearsonr
import networkx as nx
import argparse
import Make_Synonymous_Distance_Matrices_Berube_Kashtan_Using_All_Sites_Per_Pair_and_HalfGene_Record_sites_SNPs as distmatrepo
import Repo_site_catalogues_single_pair as repo


Pacific=N_Pacific+S_Pacific

def plot_singletons_doubletons_as_function_of_qcov_pid(qcutlist,pidcutlist,groupname,celllist,C1bool,returnnum):
    plt.close()
    mylists=[]
    labels=[]
    for q in qcutlist:
        for p in pidcutlist:
            myl=get_connected_components_from_all_v_all_BLAST(q,p,groupname,celllist,returnnum=returnnum,C1bool=C1bool,returnbool=True)
            mylists.append(myl)
            labels.append('Qcut %s Pidcut %s' %(str(q),str(p)))

    xlist=list(range(1,returnnum+1))
    for i in range(len(labels)):
        plt.plot(xlist,mylists[i],'o-',label=labels[i])
    plt.legend()
    plt.xlabel('Number of Cells in Connected Component')
    plt.ylabel('Number of Connected Components')
    plt.title('Number of Singletons Through %i-tons for %s\nAs a function of Coverage and Identity Cuts' %(returnnum,groupname))
    if C1bool:
        plt.show()
    if not C1bool:
        plt.savefig('%s_Singletons_as_funct_of_cutoffs.png' %(groupname))
        plt.close()
            
    
#########################################################

def count_quasi_core(groupname,ncells,orthosizes,C1bool,qcut,pidcut):
    #get number of orthogroups with e.g. more than 30,40,50 , etc perc of cells
    mycuts=[.1,.2,.3,.4,.5,.6,.7,.8,.9]
    myngroups=[]
    for i in range(len(mycuts)):
        cutsize=mycuts[i]*ncells
        mycount=0
        for item in orthosizes:
            if item>=cutsize:
                mycount+=1
        myngroups.append(mycount)

    plt.plot(mycuts,myngroups,'-o')
    plt.xlabel('Fraction of Cells in the OrthoGroup')
    plt.ylabel('Number of Orthogroups')
    plt.title('%s (%i Cells) Quasi-Core Orthogroups: Number of Orthogroups vs Coverage\n(Identity %i and Coverage Cutoff %i)' %(groupname,ncells,pidcut,qcut))
    if C1bool:
        plt.show()
    if not C1bool:
        plt.savefig('%s_Quasi_Core_Orthgroup_Number_vs_Size.png' %groupname)
        plt.close()

def check_for_same_cell(qseq,sseq,C1bool):
    idx1=qseq.find('_')
    if C1bool:
        cell1=qseq[idx1+1:]
    if not C1bool:
        cell1=qseq[:idx1]
    idx2=sseq.find('_')
    if C1bool:
        cell2=sseq[idx2+1:]
    if not C1bool:
        cell2=sseq[:idx2]
    if cell1==cell2:
        return True
    if cell1=='B245a_518D8' or cell2=='B245a_518D8':
        return True
    return False

def get_connected_components_from_all_v_all_BLAST(qcut,pidcut,groupname,celllist,returnnum=0,C1bool=0,returnbool=False):
    print(pidcut,qcut)
    if C1bool:
        infileblast='../input_data/%s_nonCore_ORFs_All_Against_All_2.out' %groupname
    if not C1bool:
        infileblast='Berube_HLII_NonCore/%s_nonCore_ORFs_All_Against_All_notask.out' %groupname
    edges=[]#in any connected component, later make sure that no cell is listed more than once
    fieldnames=['qseqid','sseqid','pident','qcovhsp','length']

    bidir_edges=[]
    with open(infileblast,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t',fieldnames=fieldnames)
        for row in csvreader:
            
            query=row['qseqid']
            pid=float(row['pident'])
            qcov=float(row['qcovhsp'])
            if pid<pidcut or qcov<qcut:
                continue
            sseq=row['sseqid']
            if check_for_same_cell(query,sseq,C1bool) and sseq!=query:
                #print(sseq,query)
                #t=input('t')
                continue#new 1-3
            if '245a_518D8' in sseq or '245a_518D8' in query:
                continue
            edges.append((query,sseq))
            #newseqs=sorted([query,sseq])
            #edges.append((newseqs[0],newseqs[1],0.5))
            
            #if (sseq,query) in edges:
                #bidir_edges.append((query,sseq))
    edges=prune_edges_to_bidirectional(edges)
    #edges=bidir_edges
    
    orthogroups=get_graph_connected_comps(edges,C1bool,groupname,qcut,pidcut)

    orthogroups=filter_orthogroups(orthogroups,C1bool,groupname,celllist,qcut,pidcut)###edit 1-5
    
    with open('../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut)),'w') as myoutfile:
        for item in orthogroups:
            myoutfile.write('%s\n' %str(item))
            
    mylist=list(range(len(celllist)+2))
    mybins=[mylist[i]-0.5 for i in range(len(mylist))]
    orthosizes=[]
    nover=0
    nover_70=0#clusters longer than 70% celllist length
    for item in orthogroups:
        orthosizes.append(len(item))
        if len(item)>100:
            print('orthgroup size %i' %len(item))
        if len(item)>len(celllist):
            nover+=1
        if len(item)>0.7*len(celllist):
            nover_70+=1
    print('%i (%f) clusters larger than celllist length; %i (%f) clusters larger than 70 percent of celllist length' %(nover,float(nover)/len(orthogroups),nover_70,float(nover_70)/len(orthogroups)))


    plt.hist(orthosizes)
    plt.xlabel('Number of Cells in Connected Component')
    plt.ylabel('Number of Connected Components')
    plt.title('%s Size of non-C1-Core ORF Putative Orthologous Gene Groups\nCoverage cut %i, Identity Cut %i' %(groupname,qcut,pidcut))
    if C1bool:
        plt.show()
        count_quasi_core(groupname,len(celllist),orthosizes,C1bool,qcut,pidcut)
    if not C1bool:
        #plt.savefig('%s_NonCore_OrtoGroups_Sizes1.png' %groupname)
        plt.close()
    
    plt.hist(orthosizes,bins=mybins)
    plt.xlabel('Number of Cells in Connected Component')
    plt.ylabel('Number of Connected Components')
    plt.title('%s Size of non-C1-Core ORF Putative Orthologous Gene Groups\nCoverage cut %i, Identity Cut %i' %(groupname,qcut,pidcut))
    plt.yscale('log')
    #plt.show()
    if C1bool:
        plt.show()
    if not C1bool:
        plt.savefig('%s_NonCore_OrtoGroups_Sizes_bins_q%s_p%s.png' %(groupname,str(qcut),str(pidcut)))
        plt.close()

    if returnbool:
        #return up to some number, say 10
        myl=[]#number of singletons, doubletons,etc up to returnnum-tons
        for i in range(1,returnnum+1):
            myl.append(orthosizes.count(i))
        print('Qcut %s Pid %s Group %s: number of clusters with 1-%i cells:' %(str(qcut),str(pidcut),groupname,returnnum))
        print(myl)
        return myl


def filter_orthogroups(orthogroups,C1bool,groupname,celllist,qcut,pidcut):
    #take out self-hits; repeated cells per orthogorup. is this helpful? how does this change results? plot sizes before and after
    oldsizes=[len(orthogroups[i]) for i in range(len(orthogroups))]
    newgroups=[]
    for item in orthogroups:
        if len(item)==1:
            newgroups.append(item)
            continue
        newitem=filter_item(item,C1bool)
        newgroups.append(newitem)
    print('done')
    newsizes=[len(newgroups[i]) for i in range(len(newgroups))]

    mybins=[0.5]
    for i in range(len(celllist)+1):
        mybins.append(i+0.5)
    print(len(oldsizes),len(newsizes),mybins)
    plt.hist([oldsizes,newsizes],label=['Old Sizes','New Sizes'],bins=mybins)
    # plt.hist(np.array([oldsizes,newsizes],dtype=object),label=['Old Sizes','New Sizes'],bins=mybins)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Number of Cells in Orthogroup')
    plt.ylabel('Number of Groups')
    if C1bool:
        plt.show()
    if not C1bool:
        plt.savefig('Orthogroups_before_after_filtering_%s_q%s_p%s.png' %(groupname,str(qcut),str(pidcut)))
        plt.close()


    plt.hist([oldsizes,newsizes],label=['Old Sizes','New Sizes'],bins=50)#mybins)
    # plt.hist(np.array([oldsizes,newsizes],dtype=object),label=['Old Sizes','New Sizes'],bins=mybins)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Number of Cells in Orthogroup')
    plt.ylabel('Number of Groups')
    if C1bool:
        plt.show()
    if not C1bool:
        plt.savefig('fewbinsOrthogroups_before_after_filtering_%s_q%s_p%s.png' %(groupname,str(qcut),str(pidcut)))
        plt.close()
    
    return newgroups



def filter_item(l,C1bool):
    print(l)
    newl=[l[i] for i in range(len(l))]
    for i in range(len(l)):
        print(i,len(l))
        for j in range(i):
            print(j)
            samebool=check_for_same_cell(l[i],l[j],C1bool)
            if samebool:
                if l[j] in newl:
                    newl.remove(l[j])
               # break
    return newl
    
def prune_edges_to_bidirectional(edges):
    newedges=[]
    print('pruning')
    i=0

    G=nx.DiGraph()
    G.add_edges_from(edges)
    for item in edges:
        num=G.number_of_edges(item[0],item[1])+G.number_of_edges(item[1],item[0])
        print(num)
        #if num<2:
            #t=input('t')
        if num>1:
            newedges.append(item)
    '''
    #edges=[(1,2),(4,5),(5,4),(3,3),(2,1)]
    
    myset=set(edges)
    for item in edges:
        
        i+=1
        if i%1000==0:
            print(i,len(edges))
        (c1,c2)=item
        if c1==c2:
            myset.remove(item)
            continue
        #if (c2,c1) in myset:#edges:
        if len(list(set([(c2,c1)]).intersection(myset)))>0:
            newedges.append(item)
            edges.remove((c2,c1))
            myset.remove((c2,c1))
            myset.remove(item)
    '''
    print('done pruning')
    return newedges

def get_graph_connected_comps(edges,C1bool,groupname,qcut,pidcut):
    G=nx.Graph()
    G.add_edges_from(edges)


    #to_remove = [(a,b) for a,b, attrs in G.edges(data=True) if attrs["weight"] <1]
    #G.remove_edges_from(to_remove)
    
    #G.add_edges_from(edges)#,weight=0.5)
    orthogroups=[]
    for clq in nx.connected_components(G):
        miniclqlist=[]
        for item in clq:
            miniclqlist.append(item)
        filter_group_for_cell_duplicates(miniclqlist)
        orthogroups.append(miniclqlist)
    groupsizes=[len(orthogroups[i]) for i in range(len(orthogroups))]
    if 1:
        mylist=list(range(max(groupsizes)+1))
        mybins=[-0.5]
        for i in range(len(mylist)):
            mybins.append(i+0.5)
        plt.hist(groupsizes,bins=mybins)
        plt.xlabel('Number of Cells in Connected Component')
        plt.ylabel('Number of Connected Components')
        plt.title('%s Size of non-C1-Core ORF Putative Orthologous Gene Groups\nCoverage cut %i, Identity Cut %i' %(groupname,qcut,pidcut))
        plt.yscale('log')
        if C1bool:
            plt.show()
        if not C1bool:
            plt.savefig('%s_NonCore_OrtoGroups_Sizes_log.png' %groupname)
            plt.close()
    return orthogroups

def filter_group_for_cell_duplicates(group):
    cells=[]
    for item in group:
        idx=item.index('_')
        cell=item[idx+1:]
        #if cell in cells:
            #print('Cell %s already in this orthologous gene group (%i ORFs)' %(cell,len(group)))
        if cell not in cells:
            cells.append(cell)

######################################################################

def singleton_ORFs_per_cell(celllist,groupname):#go through the output file and get number singleton ORFs per cell
    single_dict={}#cell:count
    for item in celllist:
        single_dict[item]=0
    with open('../input_data/%s_Orthologous_Gene_Group_Qcut80_pid80.txt' %groupname,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    for line in lines:
        mylist=ast.literal_eval(line)
        if len(mylist)>1:
            continue
        s=mylist[0]
        idx=s.find('_')
        cell=s[idx+1:]
        if cell not in celllist:
            if groupname=='C1':
                print('%s not in celllist, %i' %(cell,C1_without_ribotype2.index(cell)))
            if groupname!='C1':
                print('%s not in celllist' %cell)
            t=input('t')
        single_dict[cell]+=1
    single_nums=[single_dict[cell] for cell in single_dict]

    newsingle_dict={}
    for item in single_dict:
        newsingle_dict[C1_without_ribotype2.index(item)]=single_dict[item]
    print(newsingle_dict)#nums)
    mylist=list(range(max(single_nums)+2))
    mybins=[mylist[i]-0.5 for i in range(len(mylist))]
    plt.hist(single_nums,bins=mybins)
    plt.xlabel('Number of Singleton non-C1-Core ORFs per Cell')
    plt.ylabel('Number of C1 Cells')
    plt.title('Singleton non-C1-Core ORFs per C1 Cell')
    plt.show()


    plt.hist(single_nums)#,bins=mybins)
    plt.xlabel('Number of Singleton non-C1-Core ORFs per Cell')
    plt.ylabel('Number of C1 Cells')
    plt.title('Singleton non-C1-Core ORFs per C1 Cell')
    plt.show()
            
####################################################################################################
def rehead_fasta(celllist):
    name_num_dict=blastrepo.get_kashtan_number_name()
    for cell in celllist:
        if cell not in name_num_dict:
            continue
        print(cell)
        num=name_num_dict[cell]
        qcut=25
        hitlist=blastrepo.get_list_of_cells_ORFs_with_hit_in_C1(cell,num,qcut,excessbool=excessbool)
        #now get a dictionary of ORFs that are not in the hitlist
        cellfile='/scratch/groups/dsfisher/Prochlorococcus/Data_Kashtan/Gene_Fastas/%s.genes.fna' %str(num)
        newcellfile='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/All_NT_Close_Pair_Blasting/nonC1Core_Fastas/%s_nonC1Core_AllNT.fna' %cell
        seqdict,totlyes,totlnon=blastrepo.load_ORF_sequences_allNT(cellfile,hitlist)
       
        nonC1seqdict={}
        for item in seqdict:
            if int(item) in hitlist:
                continue
            nonC1seqdict[item]=seqdict[item]
        write_to_fasta_reheaded(newcellfile,seqdict,cell)
        
def write_to_fasta_reheaded(fastafile,seqdict,cellname):
    print("WRITING")
    with open(fastafile,'w') as outfile:
        for item in seqdict:
            outfile.write('>%s_%s\n' %(item,cellname))
            outfile.write('%s\n' %seqdict[item])

######################################
#now do for berube
'''
/scratch/groups/dsfisher/Prochlorococcus/Data_Berube/All_Berube_Genes.fasta
get HLII genes of this
then, go to /scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/Kashtan/Filtered_with_HLII_SAGs/[berube sag_MIT9301nt.out]

go through for each hlii cell and make list.  for each line in csv file, ind=qseqid.rfind('_').  if int(qseqid[ind+1:]) in hlii core genes, add to list
exclude these orfs

write a new all hlii non=hlii-core fasta file.
blast all against all

or could just go one hlii cell at a time. load its genes.fna file. remove the genes in the filtered file hitting a core 9301 gene. write. onto next cell

'''
def make_berube_cluster_size_plots(celllist,groupname,qcut,pidcut):
    #qcut,pidcut=80,80
    infile='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))#2
    #infile='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/input_data/%s_Orthologous_Gene_Group_Qcut80_pid80.txt' %groupname
    sizes=[]
    with open(infile,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    for line in lines:
        sizes.append(len(ast.literal_eval(line)))

    mylist=list(range(len(celllist)+2))
    mybins=[mylist[i]-0.5 for i in range(len(mylist))]
    

    print('%i singletons, %i doubletons, %i tripletons' %(sizes.count(1),sizes.count(2),sizes.count(3)))

    corecovs=get_covered_core()
   

    if 1:
        plt.hist([corecovs,sizes],bins=mybins,label=['HLII Core','Flexible +'],stacked=True)
        plt.legend(fontsize=15)
        plt.xlabel('Number of SAGs with flexible gene', fontsize=20)
        plt.ylabel('Number of flexible gene groups', fontsize=20)
        ax=plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.tight_layout()
        plt.show()

        plt.hist([sizes,corecovs],bins=mybins,label=['Flexible +','HLII Core'],stacked=True)
        #plt.hist([corecovs,sizes],bins=mybins,label=['HLII Core','Flexible +'],stacked=True)
        plt.legend(fontsize=15)
        plt.yscale('log')
        plt.xlabel('Number of SAGs with flexible gene', fontsize=20)
        plt.ylabel('Number of flexible gene groups', fontsize=20)
        ax=plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.tight_layout()
        plt.show()
    
    plt.hist(sizes,bins=mybins)
    for i in range(5):
        print('%i cells: %i' %(i, sizes.count(i)))
    print('4 or more cells: %i' %(sum([sizes.count(j) for j in range(4,35)])))

   
    if 1:
        plt.xlabel('Number of SAGs with flexible gene', fontsize=20)
        plt.ylabel('Number of flexible gene groups', fontsize=20)
        ax=plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.tight_layout()
    if 0:
        plt.xlabel('Number of Cells in Connected Component')
        plt.ylabel('Number of Connected Components')
        plt.title('Size of non-C1-Core ORF Putative Orthologous Gene Groups\nCoverage cut %i, Identity Cut %i' %(qcut,pidcut))
        #plt.show()
    plt.yscale('log')
        #plt.savefig('%s_NonCore_OrtoGroups_Sizes_bins_logy_q%s_p%s.png' %(groupname,str(qcut),str(pidcut)))
        #plt.close()
    if 0:#add 1/f plot
        midpts=mylist
        #get norm
        if groupname=='HLII':
            myrange=20
        norm=float(sum(sizes[:myrange]))/(float(sum([1.0/j for j in range(1,1+myrange)])))
        finverse=[float(norm)/midpts[i] for i in range(1,len(midpts))]
        finverse2=[0.2*float(norm)/midpts[i] for i in range(1,len(midpts))]
        plt.plot(midpts[1:],finverse,'o',label='1/f')
        plt.plot(midpts[1:],finverse2,'o',label='1/f')
        plt.legend()
    plt.show()
    #count_quasi_core(groupname,len(celllist),sizes,False,qcut,pidcut)
###########################
def make_Berube_allgenes_fasta(celllist,groupname):

    outfile='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/Berube_HLII_NonCore/All_NonHLIICore_ORFs_%s_2.fasta' %groupname
    
    filtered_prefix='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/Kashtan/Filtered_with_HLII_SAGs/'
    fasta_prefix='/scratch/groups/dsfisher/Prochlorococcus/Data_Berube/Fasta_files/'
    cellfile_name_dict=get_Berube_cell_file_names(celllist)

    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    
    with open(outfile,'w') as myoutfile:
        for cell in celllist:
            do_single_cell(cell,cellfile_name_dict,filtered_prefix,fasta_prefix,mycores,myoutfile)

            

def do_single_cell(cell,cellfile_name_dict,filtered_prefix,fasta_prefix,mycores,myoutfile):
    if cell not in cellfile_name_dict:
        print('Cell %s not in cellfile_name_dict' %cell)
        t=input('t')
        return
    filtered_file='%s%s_MIT9301nt.out' %(filtered_prefix,cell)
    fasta_file='%s%s' %(fasta_prefix,cellfile_name_dict[cell])

    #now, get the ORF names that hit a HLII core gene in the filtered file
    core_orfs=get_ORFs_with_core_hit(mycores,filtered_file)
    print('%s: %i core ORFs' %(cell,len(core_orfs)))
    seqdict=get_seqs_of_noncore_orfs(core_orfs,cell,fasta_file)
    if len(list(seqdict.keys()))==0:
        return
    write_seqdict(myoutfile,seqdict)

def write_seqdict(myoutfile,seqdict):
    for item in seqdict:
        myoutfile.write('>%s\n' %item)
        myoutfile.write('%s\n' %(seqdict[item]))
                        
def get_seqs_of_noncore_orfs(core_orfs,cell,fasta_file):
    #go along in fasta file. if noncore, rewrite name. store as dict
    seqdict={}
    if not exists(fasta_file):
        print('File %s does not exist' %(fasta_file))
        t=input('t')
        return {}
    with open(fasta_file,'r') as myinfile:
        lines=myinfile.read()
    cellseqs=re.split('>',lines)
    print('%s: %i total ORFs' %(cell,len(cellseqs)))
    for i in range(len(cellseqs)):
        sublist=re.split('\n',cellseqs[i])
        if len(sublist)<2:
            continue
        title=sublist[0]
        ind=title.find(' ')
        myorf=title[:ind]
        #print(core_orfs)
        #print(title)
        #print(myorf,myorf in core_orfs)
        #t=input('t')
        if myorf in core_orfs:
            continue
        seq=''.join(sublist[1:])
        seq=seq.replace("N","-")
        seqdict['%s_%s' %(cell,str(myorf))]=seq
    print('%s: returning %i non core ORFs' %(cell,len(list(seqdict.keys()))))
    return seqdict

def get_ORFs_with_core_hit(mycores,filtered_file):
    core_orfs=[]
    with open(filtered_file,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            q=row['qseqid']
            s=row['sseqid']
            ind=q.rfind('_')
            qnum=int(q[ind+1:])
            #print(qnum,qnum in mycores)
            #t=input('t')
            #if t=='y':
                #print(mycores)
            #find: not a problem with findnig core orfs
            if qnum in mycores:
                core_orfs.append(s)
    return core_orfs
        
def get_Berube_cell_file_names(celllist):
    with open('/scratch/groups/dsfisher/Prochlorococcus/Data_Berube/Cell_File_Names_List.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    namedict={}#cell:filename
    for cell in celllist:
        for line in lines:
            if cell in line:
                namedict[cell]=line
                break
    return namedict
#########################################


def examine_doubletons_etc(celllist,groupname,C1bool,median_mean='Mean'):
    qcut,pidcut=60,80
  
    infile='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))#2
    print(infile)

    if C1bool:
        ldict=get_orf_len_dict(groupname,celllist)
    oldgroupname=groupname
    groupname+=' %i cov %i pid' %(qcut,pidcut)
    '''
    want to know:
    doubletons+ from a clique?
    for each cell pair, how often doubleton? vs doubleton wg dist?
    for each triple, how often tripleton? verse tripleton wg diam?


    FOR EACH CELL PAIR, N DOUBLETOND, N IDENTICAL DOUBLETONS
    '''
   
     
    overalllist,wgmatrix,totdivs,tottriple_diams=load_WG_divs(celllist,median_mean)
    doubleton_dict,tripleton_dict={},{}
    doubleton_wg,tripleton_wg=[],[]
    doubleton_divs,tripleton_diams=[],[]
    n_doubletons=0
    n_id,n_id_clique=0,0
    len_id,len_notid=[],[]
    n_notid_clique=0
    n_sing,n_trip=0,0
    tot_lengths={}#number of cells in the cluster: [list of lengths of the orfs]
    
    with open(infile,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    print('%i orthogroups' %len(lines))

    n_flex_genes=sum([len(ast.literal_eval(lines[i])) for i in range(len(lines))])
    n_singleton=sum([1 for i in range(len(lines)) if len(ast.literal_eval(lines[i]))==1])
    print('%i total flex genes. %i singletons. %f genes are singletons, %f groups are singletons' %(n_flex_genes,n_singleton,float(n_singleton)/n_flex_genes,float(n_singleton)/len(lines)))
    
    for line in lines:
        l=ast.literal_eval(line)
        #if len l>3 continue. then get the cells. want to know their WG dist. load that in
        my_len_l=len(l)
        my_orfs_lens_this=[]
        if len(l)==1:
            n_sing+=1
        if len(l)==3:
            n_trip+=1
        if C1bool:
            for i in range(len(l)):
                my_orfs_lens_this.append(ldict[l[i]])
            

            this_orf_length=repo.avelist(my_orfs_lens_this)###############
            if my_len_l in tot_lengths:
                tot_lengths[my_len_l].append(this_orf_length)
            
            if my_len_l not in tot_lengths:
                tot_lengths[my_len_l]=[this_orf_length]
        
        if len(l)>3:
            continue
        if len(l)==2:
           
            #print(l)
            #t=input('t')
            s1=l[0]
            s2=l[1]
            #print(s1,s2)
            
            idx1=s1.find('_')
            if C1bool:
                c1=s1[idx1+1:]
            idx2=s2.find('_')
            if C1bool:
                c2=s2[idx2+1:]
            if not C1bool:
                c2=s2[:idx2]
                c1=s1[:idx1]
            if c1=='AG-347-K23' or c2=='AG-347-K23':
                continue
            n_doubletons+=1
            if c1==c2:
                print('Cell %s doubleton mistake %s %s' %(c1,s1,s2))
                #t=input('t')
                continue

            doubdiv,doubcov=0,0#look_up_ORF_pair_divergence(oldgroupname,C1bool,s1,s2)
            if doubdiv>0:
                len_notid.append(doubcov)
                if c1 in clique1_without_ribotype and c2 in clique1_without_ribotype:
                    n_notid_clique+=1
                if c1 in cl2_without_ribotype and c2 in cl2_without_ribotype:
                    n_notid_clique+=1
                if c1 in cl4_without_ribotype and c2 in cl4_without_ribotype:
                    n_notid_clique+=1
                if c1 in cl5_without_ribotype and c2 in cl5_without_ribotype:
                    n_notid_clique+=1
                if c1 in cl6_without_ribotype and c2 in cl6_without_ribotype:
                    n_notid_clique+=1
            if doubdiv==0:
                len_id.append(doubcov)
                n_id+=1
                if c1 in clique1_without_ribotype and c2 in clique1_without_ribotype:
                    n_id_clique+=1
                if c1 in cl2_without_ribotype and c2 in cl2_without_ribotype:
                    n_id_clique+=1
                if c1 in cl4_without_ribotype and c2 in cl4_without_ribotype:
                    n_id_clique+=1
                if c1 in cl5_without_ribotype and c2 in cl5_without_ribotype:
                    n_id_clique+=1
                if c1 in cl6_without_ribotype and c2 in cl6_without_ribotype:
                    n_id_clique+=1
            doubleton_divs.append(doubdiv)
            
            ind1=overalllist.index(c1)
            ind2=overalllist.index(c2)
            [dist,n]=wgmatrix[ind1][ind2]
            mypair=sorted([c1,c2])
            mypair=(mypair[0],mypair[1])
            doubleton_wg.append(dist)
            c1=mypair[0]
            c2=mypair[1]
            if mypair in doubleton_dict:
                doubleton_dict[mypair]+=1
            if mypair not in doubleton_dict:
                doubleton_dict[mypair]=1
        if len(l)==3:
            s1=l[0]
            s2=l[1]
            s3=l[2]
            #get_tripleton_orf_diams(s1,s2,s3,C1bool,groupname,tripleton_diams)
            
            idx3=s3.find('_')
            if C1bool:
                c3=s3[idx3+1:]
            idx1=s1.find('_')
            if C1bool:
                c1=s1[idx1+1:]
            idx2=s2.find('_')
            if C1bool:
                c2=s2[idx2+1:]
            if not C1bool:
                c1=s1[:idx1]
                c2=s2[:idx2]
                c3=s3[:idx3]
            if c1=='AG-347-K23' or c2=='AG-347-K23' or c3=='AG-347-K23':
                continue
            ind1=overalllist.index(c1)
            ind2=overalllist.index(c2)
            ind3=overalllist.index(c3)
            if c1==c2 or c2==c3 or c3==c1:
                print('Cells %s %s %s tripleton %s %s %s' %(c1,c2,c3,s1,s2,s3))
                t=input('t')
            [dist,n]=wgmatrix[ind1][ind2]
            [dist2,n]=wgmatrix[ind1][ind3]
            [dist3,n]=wgmatrix[ind2][ind3]
            mydiam=max([dist,dist2,dist3])
            mypair=sorted([c1,c2,c3])
            mypair=(mypair[0],mypair[1],mypair[2])
            tripleton_wg.append(mydiam)
            c1=mypair[0]
            c2=mypair[1]
            if mypair in tripleton_dict:
                tripleton_dict[mypair]+=1
            if mypair not in tripleton_dict:
                tripleton_dict[mypair]=1
    print('%s: %i singleton, %i doubleton, %i tripleton' %(groupname,n_sing,n_doubletons,n_trip))
    print('%i total doubletons, of which %i are in cliques; %i doubletons identical, of which %i are in cliques' %(len(doubleton_divs),n_notid_clique+n_id_clique,n_id,n_id_clique))

    #print('Average length of identical doubletons: %f; average length of non identical doubles: %f' %(repo.avelist(len_id),repo.avelist(len_notid)))


    if 0:#C1bool:
        for item in sorted(list(tot_lengths.keys())):
            print('%i Cells, Ave Length %f' %(item,repo.avelist(tot_lengths[item])))
            plt.bar(item,height=repo.avelist(tot_lengths[item]),color='blue')
        plt.xlabel('Number of Cells in Orthologous Cluster')
        plt.ylabel('Average Length of ORFs in Clusters')
        plt.title('%s Average Length of Orthologous Gene Groups\nAs a Function of Number of Cells' %groupname)
        plt.show()

        
        plt.hist([len_id,len_notid],bins=40,label=['Identical','Not Identical'])
        plt.xlabel('Length of Match')
        plt.ylabel('Number of doubletons')
        plt.title('%s Lengths of Doubleton Flexible ORFs' %groupname)
        plt.legend()
        plt.show()


        plt.hist([len_id,len_notid],bins=60,label=['Identical','Not Identical'],density=1)
        plt.xlabel('Length of Match')
        plt.ylabel('Number of doubletons')
        plt.title('Normalized %s Lengths of Doubleton Flexible ORFs' %groupname)
        plt.legend()
        plt.show()

        logbins=np.logspace(1,np.log10(max([max(len_id),max(len_notid)])),40)
        plt.hist([len_id,len_notid],bins=logbins,label=['Identical','Not Identical'],density=1)
        plt.xlabel('Length of Match')
        plt.ylabel('Number of doubletons')
        plt.title('Normalized %s Lengths of Doubleton Flexible ORFs' %groupname)
        plt.legend()
        plt.xscale('log')
        plt.show()
    
    #print(doubleton_dict)
    if not C1bool:
        print('Doubletons')
        group1,group2=Atlantic,Pacific
        check_same_sample_same_ocean(doubleton_dict,group1,group2)
        get_doubleton_ocean_null(group1,group2,celllist)
        print('Tripletons')
        check_same_sample_same_ocean(tripleton_dict,group1,group2)
        get_tripleton_ocean_null(group1,group2,celllist)
        
    
    plt.hist(np.array([doubleton_wg,totdivs],dtype=object),label=['Pairs with Doubleton Genes','All Cell Pairs'],bins=15,weights=[[1.0 for i in range(len(doubleton_wg))],[float(len(doubleton_wg))/len(totdivs) for i in range(len(totdivs))]])#density=True)
    plt.xlabel('Whole Genome %s Divergence' %(median_mean))
    plt.ylabel('Number of Cell Pairs, Normalized')
    plt.legend()
    plt.title('%s WG Divergences of Cells with Doubleton Genes\n %i, %i entries' %(groupname,len(doubleton_wg),len(totdivs)))
    plt.show()

  
    
    plt.hist(np.array([tripleton_wg,tottriple_diams],dtype=object) ,label=['Pairs with Tripleton Genes','All Cell Trios'],bins=15,weights=[[1.0 for i in range(len(tripleton_wg))],[float(len(tripleton_wg))/len(tottriple_diams) for i in range(len(tottriple_diams))]])#density=1)
    plt.xlabel('Maximum Whole Genome %s Divergence' %(median_mean))
    plt.legend()
    plt.ylabel('Number of Cell Pairs, Normalized')
    plt.title('%s WG Max Divergences of Cells with Tripleton Genes' %groupname)
    plt.show()

    for item in doubleton_dict:
        if doubleton_dict[item]>10:
            print(item,doubleton_dict[item])
            continue

    myl=[doubleton_dict[item] for item in doubleton_dict]
    maxl=max(myl)
    mybins=[0.5]
    for i in range(maxl+1):
        mybins.append(i+.5)
    doubleton_sum=sum([doubleton_dict[item] for item in doubleton_dict])
    print('%i doubletons total' %doubleton_sum,n_doubletons)
    plt.hist([doubleton_dict[item] for item in doubleton_dict],bins=mybins)
    plt.title('Number of doubletons per cell pair in %s if nonzero (%i pairs)' %(groupname,float(len(celllist)*(-1+len(celllist)))*.5))
    plt.xlabel('Number of doubletons per cell pair')
    plt.ylabel('Number of cell pairs')
    plt.show()


    for item in tripleton_dict:
        if tripleton_dict[item]>1:
            print(item,tripleton_dict[item])

    ntrios=float(len(celllist)*(len(celllist)-1)*(len(celllist)-2))/6
    myl=[tripleton_dict[item] for item in tripleton_dict]
    maxl=max(myl)
    mybins=[0.5]
    for i in range(maxl+1):
        mybins.append(i+.5)
    plt.hist([tripleton_dict[item] for item in tripleton_dict],bins=mybins)
    plt.title('Number of tripletons per cell trio in %s if nonzero (%i trios)' %(groupname,ntrios))
    plt.xlabel('Number of tripletons per cell trio')
    plt.ylabel('Number of cell trios')
    plt.show()

    
    plt.hist(doubleton_divs,bins=40)
    plt.xlabel('BLAST nucleotide divergence of doubletons')
    plt.ylabel('Number of doubletons')
    plt.title('%s Doubleton Divergences using non-core Genes' %(groupname))
    plt.show()


    d_nonzero=[doubleton_divs[i] for i in range(len(doubleton_divs))  if doubleton_divs[i]>0]
    logbins=np.logspace(np.log10(min(d_nonzero))-1,0,40)
    #print(logbins)
    #print('Divs')
    #print(doubleton_divs)
    #print(np.log10(min(d_nonzero))-0.5,doubleton_divs.count(0),np.log10(min(d_nonzero))-0.75)
    plt.hist(doubleton_divs,bins=logbins)
    plt.bar((min(d_nonzero))*0.25,doubleton_divs.count(0),width=min(d_nonzero)*0.15,label='Identical',color='orange')
    # plt.bar(np.log10(min(d_nonzero))-0.5,doubleton_divs.count(0),width=np.log10(min(d_nonzero))-0.75,color='orange',label='Identical')
    plt.legend()
    plt.xscale('log')
    plt.xlabel('BLAST nucleotide divergence of doubletons')
    plt.ylabel('Number of doubletons')
    plt.title('%s Doubleton Divergences using non-core Genes\n%i Total, %.2f of which Identical' %(groupname,len(doubleton_divs),float(doubleton_divs.count(0))/len(doubleton_divs)))
    plt.show()


   
    plt.hist(tripleton_diams,bins=40)#logbins)
    plt.xscale('log')
    #plt.xlabel('Maximum nucleotide divergence of tripletons')
    plt.ylabel('Number of tripletons')
    plt.title('%s Tripleton Divergence Diameters using non-core Genes' %(groupname))
    plt.show()

def check_same_sample_same_ocean(mydict,group1,group2):
    #each entry in the dict is a pair, triple, etc of cells. check their sample and if same ocean, for berube
    samesample=0
    differentsample=0
    sameocean=0
    differentocean=0
    for item in mydict:
        mylen=len(item)
        samples=[]
        oceans=[]
        for i in range(mylen):
            c=item[i]
            sample=c[3:6]
            samples.append(sample)
            if c in group1:#Atlantic:
                oceans.append('Atlantic')
            if c in group2:#N_Pacific or c in S_Pacific:
                oceans.append('Pacific')
        if samples.count(samples[0])==len(samples):
            samesample+=1
        if samples.count(samples[0])!=len(samples):
            differentsample+=1
        if  len(oceans)<mylen:
            #print(item)
            #t=input('t')
            continue
       
        if oceans.count(oceans[0])==len(oceans):
            sameocean+=1
        if oceans.count(oceans[0])!=len(oceans):
            differentocean+=1
    print('%i groups considered. %f from same sample, %f from same ocean, %i (%f) not in Atl/Pac' %(len(list(mydict.keys())),float(samesample)/(samesample+differentsample),float(sameocean)/(sameocean+differentocean),len(list(mydict.keys()))-sameocean-differentocean,(len(list(mydict.keys()))-sameocean-differentocean)/float(len(list(mydict.keys())))))
    

def get_orf_len_dict(groupname,celllist):
    print(groupname)
    if groupname=='C1':
        infile='../input_data/All_nonC1Core_ORFs.fasta'
        
    with open(infile,'r') as myinfile:
        lines=myinfile.read()
    cellseqs=re.split('>',lines)
    ldict={}
    for i in range(len(cellseqs)):
        sublist=re.split('\n',cellseqs[i])
        seq=''.join(sublist[1:])
        ldict[sublist[0]]=len(seq)
    return ldict
    
def get_doubleton_ocean_null(group1,group2,celllist):
    #wtk: how many of the doubletons do we expect to be from the same ocean?
    '''
    given: number of doubletons. fraction of cells from each ocean
    ideally: coverage of each cell
    then: we have total number of non core ORFs from PAC, Atl.  for a given doubleton, its chance from the same ocean is (number of orfs from same ocean-number of this cell's orfs)/(number of orfs from different ocean+number from same - this cells number of orfs)
    '''
    group1=list(set(group1).intersection(set(celllist)))#cells from ocean 1, ocean 2
    group2=list(set(group2).intersection(set(celllist)))
    
    with open('../input_data/HLII_Noncore_ORFs_per_Cell.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    ndict={}
    for item in lines:
        l=re.split(' ',item)
        ndict[l[-1][1:]]=int(l[-2])#just gets cell: number of noncore ORFs
    count1,count2=0,0
    for item in group1:
        count1+=ndict[item]
    for item in group2:
        count2+=ndict[item]
    #number of noncore orfs total in each ocean
    sampledict={}
    for item in ndict:
        if item not in celllist:
            continue
        samp=item[3:6]
        if samp in sampledict:
            sampledict[samp]+=ndict[item]
        if samp not in sampledict:
            sampledict[samp]=ndict[item]
    #wtk: expected number from same ocean, same sample
    #expected same ocean by chance
    totcount=sum([ndict[item] for item in celllist])

    f1,f2=float(count1)/(totcount),float(count2)/(totcount)
    p_same_ocean=f1*((float(count1)-1.0)/(totcount-1))+f2*((float(count2)-1.0)/(totcount-1))
    
    #f1,f2=float(count1)/(count1+count2),float(count2)/(count1+count2)
    #p_same_ocean=f1*(f1-1.0/(count1+count2-1))+f2*(f2-1.0/(count1+count2-1))
    print('Probability a doubleton is from the same ocean if random shuffle: %f' %(p_same_ocean))
    #now same sample
    totsample=sum([sampledict[item] for item in sampledict])
    p_same_sample=0
    for item in sampledict:
        p_same_sample+=float(sampledict[item]*(sampledict[item]-1))/(totsample*(totsample-1))
    print('Probability a doubleton is from the same sample if random shuffle: %f' %(p_same_sample))


def get_tripleton_ocean_null(group1,group2,celllist):
    group1=list(set(group1).intersection(set(celllist)))
    group2=list(set(group2).intersection(set(celllist)))
    with open('../input_data/HLII_Noncore_ORFs_per_Cell.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    ndict={}
    for item in lines:
        l=re.split(' ',item)
        ndict[l[-1][1:]]=int(l[-2])
    count1,count2=0,0
    for item in group1:
        count1+=ndict[item]
    for item in group2:
        count2+=ndict[item]
    sampledict={}
    for item in ndict:
        if item not in celllist:
            continue
        samp=item[3:6]
        if samp in sampledict:
            sampledict[samp]+=ndict[item]
        if samp not in sampledict:
            sampledict[samp]=ndict[item]
    #wtk: expected number from same ocean, same sample
    #expected same ocean by chance

    totcount=sum([ndict[item] for item in celllist])

    f1,f2=float(count1)/(totcount),float(count2)/(totcount)
    p_same_ocean=f1*(f1-1.0/(totcount-1))*(f1-2.0/(totcount-2))+f2*(f2-1.0/(totcount-1))*(f2-2.0/(totcount-2))
    
    #f1,f2=float(count1)/(count1+count2),float(count2)/(count1+count2)
    #p_same_ocean=f1*(f1-1.0/(count1+count2-1))*(f1-2.0/(count1+count2-2))+f2*(f2-1.0/(count1+count2-1))*(f2-2.0/(count1+count2-2))
    print('Probability a Tripleton is from the same ocean if random shuffle: %f' %(p_same_ocean))
    #now same sample
    totsample=sum([sampledict[item] for item in sampledict])
    p_same_sample=0
    for item in sampledict:
        p_same_sample+=float(sampledict[item]*(sampledict[item]-1))/(totsample*(totsample-1))*(float(sampledict[item]-2))/(totsample-2)
    print('Probability a Tripleton is from the same sample if random shuffle: %f' %(p_same_sample))
    
def load_WG_divs(celllist,median_mean,tfna='a'):
    overalllist=distmatrepo.overall_order_list()
    infilewg='../input_data/Cell_Pair_Matrices/Berube_Kash_WG_Matrices_%s.npy' %(tfna)
    wgarray=distmatrepo.load_arrays(infilewg)
    distmat=wgarray['%sWG_%s' %(median_mean,tfna)]
    totdivs=[]
    tottriple_diams=[]
    for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            [d,n]=distmat[idxi][idxj]
            totdivs.append(d)
            for k in range(j):
                idxk=overalllist.index(celllist[k])
                [d2,n]=distmat[idxk][idxi]
                [d3,n]=distmat[idxj][idxk]
                maxd=max([d,d2,d3])
                tottriple_diams.append(maxd)
    return overalllist,distmat,totdivs,tottriple_diams
    
###########################################################################################

def get_tripleton_orf_diams(s1,s2,s3,C1bool,groupname,tripleton_diams):
    d1,c1=look_up_ORF_pair_divergence(groupname,C1bool,s1,s2)
    d2,c2=look_up_ORF_pair_divergence(groupname,C1bool,s3,s2)
    d3,c3=look_up_ORF_pair_divergence(groupname,C1bool,s1,s3)
    if max([d1,d2,d3])==0:
        return
    tripleton_diams.append(max([d1,d2,d3]))

def look_up_ORF_pair_divergence(groupname,C1bool,orf1,orf2):
    if C1bool:
        infileblast='../input_data/%s_nonCore_ORFs_All_Against_All_2.out' %groupname
    if not C1bool:
        infileblast='Berube_HLII_NonCore/%s_nonCore_ORFs_All_Against_All_notask.out' %groupname
    fieldnames=['qseqid','sseqid','pident','qcovhsp','length']
    with open(infileblast,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t',fieldnames=fieldnames)
        for row in csvreader:
            query=row['qseqid']
            #print(orf1,orf2,query)
            #t=input('t')
            if query !=orf1 and query!=orf2:
                continue
            sseq=row['sseqid']
            if sseq!=orf1 and sseq!=orf2:
                continue
            if sseq==query:
                continue
            pid=float(row['pident'])
            qcov=float(row['qcovhsp'])
            #print(qcov)
            length=int(row['length'])
            return 1.0-0.01*pid,length#qcov*0.01
    print(orf1,orf2,' not found')
    #t=input('t')
    return 0,0


def dictionary_ORF_pair_divergence(celllist,C1bool,groupname):
    if C1bool:
        infileblast='../input_data/%s_nonCore_ORFs_All_Against_All_2.out' %groupname
    if not C1bool:
        infileblast='Berube_HLII_NonCore/%s_nonCore_ORFs_All_Against_All_notask.out' %groupname
    fieldnames=['qseqid','sseqid','pident','qcovhsp','length']
    mydict={}#want to get queries for all cellsin celllist...
    
    with open(infileblast,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t',fieldnames=fieldnames)
        for row in csvreader:
            query=row['qseqid']
            sseq=row['sseqid']
            if sseq==query:
                continue
            if not check_if_in_celllist(query,celllist) or not check_if_in_celllist(sseq,celllist):
                continue          
            pid=float(row['pident'])
            qcov=float(row['qcovhsp'])
            length=int(row['length'])
            mykey=sorted([query,sseq])
            mynewkey=(mykey[0],mykey[1])
            mydict[mynewkey]=[ 1.0-0.01*pid,length]#qcov*0.01
    return mydict

def check_if_in_celllist(seq,celllist):
    for item in celllist:
        if item in seq:
            return True
    return False

############################################################
def get_frac_identical_similar_lengths(groupname,celllist,upperlengthcut=600,lowerlengthcut=0):
   
    mygenes=[]
    with open('../input_data/MIT9301_Gene_Lengths.csv','r') as mif:
        csvreader=csv.DictReader(mif,delimiter='\t')
        for row in csvreader:
            l=int(row['Length'])
            if l>=lowerlengthcut and l<=upperlengthcut:
                mygenes.append(int(row['GeneNum']))

    titlegroupname='Berube_Kash'
    suffix='_Discrete_Including_Non_Core'
    infilematrices='../input_data/Cell_Pair_Matrices/Each_Cell_Pair_Sites/Matrices_%s_All_NT%s.npy' %(titlegroupname,suffix)
   
    genearrays=distmatrepo.load_arrays(infilematrices)
    overalllist=distmatrepo.overall_order_list()
    divlist=[]
    for item in list(genearrays.keys()):
        mygene=get_gene_from_key_to_check(item)
        if mygene not in mygenes:
            continue
        get_distances_single_mat(genearrays[item],celllist,overalllist,item,divlist)


    print('%i Genes Length %i - %i' %(len(mygenes),lowerlengthcut,upperlengthcut))
    print('%i pair divergences , of which %i (%f) identical' %(len(divlist),divlist.count(0),float(divlist.count(0))/len(divlist)))

    plt.hist(divlist,bins=40)
    plt.title('%s Pair Divergences, Genes %i-%i Length\n%i Pairs, of which %i (%f) identical' %(groupname,upperlengthcut,lowerlengthcut,len(divlist),divlist.count(0),float(divlist.count(0))/len(divlist)))
    plt.ylabel('Number of cell pairs x genes')
    plt.xlabel('Nucleotide divergence')
    plt.show()

    d_nonzero=[divlist[i] for i in range(len(divlist))  if divlist[i]>0]
    logbins=np.logspace(np.log10(min(d_nonzero))-1,0,40)
    plt.hist(divlist,bins=logbins)
    plt.bar((min(d_nonzero))*0.25,divlist.count(0),width=min(d_nonzero)*0.15,label='Identical',color='orange')
    plt.legend()
    plt.title('%s Pair Divergences, Genes %i-%i Length\n%i Pairs, of which %i (%f) identical' %(groupname,upperlengthcut,lowerlengthcut,len(divlist),divlist.count(0),float(divlist.count(0))/len(divlist)))
    plt.ylabel('Number of cell pairs x genes')
    plt.xlabel('Nucleotide divergence')
    plt.xscale('log')
    plt.show()
    
def get_distances_single_mat(mat,celllist,overalllist,mykey,divlist):
    for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            [ndiff,length]=mat[idxi][idxj]
            if length>10:
                divlist.append(float(ndiff)/length)


def get_gene_from_key_to_check(key):
    idx1=key.find('_')
    idx2=key.find('_',idx1+1)
    g=int(key[idx1+1:idx2])
    return g


####################
def look_for_cell_groups_in_orthogroups(totgroupname):
    #WTK: for a given size n cells of orthogroup, how many of the genes have the same set of n cells?
    nplotto=50
    qcut,pidcut=60,80
  
    infile='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(totgroupname,str(qcut),str(pidcut))#2
    print(infile)

    ndict={}#ncells:[list of sorted celllists]. then count occurences of celllist
    with open(infile,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    for line in lines:
        l=ast.literal_eval(line)
        myn=len(l)
        newl=[]
        for item in l:
            idx1=item.find('_')
            newl.append(item[idx1+1:])
        newl=sorted(newl)
        if myn in ndict:
           ndict[myn].append(newl)
        if myn not in ndict:
            ndict[myn]=[newl]

    ncells=sorted(list(ndict.keys()))
    
    ngrps=[len(ndict[ncells[j]]) for j in range(len(ncells))]
    fracs=[]#a list for each entry in ncells
    for i in range(len(ncells)):
        minilist=[]#to be appended to fracs
        mylist=ndict[ncells[i]]#list of cell groups
        #now count how many time each group occurs
        myset=get_uniques(mylist)#set(mylist)
        for item in myset:
            myfrac=float(mylist.count(item))/len(mylist)
            minilist.append(myfrac)
            if myfrac>0.1 and ncells[i]<nplotto and mylist.count(item)>2:
                print('%i cells Frac %f' %(ncells[i],myfrac), item)
        fracs.append(minilist)

    nplotto=len(ncells)-1
    ticklabs=[]
    for i in range(nplotto):
        mylab='%i Cells\n(%i Groups)' %(ncells[i],ngrps[i])
        ticklabs.append(mylab)
        myfracs=fracs[i]
        for f in myfracs:
            plt.plot(ncells[i],f,'o',color='blue',alpha=0.3)
    plt.ylabel('Fraction of Orthogroups containing the same cells')
    ax=plt.gca()
    ax.set_xticks(ncells[:nplotto])
    ax.set_xticklabels(ticklabs)
    plt.title('%s Flexible Gene Orthogroups: For grps of a\nn cells, how often are these n cells the same?' %totgroupname)
    plt.show()
        
def get_uniques(mainlist):
    ulist=[]
    for item in mainlist:
        if item not in ulist:
            ulist.append(item)
    return ulist

def look_for_cliques_in_orthogroups(celllist_list,groupname_list,totgroupname):
    qcut,pidcut=60,80
  
    infile='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(totgroupname,str(qcut),str(pidcut))#2
    print(infile)

    tot_cells_in_orthogroup=[]
    clique_cells_in_orthogroup={}
    for item in groupname_list:
        clique_cells_in_orthogroup[item]=[]
    with open(infile,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    for line in lines:
        l=ast.literal_eval(line)
        tot_cells_in_orthogroup.append(len(l))
        newl=[]
        for item in l:
            idx1=item.find('_')
            newl.append(item[idx1+1:])
        for i in range(len(celllist_list)):
            cl=celllist_list[i]
            gn=groupname_list[i]
            noverlap=len(list(set(cl).intersection(set(newl))))
            clique_cells_in_orthogroup[gn].append(noverlap)
    #how to plot? WTK: for each clique, distribution of number of cells per orthogroup vs the total size of the orthogroup
    #also: orthogroups that draw from multiple cliques

    for i in range(len(groupname_list)):
        
        item=groupname_list[i]
        #this might be better as a violin plot series
        plt.plot(clique_cells_in_orthogroup[item],tot_cells_in_orthogroup,'o',alpha=0.3)
        plt.ylabel('Total number of %s cells in Orthogroup' %totgroupname)
        plt.xlabel('Number of %s Cells in Orthogroup (of %i)' %(item,len(celllist_list[i])))
        plt.title('Flexible Gene Orthogroups in %s:\nLooking for groups from %s' %(totgroupname,item))
        x=np.linspace(0,max(tot_cells_in_orthogroup),100)
        plt.plot(x,x,'-')

        plt.show()


        myncells=list(range(len(celllist_list[i])+1))
        mydists=[]
        for littlen in myncells:
            littlelist=[]
            for k in range(len(tot_cells_in_orthogroup)):
                if clique_cells_in_orthogroup[item][k]==littlen:
                    littlelist.append(tot_cells_in_orthogroup[k])
            if littlelist==[]:
                littlelist=[float('nan'),float('nan')]
            mydists.append(littlelist)
            if littlen==len(celllist_list[i]) or littlen==len(celllist_list[i])-1:
                print('%s: %i genes with %i cells exclusively' %(item,littlelist.count(littlen),littlen))
                      
        
        plt.violinplot(mydists,positions=myncells,showmedians=True)
        plt.ylabel('Total number of %s cells in Orthogroup' %totgroupname)
        plt.xlabel('Number of %s Cells in Orthogroup (of %i)' %(item,len(celllist_list[i])))
        plt.title('Flexible Gene Orthogroups in %s:\nLooking for groups from %s' %(totgroupname,item))
        x=np.linspace(0,max(tot_cells_in_orthogroup),100)
        plt.plot(x,x,'-')
        plt.show()
        
        #now make a plot with the y axis the total number of other cliqe cells in the orthogroup
        othercliques=[]
        for j in range(len(tot_cells_in_orthogroup)):
            myn=0
            for otheritem in groupname_list:
                if otheritem==item:
                    continue
                myn+=clique_cells_in_orthogroup[otheritem][j]
            othercliques.append(myn)
        plt.plot(clique_cells_in_orthogroup[item],othercliques,'o',alpha=0.3)
        plt.ylabel('Number of cells from other %s Cliques in the Orthogroup' %totgroupname)
        plt.xlabel('Number of %s Cells in Orthogroup (of %i)' %(item,len(celllist_list[i])))
        plt.title('Flexible Gene Orthogroups in %s:\nSharing orthogroups between cliques and %s' %(totgroupname,item))
        plt.show()

        mydists=[]
        for littlen in myncells:
            littlelist=[]
            for k in range(len(tot_cells_in_orthogroup)):
                if clique_cells_in_orthogroup[item][k]==littlen:
                    littlelist.append(othercliques[k])
            if littlelist==[]:
                littlelist=[float('nan'),float('nan')]
            mydists.append(littlelist)
           

        plt.violinplot(mydists,positions=myncells,showmedians=True)
        plt.ylabel('Number of cells from other %s Cliques in the Orthogroup' %totgroupname)
        plt.xlabel('Number of %s Cells in Orthogroup (of %i)' %(item,len(celllist_list[i])))
        plt.title('Flexible Gene Orthogroups in %s:\nSharing orthogroups between cliques and %s' %(totgroupname,item))
        plt.show()

        t=input('t')



######################################

def plot_flexible_genes_spanning_ribotypes(totgroupname,sizecut,fraccut,celllist):
    #wtk: for flex genes in more than 1 of C123, what are the divergences within vs between? shite
    qcut,pidcut=60,80
  
    infile='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(totgroupname,str(qcut),str(pidcut))#2
    print(infile)

    mydict=dictionary_ORF_pair_divergence(celllist,True,totgroupname)#divergences

    outfile='%s_Inter_Intra_Divergences_Flexible_Genes.csv' %totgroupname
    fieldnames=['NC1','NC2','NC3','Intra','Inter']
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
        outwriter.writeheader()

        with open(infile,'r') as myinfile:
            lines=[line.rstrip() for line in myinfile if len(line)>0]
        for line in lines:
            l=ast.literal_eval(line)

            #wtk: span C123? how many cells? look at large clusters first
            lookbool,minilists=select_orthogroup(l,sizecut,fraccut)
            if not lookbool:
                continue
            print('Looking into orthogroup size %i' %len(l))
            #now get intra and interdivs and plot

            intra,inter=get_inter_intra_divergences(minilists,totgroupname,mydict)

            outwriter.writerow({'NC1':len(minilists[0]),'NC2':len(minilists[1]),'NC3':len(minilists[2]),'Intra':str(intra),'Inter':str(inter)})

            
            '''
            plt.hist([intra,inter],bins=50,label=['Intra-Ribotype','Inter-Ribotype'])
            plt.legend()
            plt.title('%s Divergences at Flexible Gene Orthogroup with %i Cells' %(totgroupname,len(l)))
            plt.ylabel('Number of Cell Pairs')
            plt.xlabel('Blast nucleotide Divergence')
            plt.show()
            t=input('t')
            '''
        #########
       
          

def get_inter_intra_divergences(minilists,totgroupname,mydict):
    intra=[]
    inter=[]
    for i in range(len(minilists)):
        print('new ribotype')
        li=minilists[i]
        for ii in range(len(li)):
            #doing all cells not
            for jj in range(ii):
                cii=li[ii]
                cjj=li[jj]
                div,length=new_look_up_ORF_pair_divergence(cii,cjj,mydict)#totgroupname,True,cii,cjj)
                if length>11:
                    intra.append(div)
        for j in range(i):
            lj=minilists[j]
            for iii in range(len(li)):
                for jjj in range(len(lj)):
                    ciii=li[iii]
                    cjjj=lj[jjj]
                    div,length=new_look_up_ORF_pair_divergence(ciii,cjjj,mydict)#totgroupname,True,ciii,cjjj)
                    if length>11:
                        inter.append(div)
    return intra,inter
                
def new_look_up_ORF_pair_divergence(orf1,orf2,mydict):
    mykey=sorted([orf1,orf2])
    mynewkey=(mykey[0],mykey[1])
    if mynewkey not in mydict:
        return 0,0
    l=mydict[mynewkey]
    return l[0],l[1]

def select_orthogroup(grouplist,sizecut,fraccut):
    #sizecut is minimum number of cells in orthogroup.  fraccut is the max fraction of cells that can be from the same orthogroup
    l1,l2,l3=C1_without_ribotype2,C2_without,C3_without
    if len(grouplist)<sizecut:# or len(grouplist)>20:###edt
        return False,[]
    minilists=[[],[],[]]
    for item in grouplist:
        idx1=item.find('_')
        mycell=item[idx1+1:]
        if mycell in l1:
            minilists[0].append(item)
        if mycell in l2:
            minilists[1].append(item)
        if mycell in l3:
            minilists[2].append(item)
    myfracs=[float(len(minilists[i]))/len(grouplist) for i in range(len(minilists))]
    maxfrac=max(myfracs)
    if maxfrac>fraccut:
        return False,[]
    if len(minilists[0])>20:
        return False, []#a C1 core gene ish
    return True,minilists

###############

def estimate_flexgenes_per_cell(celllist,groupname):
    outfile='%s_Core_Coverage_using_1581_Core_Genes.csv' %groupname
    fieldnames=['Cell','CoreCov']
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
        outwriter.writeheader()
        for c in celllist:
            flexfile='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/All_NT_Close_Pair_Blasting/nonC1Core_Fastas/%s_nonC1Core_AllNT.fna' %c
            corefile='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/All_NT_Close_Pair_Blasting/YESC1Core_Fastas/%s_YESC1Core_AllNT.fna' %c
            ncore=count_seqs(corefile)
            if ncore=='-':
                continue
            outwriter.writerow({'Cell':c,'CoreCov':float(ncore)/1581.})
        #there are supposed to be 1581 (number in reheaded_converted_core.faa
    
def count_seqs(infile):
    totseqdict={}#all nucleotides
    if not exists(infile):
        print('File %s does not exist' %(infile))
        return '-'
    with open(infile,'r') as myinfile:
        lines=myinfile.read()
    cellseqs=re.split('>',lines)
    return len(cellseqs)

def get_covered_core():
    with open('../input_data/HLII_Cells_Covered_Per_Orthogroup.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    covs=[]
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    for line in lines:
        idx1=line.find('_')
        gene=int(line[4:idx1])
        if gene not in mycores:
            continue
        idx2=line.find(':')
        cov=int(line[idx2+1:])
        covs.append(cov)
    return covs

###########################################################################################
if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run.')
    args = parser.parse_args()
    print(args)
    prog=args.prog

    if prog==0:
        excessbool=1
        #groupname='C1'
        #celllist=C1_without_ribotype2

        #groupname='C3'
        #celllist=C3_without

        celllist=C1_without_ribotype2+C2_without+C3_without
        groupname='C1_C2_C3'
        
        #rehead_fasta(celllist)
        
        pidcut=40
        qcut=40
        get_connected_components_from_all_v_all_BLAST(qcut,pidcut,groupname,celllist,C1bool=1)


        singleton_ORFs_per_cell(celllist,groupname)

        #{0: 64, 1: 7, 2: 58, 3: 0, 4: 77, 5: 103, 6: 18, 7: 80, 8: 10, 9: 70, 10: 35, 11: 22, 12: 61, 13: 88, 14: 117, 15: 58, 16: 27, 17: 42, 18: 30, 19: 34, 20: 107, 21: 14, 22: 0, 23: 143, 24: 77, 25: 30, 26: 68, 27: 83, 28: 93, 29: 82, 30: 59, 31: 19, 32: 104, 33: 941, 34: 6, 35: 49, 36: 116, 37: 23, 38: 51, 39: 9, 40: 87, 41: 36, 42: 23, 43: 33, 44: 39, 45: 17, 46: 53, 47: 44, 48: 84, 49: 2, 50: 30, 51: 3, 52: 22}

    if prog==1:
        groupname='HLII'
        celllist=hliilist2
        #groupname='Blue_Basin_NoClose'
        #celllist=BB_list
        
        pidcut,qcut=80,60
        #make_Berube_allgenes_fasta(celllist,groupname)
        #get_connected_components_from_all_v_all_BLAST(qcut,pidcut,groupname,celllist)
        make_berube_cluster_size_plots(celllist,groupname,qcut,pidcut)

    if prog==2:
        groupname='C1'#_bidir'
        C1bool=1
        celllist=C1_without_ribotype2
        #celllist.remove(C1_without_ribotype[33])
        returnnum=10
        qcutlist=[40,60,80]
        pidcutlist=[40,80]
        #plot_singletons_doubletons_as_function_of_qcov_pid(qcutlist,pidcutlist,groupname,celllist,C1bool,returnnum)

        examine_doubletons_etc(celllist,groupname,C1bool)
        #get_frac_identical_similar_lengths(groupname,celllist)

    if prog==-2:
        C1bool=0
        #celllist=BB_list
        #groupname='Blue_Basin_NoClose'
        celllist=hliilist2
        groupname='HLII'
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
        examine_doubletons_etc(celllist,groupname,C1bool)

        '''
       
        17744 orthogroups
90189 total flex genes. 14399 singletons. 0.159654 genes are singletons, 0.811486 groups are singletons
        '''
    if prog==3:
        groupname='HLII'#_bidir'
        C1bool=0
        celllist=hliilist2

        returnnum=10
        qcutlist=[40,60,80]
        pidcutlist=[40,80]
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
        plot_singletons_doubletons_as_function_of_qcov_pid(qcutlist,pidcutlist,groupname,celllist,C1bool,returnnum)
       
        #examine_doubletons_etc(celllist,groupname,C1bool)

    
        
    if prog==4:
        groupname='Blue_Basin_NoClose'#BB_bidir'
        C1bool=0
        celllist=BB_list

        returnnum=10
        qcutlist=[40,60,80]
        pidcutlist=[40,80]
        plot_singletons_doubletons_as_function_of_qcov_pid(qcutlist,pidcutlist,groupname,celllist,C1bool,returnnum)

    if prog==5:
        get_doubleton_ocean_null()

    if prog==6:
        totgroupname='C1'
        celllist_list=[clique1_without_ribotype,cl2_without_ribotype,cl4_without_ribotype,cl5_without_ribotype,cl6_without_ribotype]
        groupname_list=['Clique1','Clique2','Clique4','Clique5','Clique6']
        look_for_cliques_in_orthogroups(celllist_list,groupname_list,totgroupname)

    if prog==-6:#instead of cliques, C1-C3
        totgroupname='C1_C2_C3'
        celllist_list=[C1_without_ribotype2,C2_without,C3_without]
        groupname_list=['C1','C2','C3']
        look_for_cliques_in_orthogroups(celllist_list,groupname_list,totgroupname)

        
    if prog==7:
        #totgroupname='C1'
        #totgroupname='C3'
        totgroupname='C1_C2_C3'
        look_for_cell_groups_in_orthogroups(totgroupname)


    if prog==8:
        totgroupname='C1_C2_C3'
        sizecut=10
        fraccut=.79
        celllist=C1_without_ribotype2+C2_without+C3_without
        plot_flexible_genes_spanning_ribotypes(totgroupname,sizecut,fraccut,celllist)

    if prog==9:
        groupname="C1"
        celllist=C1_without_ribotype2
        estimate_flexgenes_per_cell(celllist,groupname)
'''
are we excluding the core genes?
hlii non core: 245356 orfs. 167 cells. 1469 ORFs/cell

all cells all genes 873996

found problem: core_orfs were strings but myorf was int
now 91388 ORFs; 547 per cell

41048 orfs for bb no close non 9301 core

c1 non core: 79583. 1530 per cell wtf
found problem: redoing for C1 with C1-allnt_with_orthogroups file
no 20367 ORFs, about 400/cell
'''



'''
c1 doubletons (80,80,bidir)
('B241_526B17', 'B241_529O19') 26 (not in close cliques)
('B241_527L16', 'B243_495N16') 43 (clique 1 close pair)
('B241_527P5', 'B241_528K19') 33 (cl1 close pair)
('B241_529J11', 'B245a_518E10') 40 (cl2 close pair)---found by kashtan
('B243_496N4', 'B243_498L10') 20 (not in close cliques)
('B243_498F21', 'B245a_521O20') 19 (O20 in cl1; F21 not in close clique)



triples new? 3-13-24 (qcov redone 60,80)

'B241_526B17', 'B241_527L15', 'B241_529J16') 2
('B241_526B22', 'B241_528N8', 'B245a_520B18') 3
('B241_526B22', 'B243_498F21', 'B245a_520B18') 6
('B241_527G5', 'B241_527N11', 'B241_528N8') 2
('B241_527L15', 'B243_496E10', 'B243_498J20') 2
('B241_527L15', 'B243_498P3', 'B245a_518K17') 3
('B241_527L15', 'B245a_518K17', 'B245a_521M10') 2
('B241_527L16', 'B241_528N20', 'B243_496E10') 2
('B241_527P5', 'B241_528K19', 'B245a_521K15') 2
('B241_527P5', 'B243_496E10', 'B243_498L10') 3
('B241_528N8', 'B241_529J16', 'B243_495N16') 2
('B241_529C4', 'B243_495K23', 'B245a_520E22') 2
('B243_495K23', 'B243_496E10', 'B245a_520E22') 2
('B243_496N4', 'B243_498F21', 'B243_498L10') 2



BB:
587 total doubletons, 
Doubletons
360 groups considered. 0.197222 from same sample, 0.520492 from same ocean, 116 (0.322222) not in Atl/Pac
Probability a doubleton is from the same ocean if random shuffle: 0.365313
Probability a doubleton is from the same sample if random shuffle: 0.227431
Tripletons
124 groups considered. 0.080645 from same sample, 0.287671 from same ocean, 51 (0.411290) not in Atl/Pac
Probability a Tripleton is from the same ocean if random shuffle: 0.168825
Probability a Tripleton is from the same sample if random shuffle: 0.065759
('AG-347-L02', 'AG-402-N17') 22
('AG-355-I04', 'AG-355-N02') 37
587 doubletons total 587
('AG-347-E23', 'AG-347-I04', 'AG-355-A09') 4
('AG-347-E23', 'AG-347-K19', 'AG-355-A09') 3
('AG-347-I06', 'AG-355-P15', 'AG-424-P16') 2
('AG-347-I21', 'AG-347-M23', 'AG-355-A09') 2
('AG-347-L02', 'AG-402-A04', 'AG-402-N17') 7
('AG-347-L02', 'AG-418-J17', 'AG-449-O05') 3
('AG-347-L02', 'AG-355-N16', 'AG-402-N17') 9
('AG-347-G20', 'AG-347-L17', 'AG-347-L20') 9
('AG-355-A09', 'AG-402-K16', 'AG-418-D13') 2
('AG-347-G20', 'AG-355-I04', 'AG-402-A04') 2
('AG-355-I04', 'AG-355-J09', 'AG-355-P18') 3
('AG-355-J09', 'AG-402-K16', 'AG-418-J17') 14
('AG-347-J06', 'AG-355-O17', 'AG-442-B03') 3
('AG-355-B23', 'AG-402-K16', 'AG-418-D13') 4
('AG-418-J17', 'AG-442-N07', 'AG-459-B06') 8
('AG-347-C10', 'AG-402-A04', 'AG-402-N23') 2


hlii

1408 total doubletons, 
Doubletons
847 groups considered. 0.168831 from same sample, 0.545455 from same ocean, 319 (0.376623) not in Atl/Pac
Probability a doubleton is from the same ocean if random shuffle: 0.325808
Probability a doubleton is from the same sample if random shuffle: 0.170578
Tripletons
420 groups considered. 0.040476 from same sample, 0.310185 from same ocean, 204 (0.485714) not in Atl/Pac
Probability a Tripleton is from the same ocean if random shuffle: 0.138615
Probability a Tripleton is from the same sample if random shuffle: 0.039429

('AG-347-L02', 'AG-402-N17') 22
('AG-347-J23', 'AG-355-A09') 32
('AG-355-I04', 'AG-355-N02') 36
('AG-402-F05', 'AG-402-L23') 21
('AG-402-N08', 'AG-418-I21') 20
('AG-355-P15', 'AG-418-I20') 13
('AG-347-J05', 'AG-347-M15') 24
('AG-412-A14', 'AG-449-D16') 17
('AG-347-C10', 'AG-670-M15') 11
('AG-347-K02', 'AG-424-P23') 16
('AG-347-M15', 'AG-459-N19') 18
('AG-418-F16', 'AG-436-E22') 11
1408 doubletons total 1408
('AG-347-E23', 'AG-347-I04', 'AG-355-A09') 4
('AG-347-E23', 'AG-347-K19', 'AG-355-A09') 3
('AG-347-G18', 'AG-402-F05', 'AG-418-D13') 2
('AG-347-G18', 'AG-355-N18', 'AG-355-N23') 3
('AG-347-I04', 'AG-424-A03', 'AG-459-N19') 2
('AG-347-I19', 'AG-347-K17', 'AG-347-L19') 7
('AG-347-I21', 'AG-347-K17', 'AG-402-K16') 2
('AG-347-J14', 'AG-418-I20', 'AG-449-D16') 2
('AG-347-J19', 'AG-355-L02', 'AG-388-A01') 8
('AG-347-J20', 'AG-347-K02', 'AG-424-P23') 5
('AG-347-K17', 'AG-418-O03', 'AG-429-C19') 3
('AG-347-K20', 'AG-402-G23', 'AG-449-K21') 3
('AG-347-L02', 'AG-402-A04', 'AG-402-N17') 7
('AG-347-L02', 'AG-355-N16', 'AG-402-N17') 5
('AG-347-G20', 'AG-347-L17', 'AG-347-L20') 9
('AG-355-B18', 'AG-355-J09', 'AG-355-N18') 2
('AG-347-J05', 'AG-355-I04', 'AG-402-O16') 2
('AG-355-I04', 'AG-402-O16', 'AG-459-D04') 4
('AG-355-I04', 'AG-355-J09', 'AG-355-P18') 2
('AG-355-K20', 'AG-402-N08', 'AG-418-I21') 3
('AG-355-L02', 'AG-355-L21', 'AG-355-P11') 2
('AG-388-F11', 'AG-402-I23', 'AG-418-P13') 2
('AG-402-F05', 'AG-402-L23', 'AG-412-J13') 15
('AG-347-C10', 'AG-402-K22', 'AG-402-L23') 2
('AG-355-K23', 'AG-402-N08', 'AG-449-J16') 2
('AG-355-B23', 'AG-402-K16', 'AG-418-D13') 4
('AG-335-I15', 'AG-418-J17', 'AG-418-M21') 2
('AG-335-I15', 'AG-347-I23', 'AG-418-M21') 2
('AG-418-J17', 'AG-442-N07', 'AG-459-B06') 8
('AG-347-J05', 'AG-347-M15', 'AG-459-N19') 12
('AG-347-K22', 'AG-418-O03', 'AG-429-C19') 5
('AG-355-M18', 'AG-429-E20', 'AG-449-D16') 2
('AG-402-C22', 'AG-424-A14', 'AG-432-K16') 3
('AG-355-M18', 'AG-355-N18', 'AG-449-D22') 2
('AG-402-F05', 'AG-412-C21', 'AG-449-D22') 2
('AG-347-C10', 'AG-347-I22', 'AG-670-M15') 3
('AG-347-E03', 'AG-347-K02', 'AG-424-P23') 2
('AG-347-G20', 'AG-402-A04', 'AG-449-J16') 2
('AG-347-J06', 'AG-355-N18', 'AG-402-K22') 2
('AG-347-O22', 'AG-418-P06', 'AG-449-C14') 8
('AG-355-N18', 'AG-355-P15', 'AG-418-B17') 3
('AG-412-C21', 'AG-418-F16', 'AG-436-E22') 2



#################


c1 60 80 after qcov correction (previously, had stored length in the qcov field)

('B241_526B17', 'B241_529O19') 26 (not in close clique) 
('B241_527L16', 'B243_495N16') 43 (cl1)
('B241_527P5', 'B241_528K19') 33 (cl2)
('B241_529J11', 'B245a_518E10') 40 (cl2)
('B243_496N4', 'B243_498L10') 20 (not in clique)
('B243_498F21', 'B245a_521O20') 19 (o20 in cl1, f21 not in clique)
1025 doubletons total 1025


('B241_526B17', 'B241_527L15', 'B241_529J16') 2 (l15 in cl5)
('B241_526B22', 'B241_528N8', 'B245a_520B18') 3 (b22 in cl6, n8 in cl4)
('B241_526B22', 'B243_498F21', 'B245a_520B18') 6 (b22 in cl6)
('B241_527G5', 'B241_527N11', 'B241_528N8') 2 (n8 in cl4)
('B241_527L15', 'B243_496E10', 'B243_498J20') 2 (e10 in cl2, l15 in cl5)
('B241_527L15', 'B243_498P3', 'B245a_518K17') 3 (p3 and l15 in cl5)
('B241_527L15', 'B245a_518K17', 'B245a_521M10') 2
('B241_527L16', 'B241_528N20', 'B243_496E10') 2
('B241_527P5', 'B241_528K19', 'B245a_521K15') 2
('B241_527P5', 'B243_496E10', 'B243_498L10') 3
('B241_528N8', 'B241_529J16', 'B243_495N16') 2
('B241_529C4', 'B243_495K23', 'B245a_520E22') 2
('B243_495K23', 'B243_496E10', 'B245a_520E22') 2
('B243_496N4', 'B243_498F21', 'B243_498L10') 2

none look to be from same clique


below: with proper qcov , 80,80
('B241_526B17', 'B241_529O19') 25
('B241_527L16', 'B243_495N16') 45
('B241_527P5', 'B241_528K19') 28
('B241_529J11', 'B245a_518E10') 40
('B243_496N4', 'B243_498L10') 18
('B243_498F21', 'B245a_521O20') 18
890 doubletons total 890
('B241_526B22', 'B243_498F21', 'B245a_520B18') 6
('B241_526D20', 'B241_527G5', 'B243_498M14') 2
('B241_527G5', 'B241_527N11', 'B241_528N8') 2
('B241_527L15', 'B243_498P3', 'B245a_518K17') 2
('B241_527L15', 'B245a_518K17', 'B245a_521M10') 3
('B241_527P5', 'B241_528K19', 'B245a_521K15') 2
('B241_528N8', 'B243_498G3', 'B243_498L10') 2
('B241_529C4', 'B243_495K23', 'B245a_520E22') 2
('B243_496N4', 'B243_498L10', 'B245a_519L21') 2
('B243_496N4', 'B243_498F21', 'B243_498L10') 2


'''
