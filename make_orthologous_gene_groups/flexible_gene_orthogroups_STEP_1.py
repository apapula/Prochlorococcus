import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac
import scipy.spatial.distance as ssd
import networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy import stats
import statistics
import math
from matplotlib import colors
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,AllKashBATS_without_ribotype, BB_closegrp_1,HLII_NoClose
import random
import argparse
import Get_R_From_Cell_Pair_Matrices as rrepo
import ast
import re
import Repo_site_catalogues_single_pair as repo
from os.path import exists


'''
This file takes in an all-v-all blast of non-core ORFs from a given cell group.  With input pident and qcovhsp cuts, it creates a graph with edges connecting ORFs that meet these cuts and takes connected components of this graph to serve as flex gene orthogroups

It then filters out duplicate entries: a given sag having multiple orfs per orthogroup

these orthogroups compare well with those yielded by IMG using gff or pfaff annotations, but goes far beyond, including all orfs instead of just annotated

Then, it creates files listing the name (identification number) of ORFs (1 file for each orthogroups).  These are used in blastdbcmd in step 2 file to create fasta files for each flex orthogroup

This file contains step A (p=0) and step B (p=1)

Internal note: made from C1_Flexible_Genes + Flexible_gene_divergences
'''
###########################################
#Step A: get connected components from all vs all blast (prog=0)
#########################################################

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
                continue
            if '245a_518D8' in sseq or '245a_518D8' in query:
                continue#this cell is erroneous
            edges.append((query,sseq))
    edges=prune_edges_to_bidirectional(edges)   
    
    orthogroups=get_graph_connected_comps(edges,C1bool,groupname,qcut,pidcut)

    orthogroups=filter_orthogroups(orthogroups,C1bool,groupname,celllist,qcut,pidcut)
    
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
    print('done pruning')
    return newedges

def filter_orthogroups(orthogroups,C1bool,groupname,celllist,qcut,pidcut):
    #take out self-hits; repeated cells per orthogroup.
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
    return newl
    


def get_graph_connected_comps(edges,C1bool,groupname,qcut,pidcut):
    G=nx.Graph()
    G.add_edges_from(edges)

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
        if cell not in cells:
            cells.append(cell)

###########################################################
#Step B  (prog=1) makes files naming all ORFs per orthogroup, to be used in step 2 file for blastdbcmd from a database
#########################################################

def make_files_for_orthogroups(celllist,groupname,sizecut,qcut,pidcut):
    group_cells,group_ids=load_ortho_group_dict_from_clustering(celllist,groupname,sizecut,qcut,pidcut)
    filenamelist=[]
    #for each item in group_ids, title the file FlexibleGene_(item+2000)_ORFnames.txt. add this name to list. then write the orfs in the file
    for item in group_ids:
        if groupname!='C1':
            mytitle='FlexibleGene_%i_ORFnames.txt' %(int(item)+2000)
        if groupname=='C1':
            mytitle='FlexibleGene_%i_ORFnames.txt' %(int(item)+3000)
        filenamelist.append(mytitle)
        with open(mytitle,'w') as outf:
            for orf in group_ids[item]:
                outf.write('%s\n' %orf)
    with open('Flexible_Gene_File_Name_List.txt','w') as outf:
        for i in range(len(filenamelist)):
            outf.write('%s\n' %(filenamelist[i]))

def load_ortho_group_dict_from_clustering(celllist,groupname,sizecut,qcut,pidcut):
    f='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    groups_cells,groups_ids={},{}

    lower_cut=4
    lineindex=0
    for line in lines:
        group=ast.literal_eval(line)
        if len(group)<lower_cut:
            continue    
        if len(group)>sizecut:
            continue
        newgroup=[]
        newid=[]
        for c in group:
            idx=c.find('_')
            if groupname=='C1':
                mynewid=c[:idx]
                cell=c[idx+1:]
            if groupname!='C1':
                cell=c[:idx]
                mynewid=c[idx+1:]
            if cell not in celllist:
                print(cell)
                print(celllist)
                t=input('t')
                continue
            newgroup.append(cell)
            newid.append(mynewid)
        groups_cells[lineindex]=newgroup
        groups_ids[lineindex]=newid
        lineindex+=1
    return groups_cells,groups_ids

######################################




    
##########################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-n','--n',type=int)
    parser.add_argument('-pidcut','--pident',type=int,default=40)
    parser.add_argument('-qcut','--qcut',type=int,default=40)
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    n=args.n
    tfna=args.tfna
    pidcut=args.pident
    qcut=args.qcut

    c1bool=0
    if n==0:
        groupname='HLII'
        celllist=hliilist2
        sizecut=34
    if n==2:
        groupname='Blue_Basin_NoClose'
        celllist=BB_list
        sizecut=14
   
    if n==1:
        celllist=C1_without_ribotype
        sizecut=14
        groupname='C1'
        c1bool=1
        
    if prog==0:
        get_connected_components_from_all_v_all_BLAST(qcut,pidcut,groupname,celllist,C1bool=c1bool)

    if prog==1:
       
        make_files_for_orthogroups(celllist,groupname,sizecut,qcut,pidcut)
      
   
