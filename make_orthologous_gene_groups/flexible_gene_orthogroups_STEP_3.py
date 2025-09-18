import numpy as np
import math
import ast
import statistics
from os.path import exists
import statistics
import argparse
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,AllKashBATS_without_ribotype, BB_closegrp_1,C3
import matplotlib.pyplot as plt
import Repo_site_catalogues_single_pair as repo
import random

'''
This file creates distance matrices for each flex gene orthogroup, for C1 (p=1) or all HLII (p=0)

As per usual, users will need to adjust paths
'''


def write_matrices(tfna,outfile,cellsubset,groupname):
    genearrays={}#gene label:array of matrix
    genearraysrand={}
    sample_genenums=list(range(731))
    if groupname=='C1':
        sample_genenums=list(range(228))
    mylist=overall_order_list()
    if cellsubset==[]:
        cellsubset=mylist
    for i in range(len(sample_genenums)):
        mygene=sample_genenums[i]+2000
        if groupname=='C1':
            mygene+=1000
        print('Gene ',mygene)
        do_single_gene(genearrays,mylist,tfna,mygene,genearraysrand,cellsubset,groupname)
    #save to array
    save_to_array(outfile,genearrays)
        
def save_to_array(outfile,genearrays):
    np.save(outfile,genearrays)
    
def do_single_gene(genearrays,overalllist,tfna,mygene,genearraysrand,cellsubset,groupname):
    dim=len(overalllist)
    initialize_arrays(dim,genearrays,genearraysrand,tfna,mygene)
    if tfna!='aa':
        infile='Flexible_Orthogroups_Aligned/FlexibleGene%i.fasta' %mygene
    if tfna=='aa':
        infile='Flexible_Orthogroups_Aligned/FlexibleGene%i.faa' %mygene
    if groupname=='C1':
        infile='C1_'+infile
   
    for i in range(len(overalllist)):
        celli=overalllist[i]
        if celli not in cellsubset:
            continue
        for j in range(i):
            cellj=overalllist[j]
            if cellj not in cellsubset:
                continue
            do_single_cell_pair(genearrays,tfna,mygene,genearraysrand,celli,cellj,i,j,infile)


        
            
        
def do_single_cell_pair(genearrays,tfna,mygene,genearraysrand,cell1,cell2,idx1,idx2,infile):
    mylist=[cell1,cell2]
    #print(mylist)
    if tfna in ['t','f','n','a','atvscg','aa']:
        tfna2=tfna
        if tfna=='aa':
            tfna2='aminoacid'
        ntdict,positionlist=repo.load_single_gene_sequences(infile,mylist,tfna2,gapsbool=1)
    if tfna=='tf':#two and fourfold combined
        ntdictt,positionlistt=repo.load_single_gene_sequences(infile,mylist,'t',gapsbool=1)
        ntdictf,positionlistf=repo.load_single_gene_sequences(infile,mylist,'f',gapsbool=1)
        if len(list(ntdictf.keys()))<2:
            return
        ntdict={}
        for item in ntdictt:
            ntdict[item]=ntdictt[item]+ntdictf[item]
        
    if len(list(ntdict.keys()))!=2:
        return
    #print(list(ntdict.keys()))
    #done loading seqs for those two cells. now, break into fractions. if frac>1, do random as well.
    ndiff,mymc=get_pair_divergence(ntdict,cell1,cell2)
    #print(mydiv)
    if mymc<10:
        return
    genearrays['Gene_%i_%s_Pairs' %(mygene,tfna)][idx1][idx2]=[ndiff,mymc]
    genearrays['Gene_%i_%s_Pairs' %(mygene,tfna)][idx2][idx1]=[ndiff,mymc]
       
def get_only_syn_sites(ntdict):
    newdict={}#weed out gaps
    s1=''
    s2=''
    cells=list(ntdict.keys())
    seq1=ntdict[cells[0]]
    seq2=ntdict[cells[1]]
    for i in range(len(seq1)):
        if seq1[i] in ['A','T','C','G'] and seq2[i] in ['A','T','C','G']:
            s1+=seq1[i]
            s2+=seq2[i]
    newdict[cells[0]]=s1
    newdict[cells[1]]=s2
    return newdict

def get_pair_divergence(ntdict,cell1,cell2):
    #ntdict can be random or not, fraction or not.
    seqi=ntdict[cell1]
    seqj=ntdict[cell2]
    ndiff=sum(1 for a, b in zip(seqi,seqj) if a != b and a!='-' and b!='-' and a!='N' and b!='N')
    nmc=sum(1 for a, b in zip(seqi,seqj) if a!='-' and b!='-' and a!='N' and b!='N')
    if nmc<10:#==0:
        return 0,0
    return ndiff,nmc

    
   


    
def initialize_arrays(dim,genearrays,genearraysrand,tfna,mygene):
    mat=np.array([[[0,0] for i in range(dim)] for j in range(dim)])
    
    genearrays['Gene_%i_%s_Pairs' %(mygene,tfna)]=mat
    
   
            
def overall_order_list():
    mylist=[]
    for i in range(len(hliilist2)):
        mylist.append(hliilist2[i])
    for i in range(len(AllKashBATS_without_ribotype)):
        mylist.append(AllKashBATS_without_ribotype[i])
    return mylist


if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB',required=False,default=-1)
   
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
   
   
    n=args.n
    tfna=args.tfna
    suffix='_Flexible2'
    
    
    groupnamet='Berube_Kash'
    if n==0:
        groupname='HLII'
        celllist=hliilist2
       
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype
    if n==2:
        groupname='BB_NoClose'
        celllist=BB_list#HLII_BB_NoCloseCells_Median04
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
    if n==3:
        celllist=clique1_without_ribotype
        groupname='Clique1'
    if n==4:
        celllist=BB_closegrp_1
        groupname='Closegrp1'
    if n==-3:
        groupname='C3'
        celllist=C3
        celllist=[celllist[i][:-3] for i in range(len(celllist))]
        
    if prog==0 and n!=-1:
        groupnamet=groupname#'HLII'
   
    if tfna=='t':
        outfile='Matrices_%s_TWOFOLDDEGENERATE%s' %(groupnamet,suffix)
    if tfna=='tf':
        
        outfile='Matrices_%s_TWO_AND_FOUR%s' %(groupnamet,suffix)
    if tfna=='aa':
        outfile='Matrices_Berube_Kash_AA_Discrete_Including_Non_Core'
    if tfna=='f':
        outfile='Matrices_%s_FOURFOLDDEGENERATE%s' %(groupnamet,suffix)
    if tfna=='n':
        outfile='Matrices_%s_NONSYN_SECONDNT%s' %(groupnamet,suffix)
    if tfna=='a':
        outfile='Matrices_%s_All_NT%s' %(groupnamet,suffix)
   
    outfile='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfile
    

    outfilewhole=outfile
   
    if tfna in ['t','f','a','aa']:
        outfilewhole+='.npy'
        infilematrices=outfilewhole
    if tfna=='tf':
        infilematrices=['%s_firsthalf.npy' %outfilewhole,'%s_secondhalf.npy' %outfilewhole]

    

    outfile='../input_data/%s' %outfile


    if prog==0:
        if n==-1:
            cellsubset=[]
        if n!=-1:
            cellsubset=celllist
        print(outfile)
        groupname='HLII'
        write_matrices(tfna,outfile,cellsubset,groupname)
        
    if prog==1:
      groupname='C1'
      outfile=outfile+'_C1'
      cellsubset=C1_without_ribotype
      write_matrices(tfna,outfile,cellsubset,groupname)
