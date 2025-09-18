import numpy as np
import ast
import csv
import math
from scipy import stats
import networkx as nx
import random
from os.path import exists
import statistics
import argparse
import matplotlib.pyplot as plt
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list
import Repo_site_catalogues_single_pair as repo

'''
This script writes an npz file recording cell SNP statistics for pairs of SNPs within the same gene, for the cell list in question, for core genes in the file CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt.  This information is used for linkage statistics in other programs.  twofold degenerate, fourfold degenerate, all snps, or second letter SNPs can be used

'''


    
def compile_Pairs_SNPs_to_outfile(outfile,celllist,downsamplesize,tfna,samplesize,coveragecut=0.02,mediancut=0):#tfna twofold fourfold nonsyn (second codon), all
    
    if mediancut:
        outfile=outfile[:-4]+'_MEDIAN_LENGTH_CUT.npz'
    print('Downsampling to ',downsamplesize)
    if tfna not in ['t','a','f','n']:
        print('Twofold Fourfold All SNPs or Nonsyn second codons?')
        return
    covcut=len(celllist)*coveragecut
    gene=[]
    site1=[]
    site2=[]
    naa=[]
    nat=[]
    nac=[]
    nag=[]
    ntt=[]
    ntc,ntg,ncc,ncg,ngg,nta,nca,nga,ngc,nct,ngt=[],[],[],[],[],[],[],[],[],[],[]
    myns=[naa,nat,nac,nag,nta,ntt,ntc,ntg,nca,nct,ncc,ncg,nga,ngt,ngc,ngg]
    
    sample_genenums=repo.make_sample_all_single_genes()#1,1 is length(number of genes), clusterbool
    if samplesize==0:
        samplesize=len(sample_genenums)
  
    with open('CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]  


    
    for i in range(len(mycores)):
       
        mygene=mycores[i]
        print('Gene %i' %mygene)
        add_stats_to_list_pair(mygene,gene,site1,site2,myns,covcut,downsamplesize,tfna,mediancut)
    repo.write_to_file_array_site_pairs(outfile,celllist,gene,site1,site2,myns)
    outcsv=outfile[:-4]+'.csv'
    repo.write_to_csv_site_pairs(outcsv,celllist,gene,site1,site2,myns)
    
        
def get_length(ntdict):
    cells=list(ntdict.keys())
    lens=[]
    for item in cells:
        l=len(ntdict[item])-ntdict[item].count('-')
        lens.append(l)
   
    #print('average len %f ' %(float(sum(lens))/len(lens)))
    return float(sum(lens))/len(lens)

def add_stats_to_list_pair(mygene,gene,site1,site2,myns,covcut,downsamplesize,tfna,mediancut):
    infile='Aligned_Gene_Groups/Gene%i_Kashtan_and_HLII.fasta' %mygene
    ntdict,positionlist=repo.load_single_gene_sequences(infile,celllist,tfna)
    if len(list(ntdict.keys()))==0:
        return
    #median length cut
    if mediancut:
        mylen=get_length(ntdict)#len(ntdict[list(ntdict.keys())[0]])
        if tfna=='t':
            medlength=45
        if tfna=='f':
            medlength=40
        if tfna=='tf':
            medlength=''
        if tfna=='watcg':
            medlength=283
        if tfna=='a':
            medlength=823
        print(mylen,medlength)
        
        if mylen<medlength:
            return
      
    mycells=list(ntdict.keys())
    if len(mycells)<covcut:
        return
    get_stats_pair(ntdict,mygene,gene,site1,site2,myns,covcut,downsamplesize,positionlist)
    

def get_stats_pair(ntdict,mygene,gene,site1,site2,myns,covcut,downsamplesize,positionlist):
    celllist=list(ntdict.keys())
    mylen=len(ntdict[celllist[0]])
    allele_possibilities=[]
    nucleotides=['A','T','C','G']
    for i in range(4):
        for j in range(4):
            allele_possibilities.append([nucleotides[i],nucleotides[j]])
            
    for i in range(len(positionlist)):
        pi=positionlist[i]
        for j in range(i):
            pj=positionlist[j]
            pair_alleles=[]
            for item in celllist:
                basei=ntdict[item][pi]
                basej=ntdict[item][pj]
                if basei == '-' or basej=='-':
                    continue
                pair_alleles.append([basei,basej])
            if len(pair_alleles)==0:
                continue
            
            nmc=len(pair_alleles)
            if downsamplesize=='':
                gene.append(mygene)
                site1.append(pi)
                site2.append(pj)
                for k in range(len(allele_possibilities)):
                    myns[k].append(pair_alleles.count(allele_possibilities[k]))
                continue
            if downsamplesize>0:
                if nmc>=downsamplesize:
                    gene.append(mygene)
                    site1.append(pi)
                    site2.append(pj)
                    #downsample, then append
                    mini_myns=[pair_alleles.count(allele_possibilities[k]) for k in range(len(allele_possibilities))]#check this
                    if nmc!=sum(mini_myns):
                        print('Problem with mutual coverage and cell pairs: mutual cov does not equal the sum over the found allele pairs %i %i' %(nmc,sum(mini_myns)))
                        t=input('t')
                    new_minimyns=repo.downsample_site_pair(mini_myns,nmc,downsamplesize,verbosebool=0)
                    for m in range(len(new_minimyns)):
                        myns[m].append(new_minimyns[m])
                        
                    

     
    

    

    
    





#######################################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB')
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    parser.add_argument('-ds','--dsbool',type=int,default=0,required=False)#to downsample or not
    parser.add_argument('-p','--prog',type=int,default=0,required=False)
    parser.add_argument('-ss','--samplesize',default=0,type=int)#number of genes to do
    args = parser.parse_args()
    print(args)
    
    n=args.n
    prog=args.prog
    
    if n==0:
        groupname='HLII'
        celllist=hliilist2
        downsamplesize=80#120
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype
        downsamplesize=math.floor(0.5*len(celllist))
    if n==2:
        groupname='BB_NoClose'
        celllist=BB_list
        downsamplesize=math.floor(0.5*len(celllist))

    if n==3:
        celllist=clique1_without_ribotype
        downsamplesize=7
        groupname='Clique1'

    samplesize=args.samplesize
    tfna=args.tfna
    suffix=''
    if tfna=='t':
        outfile='%s_TWOFOLDDEGENERATE_Intragene_Site_Pairs_NEWDS_%s.npz' %(groupname,suffix)
    if tfna=='f':
        outfile='%s_FOURFOLDDEGENERATE_Intragene_Site_Pairs_NEWDS_%s.npz' %(groupname,suffix)
    if tfna=='n':
        outfile='%s_NONSYN_SECONDNT_Intragene_Site_Pairs_NEWDS.npz' %(groupname)
    if tfna=='a':
        outfile='%s_All_Intragene_Site_Pairs_NEWDS.npz' %(groupname)
    if samplesize==0:
        outfile='Intragene_SNP_Pair_Catalogs/'+outfile
    
    if samplesize!=0:
        outfile='Intragene_SNP_Pair_Catalogs/Testing_Subsets_of_Genes/'+outfile
        outfile=outfile[:-4]+'_%i-Gene-Sample.npz' %samplesize
    outfile_dss=outfile[:-4]+'_DownSample_%iCells.npz' %downsamplesize
    
    if prog==0:
        if not args.dsbool:
            downsamplesize=''
        
        if downsamplesize!='':
            outfile=outfile_dss
        print('outfile %s' %outfile)
        compile_Pairs_SNPs_to_outfile(outfile,celllist,downsamplesize,tfna,samplesize,mediancut=1)
   
   
