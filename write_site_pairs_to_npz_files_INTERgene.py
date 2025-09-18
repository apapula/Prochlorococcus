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
This script writes an npz file recording cell SNP statistics for pairs of SNPs within two different genes, for the cell list in question, for core genes in the file CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt.  This information is used for linkage statistics in other programs.  twofold degenerate, fourfold degenerate, all snps, or second letter SNPs can be used
'''




def compile_Pairs_SNPs_to_outfile(outfile,celllist,downsamplesize,tfna,samplesize,gene_sep_bins,coveragecut=0.02):#tfna twofold fourfold nonsyn (second codon), all
    #if downsamplesize=='', don't downsample.
    print('Downsampling to ',downsamplesize)
    if tfna not in ['t','a','f','n']:
        print('Twofold Fourfold All SNPs or Nonsyn second codons?')
        return
    covcut=len(celllist)*coveragecut
    gene1=[]
    gene2=[]
    site1=[]
    site2=[]
    naa=[]
    nat=[]
    nac=[]
    nag=[]
    ntt=[]
    ntc,ntg,ncc,ncg,ngg,nta,nca,nga,ngc,nct,ngt=[],[],[],[],[],[],[],[],[],[],[]
    myns=[naa,nat,nac,nag,nta,ntt,ntc,ntg,nca,nct,ncc,ncg,nga,ngt,ngc,ngg]

    ##EDIT
    for i in range(len(gene_sep_bins)-1):
        print(gene_sep_bins[i])
        do_sep_bin(gene_sep_bins,i,gene1,gene2,site1,site2,myns,samplesize,covcut)
    
    repo.write_to_file_array_site_pairs_intergene(outfile,celllist,gene1,gene2,site1,site2,myns)
    outcsv=outfile[:-4]+'.csv'
    repo.write_to_csv_site_pairs_intergene(outcsv,celllist,gene1,gene2,site1,site2,myns)
    
        

def do_sep_bin(gene_sep_bins,iiii,gene1,gene2,site1,site2,myns,samplesize,covcut):
    seplow=gene_sep_bins[iiii]
    sephigh=gene_sep_bins[iiii+1]
    genenums1,genenums2=repo.get_longrange_gene_sample(seplow,sephigh,samplesize)#1,1 is length(number of genes), clusterbool
   
    for i in range(min(samplesize,len(genenums1))):#len(sample_genenums)):
        if i%10==0:
            print('%f (%i of %i)' %(float(i)/(samplesize),i,(samplesize)))
        mygene1=genenums1[i]
        mygene2=genenums2[i]
        print('Gene %i %i' %(mygene1,mygene2))
        add_stats_to_list_pair_intergene(mygene1,mygene2,gene1,gene2,site1,site2,myns,covcut,downsamplesize,tfna)

        
def add_stats_to_list_pair_intergene(mygene1,mygene2,gene1,gene2,site1,site2,myns,covcut,downsamplesize,tfna):
    infile1='Aligned_Gene_Groups/Gene%i_Kashtan_and_HLII.fasta' %mygene1
    ntdict1,positionlist1=repo.load_single_gene_sequences(infile1,celllist,tfna)
    if len(list(ntdict1.keys()))==0:
        return
   
    mycells1=list(ntdict1.keys())
    if len(mycells1)<covcut:
        return
    infile2='Aligned_Gene_Groups/Gene%i_Kashtan_and_HLII.fasta' %mygene2
    ntdict2,positionlist2=repo.load_single_gene_sequences(infile2,celllist,tfna)
    if len(list(ntdict2.keys()))==0:
        return
   
    mycells2=list(ntdict2.keys())
    if len(mycells2)<covcut:
        return
    ntdict1,ntdict2=get_mc_cells(ntdict1,ntdict2)
    if len(list(ntdict1.keys()))<covcut:
        return#new 9-1-23
    get_stats_pair_intergene(ntdict1,ntdict2,mygene1,mygene2,gene1,gene2,site1,site2,myns,covcut,downsamplesize,positionlist1,positionlist2)
    
#EDIT BELOW


def get_stats_pair_intergene(ntdict1,ntdict2,mygene1,mygene2,gene1,gene2,site1,site2,myns,covcut,downsamplesize,positionlist1,positionlist2):
    celllist=list(ntdict1.keys())#should be the same for ntdict 1 and 2 due to get mutucal cov
    mylen=len(ntdict1[celllist[0]])
    allele_possibilities=[]
    nucleotides=['A','T','C','G']
    for i in range(4):
        for j in range(4):
            allele_possibilities.append([nucleotides[i],nucleotides[j]])
            
    for i in range(len(positionlist1)):
        pi=positionlist1[i]
        for j in range(len(positionlist2)):
            pj=positionlist2[j]
            pair_alleles=[]
            for item in celllist:
                basei=ntdict1[item][pi]
                basej=ntdict2[item][pj]
                if basei == '-' or basej=='-':
                    continue
                pair_alleles.append([basei,basej])
            if len(pair_alleles)==0:
                continue
            
            nmc=len(pair_alleles)
            if downsamplesize=='':
                gene1.append(mygene1)
                gene2.append(mygene2)
                site1.append(pi)
                site2.append(pj)
                for k in range(len(allele_possibilities)):
                    myns[k].append(pair_alleles.count(allele_possibilities[k]))
                continue
            if downsamplesize>0:
                if nmc>=downsamplesize:
                    gene1.append(mygene1)
                    gene2.append(mygene2)
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
                        
                    

     
    

    
def get_mc_cells(ntdict1,ntdict2):
    #only want cells covered at both genes for site pairs
    new1,new2={},{}
    for item in ntdict1:
        if item in ntdict2:
            new1[item]=ntdict1[item]
            new2[item]=ntdict2[item]
    return new1,new2
    
    





#######################################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB')
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    parser.add_argument('-ds','--dsbool',type=int,default=1,required=False)#to downsample or not
    parser.add_argument('-p','--prog',type=int,default=0,required=False)
    parser.add_argument('-ss','--samplesize',default=5000,type=int)#number of genes to do
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
        outfile='%s_TWOFOLDDEGENERATE_INTERgene_Site_Pairs_NEWDS_%s.npz' %(groupname,suffix)
    if tfna=='f':
        outfile='%s_FOURFOLDDEGENERATE_INTERgene_Site_Pairs_NEWDS_%s.npz' %(groupname,suffix)
    if tfna=='n':
        outfile='%s_NONSYN_SECONDNT_INTERgene_Site_Pairs_NEWDS.npz' %(groupname)
    if tfna=='a':
        outfile='%s_All_INTERgene_Site_Pairs_NEWDS.npz' %(groupname)
    if samplesize==5000:
        outfile='Different_Gene_SNP_Pair_Catalogs/'+outfile
    
    if samplesize!=5000:
        outfile='Different_Gene_SNP_Pair_Catalogs/Testing_Subsets_of_Genes/'+outfile
        outfile=outfile[:-4]+'_%i-Gene-Sample.npz' %samplesize
    outfile_dss=outfile[:-4]+'_DownSample_%iCells___.npz' %downsamplesize
    
    if prog==0:
        if not args.dsbool:
            downsamplesize=''
        
        if downsamplesize!='':
            outfile=outfile_dss
        print('outfile %s' %outfile)
        gene_sep_bins=[1,2,3,4,5,6,7,8,9,10,25,63,158,398,1000]
        compile_Pairs_SNPs_to_outfile(outfile,celllist,downsamplesize,tfna,samplesize,gene_sep_bins)
   
