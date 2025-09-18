

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
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,HLII_NoClose
import Repo_site_catalogues_single_pair as repo


'''
This program creates a catalog of SNPs within a given cell group, for twofold degenerate, fourfold degenerate, all sites, or second codon sites.  Downsampling can be done according to average cell coverage
'''


def compile_Single_SNPs_to_outfile(outfile,celllist,downsamplesize,tfna,coveragecut=0.02):#tfna twofold fourfold nonsyn (second codon), all
    #THIS JUST USES CORE GENES
    if tfna not in ['t','a','f','n']:
        print('Twofold Fourfold All SNPs or Nonsyn second codons?')
        return
    covcut=len(celllist)*coveragecut
    gene=[]
    site=[]
    na=[]
    nt=[]
    nc=[]
    ng=[]
    #and downsample

    gene_ds=[]
    site_ds=[]
    na_ds=[]
    nc_ds=[]
    nt_ds=[]
    ng_ds=[]
    
    sample_genenums=repo.make_sample_all_single_genes()#1,1 is length(number of genes), clusterbool
    for i in range(len(sample_genenums)):
        if i%10==0:
            print('%f (%i of %i)' %(float(i)/len(sample_genenums),i,len(sample_genenums)))
        mygene=sample_genenums[i][0]
        add_stats_to_list(mygene,gene,site,na,nt,ng,nc,covcut,downsamplesize,site_ds,gene_ds,na_ds,nt_ds,nc_ds,ng_ds,tfna)
    repo.write_to_file_array_single_site(outfile,celllist,gene,site,na,nt,nc,ng)
    outcsv=outfile[:-4]+'.csv'
    repo.write_to_csv_single_site(outcsv,celllist,gene,site,na,nt,nc,ng)
    if downsamplesize!='':
        outfile_ds=outfile[:-4]+'_DownSample_%iCells.npz' %downsamplesize
        repo.write_to_file_array_single_site(outfile_ds,celllist,gene_ds,site_ds,na_ds,nt_ds,nc_ds,ng_ds)
        outcsv=outfile_ds[:-4]+'.csv'
        repo.write_to_csv_single_site(outcsv,celllist,gene_ds,site_ds,na_ds,nt_ds,nc_ds,ng_ds)
        


def add_stats_to_list(mygene,gene,site,na,nt,ng,nc,covcut,downsamplesize,site_ds,gene_ds,na_ds,nt_ds,nc_ds,ng_ds,tfna):
    infile='Aligned_Gene_Groups/Gene%i_Kashtan_and_HLII.fasta' %mygene
    ntdict,positionlist=repo.load_single_gene_sequences(infile,celllist,tfna)
    if len(list(ntdict.keys()))==0:
        return
   
    mycells=list(ntdict.keys())
    if len(mycells)<covcut:
        return
    get_stats(ntdict,site,gene,na,nt,nc,ng,mygene,downsamplesize,gene_ds,site_ds,na_ds,nt_ds,nc_ds,ng_ds)
    

def get_stats(ntdict,site,gene,na,nt,nc,ng,mygene,downsamplesize,gene_ds,site_ds,na_ds,nt_ds,nc_ds,ng_ds):
    celllist=list(ntdict.keys())
    mylen=len(ntdict[celllist[0]])
    
    for i in range(mylen):
        uniqsi={}#for alleles (other than '-'), number of cells
        #get stats for this size, then if coverage>=downsamplesize, downsample this site
        for c in celllist:
            allelec=ntdict[c][i]
            if allelec in uniqsi:
                uniqsi[allelec]+=1
            if allelec!='-' and allelec not in uniqsi:
                uniqsi[allelec]=1
        if len(uniqsi)==0:
            continue#new 7-11. necessary, was wrong before for fourfold and nonsyn, which included non fourfold and non-nonsyn sites as uncovered
        if 1:#don't want to stipulate polymorphic; want number of fixed sites len(uniqsi)>=2:
            gene.append(mygene)
            site.append(i)

            myna,mync,mynt,myng=0,0,0,0
            if 'A' in uniqsi:
                myna=uniqsi['A']
            if 'G' in uniqsi:
                myng=uniqsi['G']
            if 'C' in uniqsi:
                mync=uniqsi['C']
            if 'T' in uniqsi:
                mynt=uniqsi['T']
            na.append(myna)
            nt.append(mynt)
            nc.append(mync)
            ng.append(myng)
            nmc=myna+mynt+myng+mync
            if downsamplesize!='':
                if nmc>=downsamplesize:
                    newna,newnt,newnc,newng=repo.downsample_site(myna,mynt,mync,myng,nmc,downsamplesize)
                    gene_ds.append(mygene)
                    site_ds.append(i)
                    na_ds.append(newna)
                    nt_ds.append(newnt)
                    ng_ds.append(newng)
                    nc_ds.append(newnc)




#######################################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB')
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    parser.add_argument('-ds','--ds',default=1,type=int)
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    n=args.n
    ds=args.ds
        
    if n==0:
        groupname='HLII'
        celllist=hliilist2
        downsamplesize=120
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype
        downsamplesize=math.floor(0.7*len(celllist))
    if n==2:
        groupname='BB_NoClose'
        celllist=BB_list
        downsamplesize=48

    if n==3:
        celllist=clique1_without_ribotype
        downsamplesize=7
        groupname='Clique1'
    if n==5:
        celllist=HLII_NoClose
        groupname='HLII_NoClose'
        downsamplesize=math.floor(0.7*len(celllist))

   
    tfna=args.tfna
    suffix=''
    if tfna=='t':
        outfile='%s_TWOFOLDDEGENERATE_Single_Sites_%s.npz' %(groupname,suffix)
    if tfna=='f':
        outfile='%s_FOURFOLDDEGENERATE_Single_Sites_%s.npz' %(groupname,suffix)
    if tfna=='n':
        outfile='%s_NONSYN_SECONDNT_Single_Sites_.npz' %(groupname)
    if tfna=='a':
        outfile='%s_All_Single_Sites_.npz' %(groupname)
    outfile='Single_Site_Catalogs/%s' %outfile
    outfile_dss=outfile[:-4]+'_DownSample_%iCells.npz' %downsamplesize
        
    if prog==0:
        
        print('outfile %s' %outfile)
        compile_Single_SNPs_to_outfile(outfile,celllist,downsamplesize,tfna)
    
        
   
