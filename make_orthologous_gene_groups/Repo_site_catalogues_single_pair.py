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
import re
import Load_Sequences_Repo as load_sequences_repo

'''

na,nt,nc,ng,gene,site=load_single_sitefile(infile)
average_of_l=avelist(l)
na,nt,nc,ng=check_if_twofold_degen(na,nt,nc,ng,gene,site)
newna,newnt,newnc,newng=downsample_site(na,nt,nc,ng,downsamplesize
write_to_file_array_single_site(outfile,celllist,gene,site,na,nt,nc,ng)
sample_genenums=make_sample_all_single_genes_)
write_to_csv_single_site(outcsv,celllist,gene,site,na,nt,nc,ng)
seqdict=load_single_gene_sequences(infile,celllist,tfna,gapsbool=1,verbosebool=0)
write_to_csv_site_pairs(outcsv,celllist,gene,site1,site2,myns)naa,nat,nac,nag,ntt,ntc,ntg,ncc,ncg,ngg)
 write_to_file_array_site_pairs(outfile,celllist,gene,site1,site2,myns),naa,nat,nac,nag,ntt,ntc,ntg,ncc,ncg,ngg)
'''

def load_mult_arrays(arrayfilelist):
    genearrays={}
    for i in range(len(arrayfilelist)):
        npfile=np.load(arrayfilelist[i],allow_pickle=1)
        mykeys=npfile[()].keys()
        for item in mykeys:
            genearrays[item]=npfile[()][item]
    return genearrays

def get_array_name(tfna):
    prefix='../input_data/Cell_Pair_Matrices/Each_Cell_Pair_Sites/'
    if tfna=='a':
        return prefix+'Matrices_Berube_Kash_All_NT_Discrete_Including_Non_Core.npy'
    if tfna=='aa':
        return prefix+'Matrices_Berube_Kash_All_NT_Discrete.npy'
    if tfna=='n':
        return prefix+'Matrices_Berube_Kash_NONSYN_SECONDNT_Discrete_Including_Non_Core.npy'

def load_single_gene_sequences(infile,celllist,tfna,gapsbool=1,verbosebool=0):
    if gapsbool==0:
        print('Gaps must be taken into account')
        return
    seqdict,positionlist=load_sequences_repo.load_single_gene_sequences_total(infile,celllist,tfna,gapsbool,verbosebool)
    return seqdict,positionlist



def get_gene_from_key_to_check(key):
    idx1=key.find('_')
    idx2=key.find('_',idx1+1)
    g=int(key[idx1+1:idx2])
    return g

def return_tfna_string(tfna):
    if tfna=='n':
        tfna_string='Second Letter'
    if tfna=='t':
        tfna_string='Twofold'
    if tfna=='f':
        tfna_string='Fourfold'
    if tfna=='tf':
        tfna_string='Twofold+Fourfold'
    if tfna=='a':
        tfna_string='All SNPs'
    if tfna=='aa':
        tfna_string='Amino Acid'
    return tfna_string


def write_to_csv_single_site(outcsv,celllist,gene,site,na,nt,nc,ng):        
   outfields=['GeneNum','Site','A','T','C','G']
   with open(outcsv,'w') as outfile:
       outwriter=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
       outwriter.writeheader()
       for i in range(len(gene)):
           outwriter.writerow({'GeneNum':gene[i],'Site':site[i],'A':na[i],'T':nt[i],'C':nc[i],'G':ng[i]})

def write_to_csv_site_pairs(outcsv,celllist,gene,site1,site2,myns):        
   outfields=['GeneNum','Site1','Site2','AA','AT','AC','AG','TA','TT','TC','TG','CA','CT','CC','CG','GA','GT','GC','GG']
   print(len(myns),' should be 16')
   for i in range(len(myns)):
       print(len(myns[i]))
   print(len(gene),len(site1),len(site2))
   with open(outcsv,'w') as outfile:
       outwriter=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
       outwriter.writeheader()
       for i in range(len(gene)):
           outwriter.writerow({'GeneNum':gene[i],'Site1':site1[i],'Site2':site2[i],'AA':myns[0][i],'AT':myns[1][i],'AC':myns[2][i],'AG':myns[3][i],'TA':myns[4][i],'TT':myns[5][i],'TC':myns[6][i],'TG':myns[7][i],'CA':myns[8][i],'CT':myns[9][i],'CC':myns[10][i],'CG':myns[11][i],'GA':myns[12][i],'GT':myns[13][i],'GC':myns[14][i],'GG':myns[15][i]})



def write_to_csv_site_pairs_intergene(outcsv,celllist,gene1,gene2,site1,site2,myns):        
   outfields=['Gene1','Gene2','Site1','Site2','AA','AT','AC','AG','TA','TT','TC','TG','CA','CT','CC','CG','GA','GT','GC','GG']
   print(len(myns),' should be 16')
   for i in range(len(myns)):
       print(len(myns[i]))
   print(len(gene1),len(site1),len(site2))
   with open(outcsv,'w') as outfile:
       outwriter=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
       outwriter.writeheader()
       for i in range(len(gene1)):
           outwriter.writerow({'Gene1':gene1[i],'Gene2':gene2[i],'Site1':site1[i],'Site2':site2[i],'AA':myns[0][i],'AT':myns[1][i],'AC':myns[2][i],'AG':myns[3][i],'TA':myns[4][i],'TT':myns[5][i],'TC':myns[6][i],'TG':myns[7][i],'CA':myns[8][i],'CT':myns[9][i],'CC':myns[10][i],'CG':myns[11][i],'GA':myns[12][i],'GT':myns[13][i],'GC':myns[14][i],'GG':myns[15][i]})


def write_to_csv_site_pairs_old(outcsv,celllist,gene,site1,site2,naa,nat,nac,nag,ntt,ntc,ntg,ncc,ncg,ngg):        
   outfields=['GeneNum','Site1','Site2','AA','AT','AC','AG','TT','TC','TG','CC','GC','GG']
   with open(outcsv,'w') as outfile:
       outwriter=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
       outwriter.writeheader()
       for i in range(len(gene)):
           outwriter.writerow({'GeneNum':gene[i],'Site1':site1[i],'AA':naa[i],'AT':nat[i],'AC':nac[i],'AG':nag[i],'TT':ntt[i],'TC':ntc[i],'TG':ntg[i],'CG':ncg[i],'GG':ngg[i]})

def make_sample_all_single_genes():
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]        
        
    with open('HLIIRand_CoveredGeneNumbers.txt','r') as coveredinfile:
        coveredgenes=[int(line.rstrip()) for line in coveredinfile if len(line)>0]            
    samplestarts=[]#starting genes for each sample
    sample_genenums=[]
    for mygene in coveredgenes:
        if mygene not in mycores:
            continue
        mylist=[mygene]
        coveredbool=True
        for item in mylist:
            if item not in coveredgenes:
                coveredbool=False
                break
        if coveredbool==False:
            continue
        if not all(mylist) in mycores:
            continue
      
        sample_genenums.append(mylist)
    return sample_genenums

def make_sample_all_genes_including_noncore():
    with open('HLIIRand_CoveredGeneNumbers.txt','r') as coveredinfile:
        coveredgenes=[int(line.rstrip()) for line in coveredinfile if len(line)>0]            
    samplestarts=[]#starting genes for each sample
    sample_genenums=[]
    for mygene in coveredgenes:
        mylist=[mygene]
        coveredbool=True
        for item in mylist:
            if item not in coveredgenes:
                coveredbool=False
                break
        if coveredbool==False:
            continue      
        sample_genenums.append(mylist)
    return sample_genenums

def load_INTERgene_pair_file(infile):
    npfile=np.load(infile)
    gene1=npfile['gene1']
    site1=npfile['site1']
    gene2=npfile['gene2']
    site2=npfile['site2']
    myns=[]
    mybases=['a','t','c','g']
    for i in range(len(mybases)):
        for j in range(len(mybases)):
            mykey='n%s%s' %(mybases[i],mybases[j])
            l=npfile[mykey]
            myns.append(l)
            #print(mykey)
    return gene1,site1,gene2,site2,myns

def load_pair_file(infile):
    npfile=np.load(infile)
    gene=npfile['gene']
    site1=npfile['site1']
    site2=npfile['site2']
    myns=[]
    mybases=['a','t','c','g']
    for i in range(len(mybases)):
        for j in range(len(mybases)):
            mykey='n%s%s' %(mybases[i],mybases[j])
            l=npfile[mykey]
            myns.append(l)
            #print(mykey)
    return gene,site1,site2,myns

    
def write_to_file_array_single_site(outfile,celllist,gene,site,na,nt,nc,ng):
    gene=np.array(gene)
    site=np.array(site)
    na=np.array(na)
    nt=np.array(nt)
    nc=np.array(nc)
    ng=np.array(ng)
    np.savez(outfile,gene=gene,site=site,na=na,nt=nt,nc=nc,ng=ng)

def write_to_file_array_site_pairs(outfile,celllist,gene,site1,site2,myns):
    gene=np.array(gene)
    site1=np.array(site1)
    site2=np.array(site2)
    
    np.savez(outfile,gene=gene,site1=site1,site2=site2,naa=np.array(myns[0]),nat=np.array(myns[1]),nac=np.array(myns[2]),nag=np.array(myns[3]),nta=np.array(myns[4]),ntt=np.array(myns[5]),ntc=np.array(myns[6]),ntg=np.array(myns[7]),nca=np.array(myns[8]),nct=np.array(myns[9]),ncc=np.array(myns[10]),ncg=np.array(myns[11]),nga=np.array(myns[12]),ngt=np.array(myns[13]),ngc=np.array(myns[14]),ngg=np.array(myns[15]))

def write_to_file_array_site_pairs_intergene(outfile,celllist,gene1,gene2,site1,site2,myns):
    gene1=np.array(gene1)
    gene2=np.array(gene2)
    site1=np.array(site1)
    site2=np.array(site2)
    
    np.savez(outfile,gene1=gene1,gene2=gene2,site1=site1,site2=site2,naa=np.array(myns[0]),nat=np.array(myns[1]),nac=np.array(myns[2]),nag=np.array(myns[3]),nta=np.array(myns[4]),ntt=np.array(myns[5]),ntc=np.array(myns[6]),ntg=np.array(myns[7]),nca=np.array(myns[8]),nct=np.array(myns[9]),ncc=np.array(myns[10]),ncg=np.array(myns[11]),nga=np.array(myns[12]),ngt=np.array(myns[13]),ngc=np.array(myns[14]),ngg=np.array(myns[15]))


def write_to_file_array_site_pairs_old(outfile,celllist,gene,site1,site2,naa,nat,nac,nag,ntt,ntc,ntg,ncc,ncg,ngg):
    gene=np.array(gene)
    site1=np.array(site1)
    site2=np.array(site2)
    ntt=np.array(ntt)
    ntc=np.array(ntc)
    ntg=np.array(ntg)
    ncc=np.array(ncc)
    ncg=np.array(ncg)
    ngg=np.array(ngg)
    naa=np.array(naa)
    nat=np.array(nat)
    nac=np.array(nac)
    nag=np.array(nag)
    np.savez(outfile,gene=gene,site1=site1,site2=site2,naa=naa,nat=nat,nac=nac,nag=nag,ntt=ntt,ntc=ntc,ntg=ntg,ncc=ncc,ncg=ncg,ng=ngg)

def load_single_sitefile(infile):
    npfile=np.load(infile)
    gene=npfile['gene']
    site=npfile['site']
    na=npfile['na']
    nt=npfile['nt']
    nc=npfile['nc']
    ng=npfile['ng']
    return na,nt,nc,ng,gene,site

def avelist(l):
    return float(sum(l))/len(l)

def check_if_twofold_degen(na,nt,nc,ng,gene,site):
    #wtk how often twofold sites only have alleles at the two syn codons
    frac_nonsyn=[]#for each site, fraction of cells that do not have the synonymous allele
    #print(len(na),len(nt),len(nc),len(ng))
    nzero=0
    nag=[na[i]+ng[i] for i in range(len(ng))]
    nct=[nc[i]+nt[i] for i in range(len(nt))]
    newna=[]
    newnt=[]
    newng=[]
    newnc=[]#for only sites with syn snps
    for i in range(len(nag)):
        if nag[i]>0 and nct[i]==0:
            #frac_nonsyn.append(0)
            nzero+=1
            newna.append(na[i])            
            newng.append(ng[i])
            #don't have c or t bc this is not a ct twofold site
            continue
        if nag[i]==0 and nct[i]>0:
            #frac_nonsyn.append(0)
            newnt.append(nt[i])
            newnc.append(nc[i])
            nzero+=1
            continue
        if nag[i]>0 and nct[i]>0:
            frac_nonsyn.append(float(min(nag[i],nct[i]))/(nag[i]+nct[i]))
    print('Number of Nonsyn Sites: %i' %len(frac_nonsyn))
    return newna,newnt,newnc,newng

def downsample_site(myna,mynt,mync,myng,nmc,downsamplesize):
    mycells=np.random.choice(np.array(list(range(1,nmc+1))),size=downsamplesize,replace=False)
    newna,newnt,newnc,newng=0,0,0,0
    for item in mycells:
        if item<=myna:
            newna+=1
        if item>myna and item<=myna+mynt:
            newnt+=1
        if item>mynt+myna and item<=myna+mynt+mync:
            newnc+=1
        if item>mynt+myna+mync:
            newng+=1
    return newna,newnt,newnc,newng

def downsample_site_pair(myns,nmc,downsamplesize,verbosebool=0):   
    mylist=[]
    for i in range(len(myns)):
        for k in range(myns[i]):
            mylist.append(i)
    mycells=list(np.random.choice(np.array(mylist),size=downsamplesize,replace=False))
   
    new_myns=[0 for i in range(len(myns))]
    for i in range(len(new_myns)):
        new_myns[i]=mycells.count(i)

    if verbosebool:
        print('myns',myns)
        print('totlist', mylist)
        print('selected list',mycells)
        print('newlist',new_myns)
    return new_myns
    
    
def get_longrange_gene_sample(seplower,sepupper,ngenes):
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    if seplower+1==sepupper:
        glist1=mycores[:-sepupper]
        print(seplower,sepupper,len(glist1))
        genelist1=[]
        genelist2=[]
        for i in range(len(glist1)):
            myg1=glist1[i]
            myg2=myg1+seplower
            if myg2 in mycores:
                genelist1.append(myg1)
                genelist2.append(myg2)
        print('Seps %i - %i : %i gene pairs' %(seplower,sepupper, len(genelist1)))
        return genelist1,genelist2
    genes1_tries=random.choices(mycores,k=ngenes*4)
    genes1=[]
    genes2=[]
    myseps=random.choices(list(range(seplower,sepupper)),k=ngenes*4)
    
    for i in range(len(genes1_tries)):
        g1=genes1_tries[i]
        g2=g1+myseps[i]
        if g2 not in mycores:
            continue
        genes1.append(g1)
        genes2.append(g2)
        if len(genes1)>=ngenes:
            break
    print(seplower,sepupper,len(genes1))
    return genes1,genes2

def get_longrange_gene_sample_all(seplower,sepupper):
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    if seplower+1==sepupper:
        glist1=mycores[:-sepupper]
        print(seplower,sepupper,len(glist1))
        genelist1=[]
        genelist2=[]
        for i in range(len(glist1)):
            myg1=glist1[i]
            myg2=myg1+seplower
            if myg2 in mycores:
                genelist1.append(myg1)
                genelist2.append(myg2)
        print('Seps %i - %i : %i gene pairs' %(seplower,sepupper, len(genelist1)))
        return genelist1,genelist2

    genes1,genes2=[],[]
    for i in range(len(mycores)):
        gi=mycores[i]
        for j in range(seplower,sepupper):
            gj=gi+j
            if gj not in mycores:
                continue
            genes1.append(gi)
            genes2.append(gj)
    print('Seps %i- %i: %i gene pairs' %(seplower,sepupper,len(genes1)))
    return genes1,genes2


def get_gene_from_key_to_check(key):
    idx1=key.find('_')
    idx2=key.find('_',idx1+1)
    g=int(key[idx1+1:idx2])
    return g



fourfoldlist=['CT','GT','TC','CC','GC','AC','GG','CG']
twofold_list_ct=['TTT','TTC','TAC','TAT','CAT','CAC','AAT','AAC','GAT','GAC','AGT','AGC','TGT','TGC']
twofold_list_ag=['TTA','TTG','CAA','CAG','AAA','AAG','GAA','GAG','AGA','AGG']
twofoldlist=twofold_list_ct+twofold_list_ag
