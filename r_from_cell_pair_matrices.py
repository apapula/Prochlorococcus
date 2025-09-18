import numpy as np
import math
import ast
import statistics
from os.path import exists
import statistics
import argparse
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,AllKashBATS_without_ribotype, BB_closegrp_1, HLII_Outlier_Greater_Than_8_Percent, BB_closegrp_2, HLII_NoClose
import matplotlib.pyplot as plt
import Repo_site_catalogues_single_pair as repo
import random
from scipy.stats import linregress, pearsonr
import csv
import networkx as nx
import re
from matplotlib.backends.backend_pdf import PdfPages


'''
this script uses the pre-computed distance matrices to calculate linkage of gene pair divergences along the (highly syntenic) genome, here using the measure we call R: the pearson correlation coefficient of cell pair divergences at a pair of genes


'''



def get_r_with_intervals(infile_fullgene,celllist,sepbins,npairs,groupname,sherlockbool,discretebool,tfna):
   

    sep_midpoints,rs,r_1_up,r_2_up,r_1_down,r_2_down=get_r_and_separation_intervals(sepbins,npairs,celllist,infile_fullgene,discretebool,tfna)
    
    print('seps=',sep_midpoints)
    print('R=',rs)
   
    print('R_1sigma_up=',r_1_up)
    print('R_2sigma_up=',r_2_up)
    print('R_1sigma_down=',r_1_down)
    print('R_2sigma_down=',r_2_down)


def get_r_and_separation_intervals(sepbins,npairs,celllist,arrayfile,discretebool,tfna):
    overalllist=overall_order_list()
    if tfna in ['t','f','n','a','aa']:
        genearrays,mykeys=load_arrays(arrayfile)
    if tfna=='tf':
        genearrays=load_two_arrays(arrayfile)
    sep_midpoints=[]#if sepbins is [1,2,3,4,10], want result to be [1,2,3,6.5]. separations are [seplower, sepupper-1]
    rs=[]
    r_1_up,r_2_up,r_1_down,r_2_down=[],[],[],[]
    for i in range(len(sepbins)-1):
        midpt=0.5*(sepbins[i]+sepbins[i+1]-1)
        sep_midpoints.append(midpt)
        r,r1u,r2u,r1d,r2d=do_single_separation_bin_interval(sepbins[i],sepbins[i+1],genearrays,npairs,overalllist,celllist,discretebool)
        if r!='-':
            rs.append(r)
            r_1_up.append(r1u)
            r_1_down.append(r1d)
            r_2_down.append(r2d)
            r_2_up.append(r2u)
        #now print sep bin
        print('midpoints')
        print(sep_midpoints)
        print('rs')
        print(rs)
        print('r_1_up')
        print(r_1_up)
        print('r_2_up')
        print(r_2_up)
        print('r_1_down')
        print(r_1_down)
        print('r_2_down')
        print(r_2_down)
              
        
    return sep_midpoints,rs,r_1_up,r_2_up,r_1_down,r_2_down

def do_single_separation_bin_interval(seplow,sephigh,genearrays,npairs,overalllist,celllist,discretebool):
   
    mykeys=list(genearrays.keys())
    keys1,keys2=get_random_pair_sample(npairs,mykeys,seplow,sephigh)
    my_rs=[]
    for i in range(len(keys1)):
        r,slope,npairs=get_r_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool)
        if r!='-':
            my_rs.append(r)
        
    if len(my_rs)==0:
        return '-','-','-','-','-'
    r= repo.avelist(my_rs)
    r1u=np.quantile(my_rs,.6827)
    r1d=np.quantile(my_rs,.3173)
    r2u=np.quantile(my_rs,.9545)
    r2d=np.quantile(my_rs,.0455)
    return r, r1u,r2u,r1d,r2d
######################


    
def get_r_and_separation_for_fullgene_pairs(sepbins,npairs,celllist,arrayfile,discretebool,tfna):
    overalllist=overall_order_list()
    if tfna in ['t','f','n','a','aa']:
        genearrays,mykeys=load_arrays(arrayfile)
    if tfna=='tf':
        genearrays=load_two_arrays(arrayfile)
    sep_midpoints=[]#if sepbins is [1,2,3,4,10], want result to be [1,2,3,6.5]. separations are [seplower, sepupper-1]
    rs=[]
    shuffled_rs=[]
    for i in range(len(sepbins)-1):
        midpt=0.5*(sepbins[i]+sepbins[i+1]-1)
        sep_midpoints.append(midpt)
        r,shuff_r=do_single_separation_bin(sepbins[i],sepbins[i+1],genearrays,npairs,overalllist,celllist,discretebool)
        
            rs.append(r)
            shuffled_rs.append(shuff_r)
    #print(sep_midpoints)
    #print(rs)
    return sep_midpoints,rs,shuffled_rs

def do_single_separation_bin(seplow,sephigh,genearrays,npairs,overalllist,celllist,discretebool):
    mykeys=list(genearrays.keys())
    keys1,keys2=get_random_pair_sample(npairs,mykeys,seplow,sephigh)
    my_rs=[]
    shuffled_rs=[]
    for i in range(len(keys1)):
        r,slope,mynpairs=get_r_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool)
        if r!='-':
            my_rs.append(r)
        #now need to shuffle
        sr,sslope,_=get_r_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,shufflebool=1)
        if sr!='-':
            shuffled_rs.append(sr)
    if len(my_rs)==0 or len(shuffled_rs)==0:
        return '-','-'
    return repo.avelist(my_rs),repo.avelist(shuffled_rs)
            
def get_random_pair_sample(npairs,mykeys,seplower,sepupper):
    #mykeys are the genes included
    key1,key2=[],[]
    if npairs<500000:
        #do for all gene pairs 1,2,3,4,...10 genes away. then random sample
        genes1,genes2=repo.get_longrange_gene_sample(seplower,sepupper,npairs)
        #now convert into key names
        templatekey=mykeys[0]
        idx1=templatekey.find('_')
        idx2=templatekey.find('_',idx1+1)
        prefix=templatekey[:idx1+1]
        suffix=templatekey[idx2:]
        for i in range(len(genes1)):
            if '%s%s%s' %(prefix,str(genes1[i]),suffix) not in mykeys or '%s%s%s' %(prefix,str(genes2[i]),suffix) not in mykeys:
                continue
            key1.append('%s%s%s' %(prefix,str(genes1[i]),suffix))
            key2.append('%s%s%s' %(prefix,str(genes2[i]),suffix))
    if npairs>500000:#all pairs
        genes1,genes2=repo.get_longrange_gene_sample_all(seplower,sepupper)
        templatekey=mykeys[0]
        idx1=templatekey.find('_')
        idx2=templatekey.find('_',idx1+1)
        prefix=templatekey[:idx1+1]
        suffix=templatekey[idx2:]
        for i in range(len(genes1)):
            
            key1.append('%s%s%s' %(prefix,str(genes1[i]),suffix))
            key2.append('%s%s%s' %(prefix,str(genes2[i]),suffix))
    return key1,key2



        
def get_gene_from_key_to_check(key):
    idx1=key.find('_')
    idx2=key.find('_',idx1+1)
    g=int(key[idx1+1:idx2])
    return g


def load_two_arrays(arrayfilelist):
    genearrays={}
    npfile=np.load(arrayfilelist[0],allow_pickle=1)
    mykeys=npfile[()].keys()
    for item in mykeys:
        genearrays[item]=npfile[()][item]
    npfile2=np.load(arrayfilelist[1],allow_pickle=1)
    mykeys2=npfile2[()].keys()
    for item in mykeys2:
        genearrays[item]=npfile2[()][item]
    return genearrays


def get_r_gene_pair(mat1,mat2,overalllist,celllist,key1,key2,discretebool,shufflebool=0,sherlockbool=1,pdfbool=0,pdf='',tfna=''):
    divs1=[]
    divs2=[]
    

    for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            if not discretebool:
                d1=mat1[idxi][idxj]
                d2=mat2[idxi][idxj]
                
                if d1<1 and d2<1:
                    divs1.append(d1)
                    divs2.append(d2)
            if discretebool:
                [n1,l1]=mat1[idxi][idxj]
                [n2,l2]=mat2[idxi][idxj]
                if l1<12 or l2<12:
                    continue
                d1=float(n1)/l1
                d2=float(n2)/l2
                divs1.append(d1)
                divs2.append(d2)
    if len(divs1)<10:
        if pdfbool:
            return '-','-','-'
        if not pdfbool:
            return '-','-','-'
    #now plt and r
    if shufflebool:
        random.shuffle(divs1)
        random.shuffle(divs2)
    slope, intercept, r_value, p_value, std_err = linregress(divs1,divs2)
    
    if 0:#not sherlockbool:
        x=np.linspace(0,max(max(divs1),max(divs2)),100)
        plt.plot(x,slope*x+intercept)
        plt.plot(divs1,divs2,'.',alpha=0.3)
        plt.title('slope %f intercept %f r value %f\n%s %s\n ave div half 1 %f half 2 %f' %(slope,intercept,r_value,key1,key2,repo.avelist(divs1),repo.avelist(divs2)))
        plt.show()
        print('%s %s r value %f' %(key1,key2,r_value))
    
    if math.isnan(r_value):
        print('nan')
        r_value='-'
    if pdfbool:
        if r_value>0.6:
            g1=get_gene_from_key_to_check(key1)
            g2=get_gene_from_key_to_check(key2)
            x=np.linspace(0,max(max(divs1),max(divs2)),100)
            plt.plot(x,slope*x+intercept)
            plt.plot(divs1,divs2,'.',alpha=0.3)
            plt.title('%s Genes %i and %i %s\nslope %f intercept %f r value %f' %(groupname,g1,g2,tfna,slope,intercept,r_value))
            pdf.savefig()
            plt.close()
        return r_value,covariance(divs1,divs2),len(divs1)
    if not pdfbool:
        return r_value,covariance(divs1,divs2),len(divs1)

def covariance(x,y):
    mx=repo.avelist(x)
    my=repo.avelist(y)
    n=len(x)
    cov=0
    for i in range(n):
        cov+=(x[i]-mx)*(y[i]-my)
    cov*=(1.0/(n-1))
    return cov

def load_arrays(arrayfile):
    #load into a dict
    genearrays={}
    npfile=np.load(arrayfile,allow_pickle=1)
    mykeys=npfile[()].keys()
    for item in mykeys:
        genearrays[item]=npfile[()][item]
    return genearrays,sorted(list(genearrays.keys()))

def overall_order_list():
    mylist=[]
    for i in range(len(hliilist2)):
        mylist.append(hliilist2[i])
    for i in range(len(AllKashBATS_without_ribotype)):
        mylist.append(AllKashBATS_without_ribotype[i])
    return mylist

##########################################################################################




###########################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB')#ntitle is title of matrix file. n is the group to plot
    parser.add_argument('-shr','--sherlockbool',default=1,type=int)
    parser.add_argument('-f','--frac',default=2,type=int)
    parser.add_argument('-np','--npairs',type=int,default=800000)
    parser.add_argument('-db','--discretebool',type=int,default=1)
    parser.add_argument('-cb','--corebool',type=int,default=1)
    parser.add_argument('-sg','--startgene',default=0,type=int)
    parser.add_argument('-dl','--distancelow',type=int,required=False)
    parser.add_argument('-dh','--distancehigh',type=int,required=False)
    parser.add_argument('-w','--windowsize',type=int,required=False)
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    discretebool=args.discretebool
    startgene=args.startgene
    corebool=args.corebool
    frac=args.frac
    n=args.n
    tfna=args.tfna
    suffix=''
    
    titlegroupname='Berube_Kash'
    if n==0:
        groupname='HLII'
        celllist=hliilist2
        paircut=500
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype
    if n==2:
        groupname='BB_NoClose'
        celllist=BB_list
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
        paircut=500
    if n==3:
        celllist=clique1_without_ribotype
        groupname='Clique1'
    
  
        
    if tfna=='t':
        outfilewhole='Matrices_%s_TWOFOLDDEGENERATE%s' %(titlegroupname,suffix)
        outfilefrac='Matrices_%s_TWOFOLDDEGENERATE%s' %(titlegroupname,suffixfrac)
    if tfna=='tf':
        outfilewhole='Matrices_%s_TWO_AND_FOUR%s' %(titlegroupname,suffix)
        outfilefrac='Matrices_%s_TWO_AND_FOUR%s' %(titlegroupname,suffixfrac)
    if tfna=='f':
        outfilewhole='Matrices_%s_FOURFOLDDEGENERATE%s' %(titlegroupname,suffix)
        outfilefrac='Matrices_%s_FOURFOLDDEGENERATE%s' %(titlegroupname,suffixfrac)
    if tfna=='n':
        
        outfilewhole='Matrices_Berube_Kash_NONSYN_SECONDNT_Discrete_Including_Non_Core'
        outfilefrac=''
    if tfna=='a':
       
        outfilewhole='Matrices_Berube_Kash_All_NT_Discrete_Including_Non_Core'
        outfilefrac='Matrices_%s_All_NT%s' %(titlegroupname,suffixfrac)
    if tfna=='aa':
        outfilefrac=''
        outfilewhole='Matrices_Berube_Kash_AA_Discrete_Including_Non_Core'
    outfilewhole='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfilewhole
    outfilefrac='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfilefrac
    outfilefracrandom=outfilefrac+'_RANDOM'
    outfilewhole+='.npy'
    outfilefrac+='.npy'
    outfilefracrandom+='.npy'

    if discretebool:
        if tfna=='tf':
            titlegroupname='Berube_Kash'
            suffix='_Discrete'
            outfilewhole='Matrices_%s_TWO_AND_FOUR%s' %(titlegroupname,suffix)
            outfilewhole='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfilewhole
            infilematrices=['%s_firsthalf.npy' %outfilewhole,'%s_secondhalf.npy' %outfilewhole]
            outfilewhole=infilematrices

    if args.sherlockbool==0:
        outfilewhole='../input_data/'+outfilewhole
        outfilefrac='../input_data/'+outfilefrac
        outfilefracrandom='../input_data/'+outfilefracrandom
    npairs=args.npairs

    sepbins=[1,2,3,4,5,6,7,8,9,10,25,63,158,398,1000]
    if prog==0:
        print(outfilewhole)
        print(outfilefrac)
        print(outfilefracrandom)
        plot_r_and_random(outfilewhole,outfilefrac,outfilefracrandom,frac,celllist,sepbins,npairs,groupname,args.sherlockbool,discretebool,tfna)
   
