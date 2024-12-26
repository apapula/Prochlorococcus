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
needs to work for gene fractions and whole genes'

first, load the keys of the dictionary

for gene halfs: 
plot half 1 vs half 2 for all genes. get a value. compare to random.

for whole genes: 
pick a set of N unique pairs
plot gene x vs gene y div for MC cell pairs
get R
add this to a list with seps, a list with Rs
use scipy.statistics to bin this later, and plot
To get a random shuffle, have for each gene in a pair: div_gene_x_list, div_gene_y_list in same cell pair ordering. shuffle 1 or both of these to get a random R

compare random shuffling of half genes to random site selection. random site selection should still be much higher r than shuffled

will want to just print the halfgene values and the fullgene values, plot locally
'''

bigcluster_BB_6_100= [11, 12, 14, 21, 23, 25, 35, 37, 70, 71, 77, 78, 103, 108, 111, 125, 143, 147, 148, 149, 150, 167, 172, 261, 267, 268, 269, 270, 274, 285, 287, 289, 290, 291, 292, 293, 296, 314, 315, 328, 334, 335, 337, 348, 351, 352, 354, 412, 413, 451, 482, 490, 491, 503, 504, 510, 511, 512, 514, 534, 536, 537, 538, 541, 568, 587, 588, 590, 595, 604, 622, 683, 804, 807, 808, 809, 810, 814, 815, 817, 825, 898, 904, 905, 906, 907, 908, 910, 913, 914, 945, 949, 950, 1025, 1041, 1042, 1044, 1046, 1047, 1048, 1050, 1055, 1058, 1067, 1072, 1091, 1154, 1173, 1174, 1175, 1176, 1183, 1185, 1188, 1189, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1199, 1219, 1220, 1221, 1224, 1271, 1276, 1287, 1303, 1347, 1373, 1379, 1380, 1381, 1382, 1383, 1386, 1387, 1388, 1389, 1469, 1475, 1481, 1502, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1625, 1627, 1628, 1629, 1630, 1631, 1633, 1673, 1674, 1675, 1676, 1677, 1678, 1680, 1681, 1682, 1683, 1690, 1692, 1693, 1694, 1695, 1696, 1697, 1698, 1699, 1700, 1728, 1743, 1765, 1766, 1767, 1792, 1796, 1797, 1800, 1801, 1802, 1803, 1804, 1817, 1820, 1821, 1822, 1823, 1825, 1826, 1827, 1828, 1829, 1830, 1831, 1838, 1878, 1903, 1904, 1905, 1907]


def get_r_with_intervals(infile_fullgene,celllist,sepbins,npairs,groupname,sherlockbool,discretebool,tfna):
   

    sep_midpoints,rs,r_1_up,r_2_up,r_1_down,r_2_down=get_r_and_separation_intervals(sepbins,npairs,celllist,infile_fullgene,discretebool,tfna)
    
    print('seps=',sep_midpoints)
    print('R=',rs)
    #get: 68.27, 95.45 % of data less than this value for upper sigmas; for lower, 31.73% less than this, 4.55% I guess? INCREASE NUMBER OF PAIRS
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

def plot_r_and_random(infile_fullgene,infilefrac,infilefracrandom,frac,celllist,sepbins,npairs,groupname,sherlockbool,discretebool,tfna):
    #print('fraction point')
    #genehalfr,genehalfr_random=get_gene_fraction_point(frac,celllist,infilefrac,infilefracrandom,sherlockbool)
   
    #for each sep in sepbins, get the sample of genes to use. convert to key names.  Then for each pair in the keypairs, get r. average for the sep bin. plot. ALSO DO THE RANDOM ONE!

    sep_midpoints,rs,shuffled_rs=get_r_and_separation_for_fullgene_pairs(sepbins,npairs,celllist,infile_fullgene,discretebool,tfna)
    
    if sherlockbool==1:
        #print('Gene Halves: %f Data, Random %f' %(genehalfr,genehalfr_random))
        #print('half_gene=[ %f ]' %genehalfr)
        #print('half_gene_random_sites=[%f]' %genehalfr_random)
        print('seps=',sep_midpoints)
        print('R=',rs)
        print('R_random=',shuffled_rs)
    if sherlockbool==0:
        plt.plot(0.5,genehalfr,'o',color='b',label='Data')
        plt.plot(0.5,genehalfr_random,'o',color='red',label='Random Sites')
        plt.plot(sep_midpoints,rs,'o',color='b')
        plt.plot(sep_midpoints,shuffled_rs,'o',color='orange',label='Random Shuffle')
        plt.legend()
        plt.xlabel('Separation (Number of Genes)')
        plt.ylabel('Average R for Gene Pairs')
        plt.title('Average Correlation between cell pair divergences at gene pairs\n%s (%i Gene Pairs/Wide Bin)' %(groupname,npairs))
        plt.show()
    
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
        #print(i,midpt,r,shuff_r)
        #t=input('t')
        if r!='-':
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


def get_gene_fraction_point_test(frac,celllist,infile):
    #only works for frac=2
    sherlockbool=1
    overalllist=overall_order_list()
    
    genearrays,arraykeys=load_arrays(infile)
    random_r,random_slope,r,slope=[],[],[],[]
    
    mynpairs=int(0.5*len(arraykeys))
    for i in range(mynpairs):
        key1=arraykeys[i*2]
        key2=arraykeys[i*2+1]
        #make sure the same gene
        if get_gene_from_key_to_check(key1)!=get_gene_from_key_to_check(key2):
            print('Not same Gene: %s %s' %(key1,key2))
            return
        #print('Frac Gene ',key1,key2)
        myr,myslope,_=get_r_gene_pair(genearrays[key1],genearrays[key2],overalllist,celllist,key1,key2,discretebool,sherlockbool=sherlockbool)
        
            
        if myr!='-':
            r.append(myr)
            slope.append(myslope)
        #now for random
       
    print('Fraction Gene halves R %f' %(repo.avelist(r)))#,repo.avelist(random_r)))   
    return repo.avelist(r)#,repo.avelist(random_r)

def get_gene_fraction_point(frac,celllist,infile,infilerandom):
    sherlockbool=1
    #only works for frac=2
    overalllist=overall_order_list()
    
    genearrays,arraykeys=load_arrays(infile)
    print(arraykeys)
    genearraysrand,arraykeysrand=load_arrays(infilerandom)
    random_r,random_slope,r,slope=[],[],[],[]
    
    mynpairs=int(0.5*len(arraykeys))
    for i in range(mynpairs):
        key1=arraykeys[i*2]
        key2=arraykeys[i*2+1]
        print(key1,key2)
        #make sure the same gene
        if get_gene_from_key_to_check(key1)!=get_gene_from_key_to_check(key2):
            print('Not same Gene: %s %s' %(key1,key2))
            return
        #print('Frac Gene ',key1,key2)
        myr,myslope,_=get_r_gene_pair(genearrays[key1],genearrays[key2],overalllist,celllist,key1,key2,sherlockbool=sherlockbool)
        
            
        if myr!='-':
            r.append(myr)
            slope.append(myslope)
        #now for random
        key1=arraykeysrand[i*2]
        key2=arraykeysrand[i*2+1]
        #make sure the same gene
        if get_gene_from_key_to_check(key1)!=get_gene_from_key_to_check(key2):
            print('Not same Gene: %s %s' %(key1,key2))
            return
        #print('Shuffled')
        myr,myslope,_=get_r_gene_pair(genearraysrand[key1],genearraysrand[key2],overalllist,celllist,key1,key2,sherlockbool=sherlockbool)
        if myr!='-':
            random_r.append(myr)
            random_slope.append(myslope)
    print('Fraction Gene Ordered halves R %f ; random R %f' %(repo.avelist(r),repo.avelist(random_r)))
    print('Fraction Gene Ordered halves Covariance %f ; random Covariance %f' %(repo.avelist(slope),repo.avelist(random_slope)))
    return repo.avelist(r),repo.avelist(random_r)
        
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
    '''
    print(len(overalllist))
    print(np.shape(mat1))
    print(mat1[261])
    print(mat1[190])
    print(mat1[200])
    print(mat1[201])
    t=input('t')
    '''

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
    #myr,_=pearsonr(divs1,divs2)
    #print('r %f, covariance %f' %(myr,covariance(divs1,divs2)))
    if 0:#not sherlockbool:
        x=np.linspace(0,max(max(divs1),max(divs2)),100)
        plt.plot(x,slope*x+intercept)
        plt.plot(divs1,divs2,'.',alpha=0.3)
        plt.title('slope %f intercept %f r value %f\n%s %s\n ave div half 1 %f half 2 %f' %(slope,intercept,r_value,key1,key2,repo.avelist(divs1),repo.avelist(divs2)))
        plt.show()
        print('%s %s r value %f' %(key1,key2,r_value))
    #print('R %f Slope %f' %(r_value,slope))
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

def make_R_histogram(infilefullgene,celllist,groupname,sherlockbool,npairs,tfna,distancelow=100,distancehigh=1000):
    tfnastring=repo.return_tfna_string(tfna)
    overalllist=overall_order_list()
    r,rshuf=[],[]
    if tfna!='tf':
        genearrays,mykeys=load_arrays(infilefullgene)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
        mykeys=list(genearrays.keys())
    keys1,keys2=get_random_pair_sample(npairs,mykeys,distancelow,distancehigh)
    
    for  i in range(len(keys1)):
        if i%100==0:
            print('%f (%i of %i) ' %(float(i)/len(keys1),i,len(keys1)))
        myr,slope,_=get_r_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,sherlockbool=sherlockbool)
        if myr!='-':
            r.append(myr)
        #now need to shuffle
        mysr,sslope,_=get_r_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,shufflebool=1)
        if mysr!='-':
            rshuf.append(mysr)
            
    '''
    plt.hist([r,rshuf],histtype='step',bins=50,label=['Data','Shuffle'])
    
    plt.legend()
    plt.ylabel('Number of Gene Pairs')
    plt.xlabel('Pearson R')
    ax=plt.gca()
    #ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    plt.title('%s R for %s Cell Pair Divergences\nGenes Separated %i-%i Genes, %i Pairs' %(groupname,tfnastring,distancelow,distancehigh,npairs))
    plt.savefig('%s_%s_R_Histograms_AllSites_%ipairs.png' %(groupname,tfna,npairs))
    plt.close()
    '''

    myn,newbins,mypatches=plt.hist([r,rshuf],histtype='step',bins=50,label=['Data','Shuffle'])
    
    plt.legend()
    plt.ylabel('Number of Gene Pairs')
    plt.xlabel('Pearson R')
    ax=plt.gca()
    plt.yscale('log')
    #ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    plt.title('%s R for %s Cell Pair Divergences\nGenes Separated %i-%i Genes, %i Pairs' %(groupname,tfnastring,distancelow,distancehigh,npairs))
    plt.savefig('%s_%s_R_Histograms_AllSites_%ipairs_log.png' %(groupname,tfna,npairs))
    plt.close()

    print('Counts (data, random)')
    print(myn)
    print('Bins')
    print(list(newbins))



def R_all_Pairs(infilefullgene,celllist,groupname,sherlockbool,tfna,startgene=0):
    print(infilefullgene)
    tfnastring=repo.return_tfna_string(tfna)
    overalllist=overall_order_list()
    pdf=PdfPages('%s_R_Greater_Than_0.6_%s_withNum.pdf' %(groupname,tfna))
    r,rshuf=[],[]
    if tfna!='tf':
        genearrays,mykeys=load_arrays(infilefullgene)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
        mykeys=list(genearrays.keys())

    
    sample_genenums=repo.make_sample_all_single_genes()
    #sample_genenums=bigcluster_BB_6_100
    if 1:#all pairs
        dl=100
        dh=1000
        print('R %s %s all pairs, dl %i dh %i' %(groupname,tfna,dl,dh))
        rs=[]
        for i in range(len(sample_genenums)):
            if i%100==0:
                print('%i of %i' %(i,len(sample_genenums)))
            genei=sample_genenums[i][0]
            if tfna!='aa':
                keyi='Gene_%i_%s_Pairs' %(genei,tfna)
            if tfna=='aa':
                keyi='Gene_%i_AA' %(genei)
            for j in range(i):
                genej=sample_genenums[j][0]
                if abs(genei-genej)<100:
                    continue
                if tfna!='aa':
                    keyj='Gene_%i_%s_Pairs' %(genej,tfna)
                if tfna=='aa':
                    keyj='Gene_%i_AA' %(genej)
                myr,slope,numcells=get_r_gene_pair(genearrays[keyi],genearrays[keyj],overalllist,celllist,keyi,keyj,discretebool,sherlockbool=sherlockbool,pdfbool=1,pdf=pdf,tfna=tfna)
                if myr=='-':
                    continue
                rs.append(myr)
        myn,newbins,mypatches=plt.hist(rs,histtype='step',bins=50)
        plt.close()

        print('Counts (data only)')
        print(myn)
        print('Bins')
        print(list(newbins))
                
    if 0:
        outfile='%s_R_Greater_Than_0.6_%s_withNum.csv' %(groupname,tfna)
        fieldnames=['Gene1','Gene2','R','NPairs']
        with open(outfile,'w') as myoutfile:
            outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
            outwriter.writeheader()

            for i in range(len(sample_genenums)):
                if i%100==0:
                    print('%i of %i' %(i,len(sample_genenums)))
                genei=sample_genenums[i]#[0]
                if genei<startgene:
                    continue
                if tfna!='aa':
                    keyi='Gene_%i_%s_Pairs' %(genei,tfna)
                if tfna=='aa':
                    keyi='Gene_%i_AA' %(genei)
                for j in range(i):
                    genej=sample_genenums[j]#[0]
                    if tfna!='aa':
                        keyj='Gene_%i_%s_Pairs' %(genej,tfna)
                    if tfna=='aa':
                        keyj='Gene_%i_AA' %(genej)
                    myr,slope,numcells=get_r_gene_pair(genearrays[keyi],genearrays[keyj],overalllist,celllist,keyi,keyj,discretebool,sherlockbool=sherlockbool,pdfbool=1,pdf=pdf,tfna=tfna)
                    if myr=='-':
                        continue
                    if myr>0.60:
                        outwriter.writerow({'Gene1':genei,'Gene2':genej,'R':myr,'NPairs':numcells})

        pdf.close()
    
def process_R_allpairs(groupname,tfna,rcut,paircut,corebool):
    infile='../input_data/%s_R_Greater_Than_0.6_%s_withNum.csv' %(groupname,tfna)
    #infile='../input_data/%s_R_Greater_Than_0.6_%s.csv' %(groupname,tfna)
    minsep=10

    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]  

    
    '''
    wtk: how many genes? what are the highest R? ideally would have liked to know how many covered
    create a graph. get cliques/connected components
    '''
    edges=[]
    Rvals=[]
    rdict={}
    separated_rvals100,separated_rvals10=[],[]
    with open(infile,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            
            r=float(row['R'])
            g1=int(row['Gene1'])
            g2=int(row['Gene2'])
            if corebool:
                if g1 not in mycores or g2 not in mycores:
                    continue
            Rvals.append(r)
            if abs(g1-g2)>=100:
                separated_rvals100.append(r)
            if abs(g1-g2)>=10:
                separated_rvals10.append(r)
            if r>=rcut:
                if int(row['NPairs'])<paircut:
                    continue
                if abs(g1-g2)<minsep:
                    continue
                edges.append((g1,g2))
                if g1 in rdict:
                    rdict[g1]+=1
                if g2 in rdict:
                    rdict[g2]+=1
                if g1 not in rdict:
                    rdict[g1]=1
                if g2 not in rdict:
                    rdict[g2]=1
    print(rdict)
    t=input('t')
    
    if corebool:
        groupname+=' Core'
    if 1:
        ntotpair=782000
        rvalabove=[item for item in Rvals if item>=.7]
        plt.hist(Rvals,bins=40)
        plt.xlabel('R value of Gene Pair')
        plt.ylabel('Number of Gene Pairs')
        plt.title('%s %s R values for Highly Linked Gene Pairs\n%i high R pairs of %i total pairs (%f)\n%i (%f) above .7' %(groupname,tfna,len(Rvals),ntotpair,float(len(Rvals))/ntotpair,len(rvalabove),float(len(rvalabove))/ntotpair))
        plt.show()

        print('r >.6 >9 sep %i > 99 sep %i' %(len(separated_rvals10),len(separated_rvals100)))
        plt.hist([Rvals,separated_rvals10,separated_rvals100],label=['All Gene Pairs','Separated by >9 Genes','Separated >99 Genes'],bins=40)
        plt.xlabel('R value of Gene Pair')
        plt.ylabel('Number of Gene Pairs')
        plt.title('%s %s R values for Highly Linked Gene Pairs' %(groupname,tfna))
        plt.legend()
        plt.show()

   

    rlist=[rdict[item] for item in rdict]
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    for item in mycores:
        if item not in rdict:
            rlist.append(0)
    mybins=[-0.5]
    for i in range(max(rlist)+2):
        mybins.append(i+0.5)
    plt.hist(rlist,bins=mybins)
    plt.xlabel('Number of Gene Pairs per Gene with R>%f' %rcut)
    plt.ylabel('Number of Genes')
    plt.title('%s: For each gene, the number of other genes\nwith which it has R> %f\nMin Gene Separation %i' %(groupname,rcut,minsep))
    plt.close()#show()
    
    #get_max_clique(edges,rcut,groupname,minsep)
    ntrials=1000
    nedges=len([item for item in separated_rvals10 if item>=.7])
    ngenepairs=781000
    get_graph_connected_comps(edges,rcut,groupname,minsep,ntrials,ngenepairs,nedges)


  

def get_max_clique(edges,rcut,groupname,minsep):
    cliquelengths=[]
    G=nx.Graph()
    G.add_edges_from(edges)
    max_clique=[]
    for clq in nx.clique.find_cliques(G):
        #print(clq)
        cliquelengths.append(len(clq))
        #get_group_annotation(clq)
    plt.hist(cliquelengths,bins=list(range(max(cliquelengths)+1)))
    plt.xlabel('Number of Genes in Clique')
    plt.title('%s Size of Cliques of Genes with R > %f\nMin Gene Separation %i' %(groupname,rcut,minsep))
    plt.ylabel('Number of Cliques')
    #ax=plt.gca()
    #ax.tick_params(axis='both', which='major', labelsize=25)
    plt.show()

def get_connected_component_null(ngenepairs,nedges,ntrials):
    #nedges is the number of gene pairs greater than 10 gene sep with connection
    ngenes=int(.5*(1+math.sqrt(1+8*ngenepairs)))
    print('%i pairs %i genes' %(ngenepairs,ngenes))
    p=float(nedges)/ngenepairs
    print('should have %i edges' %nedges)
    print('prob %f (noting that genes are separated by at least 10 genes)' %p)
    groupsizes=[]
    for i in range(ntrials):
        #G1=nx.erdos_renyi_graph(ngenes,p)
        G1=nx.gnm_random_graph(ngenes,nedges,seed=100)
        #nx.draw(G1, with_labels=True)
        #plt.show()
        print('g1 has %i edges' %(G1.number_of_edges()))
        #now get connected comps
        orthogroups=[]
        for clq in nx.connected_components(G1):
            miniclqlist=[]
            for item in clq:
                miniclqlist.append(item)
            #print(miniclqlist)
            miniclqlist=sorted(miniclqlist)
            if len(miniclqlist)<2:
                continue
            orthogroups.append(miniclqlist)
            groupsizes.append(len(miniclqlist))
    return groupsizes

def get_graph_connected_comps(edges,rcut,groupname,minsep,ntrials,ngenepairs,nedges):
    G=nx.Graph()
    G.add_edges_from(edges)
    nx.draw(G, with_labels=True)
    plt.show()

    orthogroups=[]
    for clq in nx.connected_components(G):
        miniclqlist=[]
        for item in clq:
            miniclqlist.append(item)
        #print(miniclqlist)
        miniclqlist=sorted(miniclqlist)
        #get_group_annotation(miniclqlist)
        orthogroups.append(miniclqlist)
    groupsizes=[len(orthogroups[i]) for i in range(len(orthogroups))]
    print(sorted(groupsizes))
    if 0:
        mylist=list(range(max(groupsizes)+1))
        mybins=[-0.5]
        for i in range(len(mylist)):
            mybins.append(i+0.5)
        plt.hist(groupsizes,bins=mybins)
        plt.xlabel('Number of Genes in Connected Component')
        plt.ylabel('Number of Connected Components')
        plt.title('%s Size of Connected Components of Genes with R > %f\nMin Gene Separation %i' %(groupname,rcut,minsep))
        #plt.yscale('log')
        plt.show()
    if 1:
        simgroupsizes=get_connected_component_null(ngenepairs,nedges,ntrials)
        mylist=list(range(max(groupsizes)+1))
        mybins=[-0.5]
        for i in range(len(mylist)):
            mybins.append(i+0.5)
        plt.hist([groupsizes,simgroupsizes],weights=[[1 for i in range(len(groupsizes))],[1.0/ntrials for i in range(len(simgroupsizes))]],bins=mybins,label=['Data','Simulation'])
        plt.xlabel('Number of Genes in Connected Component')
        plt.legend()
        plt.ylabel('Number of Connected Components')
        plt.title('%s Size of Connected Components of Genes with R > %f\nMin Gene Separation %i' %(groupname,rcut,minsep))
        #plt.yscale('log')
        plt.show()

        plt.hist([groupsizes,simgroupsizes],weights=[[1 for i in range(len(groupsizes))],[1.0/ntrials for i in range(len(simgroupsizes))]],bins=mybins,label=['Data','Simulation'])
        plt.xlabel('Number of Genes in Connected Component')
        plt.legend()
        plt.ylabel('Number of Connected Components')
        plt.title('%s Size of Connected Components of Genes with R > %f\nMin Gene Separation %i' %(groupname,rcut,minsep))
        plt.yscale('log')
        plt.show()


        
def get_group_annotation(genelist):
    annotations=[]
    with open('../input_data/MIT9301_Annotations.csv','r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            g=int(row['GeneNum'])
            if g in genelist:
                annotations.append(row['Annotation'])
    print('Genes ',genelist)
    print('Annotations ',annotations)
    t=input('t')
#####################################
def get_gene_annotation():
    with open('../input_data/MIT9301_NT.fna','r') as myinfile:
        lines=myinfile.read()
    seqs=re.split('>',lines)
    glist,annot=[],[]
    for i in range(len(seqs)):
        sublist=re.split('\n',seqs[i])
        #focus on sublist[0]
        s=sublist[0]
        #print(s)
        if len(s)==0:
            continue
        idx1=s.find(' ')
        idx2=s.rfind('_',0,idx1)
        #prin(idx1,idx2)
        g=int(s[idx2+1:idx1])
        glist.append(g)
        idx3=s.find('protein')
        idx4=s.find(']',idx3)
        annotation=s[idx3+8:idx4]
        annot.append(annotation)
    with open('../input_data/MIT9301_Annotations.csv','w') as myoutfile:
        fieldnames=['GeneNum','Annotation']
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
        outwriter.writeheader()
        for i in range(len(glist)):
            outwriter.writerow({'GeneNum':glist[i],'Annotation':annot[i]})

def get_core_annotations():
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    glist=[]
    annot=[]
    with open('../input_data/MIT9301_Annotations.csv','r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            g=int(row['GeneNum'])
            if g in mycores:
                annot.append(row['Annotation'])
                glist.append(g)
    with open('../input_data/MIT9301_CORE_Annotations.csv','w') as myoutfile:
        fieldnames=['GeneNum','Annotation']
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
        outwriter.writeheader()
        for i in range(len(glist)):
            outwriter.writerow({'GeneNum':glist[i],'Annotation':annot[i]})
##########################################################################################


def R_sliding_window(celllist,groupname,windowsize,infilefullgene,tfna):
    print('%s R sliding window size %i %s' %(groupname,windowsize,tfna))
    windowstart=[]
    aver=[]
    avencells=[]
    overalllist=overall_order_list()
   
    if tfna!='tf':
        genearrays,mykeys=load_arrays(infilefullgene)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
        mykeys=list(genearrays.keys())

    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    genelists=[]
    for i in range(len(mycores)):
        print('%i of %i (%f)' %(i,len(mycores),float(i)/len(mycores)))
        l=[]
        for j in range(windowsize):
            l.append(mycores[(i+j)%len(mycores)])
        genelists.append(l)
        #now get R
        rlist=[]
        ncellslist=[]
        for j in range(len(l)):#genelists)):
            for k in range(j):
                keyi='Gene_%i_%s_Pairs' %(l[j],tfna)
                keyj='Gene_%i_%s_Pairs' %(l[k],tfna)
                myr,slope,numcells=get_r_gene_pair(genearrays[keyi],genearrays[keyj],overalllist,celllist,keyi,keyj,discretebool,sherlockbool=1,pdfbool=0,pdf='',tfna=tfna)
                if myr=='-':
                    continue
                rlist.append(myr)
                ncellslist.append(numcells)
        aver.append(repo.avelist(rlist))
        avencells.append(repo.avelist(ncellslist))
        windowstart.append(l[0])
    #now write

    outfields=['GeneStart','Ave_R','Ave_NCells','Genes']
    with open('%s_R_Sliding_Window_%i_%s.csv' %(groupname,windowsize,tfna),'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,fieldnames=outfields,delimiter='\t')
        outwriter.writeheader()
        for i in range(len(genelists)):
            outwriter.writerow({'GeneStart':windowstart[i],'Ave_R':aver[i],'Ave_NCells':avencells[i],'Genes':genelists[i]})

def get_group_annotations():
    annotations={}
    with open('../input_data/MIT9301_Annotations.csv','r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            g=int(row['GeneNum'])
            annotations[g]=row['Annotation']
    return annotations



def plot_sliding_window(groupname,tfna,windowsize):
    aver=[]
    genestart=[]
    ncells=[]

    annotations=get_group_annotations()
    
    with open('../input_data/%s_R_Sliding_Window_%i_%s.csv' %(groupname,windowsize,tfna),'r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            aver.append(float(row['Ave_R']))
            genestart.append(int(row['GeneStart']))
            ncells.append(float(row['Ave_NCells']))
    #now plot

    #coverage_plot(genestart,aver,ncells,4000)#cutoff)
    coverage_plot(genestart,aver,ncells,5000)
    return
    
    plt.plot(genestart,aver,'.-')#,alpha=0.3)
    plt.xlabel('Position of starting gene on MIT9301 Reference',fontsize=20)
    plt.ylabel('Average R of %i Genes' %(windowsize),fontsize=20)
   
    xtickpos=[0,500,1000,1500,2000]
    plt.xticks(ticks=xtickpos)
     
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("R Averaged over %i-Gene windows\nAlong MIT9301 Reference" %(windowsize),fontsize=20)
    #plt.title('%s R %s %i-Gene windows\nAlong MIT9301 Reference' %(groupname,tfna,windowsize),fontsize=20)
    plt.tight_layout()
    plt.show()

    fig=plt.figure()
    ax1=plt.subplot(211)
    ax2=plt.subplot(212,sharex=ax1)
    ax1.plot(genestart,aver,'.-')
    ax2.plot(genestart,ncells,'.-')
    #ax2.sharex(ax1)
    plt.xlabel('Position of starting gene on MIT9301 Reference',fontsize=20)
    ax1.set_ylabel('Average R of %i Genes' %(windowsize),fontsize=20)
    ax2.set_ylabel('Average number of cell pairs covered')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.suptitle('%s R %s %i-Gene windows\nAlong MIT9301 Reference' %(groupname,tfna,windowsize),fontsize=20)
    plt.tight_layout()
    plt.show()

    plt.hist(aver,bins=100)
    plt.xlabel('Average R for %i-gene window' %windowsize,fontsize=20)
    plt.ylabel('Number of %i-Gene windows' %windowsize,fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('Distribution of average R over %i-gene\nwindows for %s %s' %(windowsize,groupname,tfna),fontsize=20)
    plt.tight_layout()
    plt.show()

    
def coverage_plot(genestart,aver,ncells,cutoff):
    newstarts=[genestart[i] for i in range(len(genestart)) if ncells[i]>cutoff]
    newave=[aver[i] for i in range(len(ncells)) if ncells[i]>cutoff]
    badstarts=[genestart[i] for i in range(len(genestart)) if ncells[i]<=cutoff]
    goodstarts=[genestart[i] for i in range(len(genestart)) if ncells[i]>cutoff]
    newbadstarts=[i for i in range(1,1908) if i not in goodstarts and i in genestart]
    
    plt.plot(newstarts,newave,'.-')#,alpha=0.3)
    #plt.bar(newbadstarts, height=[max(newave) for i in range(len(newbadstarts))],width=1,color='r',alpha=0.3)
    c=0
    for i in range(len(genestart)-1):
        if ncells[i]<cutoff:
            mywidth=genestart[i+1]-genestart[i]
            if c==0:
                plt.bar(genestart[i],height=max(newave),width=mywidth,align='edge',color='k',alpha=0.3,label='Low Cell\nPair Coverage')
                c=1
            if c>0:
                plt.bar(genestart[i],height=max(newave),width=mywidth,align='edge',color='k',alpha=0.3)
    plt.xlabel('Position of starting gene on MIT9301 Reference',fontsize=20)
    plt.ylabel('Average R of %i Genes' %(windowsize),fontsize=20)
    plt.legend(fontsize=15)
    xtickpos=[0,500,1000,1500,2000]
    plt.xticks(ticks=xtickpos)
     
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("R Averaged over %i-Gene windows\nAlong MIT9301 Reference" %(windowsize),fontsize=20)
    #plt.title('%s R %s %i-Gene windows\nAlong MIT9301 Reference' %(groupname,tfna,windowsize),fontsize=20)
    plt.tight_layout()
    plt.show()
###########################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-nt','--ntitle',type=int,help='0: HLII; 1: C1; 2: BB',required=False,default=-1)
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
    #suffix='_CellPairs_MC_No6'#'ConsensusCodon'
    suffix='_Discrete_Including_Non_Core'
    suffixfrac=suffix+'_Frac%i_HalfSites' %frac
    #if args.sherlockbool==0:
        #suffixfrac=suffix+'_Frac%i' %frac
    suffixfrac+='_Sample'
    
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
        celllist=BB_list#HLII_BB_NoCloseCells_Median04
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
        paircut=500
    if n==3:
        celllist=clique1_without_ribotype
        groupname='Clique1'
    if n==4:
        celllist=BB_closegrp_1
        groupname='Closegrp13'

    if n==-1:
        celllist=HLII_Outlier_Greater_Than_8_Percent
        groupname='Outlier_Cells'

     
    if n==5:
        celllist=BB_closegrp_2
        groupname='Closegrp7'

    if n==6:
        celllist=HLII_NoClose
        groupname='HLII_NoClose'
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
        paircut=500
    if n==7:
        #BB with 6 outliers removed
        paircut=0
        celllist=BB_list
        groupname='BB_no_6_outliers'
        if 'AG-347-K23' in celllist:
            celllist.remove('AG-347-K23')
        if 'AG-347-I04' in celllist:
            celllist.remove('AG-347-I04')
        if 'AG-347-M18' in celllist:
            celllist.remove('AG-347-M18')
        if 'AG-442-B03' in celllist:
            celllist.remove('AG-442-B03')
        if 'AG-424-P18' in celllist:
            celllist.remove('AG-424-P18')
        if 'AG-418-I21' in celllist:
            celllist.remove('AG-418-I21')
        if 'AG-347-I21' in celllist:
            celllist.remove('AG-347-I21')
            
    ntitle=args.ntitle
    if ntitle==0:
        titlegroupname='HLII'#if titlegroupname==HLII, the matrices have the same shape, but the Kashtan values are all uncovered
        
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
        #outfilewhole='Matrices_%s_NONSYN_SECONDNT%s' %(titlegroupname,suffix)
        #outfilefrac='Matrices_%s_NONSYN_SECONDNT%s' %(titlegroupname,suffixfrac)
        outfilewhole='Matrices_Berube_Kash_NONSYN_SECONDNT_Discrete_Including_Non_Core'
        outfilefrac=''
    if tfna=='a':
        #outfilewhole='Matrices_%s_All_NT%s' %(titlegroupname,suffix)
        #outfilewhole='Matrices_Berube_Kash_All_NT_CheckingMC20_No6fold.npz'
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
    if prog==1:#histogram
        infilefullgene=outfilewhole

        distancelow=args.distancelow#1#100
        distancehigh=args.distancehigh#2#1000
        
        print('R: distance low %i, distance high %i' %(distancelow,distancehigh))
        make_R_histogram(infilefullgene,celllist,groupname,args.sherlockbool,npairs,tfna,distancelow=distancelow,distancehigh=distancehigh)

    if prog==2:
        infilefullgene=outfilewhole
        R_all_Pairs(infilefullgene,celllist,groupname,args.sherlockbool,tfna,startgene=startgene)

    if prog==3:
        rcut=0.7
        #process_R_allpairs(groupname,tfna,rcut,paircut,corebool)
        #get_gene_annotation()
        get_core_annotations()
    if prog==4:
        get_r_with_intervals(outfilewhole,celllist,sepbins,npairs,groupname,args.sherlockbool,discretebool,tfna)

    if prog==5:#fraction pt
        if tfna=='aa':
            infilerandom='Cell_Pair_Matrices/Each_Cell_Pair_Sites/Matrices_Berube_Kash_AA_Discrete_Including_Non_Core_Frac2_QUESTION_FRANDOM.npy'
            infile='Cell_Pair_Matrices/Each_Cell_Pair_Sites/Matrices_Berube_Kash_AA_Discrete_Including_Non_Core_Frac2_QUESTION.npy'

            print('AA data')
            get_gene_fraction_point_test(frac,celllist,infile)
            print('AA random')
            get_gene_fraction_point_test(frac,celllist,infilerandom)
            #get_gene_fraction_point(frac,celllist,infile,infilerandom)


    if prog==6:
        windowsize=args.windowsize
        #R_sliding_window(celllist,groupname,windowsize,outfilewhole,tfna)
        plot_sliding_window(groupname,tfna,windowsize)
'''
-rw-r--r-- 1 apaps dsfisher          0 May 22 13:25 HLII_R_Greater_Than_0.6_aa_withNum.csv
-rw-r--r-- 1 apaps dsfisher          0 May 22 13:25 BB_NoClose_RPRIME_Greater_Than_0.6_aa_withNum.pdf
-rw-r--r-- 1 apaps dsfisher          0 May 22 13:25 HLII_RPRIME_Greater_Than_0.6_aa_withNum.pdf
-rw-r--r-- 1 apaps dsfisher          0 May 22 13:25 HLII_NoClose_RPRIME_Greater_Than_0.6_aa_withNum.pdf
-rw-r--r-- 1 apaps dsfisher     213542 May 22 13:23 HLII_NoClose_R_Greater_Than_0.6_aa_withNum.csv
-rw-r--r-- 1 apaps dsfisher     662796 Apr 23 19:33 HLII_NoClose_R_Greater_Than_0.6_a_withNum.csv
-rw-r--r-- 1 apaps dsfisher    3138860 Apr 23 19:32 BB_no_6_outliers_R_Greater_Than_0.6_a_withNum.pdf
-rw-r--r-- 1 apaps dsfisher          0 Apr 23 19:32 HLII_R_Greater_Than_0.6_a_withNum.csv
-rw-r--r-- 1 apaps dsfisher  198565814 Apr 23 19:31 HLII_NoClose_R_Greater_Than_0.6_aa_withNum.pdf
-rw-r--r-- 1 apaps dsfisher      12827 Apr 23 19:31 BB_NoClose_R_Greater_Than_0.6_aa_withNum.csv
-rw-r--r-- 1 apaps dsfisher      87365 Apr 23 19:31 BB_NoClose_R_Greater_Than_0.6_a_withNum.csv
-rw-r--r-- 1 apaps dsfisher    5091148 Apr 23 18:55 BB_NoClose_R_Greater_Than_0.6_aa_withNum.pdf
-rw-r--r-- 1 apaps dsfisher  690887421 Apr 23 18:51 HLII_NoClose_R_Greater_Than_0.6_a_withNum.pdf
-rw-r--r-- 1 apaps dsfisher 1279495293 Apr 23 18:51 HLII_R_Greater_Than_0.6_a_withNum.pdf
-rw-r--r-- 1 apaps dsfisher       6265 Apr 23 18:50 BB_no_6_outliers_R_Greater_Than_0.6_a_withNum.csv
-rw-r--r-- 1 apaps dsfisher   26822968 Apr 23 18:50 BB_NoClose_R_Greater_Than_0.6_a_withNum.pdf
'''



'''
sept 19 from slurm

Namespace(corebool=1, discretebool=1, frac=2, n=0, npairs=800000, ntitle=-1, prog=1, sherlockbool=1, startgene=0, tfna='tf')
Seps 100- 1000: 503099 gene pairs


Counts (data, random)
[array([2.0000e+00, 2.0000e+00, 3.0000e+00, 0.0000e+00, 9.0000e+00,
       7.0000e+00, 1.1000e+01, 1.7000e+01, 2.7000e+01, 5.2000e+01,
       7.2000e+01, 1.4200e+02, 2.6300e+02, 5.6500e+02, 1.1980e+03,
       2.2040e+03, 4.2560e+03, 7.6090e+03, 1.2187e+04, 1.6786e+04,
       2.1431e+04, 2.5592e+04, 2.9428e+04, 3.2445e+04, 3.4697e+04,
       3.5834e+04, 3.6217e+04, 3.5209e+04, 3.3107e+04, 3.0724e+04,
       2.7729e+04, 2.4111e+04, 2.0615e+04, 1.7079e+04, 1.3803e+04,
       1.0845e+04, 8.4940e+03, 6.2910e+03, 4.6390e+03, 3.1420e+03,
       2.0510e+03, 1.3190e+03, 8.1200e+02, 4.9600e+02, 2.8100e+02,
       1.6200e+02, 7.6000e+01, 4.6000e+01, 1.0000e+01, 9.0000e+00]), array([3.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00, 1.00000e+00,
       4.00000e+00, 9.00000e+00, 1.20000e+01, 2.00000e+01, 3.70000e+01,
       5.10000e+01, 7.60000e+01, 1.07000e+02, 1.67000e+02, 3.44000e+02,
       2.08800e+03, 7.85800e+04, 3.19703e+05, 9.72120e+04, 2.82900e+03,
       3.36000e+02, 2.03000e+02, 1.16000e+02, 8.20000e+01, 4.70000e+01,
       3.30000e+01, 2.40000e+01, 9.00000e+00, 4.00000e+00, 2.00000e+00,
       2.00000e+00, 1.00000e+00, 0.00000e+00, 1.00000e+00, 0.00000e+00,
       0.00000e+00, 1.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00])]


Namespace(corebool=1, discretebool=1, frac=2, n=0, npairs=200000, ntitle=-1, prog=1, sherlockbool=1, startgene=0, tfna='tf')
100 1000 200000


Counts (data, random)
[array([0.0000e+00, 0.0000e+00, 0.0000e+00, 1.0000e+00, 0.0000e+00,
       4.0000e+00, 5.0000e+00, 5.0000e+00, 5.0000e+00, 1.1000e+01,
       1.4000e+01, 2.5000e+01, 3.8000e+01, 8.7000e+01, 1.6100e+02,
       3.8100e+02, 7.3000e+02, 1.2260e+03, 2.3030e+03, 3.9740e+03,
       6.0350e+03, 7.6160e+03, 9.4980e+03, 1.1139e+04, 1.2488e+04,
       1.3565e+04, 1.4356e+04, 1.4443e+04, 1.4272e+04, 1.3659e+04,
       1.2707e+04, 1.1643e+04, 1.0040e+04, 8.7530e+03, 7.3400e+03,
       5.9220e+03, 4.7990e+03, 3.5960e+03, 2.7640e+03, 2.1230e+03,
       1.3980e+03, 9.5100e+02, 6.6100e+02, 3.8500e+02, 2.2100e+02,
       1.2000e+02, 7.0000e+01, 3.8000e+01, 2.1000e+01, 7.0000e+00]), array([1.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
       1.0000e+00, 1.0000e+00, 3.0000e+00, 5.0000e+00, 7.0000e+00,
       7.0000e+00, 1.2000e+01, 3.0000e+01, 4.1000e+01, 5.8000e+01,
       1.2100e+02, 2.5100e+02, 7.8040e+03, 9.6712e+04, 8.8155e+04,
       5.9710e+03, 1.9200e+02, 8.2000e+01, 5.1000e+01, 3.6000e+01,
       1.9000e+01, 1.9000e+01, 9.0000e+00, 4.0000e+00, 3.0000e+00,
       1.0000e+00, 2.0000e+00, 1.0000e+00, 0.0000e+00, 0.0000e+00,
       0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 1.0000e+00,
       0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
       0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00])]





Namespace(corebool=1, discretebool=1, frac=2, n=0, npairs=200000, ntitle=-1, prog=1, sherlockbool=1, startgene=0, tfna='tf')
frac isn't relevant
unsure if R or RPRIE

array([3.0000e+00, 2.0000e+00, 0.0000e+00, 0.0000e+00, 2.0000e+00,
       2.0000e+00, 6.0000e+00, 8.0000e+00, 1.3000e+01, 1.5000e+01,
       2.8000e+01, 5.6000e+01, 1.0700e+02, 2.3400e+02, 4.8000e+02,
       8.5300e+02, 1.6580e+03, 3.0270e+03, 4.9010e+03, 6.6720e+03,
       8.4800e+03, 1.0258e+04, 1.1710e+04, 1.2845e+04, 1.4046e+04,
       1.4160e+04, 1.4510e+04, 1.4051e+04, 1.3081e+04, 1.2074e+04,
       1.1077e+04, 9.5010e+03, 8.4010e+03, 6.7010e+03, 5.4810e+03,
       4.2440e+03, 3.3280e+03, 2.4550e+03, 1.7890e+03, 1.2080e+03,
       8.0300e+02, 5.1200e+02, 3.5400e+02, 2.2100e+02, 9.4000e+01,
       6.3000e+01, 3.3000e+01, 2.9000e+01, 2.0000e+00, 6.0000e+00]), array([0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       1.00000e+00, 2.00000e+00, 6.00000e+00, 7.00000e+00, 7.00000e+00,
       2.20000e+01, 3.70000e+01, 5.00000e+01, 7.40000e+01, 1.24000e+02,
       8.31000e+02, 3.14120e+04, 1.26594e+05, 3.89430e+04, 1.14000e+03,
       1.37000e+02, 8.50000e+01, 3.50000e+01, 3.10000e+01, 2.20000e+01,
       9.00000e+00, 8.00000e+00, 1.00000e+00, 3.00000e+00, 3.00000e+00,
       0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00])]
Bins
[-0.50673281 -0.47784431 -0.44895581 -0.4200673  -0.3911788  -0.3622903
 -0.33340179 -0.30451329 -0.27562478 -0.24673628 -0.21784778 -0.18895927
 -0.16007077 -0.13118227 -0.10229376 -0.07340526 -0.04451675 -0.01562825
  0.01326025  0.04214876  0.07103726  0.09992576  0.12881427  0.15770277
  0.18659128  0.21547978  0.24436828  0.27325679  0.30214529  0.33103379
  0.3599223   0.3888108   0.41769931  0.44658781  0.47547631  0.50436482
  0.53325332  0.56214182  0.59103033  0.61991883  0.64880733  0.67769584
  0.70658434  0.73547285  0.76436135  0.79324985  0.82213836  0.85102686
  0.87991536  0.90880387  0.93769237]



Namespace(corebool=1, discretebool=1, frac=2, n=0, npairs=20000, ntitle=-1, prog=1, sherlockbool
=1, startgene=0, tfna='tf')


not log PRIME Counts (data, random)
[array([1.000e+00, 0.000e+00, 4.000e+00, 1.000e+00, 1.200e+01, 1.600e+01,
       3.000e+01, 5.500e+01, 9.900e+01, 1.610e+02, 2.380e+02, 3.650e+02,
       5.710e+02, 7.980e+02, 1.003e+03, 1.260e+03, 1.595e+03, 1.727e+03,
       1.794e+03, 1.688e+03, 1.699e+03, 1.489e+03, 1.292e+03, 1.016e+03,
       8.460e+02, 6.310e+02, 4.720e+02, 3.350e+02, 2.500e+02, 1.480e+02,
       1.150e+02, 8.100e+01, 3.700e+01, 3.300e+01, 2.500e+01, 1.600e+01,
       1.000e+01, 1.200e+01, 1.000e+01, 7.000e+00, 3.000e+00, 4.000e+00,
       1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 0.000e+00, 0.000e+00,
       0.000e+00, 1.000e+00]), array([0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
       0.0000e+00, 2.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
       1.0000e+00, 0.0000e+00, 3.0000e+00, 4.0000e+00, 7.0000e+00,
       7.0000e+00, 2.3000e+01, 3.4000e+02, 4.7590e+03, 1.0918e+04,
       3.6070e+03, 2.4000e+02, 1.6000e+01, 1.0000e+01, 4.0000e+00,
       3.0000e+00, 3.0000e+00, 3.0000e+00, 1.0000e+00, 0.0000e+00,
       1.0000e+00, 1.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
       0.0000e+00, 0.0000e+00, 0.0000e+00, 1.0000e+00, 0.0000e+00,
       0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
       0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00])]
Bins
[-0.45850159 -0.43490223 -0.41130287 -0.38770351 -0.36410414 -0.34050478
 -0.31690542 -0.29330606 -0.2697067  -0.24610734 -0.22250798 -0.19890861
 -0.17530925 -0.15170989 -0.12811053 -0.10451117 -0.08091181 -0.05731244
 -0.03371308 -0.01011372  0.01348564  0.037085    0.06068436  0.08428373
  0.10788309  0.13148245  0.15508181  0.17868117  0.20228053  0.22587989
  0.24947926  0.27307862  0.29667798  0.32027734  0.3438767   0.36747606
  0.39107543  0.41467479  0.43827415  0.46187351  0.48547287  0.50907223
  0.5326716   0.55627096  0.57987032  0.60346968  0.62706904  0.6506684
  0.67426777  0.69786713  0.72146649]


'''


'''
AA fraction pt, doesn't seem right

Fraction Gene halves R 0.036237
AA random
Fraction Gene halves R -0.016793
'''
