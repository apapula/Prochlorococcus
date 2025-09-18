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
this script uses the pre-computed distance matrices to calculate linkage of gene pair divergences along the (highly syntenic) genome, here using the measure we call R tilde (or r prime): the pearson correlation coefficient of cell pair divergences at a pair of genes, with each cell pair gene divergence normalized by that cell pair's mean divergence at mutually-covered core genes

Also included in this script are other ways of plotting this metric R tilde/R prime, such as a plot of the average gene pair correlation within a sliding window of variable size
'''

bigcluster_BB_6_100= [11, 12, 14, 21, 23, 25, 35, 37, 70, 71, 77, 78, 103, 108, 111, 125, 143, 147, 148, 149, 150, 167, 172, 261, 267, 268, 269, 270, 274, 285, 287, 289, 290, 291, 292, 293, 296, 314, 315, 328, 334, 335, 337, 348, 351, 352, 354, 412, 413, 451, 482, 490, 491, 503, 504, 510, 511, 512, 514, 534, 536, 537, 538, 541, 568, 587, 588, 590, 595, 604, 622, 683, 804, 807, 808, 809, 810, 814, 815, 817, 825, 898, 904, 905, 906, 907, 908, 910, 913, 914, 945, 949, 950, 1025, 1041, 1042, 1044, 1046, 1047, 1048, 1050, 1055, 1058, 1067, 1072, 1091, 1154, 1173, 1174, 1175, 1176, 1183, 1185, 1188, 1189, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1199, 1219, 1220, 1221, 1224, 1271, 1276, 1287, 1303, 1347, 1373, 1379, 1380, 1381, 1382, 1383, 1386, 1387, 1388, 1389, 1469, 1475, 1481, 1502, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1625, 1627, 1628, 1629, 1630, 1631, 1633, 1673, 1674, 1675, 1676, 1677, 1678, 1680, 1681, 1682, 1683, 1690, 1692, 1693, 1694, 1695, 1696, 1697, 1698, 1699, 1700, 1728, 1743, 1765, 1766, 1767, 1792, 1796, 1797, 1800, 1801, 1802, 1803, 1804, 1817, 1820, 1821, 1822, 1823, 1825, 1826, 1827, 1828, 1829, 1830, 1831, 1838, 1878, 1903, 1904, 1905, 1907]


def get_rPRIME_with_intervals(infile_fullgene,celllist,sepbins,npairs,groupname,sherlockbool,discretebool,tfna,infilewg):
    wgarray,keylist=load_arrays(infilewg)
    wgmat=wgarray['MeanWG_%s' %tfna]

    sep_midpoints,rs,r_1_up,r_2_up,r_1_down,r_2_down=get_rPRIME_and_separation_intervals(sepbins,npairs,celllist,infile_fullgene,discretebool,tfna,wgmat)
    
    print('seps=',sep_midpoints)
    print('RPRIMEnew=',rs)
    
    print('RPRIMEnew_1sigma_up=',r_1_up)
    print('RPRIMEnew_2sigma_up=',r_2_up)
    print('RPRIMEnew_1sigma_down=',r_1_down)
    print('RPRIMEnew_2sigma_down=',r_2_down)
    

def plot_rPRIME_and_random(infile_fullgene,infilefrac,infilefracrandom,frac,celllist,sepbins,npairs,groupname,sherlockbool,discretebool,tfna,infilewg):
    
   
    #for each sep in sepbins, get the sample of genes to use. convert to key names.  Then for each pair in the keypairs, get r. average for the sep bin. plot. ALSO DO THE RANDOM ONE!

    wgarray,keylist=load_arrays(infilewg)
    wgmat=wgarray['MeanWG_%s' %tfna]
    sep_midpoints,rs,shuffled_rs=get_rPRIME_and_separation_for_fullgene_pairs(sepbins,npairs,celllist,infile_fullgene,discretebool,tfna,wgmat)
    
    if sherlockbool==1:
        print('seps=',sep_midpoints)
        print('RPRIMEnew=',rs)
        print('RPRIMEnew_random=',shuffled_rs)
    if sherlockbool==0:
        plt.plot(0.5,genehalfr,'o',color='b',label='Data')
        plt.plot(0.5,genehalfr_random,'o',color='red',label='Random Sites')
        plt.plot(sep_midpoints,rs,'o',color='b')
        plt.plot(sep_midpoints,shuffled_rs,'o',color='orange',label='Random Shuffle')
        plt.legend()
        plt.xlabel('Separation (Number of Genes)')
        plt.ylabel('Average R PRIME for Gene Pairs')
        plt.title('Average Correlation between cell pair divergences at gene pairs\n%s (%i Gene Pairs/Wide Bin)' %(groupname,npairs))
        plt.show()


def get_rPRIME_and_separation_intervals(sepbins,npairs,celllist,arrayfile,discretebool,tfna,wgmat):
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
        r,r1u,r2u,r1d,r2d=do_single_separation_bin_PRIME_interval(sepbins[i],sepbins[i+1],genearrays,npairs,overalllist,celllist,discretebool,wgmat)
        if r!='-':
            rs.append(r)
            r_1_up.append(r1u)
            r_1_down.append(r1d)
            r_2_down.append(r2d)
            r_2_up.append(r2u)
   
    return sep_midpoints,rs,r_1_up,r_2_up,r_1_down,r_2_down

def get_rPRIME_and_separation_for_fullgene_pairs(sepbins,npairs,celllist,arrayfile,discretebool,tfna,wgmat):
    
    overalllist=overall_order_list()
    if tfna in ['t','f','n','a','aa']:
        genearrays,mykeys=load_arrays(arrayfile)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
    sep_midpoints=[]#if sepbins is [1,2,3,4,10], want result to be [1,2,3,6.5]. separations are [seplower, sepupper-1]
    rs=[]
    shuffled_rs=[]
    for i in range(len(sepbins)-1):
        midpt=0.5*(sepbins[i]+sepbins[i+1]-1)
        sep_midpoints.append(midpt)
        r,shuff_r=do_single_separation_bin_PRIME(sepbins[i],sepbins[i+1],genearrays,npairs,overalllist,celllist,discretebool,wgmat)
        if r!='-':
            rs.append(r)
            shuffled_rs.append(shuff_r)
    #print(sep_midpoints)
    #print(rs)
    return sep_midpoints,rs,shuffled_rs



def do_single_separation_bin_PRIME_interval(seplow,sephigh,genearrays,npairs,overalllist,celllist,discretebool,wgmat):
   
    mykeys=list(genearrays.keys())
    keys1,keys2=get_random_pair_sample(npairs,mykeys,seplow,sephigh)
    my_rs=[]
    for i in range(len(keys1)):
        r,slope=get_rPRIME_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,wgmat)
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

def do_single_separation_bin_PRIME(seplow,sephigh,genearrays,npairs,overalllist,celllist,discretebool,wgmat):
    mykeys=list(genearrays.keys())
    keys1,keys2=get_random_pair_sample(npairs,mykeys,seplow,sephigh)
    my_rs=[]
    shuffled_rs=[]
    for i in range(len(keys1)):
        r,slope=get_rPRIME_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,wgmat)
        if r!='-':
            my_rs.append(r)
        #now need to shuffle
        sr,sslope=get_rPRIME_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,wgmat,shufflebool=1)
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

def get_all_pairs(mykeys,distancelow,distancehigh):
    k1,k2=[],[]
    for i in range(len(mykeys)):
        numi=get_gene_from_key_to_check(mykeys[i])
        for j in range(i):
            numj=get_gene_from_key_to_check(mykeys[j])
            num1=max(numi,numj)
            num2=min(numi,numj)
            mydist=num1-num2
            if mydist>1907./2:
                mydist=1907-mydist
            if mydist>=distancelow and mydist<=distancehigh:
                k1.append(mykeys[i])
                k2.append(mykeys[j])
                print(num1,num2,mydist)
    return k1,k2


    
        
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


def get_rPRIME_gene_pair(mat1,mat2,overalllist,celllist,key1,key2,discretebool,wgmat,shufflebool=0,sherlockbool=1,pdfbool=0,pdf='',tfna=''):
    divs1=[]
    divs2=[]
    

    for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            [wg,ngenes]=wgmat[idxi][idxj]
            if not discretebool:
                d1=mat1[idxi][idxj]
                d2=mat2[idxi][idxj]
                
                if d1<1 and d2<1:
                    divs1.append(d1/float(wg))
                    divs2.append(d2/float(wg))
            if discretebool:
                [n1,l1]=mat1[idxi][idxj]
                [n2,l2]=mat2[idxi][idxj]
                if l1<12 or l2<12:
                    continue
                d1=float(n1)/l1
                d2=float(n2)/l2
                divs1.append(d1/float(wg))
                divs2.append(d2/float(wg))
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

def make_RPRIME_histogram(infilefullgene,celllist,groupname,sherlockbool,npairs,tfna,infilewg,distancelow=100,distancehigh=1000):
    wgarray,keylist=load_arrays(infilewg)
    wgmat=wgarray['MeanWG_%s' %tfna]
    tfnastring=repo.return_tfna_string(tfna)
    overalllist=overall_order_list()
    r,rshuf=[],[]
    if tfna!='tf':
        genearrays,mykeys=load_arrays(infilefullgene)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
        mykeys=list(genearrays.keys())
    if npairs<50000:
        keys1,keys2=get_random_pair_sample(npairs,mykeys,distancelow,distancehigh)
    if npairs>=50000:
        keys1,keys2=get_all_pairs(mykeys,distancelow,distancehigh)
    
    for  i in range(len(keys1)):
        if i%100==0:
            print('%f (%i of %i) ' %(float(i)/len(keys1),i,len(keys1)))
        myr,slope=get_rPRIME_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,wgmat,sherlockbool=sherlockbool)
        if myr!='-':
            r.append(myr)
        #now need to shuffle
        mysr,sslope=get_rPRIME_gene_pair(genearrays[keys1[i]],genearrays[keys2[i]],overalllist,celllist,keys1[i],keys2[i],discretebool,wgmat,shufflebool=1)
        if mysr!='-':
            rshuf.append(mysr)
            
    
    myn,mybins,mypatches=plt.hist([r,rshuf],histtype='step',bins=50,label=['Data','Shuffle'])
    
    plt.legend()
    plt.ylabel('Number of Gene Pairs')
    plt.xlabel('Pearson R PRIME')
    ax=plt.gca()
    #ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    plt.title('%s R PRIME for %s Cell Pair Divergences\nGenes Separated %i-%i Genes, %i Pairs' %(groupname,tfnastring,distancelow,distancehigh,npairs))
    plt.savefig('%s_%s_RPRIME_Histograms_AllSites_%ipairs.png' %(groupname,tfna,npairs))
    plt.close()

    print('not log PRIME Counts (data, random)')
    print(myn)
    print('Bins')
    print(mybins)
    

    plt.hist([r,rshuf],histtype='step',bins=50,label=['Data','Shuffle'])
    
    plt.legend()
    plt.ylabel('Number of Gene Pairs')
    plt.xlabel('Pearson R PRIME')
    ax=plt.gca()
    plt.yscale('log')
    #ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    plt.title('%s R PRIME for %s Cell Pair Divergences\nGenes Separated %i-%i Genes, %i Pairs' %(groupname,tfnastring,distancelow,distancehigh,npairs))
    plt.savefig('%s_%s_RPRIME_Histograms_AllSites_%ipairs_log.png' %(groupname,tfna,npairs))
    plt.close()



def RPRIME_all_Pairs(infilefullgene,celllist,groupname,sherlockbool,tfna,infilewg,startgene=0):
    wgarray,keylist=load_arrays(infilewg)
    wgmat=wgarray['MeanWG_%s' %tfna]

    print(infilefullgene)
    
    tfnastring=repo.return_tfna_string(tfna)
    overalllist=overall_order_list()
    pdf=PdfPages('%s_RPRIME_Greater_Than_0.6_%s_withNum.pdf' %(groupname,tfna))
    r,rshuf=[],[]
    if tfna!='tf':
        genearrays,mykeys=load_arrays(infilefullgene)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
        mykeys=list(genearrays.keys())

    
    sample_genenums=repo.make_sample_all_single_genes()
    #sample_genenums=bigcluster_BB_6_100

    ngenepairs=0
    
    if 1:
        rps=[]
        for i in range(len(sample_genenums)):
            if i%100==0:
                print('%i of %i' %(i,len(sample_genenums)))
           
            genei=sample_genenums[i][0]#[i] if using bigcluster_BB_6_100
            print(genei)
            if genei<startgene:
                continue
            if tfna!='aa':
                keyi='Gene_%i_%s_Pairs' %(genei,tfna)
            if tfna=='aa':
                keyi='Gene_%i_AA' %(genei)
            for j in range(i):
                genej=sample_genenums[j][0]#j]#[0]
                if tfna!='aa':
                    keyj='Gene_%i_%s_Pairs' %(genej,tfna)
                if tfna=='aa':
                    keyj='Gene_%i_AA' %(genej)
                myr,slope,numcells=get_rPRIME_gene_pair(genearrays[keyi],genearrays[keyj],overalllist,celllist,keyi,keyj,discretebool,wgmat,sherlockbool=sherlockbool,pdfbool=1,pdf=pdf,tfna=tfna)
                if myr=='-':
                    continue
                rps.append(myr)
        myn,newbins,mypatches=plt.hist(rps,histtype='step',bins=50)
        plt.close()

        print('PRIME Counts (data only)')
        print(myn)
        print('Bins')
        print(list(newbins))
        
    if 0:
        outfile='%s_RPRIME_Greater_Than_0.6_%s_withNum.csv' %(groupname,tfna)
        fieldnames=['Gene1','Gene2','RPRIME','NPairs']
        with open(outfile,'w') as myoutfile:
            outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
            outwriter.writeheader()

            for i in range(len(sample_genenums)):
                if i%100==0:
                    print('%i of %i' %(i,len(sample_genenums)))

                genei=sample_genenums[i][0]#[i] if using bigcluster_BB_6_100
                print(genei)
                if genei<startgene:
                    continue
                if tfna!='aa':
                    keyi='Gene_%i_%s_Pairs' %(genei,tfna)
                if tfna=='aa':
                    keyi='Gene_%i_AA' %(genei)
                for j in range(i):
                    genej=sample_genenums[j][0]#j]#[0]
                    if tfna!='aa':
                        keyj='Gene_%i_%s_Pairs' %(genej,tfna)
                    if tfna=='aa':
                        keyj='Gene_%i_AA' %(genej)
                    myr,slope,numcells=get_rPRIME_gene_pair(genearrays[keyi],genearrays[keyj],overalllist,celllist,keyi,keyj,discretebool,wgmat,sherlockbool=sherlockbool,pdfbool=1,pdf=pdf,tfna=tfna)
                    if myr=='-':
                        continue
                    ngenepairs+=1
                    if myr>0.60:
                        outwriter.writerow({'Gene1':genei,'Gene2':genej,'RPRIME':myr,'NPairs':numcells})

        print('R PRIME %s %s %i total gene pairs' %(groupname,tfna,ngenepairs))
        pdf.close()
    

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
        #print(idx1,idx2)
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
##########################################################################################

def gene_cellpair_variance_distribution(groupname,celllist,tfna,infilewg,arrayfile):
    
    wgarray,keylist=load_arrays(infilewg)
    wgmat=wgarray['MeanWG_%s' %tfna]
    variances,variancesprime=[],[]
    overalllist=overall_order_list()
    if tfna in ['t','f','n','a','aa']:
        genearrays,mykeys=load_arrays(arrayfile)
    if tfna=='tf':
        genearrays=load_two_arrays(infilematrices)
    for item in list(genearrays.keys()):
        print(item)
        mat=genearrays[item]
        divs=[]
        divsprime=[]
        for i in range(len(celllist)):
            idxi=overalllist.index(celllist[i])
            for j in range(i):
                idxj=overalllist.index(celllist[j])
                [wg,ngenes]=wgmat[idxi][idxj]
                [n1,l1]=mat[idxi][idxj]
                if l1<12:
                    continue
                d1=float(n1)/l1
                d2=d1/float(wgmat[idxi][idxj][0])
                divs.append(d1)
                divsprime.append(d2)
        #now get variances
        if len(divs)<10:
            continue
        variances.append(statistics.variance(divs))
        variancesprime.append(statistics.variance(divsprime))
    plt.hist(variances,bins=200)
    plt.xlabel('Cell Pair Divergence Variance per gene')
    plt.ylabel('Number of Core genes')
    plt.title('%s %s Cell Pair Divergence Variances per Gene\nMean %f Variance %f' %(groupname,tfna,repo.avelist(variances),statistics.variance(variances)))
    plt.show()

    plt.hist(variancesprime,bins=200)
    plt.xlabel('Normalized R PRIME Cell Pair Divergence Variance per gene')
    plt.ylabel('Number of Core genes')
    plt.title('%s %s R PRIME Normalized Cell Pair Divergence Variances per Gene\nMean %f Median %f Variance %f' %(groupname,tfna,repo.avelist(variancesprime),statistics.median(variancesprime),statistics.variance(variancesprime)))
    plt.show()


    plt.hist([item for item in variancesprime if item <1.8],bins=40)
    plt.xlabel('Normalized R PRIME Cell Pair Divergence Variance per gene')
    plt.ylabel('Number of Core genes')
    plt.title('%s %s R PRIME Normalized Cell Pair Divergence Variances per Gene\nZoomed' %(groupname,tfna))
    plt.show()
#############


def RPRIME_sliding_window(celllist,groupname,windowsize,infilefullgene,tfna,infilewg):
    print('%s RPRIME sliding window size %i %s' %(groupname,windowsize,tfna))
    wgarray,keylist=load_arrays(infilewg)
    wgmat=wgarray['MeanWG_%s' %tfna]
    
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
                myr,slope,numcells=get_rPRIME_gene_pair(genearrays[keyi],genearrays[keyj],overalllist,celllist,keyi,keyj,discretebool,wgmat,sherlockbool=1,pdfbool=0,pdf='',tfna=tfna)
                if myr=='-':
                    continue
                rlist.append(myr)
                ncellslist.append(numcells)
        aver.append(repo.avelist(rlist))
        avencells.append(repo.avelist(ncellslist))
        windowstart.append(l[0])
    #now write

    outfields=['GeneStart','Ave_RPRIME','Ave_NCells','Genes']
    with open('%s_RPRIME_Sliding_Window_%i_%s.csv' %(groupname,windowsize,tfna),'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,fieldnames=outfields,delimiter='\t')
        outwriter.writeheader()
        for i in range(len(genelists)):
            outwriter.writerow({'GeneStart':windowstart[i],'Ave_RPRIME':aver[i],'Ave_NCells':avencells[i],'Genes':genelists[i]})

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

    
    with open('../input_data/%s_RPRIME_Sliding_Window_%i_%s.csv' %(groupname,windowsize,tfna),'r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            aver.append(float(row['Ave_RPRIME']))
            genestart.append(int(row['GeneStart']))
            ncells.append(float(row['Ave_NCells']))
            if 1:
                if aver[-1]>0.35:
                    if len(aver)>1 and aver[-2]>0.35:##added second to reduce clutter
                        continue
                    print('Starting gene %i aver %f' %(genestart[-1],aver[-1]))
                    glist=ast.literal_eval(row['Genes'])
                    l=[annotations[glist[k]] for k in range(len(glist))]
                    print(l)
                    
    coverage_plot(genestart,aver,ncells,5000)
    return

    #now plot
    plt.plot(genestart,aver,'.-')#,alpha=0.3)
    plt.xlabel('Position of Starting Gene on MIT9301 Reference',fontsize=20)
   
    xtickpos=[0,500,1000,1500,2000]
    plt.xticks(ticks=xtickpos)
     
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("$\\tilde{R}$ Averaged over %i-Gene Windows\nAlong MIT9301 Reference" %(windowsize),fontsize=20)
    # plt.title('%s R PRIME %s %i-Gene windows\nAlong MIT9301 Reference' %(groupname,tfna,windowsize),fontsize=20)
    # plt.ylabel('Average R PRIME of %i Genes' %(windowsize),fontsize=20)
    plt.ylabel("Average R' of %i Genes" %(windowsize),fontsize=20)
    plt.tight_layout()
    plt.show()

    return


    fig=plt.figure()
    ax1=plt.subplot(211)
    ax2=plt.subplot(212,sharex=ax1)
    ax1.plot(genestart,aver,'.-')
    ax2.plot(genestart,ncells,'.-')
    #ax2.sharex(ax1)
    plt.xlabel('Position of Starting Gene on MIT9301 Reference',fontsize=20)
    ax1.set_ylabel('Average R PRIME of %i Genes' %(windowsize),fontsize=20)
    ax2.set_ylabel('Average number of cell pairs covered')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.suptitle('%s RPRIME  %s %i-Gene windows\nAlong MIT9301 Reference' %(groupname,tfna,windowsize),fontsize=20)
    plt.tight_layout()
    plt.show()
    
    plt.hist(aver,bins=100)
    plt.xlabel('Average R PRIME for %i-gene window' %windowsize,fontsize=20)
    plt.ylabel('Number of %i-Gene windows' %windowsize,fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('Distribution of average R PRIME over %i-gene\nwindows for %s %s' %(windowsize,groupname,tfna),fontsize=20)
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
    plt.xlabel('Position of Starting Gene on MIT9301 Reference',fontsize=20)
    plt.ylabel("Average $\\tilde{R}$ of %i Genes" %(windowsize),fontsize=20)
    plt.legend(fontsize=15)
    xtickpos=[0,500,1000,1500,2000]
    plt.xticks(ticks=xtickpos)
     
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("$\\tilde{R}$ Averaged Over %i-Gene windows\nAlong MIT9301 Reference" %(windowsize),fontsize=20)
    #plt.title('%s R %s %i-Gene windows\nAlong MIT9301 Reference' %(groupname,tfna,windowsize),fontsize=20)
    plt.tight_layout()
    plt.show()
######################################


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
    suffixfrac=''
    
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

    if args.sherlockbool==0:
        outfilewhole='../input_data/'+outfilewhole
        outfilefrac='../input_data/'+outfilefrac
        outfilefracrandom='../input_data/'+outfilefracrandom
    if discretebool:
        if tfna=='tf':
            titlegroupname='Berube_Kash'
            suffix='_Discrete'
            outfilewhole='Matrices_%s_TWO_AND_FOUR%s' %(titlegroupname,suffix)
            outfilewhole='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfilewhole
            if not args.sherlockbool:
                outfilewhole='../input_data/'+outfilewhole
            infilematrices=['%s_firsthalf.npy' %outfilewhole,'%s_secondhalf.npy' %outfilewhole]
            outfilewhole=infilematrices

   
    npairs=args.npairs

    sepbins=[1,2,3,4,5,6,7,8,9,10,25,63,158,398,1000]
    if prog==0:
        infilewg='Cell_Pair_Matrices/Berube_Kash_WG_Matrices_%s.npy' %(tfna)
        print(outfilewhole)
       
        get_rPRIME_with_intervals(outfilewhole,celllist,sepbins,npairs,groupname,args.sherlockbool,discretebool,tfna,infilewg)
   

    if prog==1:
        infilewg='Cell_Pair_Matrices/Berube_Kash_WG_Matrices_%s.npy' %(tfna)
        infilefullgene=outfilewhole
        RPRIME_all_Pairs(infilefullgene,celllist,groupname,args.sherlockbool,tfna,infilewg,startgene=startgene)

  

    if prog==2:
        windowsize=args.windowsize
        infilewg='Cell_Pair_Matrices/Berube_Kash_WG_Matrices_%s.npy' %(tfna)
        plot_sliding_window(groupname,tfna,windowsize)
