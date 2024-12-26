import csv
import numpy as np
import ast
import csv
import math
from scipy import stats
from os.path import exists
import statistics
import argparse
import matplotlib.pyplot as plt
import re
import Load_Sequences_Repo as load_sequences_repo
from difflib import SequenceMatcher
from Bio import Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.cluster.hierarchy as hac
import Kashtan_Alignments_Pairwise_Event_Lengths as kashrepo
import Kashtan_Alignments_Record_Gaps_Between_Highdiv_Events as kashclusterrepo
from Cell_Lists import clique1_without_ribotype,C1_without_ribotype,clique1,cl2,cl3,cl4,cl5,cl6,hliilist2,BB_list
from scipy.stats import linregress, pearsonr
import Comparing_Close_Cliques_to_Poisson as poissonrepo
import random
import matplotlib.patches as mpatches
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

'''
font = {'family' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)
'''

'''
for the cell group, get the lambda from close cliques poisson

then, go through that groups flexible orthogroups. get n genes for each cell and n genes for each pair
C1_Orthologous_Gene_Group_Qcut40_pid40.txt 60 qcut 80 pid default?
lique1_Pairs_Blasting_withnum.out
check these give same numner. below: they do not, the guess is that the latter was made using hlii core?

20367 nonc1 core orfs total? All_nonC1Core_ORFs.fasta
find 19261 total orfs?! in C1_orthologou_gee_group_60q 80p
C1 core genes: 1454 according to pid 80 synon; 1458, 1580, 1582 using pid 40. so this checks out? could redo the clique1 stats using c1 core
'''

Atlantic=['355','363','388','412','429','432','436']
N_Pacific=['347','402']
S_Pacific= ['335','459','469','449','455']
Pacific=N_Pacific+S_Pacific

atlanticlist=['355','363','388','412','429','432','436']+['409','418','420','424','442','444']
pacificlist=['347','402','335','459','469','449','455']+['311','315','316','321','323','331','341','345','450','463','670','673','676','679','683','686']


def plot_expected_number_flex_genes():
    celllist=C1_without_ribotype
    groupname='C1'
    infile='../input_data/C1_Core_Coverage_using_1581_Core_Genes.csv'
    nflex=[]
    qcut,pidcut=80,80
    countdictsingle=count_n_genes_per_cell_from_connected_comps(groupname,qcut,pidcut,celllist)
    covdict={}
    with open(infile,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            covdict[row['Cell']]=float(row['CoreCov'])
            c=row['Cell']
            if c not in countdictsingle:
                continue
            nflex.append(float(countdictsingle[c])/float(row['CoreCov']))
    plt.hist(nflex,bins=30)
    plt.title('Estimated number of flexible genes per cell\nif fully covered (Mean %.3f Median %.1f)' %(float(sum(nflex))/len(nflex),statistics.median(nflex)))
    plt.ylabel('Number of C1 Cells')
    plt.xlabel('Number of non-C1 Core Genes')
    plt.show()

def plot_LOG_flex_genes_shared_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut):
    countdictpair=get_number_flex_genes_shared(groupname,celllist,qcut,pidcut)
    print('done1')
    lambdadict=get_poisson(celllist,groupname,pcut,chunklength)
    print(lambdadict)
    print('done2')
    countdictsingle=count_n_genes_per_cell_from_connected_comps(groupname,qcut,pidcut,celllist)

    if 1:#count vs lambda
        ll,lcp=[],[]
        for item in lambdadict:
            if lambdadict[item]=='-':
                continue
            c1=item[0]
            c2=item[1]
            mynorm=countdictsingle[c1]*float(countdictsingle[c2])
            if mynorm==0 or countdictpair[item]==0:
                print('%s %s have %i %i flex genes' %(c1,c2,countdictsingle[c1],countdictsingle[c2]))
                continue
            ll.append(lambdadict[item])
            #lcp.append(np.log(countdictpair[item]))
            lcp.append(countdictpair[item])
            if lcp[-1]>500:
                print('Pair %s %i flex genes shared?' %(item,lcp[-1]))
      
            
        if 1:
            plt.plot(ll,lcp,'o',alpha=0.3)
            plt.xlabel('Asexual Divergence',fontsize=20)
            #add_regression_line(ll,lcp)
            #add_regression_line_fixed_intercept(ll,lcp,283)
            plt.ylabel('No. of Flex. Genes Shared',fontsize=20)
            #plt.title('%s Pairs: Number of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
            ax=plt.gca()
            ax.set_xticks([0,0.002,0.004,.006])
            plt.yscale('log')
            ax.tick_params(axis='both', which='major', labelsize=25)
            plt.tight_layout()
            plt.show()

        



        
    if 1:#the expected number of shared genes is n_flex_1*n_flex_2/n_flex_fullgenome roughly
        countlistsingle=[countdictsingle[item] for item in countdictsingle if countdictsingle[item]>0]
        
        nflexfullgenome=float(sum(countlistsingle))/(len(countlistsingle)*.7)#299#from pid 80 synon 509. 299,299,443 from others. looks like pid 40 qcov 40.csv was used
        print('Using %f as expected number of flex genes per cell if complete coverage' %nflexfullgenome)
        ll,lcpn,lfpn=[],[],[]#lambda, count pair normalized, frac pair normalized
        for item in lambdadict:
            if lambdadict[item]=='-' or countdictpair[item]==0:
                continue
            c1=item[0]
            c2=item[1]
            mynorm=countdictsingle[c1]*float(countdictsingle[c2])
            if mynorm==0:
                #print('%s %s have %i %i flex genes' %(c1,c2,countdictsingle[c1],countdictsingle[c2]))
                continue
            ll.append(lambdadict[item])
            lfpn.append(countdictpair[item]*nflexfullgenome/float(mynorm))
            # lfpn.append(np.log(countdictpair[item]*nflexfullgenome/float(mynorm)))
            lcpn.append(-countdictpair[item]+float(mynorm)/nflexfullgenome)#nflexfullgenome*(1.0-lfpn[-1]))
            #n-overturned ~ nexpected-n_shared, n_expected=n1*n2/nflexfullgenome
        


        

        plt.plot(ll,lfpn,'o',alpha=0.3)
        plt.xlabel('Asexual Divergence',fontsize=20)
        plt.ylabel('Estimated Fraction of\nFlex. Genes Shared',fontsize=20)
        #plt.title('%s Pairs: Estimated Fraction of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        ax=plt.gca()
        ax.set_xticks([0,0.002,0.004,.006])
        #ax.set_yticks([0,5,10,15,20])
        plt.yscale('log')
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.tight_layout()
        plt.show()

        plt.plot(ll,lfpn,'o',alpha=0.3)
        #add_regression_line(ll,lfpn)
        #add_regression_line_fixed_intercept(ll,lfpn,1)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Estimated Fraction of Flexible Genes Shared')
        plt.title('%s Pairs: Estimated Fraction of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()

        
       
        
##############

def plot_flex_genes_shared_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut):
    countdictpair=get_number_flex_genes_shared(groupname,celllist,qcut,pidcut)
    print('done1')
    lambdadict=get_poisson(celllist,groupname,pcut,chunklength)
    print(lambdadict)
    print('done2')
    countdictsingle=count_n_genes_per_cell_from_connected_comps(groupname,qcut,pidcut,celllist)

    if 1:#count vs lambda
        ll,lcp=[],[]
        for item in lambdadict:
            if lambdadict[item]=='-':
                continue
            c1=item[0]
            c2=item[1]
            mynorm=countdictsingle[c1]*float(countdictsingle[c2])
            if mynorm==0:
                print('%s %s have %i %i flex genes' %(c1,c2,countdictsingle[c1],countdictsingle[c2]))
                continue
            ll.append(lambdadict[item])
            lcp.append(countdictpair[item])
            if lcp[-1]>500:
                print('Pair %s %i flex genes shared?' %(item,lcp[-1]))
        if 0:#same as below but without expo fit
            plt.plot(ll,lcp,'o',alpha=0.3)
            plt.xlabel('Asexual Divergence')
            plt.ylabel('Number of C1 Flexible Genes Shared')
            plt.title('%s Pairs: Number of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
            plt.show()
            
        if 1:
            plt.plot(ll,lcp,'o',alpha=0.3)
            plt.xlabel('Asexual Divergence')
            add_exponential_with_fixed_nflex(ll,lcp)
            add_exponential(ll,lcp)#add_regression_line(ll,lcp)
            
            plt.ylabel('Number of C1 Flexible Genes Shared')
            plt.title('%s Pairs: Number of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
            plt.show()

        if 1:#averaging above
            mybins=[0,.00005,.0001,.00015,.00018,.0002,.000225,.00025,.000275,.0003,.000325,.00035,.0004]
            bin_means1,bin_edges1,binnumber1=stats.binned_statistic(ll,lcp,statistic='mean',bins=mybins)         
            binmidpoints1=0.5*(bin_edges1[1:]+bin_edges1[:-1])
            plt.plot(binmidpoints1,bin_means1,'o',alpha=0.3)
            plt.xlabel('Asexual Divergence')
            plt.ylabel('Mean Number of C1 Flexible Genes Shared')
            plt.title('%s Pairs: Mean Number of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
            plt.show()



        
    if 1:#the expected number of shared genes is n_flex_1*n_flex_2/n_flex_fullgenome roughly
        countlistsingle=[countdictsingle[item] for item in countdictsingle if countdictsingle[item]>0]
        
        nflexfullgenome=float(sum(countlistsingle))/(len(countlistsingle)*.7)#299#from pid 80 synon 509. 299,299,443 from others. looks like pid 40 qcov 40.csv was used
        print('Using %f as expected number of flex genes per cell if complete coverage' %nflexfullgenome)
        ll,lcpn,lfpn=[],[],[]#lambda, count pair normalized, frac pair normalized
        for item in lambdadict:
            if lambdadict[item]=='-':
                continue
            c1=item[0]
            c2=item[1]
            mynorm=countdictsingle[c1]*float(countdictsingle[c2])
            if mynorm==0:
                #print('%s %s have %i %i flex genes' %(c1,c2,countdictsingle[c1],countdictsingle[c2]))
                continue
            ll.append(lambdadict[item])
            lfpn.append(countdictpair[item]*nflexfullgenome/float(mynorm))
            lcpn.append(-countdictpair[item]+float(mynorm)/nflexfullgenome)#nflexfullgenome*(1.0-lfpn[-1]))
            #n-overturned ~ nexpected-n_shared, n_expected=n1*n2/nflexfullgenome
        plt.plot(ll,lcpn,'o',alpha=0.3)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Estimated Number of Flexible Genes Overturned')
        plt.title('%s Pairs: Estimated Number of Flexible Genes Overturned vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()

        plt.plot(ll,lcpn,'o',alpha=0.3)
        add_exponential(ll,lcpn)#add_regression_line(ll,lcpn)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Estimated Number of Flexible Genes Overturned')
        plt.title('%s Pairs: Estimated Number of Flexible Genes Overturned vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()


        mybins=[0,.00005,.0001,.00015,.00018,.0002,.000225,.00025,.000275,.0003,.000325,.00035,.0004]
        bin_means1,bin_edges1,binnumber1=stats.binned_statistic(ll,lcpn,statistic='mean',bins=mybins)         
        binmidpoints1=0.5*(bin_edges1[1:]+bin_edges1[:-1])
        plt.plot(binmidpoints1,bin_means1,'o',alpha=0.3)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Mean Estimated Number of Flexible Genes Overturned')
        plt.title('%s Pairs: Mean Estimated Number of Flexible Genes Overturned vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()

        plt.plot(ll,lfpn,'o',alpha=0.3)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Estimated Fraction of Flexible Genes Shared')
        plt.title('%s Pairs: Estimated Fraction of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()

        plt.plot(ll,lfpn,'o',alpha=0.3)
        add_exponential_with_fixed_fracflex(ll,lfpn)
        add_exponential(ll,lfpn)#add_regression_line(ll,lfpn)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Estimated Fraction of Flexible Genes Shared')
        plt.title('%s Pairs: Estimated Fraction of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()

        
        bin_means2,bin_edges2,binnumber2=stats.binned_statistic(ll,lfpn,statistic='mean',bins=mybins)         
        binmidpoints2=0.5*(bin_edges2[1:]+bin_edges2[:-1])
        plt.plot(binmidpoints2,bin_means2,'o',alpha=0.3)
        plt.xlabel('Asexual Divergence')
        plt.ylabel('Mean Estimated Fraction of Flexible Genes Shared')
        plt.title('%s Pairs: Mean Estimated Fraction of Flexible Genes Shared vs Asexual Divergence\nPidcut %s Qcut %s Poisson Cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.show()

def add_exponential(x,y):
    #logy=[np.log(y[i]) for i in range(len(y)) if y[i]>0]
    #myx=[x[i] for i in range(len(x)) if y[i]>0]
    cutx=[x[i] for i in range(len(x)) if x[i]<=.004]
    cuty=[y[i] for i in range(len(y)) if x[i]<=.004]
    myx=np.linspace(0,max(x),100)
    popt,pcov=curve_fit(myexpo,cutx,cuty)#myx,y)
    plt.plot(myx,myexpo(myx,*popt),label='%5.3f * exp(-%5.3f x)' % tuple(popt))
    plt.legend()

def add_exponential_with_fixed_nflex(x,y):
    cutx=[x[i] for i in range(len(x)) if x[i]<=.004]
    cuty=[y[i] for i in range(len(y)) if x[i]<=.004]
    myx=np.linspace(0,max(x),100)
    popt,pcov=curve_fit(fixedexpo,cutx,cuty)#myx,y)
    plt.plot(myx,fixedexpo(myx,*popt),label='283 * exp(-%5.3f x)' % tuple(popt))
    #plt.legend()

def add_exponential_with_fixed_fracflex(x,y):
    cutx=[x[i] for i in range(len(x)) if x[i]<=.004]
    cuty=[y[i] for i in range(len(y)) if x[i]<=.004]
    myx=np.linspace(0,max(x),100)
    popt,pcov=curve_fit(fixedexpo1,cutx,cuty)#myx,y)
    plt.plot(myx,fixedexpo1(myx,*popt),label=' exp(-%5.3f x)' % tuple(popt))
    #plt.legend()
    
def fixedexpo1(x,b):
    return 1.*np.exp(-b*x)
    
def fixedexpo(x,b):
    return 283.*np.exp(-b*x)

def myexpo(x,a,b):
    return a*np.exp(-b*x)

def add_regression_line(x,y):
    newx,newy=[],[]
    for i in range(len(x)):
        if x[i]>0.004:
            continue
        newx.append(x[i])
        newy.append(y[i])
    slope, intercept, r_value, p_value, std_err = linregress(newx,newy)
    x=np.linspace(0,max(x),100)
    plt.plot(x,x*slope+intercept,label='y=%.3f x + %.3f' %(slope,intercept))
    plt.legend()

#
def add_regression_line_fixed_intercept(x,y,intercept):
    #log
    newx,newy=[],[]
    for i in range(len(x)):
        if x[i]>0.004:
            continue
        newx.append(x[i])
        newy.append(y[i]-np.log(intercept))
    newx=np.array(newx).reshape(-1,1)
    newy=np.array(newy)
    #slope, intercept, r_value, p_value, std_err = linregress(newx,newy)
    regression = LinearRegression(fit_intercept=False)
    regression.fit(newx,newy)
    slope=regression.coef_[0]
    x=np.linspace(0,max(x),100)
    plt.plot(x,x*slope+np.log(intercept),label='y=%.3f x + %.3f' %(slope,np.log(intercept)))
    plt.legend()
     
def get_number_flex_genes_shared(groupname,celllist,qcut,pidcut):
    countdict={}
    for i in range(len(celllist)):
        for j in range(i):
            ci=celllist[i]
            cj=celllist[j]
            mypair=sorted([celllist[i],celllist[j]])
            mypair=(mypair[0],mypair[1])
            countdict[mypair]=0
                          
        
    f='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    if groupname=='C1':
        groupcut=10
    if groupname=='HLII':
        groupcut=34
    if groupname=='Blue_Basin_NoClose':
        groupcut=14
    for line in lines:
        group=ast.literal_eval(line)
        if len(group)==1:
            continue
        #print(group)
    
        if len(group)>groupcut:###new april 18
            continue
        newgroup=[]
        #now find if each pair is in :(
        for c in group:
            idx=c.find('_')
            if groupname=='C1':
                cell=c[idx+1:]
            if groupname!='C1':
                cell=c[:idx]
            if cell not in celllist:
                print(cell)
                print(celllist)
                t=input('t')
                continue
            newgroup.append(cell)
        #now loop through pairs
        for i in range(len(newgroup)):
            for j in range(i):
                ci=newgroup[i]
                cj=newgroup[j]
                mypair=sorted([ci,cj])
                mypair=(mypair[0],mypair[1])
                countdict[mypair]+=1
                #print(mypair,countdict[mypair])
                #t=input('t')
    return countdict

#######################################
def check_number_flex_genes(qcut,pidcut):
    f1='../input_data/Clique1_Pairs_Blasting_withnum.out'
    f2='../input_data/C1_Orthologous_Gene_Group_Qcut60_pid80.txt'
    celllist=clique1_without_ribotype
    countdict2=count_n_genes_per_cell_from_connected_comps('C1',qcut,pidcut,celllist)
    countdict1={}
    with open(f1,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            c1=row['Cell1']
            flex1=int(row['ORFs1'])-int(row['CoreHits1'])
            if c1 in countdict1:
                if flex1!=countdict1[c1]:
                    print('Problem 1!')
                    t=input('t')
            countdict1[c1]=flex1
            c2=row['Cell2']
            flex2=int(row['ORFs2'])-int(row['CoreHits2'])
            if c2 in countdict1:
                if flex2!=countdict1[c2]:
                    print('Problem 1!')
                    t=input('t')
            countdict1[c2]=flex2
    for item in countdict1:
        print('%s: %i vs %i Flex ORFs' %(item,countdict1[item],countdict2[item]))
        


            
def count_n_genes_per_cell_from_connected_comps(groupname,qcut,pidcut,celllist):
    countdict={}
    for item in celllist:
        countdict[item]=0
    f='../input_data/C1_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    counter=0
    totlen=0
    for line in lines:
        counter+=1
        group=ast.literal_eval(line)
        if len(group)>10:
            continue###april 18
        totlen+=len(group)
        for c in group:
            idx=c.find('_')
            cell=c[idx+1:]
            if cell in countdict:
                countdict[cell]+=1
    print('Counts of flexible genes per cell %s %s %s ' %(groupname,str(qcut),str(pidcut)))
    print('%i total orfs' %(totlen))
    return countdict




def get_poisson(celllist,groupname,pcut,chunklength):#0.01 pcut
    if groupname=='C1' or '7-Cell C1 Clique':
        celllist=[item+'_C1' for item in celllist]

    lambdas=[]
    lambdadict={}
    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
    for i in range(len(celllist)):
        print('%i of %i' %(i,len(celllist)))
        for j in range(i):
            ci=celllist[i]
            cj=celllist[j]
            mypair=sorted([celllist[i],celllist[j]])
            mypair=(mypair[0][:-3],mypair[1][:-3])
            mylambda=poissonrepo.prob_comparer_get_poisson_single_pair_for_flexgenes(i,j,celllist,chunk_seqdict,chunklength,pcut)
            if mylambda=='-':
                lambdas.append('-')#keep the proper ordering here to compare to flex genes, weed out pts later
                lambdadict[mypair]='-'
                continue
            lambdas.append(mylambda/float(chunklength))
            lambdadict[mypair]=lambdas[-1]
    return lambdadict




#########################################################


def plot_n_2_3_4_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut):
    lambdadict=get_poisson(celllist,groupname,pcut,chunklength)
    countdictsingle=count_n_genes_per_cell_from_connected_comps(groupname,qcut,pidcut,celllist)
    #get doubletone, trips, quads per pair
    nlist=[2,3,4]
    countdict=get_2_3_4_per_pair(celllist,qcut,pidcut,nlist)

    #now plot
    if 1:
        for i in range(len(nlist)):
            divlist,ngenelist,fraclist=[],[],[]#fraclist is number of e.g. doubletons found/(ngene1*ngene2) not a real fraction but normalized by something, not meaningful units
            for item in lambdadict:
                if lambdadict[item]=='-':
                    continue
                c1=item[0]
                c2=item[1]
                mynorm=countdictsingle[c1]*float(countdictsingle[c2])
                if mynorm==0:
                    continue
                divlist.append(lambdadict[item])
                ngenelist.append(countdict[item][i])
                fraclist.append(float(ngenelist[-1])/mynorm)
            plt.plot(divlist,ngenelist,'o',alpha=0.3)
            plt.title('%s number of %i-ton flexible orthogroups per pair\nvs estimated asexual divergence of pair\nqcut %s pidcut %s Poisson cut %s' %(groupname,nlist[i],str(qcut),str(pidcut),str(pcut)))
            add_regression_line(divlist,ngenelist)
            
            plt.xlabel('Asexual Divergence')
            plt.ylabel('Number of %i-ton Flexible Genes Per Cell Pair' %nlist[i])
            plt.show()


            plt.plot(divlist,fraclist,'o',alpha=0.3)
            plt.title('%s Normalized (by number of flex genes of each cell) number\nof %i-ton flex genes per pair vs asex div\nqcut %s pidcut %s Poisson cut %s' %(groupname,nlist[i],str(qcut),str(pidcut),str(pcut)))
            add_regression_line(divlist,fraclist)
            plt.xlabel('Asexual Divergence')
            plt.ylabel('Normalized Number of %i-ton Flexible Genes Per Cell Pair' %nlist[i])
            plt.show()

def avelist(l):
    return float(sum(l))/len(l)
                 
def get_2_3_4_per_pair(celllist,qcut,pidcut,nlist):
    countdict={}
    for i in range(len(celllist)):
        for j in range(i):
            ci=celllist[i]
            cj=celllist[j]
            mypair=sorted([celllist[i],celllist[j]])
            mypair=(mypair[0],mypair[1])
            countdict[mypair]=[0 for i in range(len(nlist))]
                          
        
    f='../input_data/C1_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]

    for line in lines:
        group=ast.literal_eval(line)
        if len(group) not in nlist:###new april 18
            continue
        myidx=nlist.index(len(group))
        newgroup=[]
        #now find if each pair is in :(
        for c in group:
            idx=c.find('_')
            cell=c[idx+1:]
            if cell not in celllist:
                continue
            newgroup.append(cell)
        #now loop through pairs
        for i in range(len(newgroup)):
            for j in range(i):
                ci=newgroup[i]
                cj=newgroup[j]
                mypair=sorted([ci,cj])
                mypair=(mypair[0],mypair[1])
                countdict[mypair][myidx]+=1
    return countdict


#########################################

def binned_flex_gene_shared_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut):
    countdictpair=get_number_flex_genes_shared(groupname,celllist,qcut,pidcut)
    print('done1')
    lambdadict=get_poisson(celllist,groupname,pcut,chunklength)
    print('done2')
    countdictsingle=count_n_genes_per_cell_from_connected_comps(groupname,qcut,pidcut,celllist)

    #wtk: are the outlier-ish points that have a lot of shared genes all from the same cell?
    if 1:#count vs lambda
        countlistsingle=[countdictsingle[item] for item in countdictsingle if countdictsingle[item]>0]
        nflexfullgenome=float(sum(countlistsingle))/(len(countlistsingle)*.7)
        #print(lambdadict)
        ll,lcp=[],[]
        mydists=[[],[],[]]
        myfracdists=[[],[],[]]
        mycutdists=[[],[],[]]#exclude pairs with nore than 40 genes shared
        divbreaks=[.001,.0022,.004]
        #divbreaks=[.0005,.0008,.001]#.001,.0022,.004]
        for item in lambdadict:
            if lambdadict[item]=='-':
                continue
            c1=item[0]
            c2=item[1]
            mynorm=countdictsingle[c1]*float(countdictsingle[c2])
            if mynorm==0:
                #print('%s %s have %i %i flex genes' %(c1,c2,countdictsingle[c1],countdictsingle[c2]))
                continue
            ll.append(lambdadict[item])
            lcp.append(countdictpair[item])
            myfrac=float(lcp[-1])*nflexfullgenome/mynorm
            if ll[-1]<divbreaks[0]:
                mydists[0].append(lcp[-1])
                myfracdists[0].append(myfrac)
                if lcp[-1]<40:
                    mycutdists[0].append(lcp[-1])
            if ll[-1]<divbreaks[1] and ll[-1]>=divbreaks[0]:
                mydists[1].append(lcp[-1])
                myfracdists[1].append(myfrac)
                mycutdists[1].append(lcp[-1])
            if ll[-1]>=divbreaks[1] and ll[-1]<divbreaks[2]:
                mydists[2].append(lcp[-1])
                myfracdists[2].append(myfrac)
                mycutdists[2].append(lcp[-1])
                
            if lcp[-1]>20:
                print('Pair %s %i flex genes shared div %f?' %(item,lcp[-1],ll[-1]))
        plt.violinplot(mydists,showmeans=True,positions=list(range(len(divbreaks))))
        plt.title('%s Flexible genes per pair\nqcut %s pidcut %s Poison cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.xlabel('Estimated Asexual Divergence')
        plt.ylabel('Number of Flexible Genes Shared')
        ax=plt.gca()
        ax.set_xticks(list(range(len(divbreaks))))
        labs=['<%s\n(%i pairs)' %(str(divbreaks[0]),len(mydists[0])),'>=%s, <%s\n(%i pairs)' %(str(divbreaks[0]),str(divbreaks[1]),len(mydists[1])),'>=%s, <%s\n(%i pairs)' %(str(divbreaks[1]),str(divbreaks[2]),len(mydists[2]))]
        labs2=['<%s\n(%i pairs)\nMean %.1f Median %.1f' %(str(divbreaks[0]),len(mydists[0]),avelist(mydists[0]),statistics.median(mydists[0])),'>=%s, <%s\n(%i pairs)\nMean %.1f Median %.1f' %(str(divbreaks[0]),str(divbreaks[1]),len(mydists[1]),avelist(mydists[1]),statistics.median(mydists[1])),'>=%s, <%s\n(%i pairs)\nMean %.1f Median %.1f' %(str(divbreaks[1]),str(divbreaks[2]),len(mydists[2]),avelist(mydists[2]),statistics.median(mydists[2]))]
        labs3=['<%s\n(%i pairs)\nMean %.3f Median %.3f' %(str(divbreaks[0]),len(mydists[0]),avelist(myfracdists[0]),statistics.median(myfracdists[0])),'>=%s, <%s\n(%i pairs)\nMean %.3f Median %.3f' %(str(divbreaks[0]),str(divbreaks[1]),len(mydists[1]),avelist(myfracdists[1]),statistics.median(myfracdists[1])),'>=%s, <%s\n(%i pairs)\nMean %.3f Median %.3f' %(str(divbreaks[1]),str(divbreaks[2]),len(mydists[2]),avelist(myfracdists[2]),statistics.median(myfracdists[2]))]
        
        ax.set_xticklabels(labs2)
        plt.show()


        plt.violinplot(myfracdists,showmeans=True,positions=list(range(len(divbreaks))))
        plt.title('%s Normalized Flexible genes per pair\nqcut %s pidcut %s Poison cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.xlabel('Estimated Asexual Divergence')
        plt.ylabel('Number of Flexible Genes Shared/Number of flex genes in each cell')
        ax=plt.gca()
        ax.set_xticks(list(range(len(divbreaks))))
        labs=['<%s' %str(divbreaks[0]),'>=%s, <%s' %(str(divbreaks[0]),str(divbreaks[1])),'>=%s, <%s' %(str(divbreaks[1]),str(divbreaks[2]))]
        ax.set_xticklabels(labs3)
        plt.show()


        
        plt.violinplot(mycutdists,showmeans=True,positions=list(range(len(divbreaks))))
        plt.title('%s CUT at 40 genes Flexible genes per pair\nqcut %s pidcut %s Poison cut %s' %(groupname,str(qcut),str(pidcut),str(pcut)))
        plt.xlabel('Estimated Asexual Divergence')
        plt.ylabel('Number of Flexible Genes Shared')
        ax=plt.gca()
        ax.set_xticks(list(range(len(divbreaks))))
        ax.set_xticklabels(labs)
        plt.show()
#############################################################

def check_oceans(celllist,groupname,qcut,pidcut):
    f='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]

    idlist=[]
    if groupname=='HLII':
        cutsize=75
    if groupname=='Blue_Basin_NoClose':
        cutsize=14
    '''
    for each properly sized orthogroup, get fraction of cells from Atl. plot: ncells on x axis, frac Atl on y axis, one pt per group
    '''
    ncells,ncellsall,fracatl,fracatlall=[],[],[],[]#fracatl uses the proper "Atlantic" list that are labeled atlantic by berube; fracatlall uses all samples that look to be in atlantic ish region

    ngroups=0
    for line in lines:
        group=ast.literal_eval(line)
        if len(group)>cutsize:###new april 18
            continue
       
        thisncells=0
        thisnatl=0
        thisncellsall,thisnatlall=0,0
        #now find if each pair is in :(
        for c in group:
            idx=c.find('_')
            cell=c[:idx]#idx+1:]
            if cell not in celllist:
                continue
            sample=cell[3:6]
            
            if sample in Atlantic or sample in Pacific:
                thisncells+=1
                if sample in Atlantic:
                    thisnatl+=1
            if sample in atlanticlist or sample in pacificlist:
                thisncellsall+=1
                if sample in atlanticlist:
                    thisnatlall+=1
        if thisncellsall>=6:
            ngroups+=1       
        if thisncells>=2:
            ncells.append(thisncells)
            fracatl.append(float(thisnatl)/thisncells)
        if thisncellsall>=2:
            analyze_and_print_group(group,thisncellsall,thisnatlall,groupname,idlist)
            ncellsall.append(thisncellsall)
            fracatlall.append(float(thisnatlall)/thisncellsall)
    with open('../input_data/%s_All_Atl_Pac_Flex_ORF_nums.txt' %groupname,'w') as myoutfile:
        for item in idlist:
            myoutfile.write('%s\n' %item)

    print('%i groups with 6 to %i cells. %i ids' %(ngroups,cutsize,len(idlist)))
    if 1:
        plt.plot(ncells,fracatl,'o',alpha=0.3)
        plt.xlabel('Number of cells in the orthogroup')
        plt.ylabel('Fraction of cells from the Atlantic')
        plt.title('%s Searching for ocean-specific flexible genes:\nFraction of Cells in Orthogroup from Atlantic vs number of cells in orthogroup\nCell Subset' %groupname)
        plt.show()

        plt.plot(ncellsall,fracatlall,'o',alpha=0.3)
        plt.xlabel('Number of cells in the orthogroup')
        plt.ylabel('Fraction of cells from the Atlantic')
        plt.title('%s Searching for ocean-specific flexible genes:\nFraction of Cells in Orthogroup from Atlantic vs number of cells in orthogroup' %groupname)
        plt.show()
        
    if 1:
        nulls=get_ocean_null(ncells,fracatl)
        plt.hist([fracatl,nulls],bins=cutsize,label=["Data","Null"])
        plt.xlabel('Fraction of Cells from the Atlatic')
        plt.ylabel('Number of flexible gene orthogroups with <=%i cells' %cutsize)
        plt.title('%s  Searching for ocean-specific flexible genes:\nFraction of Cells in Orthogroup from Atlantic\nCell Subset' %groupname)
        plt.legend()
        plt.show()


        nulls=get_ocean_null(ncellsall,fracatlall)
        plt.hist([fracatlall,nulls],bins=cutsize,label=['Data','Null'])
        plt.xlabel('Fraction of Cells from the Atlatic')
        plt.ylabel('Number of flexible gene orthogroups with <=%i cells' %cutsize)
        plt.title('%s  Searching for ocean-specific flexible genes:\nFraction of Cells in Orthogroup from Atlantic' %groupname)
        plt.legend()
        plt.show()

    if 1:
        violinplotter(ncells,fracatl,cutsize,groupname)
        violinplotter(ncellsall,fracatlall,cutsize,groupname)


def analyze_and_print_group(group,thisncells,thisnatl,groupname,idlist):
    if groupname=='Blue_Basin_NoClose':
        ncut=5
        fcuta=.8
        fcutp=.1
    if groupname=='HLII':
        ncut=6
        fcuta=.8#.8
        fcutp=.2#.1
    if thisncells<ncut:
        return
    fatl=float(thisnatl)/thisncells
    if fatl<fcuta and fatl>fcutp:
        return
    orfs=[]
    for item in group:
        idx=item.find('_')
        orfs.append(item[idx+1:])
        #if orfs[-1] in ['2717327642','2667801843','2667788573','2667781951','2667774782','2667765717']:#phosphonate etc #iron ['2667765711','2717253521','2717334565']:
            #print('found orf, fatl %f' %fatl)
            #t=input('t')
    idlist.append(orfs[0])
    if fatl>=fcuta:
        print('Atlantic %f' %fatl)
        print(orfs)
        
    if fatl<=fcutp:
        #print('Pacific %f' %fatl)
        print(orfs[0])
        #t=input('t')
    #for item in orfs:
        #print(int(item))
    #idlist.append(orfs[0])
    
    #t=input('t')
    
def find_index(mylist,mynum):
    #mylist is first right bound through last left bound; no 0 oand no max
    if mynum<mylist[0]:
        return 0
    for i in range(len(mylist)-1):
        if mynum<mylist[i+1] and mynum>=mylist[i]:
            return i+1
    return len(mylist)

def add_label(violin, label,labels):
        color = violin["bodies"][0].get_facecolor().flatten()
        labels.append((mpatches.Patch(color=color), label))

def violinplotter(ncells,fracatl,cutsize,groupname):

    #plt.hist(ncells,bins=max(ncells)+1)
    #plt.show()
    if groupname=='Blue_Basin_NoClose':
        mylist=[2.5,3.5,4.5,5.5,8.5,11.5]
    if groupname=='HLII':
        mylist=[2.5,3.5,4.5,5.5,8.5,11.5,15.5,20.5,25.5]
    
    mydists=[[] for i in range(len(mylist)+1)]
    for i in range(len(ncells)):
        myn=ncells[i]
        myidx=find_index(mylist,myn)
        mydists[myidx].append(fracatl[i])

    labs=[]
    for k in range(len(mydists)):
        if k==0:
            mylab='%i cells\n(%i og.s)' %(math.floor(mylist[k]),len(mydists[k]))
        if k>0 and k<len(mylist):
            if math.ceil(mylist[k-1])!=math.floor(mylist[k]):
                mylab='%i-%i cells\n(%i og.s)' %(math.ceil(mylist[k-1]),math.floor(mylist[k]),len(mydists[k]))
            if math.ceil(mylist[k-1])==math.floor(mylist[k]):
                mylab='%i cells\n(%i og.s)' %(math.floor(mylist[k]),len(mydists[k]))
                
        if k==len(mylist):
            mylab='%i-%i cells\n(%i og.s)' %(math.ceil(mylist[k-1]),max(ncells),len(mydists[k]))
        labs.append(mylab)
    '''
    for k in range(len(mydists)):
        if k==0:
            mylab='%i cells\n(%i orthogroups)' %(math.floor(mylist[k]),len(mydists[k]))
        if k>0 and k<len(mylist):
            if math.ceil(mylist[k-1])!=math.floor(mylist[k]):
                mylab='%i-%i cells\n(%i orthogroups)' %(math.ceil(mylist[k-1]),math.floor(mylist[k]),len(mydists[k]))
            if math.ceil(mylist[k-1])==math.floor(mylist[k]):
                mylab='%i cells\n(%i orthogroups)' %(math.floor(mylist[k]),len(mydists[k]))
                
        if k==len(mylist):
            mylab='%i-%i cells\n(%i orthogroups)' %(math.ceil(mylist[k-1]),max(ncells),len(mydists[k]))
        labs.append(mylab)
    '''
    #now get null
    ntot=sum(ncells)
    listatl=[ncells[i]*fracatl[i] for i in range(len(ncells))]
    natl=sum(listatl)
    npac=ntot-natl
    mynulls=[[] for i in range(len(mylist)+1)]
    for i in range(len(ncells)):
        nc=ncells[i]
        #get expected fraction from random
        myrand=random.choices([0,1], weights=[npac,natl],  k=nc)
        myfrac=float(sum(myrand))/nc
        myidx=find_index(mylist,nc)
        mynulls[myidx].append(myfrac)
    print([len(mydists[i]) for i in range(len(mydists))])
    positions=list(range(len(mydists)))
    positions=[2*positions[i] for i in range(len(positions))]
    positions_null=[0.5+positions[i] for i in range(len(positions))]

    labels = []
    
    add_label(plt.violinplot(mydists,positions=positions,showmeans=True),'Data',labels)
    add_label(plt.violinplot(mynulls,positions=positions_null,showmeans=True),'Null',labels)
    plt.xlabel('Number of Cells in Orthogroup',fontsize=30)
    plt.ylabel('Fraction of Cells from Atlantic',fontsize=30)
    plt.legend(*zip(*labels),fontsize=20)
    ax=plt.gca()
    ax.set_xticks(positions)
    ax.set_xticklabels(labs)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #plt.title('%s Searching for ocean-specific flexible genes' %groupname)
    plt.show()


def old_violinplotter(ncells,fracatl,cutsize,groupname):

    plt.hist(ncells,bins=max(ncells)+1)
    plt.show()
    
    mydists=[[] for i in range(2,cutsize+1)]
    for i in range(len(ncells)):
        myn=ncells[i]
        myidx=myn-2
        mydists[myidx].append(fracatl[i])
    #now get null
    ntot=sum(ncells)
    listatl=[ncells[i]*fracatl[i] for i in range(len(ncells))]
    natl=sum(listatl)
    npac=ntot-natl
    mynulls=[[] for i in range(2,cutsize+1)]
    for i in range(len(ncells)):
        nc=ncells[i]
        #get expected fraction from random
        myrand=random.choices([0,1], weights=[npac,natl],  k=nc)
        myfrac=float(sum(myrand))/nc
        mynulls[nc-2].append(myfrac)
   
    positions=list(range(2,cutsize+1))
    positions_null=[0.5+positions[i] for i in range(len(positions))]
    plt.violinplot(mydists,positions=positions,showmeans=True)#,label='Data')
    plt.violinplot(mynulls,positions=positions_null,showmeans=True)
    plt.xlabel('Number of Cells in Orthogroup')
    plt.ylabel('Fraction of Cells from Atlantic')
    plt.legend()
    plt.title('%s Searching for ocean-specific flexible genes' %groupname)
    plt.show()
    
def get_ocean_null(ncells,fracatl):
    #want: number of flexible genes in group1 and group 2. then draw from all items in ncells
    ntot=sum(ncells)
    listatl=[ncells[i]*fracatl[i] for i in range(len(ncells))]
    natl=sum(listatl)
    npac=ntot-natl
    nulls=[]
    for i in range(len(ncells)):
        nc=ncells[i]
        #get expected fraction from random
        myrand=random.choices([0,1], weights=[npac,natl],  k=nc)
        nulls.append(float(sum(myrand))/nc)
    return nulls


def find_annotations():
    #use .gff files. sherlock only
    fieldnames=['Contig','IMG','CDS','start','stop','dot','strand','zero','ID']
    #or maybe we can just make a list of the cogs and grep.
#########################
if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run.')
    parser.add_argument('-n','--n',type=int)
    args = parser.parse_args()
    print(args)
    prog=args.prog
    n=args.n

    if n==3:
        celllist=clique1_without_ribotype
        groupname='7-Cell C1 Clique'
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype


    if n==0:
        groupname='HLII'
        celllist=hliilist2

    if n==2:
        groupname='Blue_Basin_NoClose'
        celllist=BB_list
    if 'B241_526K3' in celllist:
        celllist.remove('B241_526K3')
    if 'B243_495N4' in celllist:
        celllist.remove('B243_495N4')
        
    if prog==0:
        check_number_flex_genes(qcut,pidcut)

    if prog==1:
        chunklength=1000
        pcut=0.01
        qcut=80#have tried 60 qcut, 80 pidcut
        pidcut=80
        plot_LOG_flex_genes_shared_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut)


    if prog==2:
        chunklength=1000
        pcut=0.01
        qcut=80#have tried 60 qcut, 80 pidcut. now qcut 80, pidcut 80
        pidcut=80
        #plot_n_2_3_4_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut)
        binned_flex_gene_shared_v_lambda(celllist,groupname,qcut,pidcut,chunklength,pcut)

    if prog==3:
        qcut=60
        pidcut=80
        check_oceans(celllist,groupname,qcut,pidcut)
    if prog==4:
        plot_expected_number_flex_genes()
'''





air ('B241_527G5', 'B241_527N11') 25 flex genes shared div 0.001963?
Pair ('B241_527L16', 'B241_528K19') 25 flex genes shared div 0.000794?
Pair ('B241_527P5', 'B241_528K19') 65 flex genes shared div 0.000129?
Pair ('B241_527N11', 'B241_529C4') 29 flex genes shared div 0.002976?
Pair ('B241_526B17', 'B241_529O19') 47 flex genes shared div 0.000795?
Pair ('B241_527L16', 'B243_495N16') 73 flex genes shared div 0.000042?
Pair ('B241_527P5', 'B243_495N16') 26 flex genes shared div 0.000728?
Pair ('B241_528K19', 'B243_495N16') 38 flex genes shared div 0.000737?
Pair ('B243_495K23', 'B243_496E10') 22 flex genes shared div 0.000729?
Pair ('B241_527G5', 'B243_498L10') 21 flex genes shared div 0.001909?
Pair ('B243_496N4', 'B243_498L10') 47 flex genes shared div 0.000476?
Pair ('B243_498F21', 'B243_498L10') 27 flex genes shared div 0.002292?
Pair ('B241_529J11', 'B245a_518E10') 42 flex genes shared div 0.000044?
Pair ('B241_527L16', 'B245a_519O11') 23 flex genes shared div 0.000541?
Pair ('B241_528K19', 'B245a_519O11') 27 flex genes shared div 0.000747?
Pair ('B243_495N16', 'B245a_519O11') 30 flex genes shared div 0.000514?
Pair ('B243_495K23', 'B245a_520E22') 23 flex genes shared div 0.000852?
Pair ('B241_527L16', 'B245a_521B10') 23 flex genes shared div 0.000851?
Pair ('B241_528K19', 'B245a_521B10') 32 flex genes shared div 0.000597?
Pair ('B243_495N16', 'B245a_521B10') 26 flex genes shared div 0.000761?
Pair ('B245a_519O11', 'B245a_521B10') 30 flex genes shared div 0.000729?
Pair ('B241_527L16', 'B245a_521O20') 24 flex genes shared div 0.000849?
Pair ('B241_528K19', 'B245a_521O20') 31 flex genes shared div 0.000803?
Pair ('B243_495N16', 'B245a_521O20') 27 flex genes shared div 0.000812?
Pair ('B243_498F21', 'B245a_521O20') 26 flex genes shared div 0.002430?
Pair ('B245a_519O11', 'B245a_521O20') 27 flex genes shared div 0.000818?
Pair ('B245a_521B10', 'B245a_521O20') 27 flex genes shared div 0.000813?


NOT all same cell for pairs with >40 genes shared



BB cells with atlantic-specific flex genes: NOT all the same cells. could investigate: linkage of cells? (coverage will alter)

 1 AG-347-E23
      1 AG-347-I04
      4 AG-355-A09
     13 AG-355-A18
     11 AG-355-B23
      5 AG-355-G23
     11 AG-355-I04
      1 AG-355-I20
      9 AG-355-J04
     12 AG-355-J09
     11 AG-355-K03
      4 AG-355-K13
     17 AG-355-L02
     17 AG-355-M02
     19 AG-355-N02
      8 AG-355-N16
     16 AG-355-N18
     22 AG-355-O17
     21 AG-355-P15
      7 AG-355-P16
     14 AG-355-P18
      2 AG-402-I23
      1 AG-402-N17
     10 AG-418-D13
      3 AG-418-I21
      3 AG-418-O03
      6 AG-418-P06
      6 AG-418-P13
     21 AG-424-E18
     17 AG-424-P16
      7 AG-424-P18
      2 AG-459-N19

'''


'''
pacific .2 .8 orfs

Ga0171825_127	img_core_v400	CDS	11970	12647	.	-	0	ID=2717240750;locus_tag=Ga0171825_12710;product=hypothetical protein
Ga0171826_111	img_core_v400	CDS	11264	11380	.	-	0	ID=2717242073;locus_tag=Ga0171826_11123;product=hypothetical protein
Ga0115692_104	img_core_v400	CDS	50531	50761	.	+	0	ID=2667835910;locus_tag=Ga0115692_10468;product=hypothetical protein
Ga0115694_113	img_core_v400	CDS	8842	9003	.	-	0	ID=2667840163;locus_tag=Ga0115694_11317;product=hypothetical protein
Ga0115696_107	img_core_v400	CDS	49879	50544	.	-	0	ID=2667843100;locus_tag=Ga0115696_10755;product=PKHD-type hydroxylase
Ga0115698_105	img_core_v400	CDS	22745	23293	.	+	0	ID=2667845568;locus_tag=Ga0115698_10527;product=hypothetical protein
Ga0115698_105	img_core_v400	CDS	23404	24084	.	+	0	ID=2667845569;locus_tag=Ga0115698_10528;product=hypothetical protein
Ga0115706_102	img_core_v400	CDS	126793	127419	.	-	0	ID=2667826793;locus_tag=Ga0115706_102171;product=precorrin-6Y C5,15-methyltransferase (decarboxylating)
Ga0115710_103	img_core_v400	CDS	24930	26738	.	+	0	ID=2667690765;locus_tag=Ga0115710_10354;product=Predicted O-linked N-acetylglucosamine transferase, SPINDLY family
Ga0115711_114	img_core_v400	CDS	5270	5407	.	-	0	ID=2667693489;locus_tag=Ga0115711_1143;product=hypothetical protein
Ga0115711_114	img_core_v400	CDS	5533	5757	.	+	0	ID=2667693491;locus_tag=Ga0115711_1145;product=hypothetical protein
Ga0115712_108	img_core_v400	CDS	7614	8606	.	-	0	ID=2667693978;locus_tag=Ga0115712_1087;product=DNA gyrase subunit A
Ga0115714_102	img_core_v400	CDS	143758	144732	.	+	0	ID=2667695961;locus_tag=Ga0115714_102165;product=Sulfotransferase family protein
Ga0115714_135	img_core_v400	CDS	1	168	.	+	0	ID=2667697303;locus_tag=Ga0115714_1351;product=hypothetical protein
Ga0171833_118	img_core_v400	CDS	4084	4890	.	+	0	ID=2717150940;locus_tag=Ga0171833_1187;product=Fibronectin-binding protein A N-terminus (FbpA)
Ga0115716_101	img_core_v400	RNA	192203	192265	.	-	0	ID=2667699332;locus_tag=Ga0115716_101237
Ga0115719_103	img_core_v400	CDS	23663	23914	.	-	0	ID=2667704738;locus_tag=Ga0115719_10330;product=hypothetical protein
Ga0115720_128	img_core_v400	CDS	4256	5374	.	+	0	ID=2667707618;locus_tag=Ga0115720_1285;product=GDPmannose 4,6-dehydratase
Ga0115721_107	img_core_v400	CDS	19223	19603	.	+	0	ID=2667708678;locus_tag=Ga0115721_10735;product=hypothetical protein
Ga0171859_103	img_core_v400	CDS	93642	93752	.	-	0	ID=2717278993;locus_tag=Ga0171859_103104;product=hypothetical protein
Ga0115745_103	img_core_v400	CDS	162245	162373	.	-	0	ID=2667742509;locus_tag=Ga0115745_103196;product=hypothetical protein
Ga0171867_103	img_core_v400	CDS	20920	21147	.	+	0	ID=2717664108;locus_tag=Ga0171867_10330;product=hypothetical protein
Ga0115753_104	img_core_v400	CDS	48500	48961	.	-	0	ID=2667756691;locus_tag=Ga0115753_10459;product=NADH dehydrogenase
Ga0171988_105	img_core_v400	CDS	19044	19229	.	-	0	ID=2717731135;locus_tag=Ga0171988_10538;product=hypothetical protein
Ga0172004_140	img_core_v400	CDS	2206	2529	.	-	0	ID=2717741783;locus_tag=Ga0172004_1405;product= transcriptional regulator, ArsR family
Ga0172036_171	img_core_v400	CDS	3	662	.	+	0	ID=2717403657;locus_tag=Ga0172036_1711;product=Major Facilitator Superfamily protein
'''
