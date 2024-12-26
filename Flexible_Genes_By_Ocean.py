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
import Plotting_Repo as plotrepo
import random
import argparse
import Make_Synonymous_Distance_Matrices_Berube_Kashtan_Using_All_Sites_Per_Pair_and_HalfGene_Record_sites_SNPs as distmatrepo
import Flexible_Gene_Shared as flexrepo
import Official_Geography_Plot as oceanrepo

HLII_BB_NoCloseCells_Median04=['AG-347-E23', 'AG-347-I04', 'AG-347-I06', 'AG-347-I19', 'AG-347-I21', 'AG-347-J14', 'AG-347-K16', 'AG-347-K19', 'AG-347-K20', 'AG-347-L02', 'AG-347-L17', 'AG-347-L20', 'AG-355-A09', 'AG-355-I04', 'AG-355-J09', 'AG-355-M02', 'AG-355-N02', 'AG-355-N16', 'AG-355-O17', 'AG-355-P16', 'AG-402-I23', 'AG-402-L23', 'AG-418-D13', 'AG-424-E18', 'AG-442-B03', 'AG-442-N07', 'AG-355-N21', 'AG-402-O23', 'AG-347-C10', 'AG-347-G20', 'AG-347-G22', 'AG-347-J06', 'AG-347-K10', 'AG-347-K22', 'AG-347-K23', 'AG-347-L19', 'AG-347-L21', 'AG-347-M18', 'AG-347-M23', 'AG-355-A18', 'AG-355-B23', 'AG-355-G23', 'AG-355-I20', 'AG-355-J04', 'AG-355-K03', 'AG-355-K13', 'AG-355-L02', 'AG-355-N18', 'AG-355-P15', 'AG-355-P18', 'AG-402-A04', 'AG-402-G23', 'AG-402-K16', 'AG-402-N17', 'AG-402-N23', 'AG-418-I21', 'AG-418-J17', 'AG-418-O03', 'AG-418-P06', 'AG-418-P13', 'AG-424-P16', 'AG-424-P18', 'AG-436-D21', 'AG-449-O05', 'AG-455-E15', 'AG-459-B06', 'AG-459-D04', 'AG-459-N19', 'AG-459-O09']#about 70 cells

HLII_Outliers_NoCloseCells_Median04=['AG-347-J19', 'AG-347-J20', 'AG-347-J21', 'AG-347-K17', 'AG-347-K18', 'AG-355-B18', 'AG-355-P07', 'AG-388-A01', 'AG-402-K21', 'AG-402-M15', 'AG-402-N08', 'AG-402-O16', 'AG-412-J13', 'AG-418-I20', 'AG-418-M21', 'AG-335-I15', 'AG-335-J02', 'AG-347-B08', 'AG-347-J05', 'AG-347-K21', 'AG-355-J21', 'AG-424-L22', 'AG-432-G10', 'AG-459-P20', 'AG-469-M13', 'AG-670-J21', 'AG-347-B23', 'AG-347-I15', 'AG-347-K02', 'AG-347-K15', 'AG-347-M15', 'AG-347-N23', 'AG-355-M18', 'AG-388-F11', 'AG-402-C22', 'AG-402-K22', 'AG-412-A14', 'AG-412-C21', 'AG-418-C17', 'AG-418-F08', 'AG-418-G18', 'AG-418-G23', 'AG-418-K17', 'AG-424-M03', 'AG-424-P23', 'AG-436-A04', 'AG-436-E22', 'AG-449-G23', 'AG-449-J16', 'AG-459-A02', 'AG-459-E08']



Atlantic=[item for item in hliilist2 if item[3:6] in ['355','363','388','412','429','432','436']]
N_Pacific=[item for item in hliilist2 if item[3:6] in ['347','402']]
S_Pacific=[item for item in hliilist2 if item[3:6] in ['335','459','469','449','455']]
Pacific=N_Pacific+S_Pacific

Atl_NoClose=[]
SPac_NoClose=[]
NPac_NoClose=[]

Atl_JustClose,SPac_JustClose,NPac_JustClose=[],[],[]
for item in Atlantic:
    if item in HLII_BB_NoCloseCells_Median04 or item in HLII_Outliers_NoCloseCells_Median04:
        Atl_NoClose.append(item)
for item in S_Pacific:
    if item in HLII_BB_NoCloseCells_Median04 or item in HLII_Outliers_NoCloseCells_Median04:
        SPac_NoClose.append(item)
for item in N_Pacific:
    if item in HLII_BB_NoCloseCells_Median04 or item in HLII_Outliers_NoCloseCells_Median04:
        NPac_NoClose.append(item)

for item in Atlantic:
    if item not in HLII_BB_NoCloseCells_Median04 and item not in HLII_Outliers_NoCloseCells_Median04:
        Atl_JustClose.append(item)
for item in S_Pacific:
    if item not in HLII_BB_NoCloseCells_Median04 and item not in HLII_Outliers_NoCloseCells_Median04:
        SPac_JustClose.append(item)
for item in N_Pacific:
    if item not in HLII_BB_NoCloseCells_Median04 and item not in HLII_Outliers_NoCloseCells_Median04:
        NPac_JustClose.append(item)

        

    

'''

number of flex genes shared by a pair vs frac time pair is from the same ocean
'''





def prob_same_ocean_vs_shared(celllist,groupname,qcut,pidcut,group1,group2,totgroup):#group1,group2 are atlantic, pacific
    countdictpair=flexrepo.get_number_flex_genes_shared(groupname,celllist,qcut,pidcut)
    #print(countdictpair)
    nshared,sameocean=[],[]
    for item in countdictpair:
        c1=item[0]
        c2=item[1]
        if c1 not in group1 and c1 not in group2:
            continue
        if c2 not in group1 and c2 not in group2:
            continue
        ns=countdictpair[item]
        nshared.append(ns)
        if c1 in group1 and c2 in group1:
            sameocean.append(1)
            continue
        if c2 in group2 and c1 in group2:
            sameocean.append(1)
            continue
        sameocean.append(0)
    #now binned stats

    #print(nshared)
    mybins1=[1.5]
    for i in range(10):
        mybins1.append(2*i+2.5)
    mybins1.append(30)
    mybins1.append(40)
    mybins1.append(50)
    mybins1.append(60)
    mybins1.append(70)


    
   
    print('%i cell pairs? %i expected' %(len(nshared),0.5*(len(celllist))*(-1+len(celllist))))
    sharedbins=[-.5]
    for i in range(max(nshared)):
        sharedbins.append(i+.5)
    plt.hist(nshared,bins=sharedbins)
    plt.ylabel('Number of SAG pairs',fontsize=20)
    plt.xlabel('Number of shared flexible genes',fontsize=20)
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=25)r2_null(celllist,groupname,sizecut,qcut,pidcut)
    if group==1:
        plt.title('Number of flexible genes shared per SAG pair\nExcluding closely-related cell groups',fontsize=20)
    if group==0:
        plt.title('Number of flexible genes shared per SAG pair\nAll quasi-random HLII SAGs',fontsize=20)
    if group==3:
        plt.title('Number of flexible genes shared per SAG pair\nUsing only closely-related cell groups',fontsize=20)
    plt.tight_layout()
    plt.show()

    plt.hist(nshared,bins=sharedbins)
    plt.ylabel('Number of SAG pairs',fontsize=20)
    plt.xlabel('Number of shared flexible genes',fontsize=20)
    plt.yscale('log')
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=25)
    if group==1:
        plt.title('Number of flexible genes shared per SAG pair\nExcluding closely-related cell groups',fontsize=20)
    if group==0:
        plt.title('Number of flexible genes shared per SAG pair\nAll quasi-random HLII SAGs',fontsize=20)
    if group==3:
        plt.title('Number of flexible genes shared per SAG pair\nUsing only closely-related cell groups',fontsize=20)
    plt.tight_layout()
    plt.show()


    #return
    
    mybins1=20

    
    bin_means1,bin_edges1,binnumber1=stats.binned_statistic(nshared,sameocean,statistic='mean',bins=mybins1)
    ngenes_per_bin1,bin_edges1,binnumber1=stats.binned_statistic(nshared,sameocean,statistic='count',bins=mybins1)
    binmidpoints1=0.5*(bin_edges1[1:]+bin_edges1[:-1])
    ax=plt.gca()
    ax.plot(binmidpoints1,bin_means1,'o',color='blue')

    p1,var1=oceanrepo.get_null_closest_relative_distribution_simple([group1,group2],ngenes_per_bin1)
    ax.axhline(y=p1,color='cyan',label='Well-Mixed Null')
    #ax.axhspan(p1-math.sqrt(var1),p1+math.sqrt(var1),alpha=0.2,color='cyan')
    myy11=[p1-math.sqrt(var1[j]) for j in range(len(var1))]
    myy21=[p1+math.sqrt(var1[j]) for j in range(len(var1))]
    plt.fill_between(binmidpoints1,myy21,y2=myy11,color='cyan',alpha=0.2)

   
    ax.set_ylim(bottom=0)#[0,1.2])
    ax.set_xlabel('Number of Shared Covered Flexible Genes',fontsize=20)
    ax.set_ylabel('Fraction of Cells from Same Ocean',fontsize=20)
     
    ax.legend(fontsize=15)
    #ax.set_xticks([0,0.005,.01,.015,.02,.025])
    #ax.set_xticklabels(['0%','0.5%','1%','1.5%','2%','2.5%'])
    ax.set_yticks([0,0.25,0.5,.75,1])
    ax.set_yticklabels(['0','0.25','0.50','0.75','1.00'])
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.title('Fraction of Cell Pairs from the Same Ocean\nvs Number of Shared Flexible Genes' ,fontsize=20)
    plt.tight_layout()
    plt.show()




if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-n','--n',type=int)
    parser.add_argument('-g','--group',type=int,default=0)
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    group=args.group
    n=args.n

    if n==0:
        groupname='HLII'
        celllist=hliilist2
    if n==2:
        groupname='Blue_Basin_NoClose'
        celllist=BB_list

    if group==0:#atlantic vs pac includign close, then witout close
        group1=Atlantic
        group2=Pacific
    if group==1:
        group1=Atl_NoClose
        group2=SPac_NoClose+NPac_NoClose
    if group==3:
        group1=Atl_JustClose
        group2=SPac_JustClose+NPac_JustClose

    group1=list(set(group1).intersection(set(celllist)))
    group2=list(set(group2).intersection(set(celllist)))


    if prog==0:
        qcut=80
        pidcut=80
        prob_same_ocean_vs_shared(celllist,groupname,qcut,pidcut,group1,group2,group)#group1,group2 are atlantic, pacific
