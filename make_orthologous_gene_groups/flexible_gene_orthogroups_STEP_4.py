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
import random
import argparse
import Get_R_From_Cell_Pair_Matrices as rrepo
import ast
import re
import Repo_site_catalogues_single_pair as repo
from os.path import exists


'''
This file finds divergences within flex orthogroups using fastas created in step 3 (p=0) in step a

in step b (optional), it returns annotations (p=1)

Internal note: made from Flexible_gene_divergences
'''




################################################################
#step A: get flex gene divergences (p=0)
################################################################


def get_relative_div_dists(arrayfile,infilewg,tfna,celllist,groupname,meanmedian):
    overalllist=rrepo.overall_order_list()
    genearrays,mykeys=rrepo.load_arrays(arrayfile)

    dict_annotation={}
    if groupname!='C1':
        dict_annotation=load_annotation_dict()
    

    wgarray,keylist=rrepo.load_arrays(infilewg)
    wgmat=wgarray['%sWG_%s' %(meanmedian,tfna)]
    
    meandivs,mediandivs,maxdivs=[],[],[]
    meanrats,medianrats,maxrats=[],[],[]#rat is ratio of cell pair div at gene to cell pair WG div, do with mean or median

    diamratio=[]
    wgdiam,genediam=[],[]
    #do single gene
    ncellscovered=[]

    alldivs,allwg=[],[]

    for item in list(genearrays.keys()):
        print(item)
        genenum=repo.get_gene_from_key_to_check(item)
        mat=genearrays[item]
        do_single_gene(wgmat,celllist,mat,overalllist,genenum,tfna,meandivs,mediandivs,maxdivs,meanrats,medianrats,maxrats,meanmedian,diamratio,ncellscovered,alldivs,allwg,wgdiam,genediam,dict_annotation,arrayfile,plotbool=0)
    plot_summary_statistics(meandivs,mediandivs,maxdivs,meanrats,medianrats,maxrats,tfna,groupname,meanmedian,diamratio,ncellscovered)

    if 1:#plot plt and hist2d of all cell pairs at all genes, WG dist and flex dist
        plt.plot(allwg,alldivs,'o',alpha=0.3)
        plt.xlabel('Whole genome divergence (%s)' %tfna)
        plt.ylabel('Flexible Gene Divergence (%s)' %tfna)
        plt.title('%s Divergences at flexible genes vs WG' %(groupname))
        x=np.linspace(0,max([max(alldivs),max(allwg)]),100)
        plt.plot(x,x,'-')
        plt.axis('square')
        plt.show()

        #now as 2dhist

        xlabel='Whole Genome Divergence'
        ylabel='Flexible Gene Divergence'
        mybins=100#[i*.0025 for i in range(math.ceil(400*max([max(allwg),max(alldivs)])))]
        plt.hist2d(allwg,alldivs, bins=[mybins,mybins], norm=colors.LogNorm())
        #plt.xscale('log')
        #plt.yscale('log')
        cbar=plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('Cell Pairs x Genes', rotation=270,fontsize=15)

        #plt.xlabel('Log %s' %xlabel,fontsize=20)
        #plt.ylabel('Log %s' %ylabel,fontsize=20)
        ax=plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=20)
        #plt.title('Log %s vs Log %s' %(ylabel,xlabel),fontsize=20)
        plt.xlabel(xlabel,fontsize=20)
        plt.ylabel(ylabel,fontsize=20)
        #plt.title('%s %s vs %s' %(groupname,xlabel,ylabel))
        plt.tight_layout()
        plt.show()


    if 1:#plot 2dhist of gene diameter vs WG diameter
        

        xlabel='Whole Genome Diameter'
        ylabel='Flexible Gene Diameter'
        #mybins=[i*.0025 for i in range(32)]
        mybins=20#[i*.025 for i in range(math.ceil(40*max([max(wgdiam),max(genediam)])))]
        plt.hist2d(wgdiam,genediam, bins=[mybins,mybins], norm=colors.LogNorm())
        #plt.xscale('log')
        #plt.yscale('log')
        cbar=plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('Number of Flexible Genes', rotation=270,fontsize=15)

        #plt.xlabel('Log %s' %xlabel,fontsize=20)
        #plt.ylabel('Log %s' %ylabel,fontsize=20)
        ax=plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=20)
        #plt.title('Log %s vs Log %s' %(ylabel,xlabel),fontsize=20)
        plt.xlabel(xlabel,fontsize=20)
        plt.ylabel(ylabel,fontsize=20)
        #plt.title('%s %s vs %s' %(groupname,xlabel,ylabel))
        plt.tight_layout()
        plt.show()

        
    if 0:#plot estimated number of cells covered per gene (not always integers becuase, for example, two cells could be covered but their mutual coverage is so low that they don't make it into the gene divergences list)
        print(sorted(ncellscovered))
        mybins=[-.5]
        for i in range(math.ceil(max(ncellscovered))+1):
            mybins.append(i+.5)
        plt.hist(ncellscovered,bins=mybins)
        #plt.xlabel('Number of Cells Covered')
        plt.ylabel('Number of Flexible Gene Orthogroups')
        plt.title('%s Number of cells covered at\nflexible gene orthogroups' %groupname)
        plt.show()
    
def plot_summary_statistics(meandivs,mediandivs,maxdivs,meanrats,medianrats,maxrats,tfna,groupname,meanmedian,diamratio,ncellscovered):
    plot_hist(groupname,tfna,meandivs,'Mean Flexible Gene Divergence')
    plot_hist(groupname,tfna,mediandivs,'Median Flexible Gene Divergence')
    plot_hist(groupname,tfna,maxdivs,'Max Flexible Gene Divergence')
    plot_hist(groupname,tfna,meanrats,'Mean Gene/WG (%s) Ratio' %meanmedian)
    plot_hist(groupname,tfna,medianrats,'Median Gene/WG (%s) Ratio' %meanmedian)
    plot_hist(groupname,tfna,maxrats,'Max Gene/WG (%s) Ratio' %meanmedian)
    plot_hist(groupname,tfna,diamratio,'Diameter Ratio (Gene/WG (%s))' %meanmedian)
    
    
def plot_hist(groupname,tfna,mylist,listname):
    plt.hist(mylist,bins=200)
    plt.xlabel('%s' %listname)
    plt.ylabel('Number of Cell Pairs')
    plt.title('%s Flexible Genes (%s Divergence)\n%s' %(groupname,tfna,listname))
    plt.show()

def load_annotation_dict():
    dict_annotation={}
    with open('../input_data/HLII_Flexible_Gene_Annotations.csv','r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            g=int(row['GeneNum'])
            annot=row['Annotation']
            dict_annotation[g]=annot
    return dict_annotation
        

    

################################################################
#step B: get annotations (p=1) (optional)
################################################################

def get_annotations_hlii():
    #look at genenum fasta file.  get headings.  get ORF id out of heading. then get gff annotation.
    #goal: file with genenum: id
    outfile='HLII_Flexible_Gene_Annotations.csv'
    outfields=['GeneNum','Annotation']
    genenums,annotations=[],[]
    for i in range(731):
        gnum=2000+i
        print(gnum)
        #load file. get header
        myfile='Flexible_Orthogroups_Aligned/FlexibleGene%i.fasta' %gnum
        
        orfid=get_orfid(myfile)
        genenums.append(gnum)
        annotations.append(get_gff_annotation(orfid))
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,fieldnames=outfields,delimiter='\t')
        outwriter.writeheader()
        for i in range(len(genenums)):
            outwriter.writerow({'GeneNum':genenums[i],'Annotation':str(annotations[i])})
                           
def get_gff_annotation(myidlist):
    with open('../input_data/All_HLII_Gff.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    notelist=[]
    if len(myidlist)>=1:#for myid in myidlist[1:3]:
        myid=myidlist[1]
        print(myid)
        for line in lines:
            if re.search(myid,line):
                #print(line)
                if 'product' not in line:
                    notelist.append('Unannotated')
                    continue
                idx=line.index('product')
                
                print(line[idx+8:])
                notelist.append(line[idx+8:])
    return notelist

def get_orfid(infile):
    if not exists(infile):
        print('File %s does not exist' %(infile))
        return []
    with open(infile,'r') as myinfile:
        lines=myinfile.read()
    cellseqs=re.split('>',lines)
    headers=[]
    for i in range(len(cellseqs)):
        sublist=re.split('\n',cellseqs[i])
        headers.append(sublist[0])
        
    #now get orfs
    orfs=[]
    
    for i in range(len(headers)):
        item=headers[i]
        idx=item.find(' ')
        orfs.append(item[:idx])
    return orfs
        



###########################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-n','--n',type=int)
    parser.add_argument('-tfna','--tfna',type=str)
  
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    n=args.n
    tfna=args.tfna
   

    c1bool=0
    if n==0:
        groupname='HLII'
        celllist=hliilist2
        sizecut=34
    if n==2:
        groupname='Blue_Basin_NoClose'
        celllist=BB_list
        sizecut=14
   
    if n==1:
        celllist=C1_without_ribotype
        sizecut=14
        groupname='C1'
        c1bool=1
        
   
      
    if prog==0:
        suffix='_Flexible'
        groupnamet='Berube_Kash'
        meanmedian='Mean'
        infilewg='../input_data/Cell_Pair_Matrices/Berube_Kash_WG_Matrices_%s.npy' %(tfna)
         
        if tfna=='tf':
            arrayfile='../input_data/Cell_Pair_Matrices/Each_Cell_Pair_Sites/Matrices_%s_TWO_AND_FOUR%s.npy' %(groupnamet,suffix)
        if tfna=='a':
            arrayfile='../input_data/Cell_Pair_Matrices/Each_Cell_Pair_Sites/Matrices_%s_All_NT%s.npy' %(groupnamet,suffix)

        if groupname=='C1':
            arrayfile=arrayfile[:-4]+'2_C1.npy'
        get_relative_div_dists(arrayfile,infilewg,tfna,celllist,groupname,meanmedian)

    if prog==1:
        get_annotations_hlii()
