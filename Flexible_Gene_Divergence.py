

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
#import Flexible_Gene_Shared as flexrepo
import Get_R_From_Cell_Pair_Matrices as rrepo
import ast
import re
import Repo_site_catalogues_single_pair as repo
from os.path import exists
import Plot_Single_Gene as plotrepo2

'''
will want to: make divergence matrices. need to know which SAG is which sequence.  Are the titles correctly written in the same format as used to write matrices?

then for each gene: get the covered SAGs.  plot: y axis: div at this gene; x axis: WG mean/median div. will want to use synonymous

summary statistic for a given flex gene: 
a. mean, median cell pair divs, and max
b. mean, median ratio of cell pair div to WG div, and max

Question: multiple transfers into HLII or not?
1. if multiple transfers in from a more diverse pop than HLII, then diversity will be greater than HLII diversity
2. 
'''



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
        
###########

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
    #plot_summary_statistics(meandivs,mediandivs,maxdivs,meanrats,medianrats,maxrats,tfna,groupname,meanmedian,diamratio,ncellscovered)

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
        plt.xlabel('Number of Cells Covered')
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
        
def do_single_gene(wgmat,celllist,mat,overalllist,genenum,tfna,meandivs,mediandivs,maxdivs,meanrats,medianrats,maxrats,meanmedian,diamratio,ncellscovered,alldivs,allwg,wgdiam,genediam,dict_annotation,infilematrices,plotbool=1):
    #want a list of cell pair divs and relative cell pair divs. option to visualize
    relativedivs,genedivs,wgdivs=[],[],[]
    mynsites=[]
    for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            [nsnps,nmc]=mat[idxi][idxj]
            if nmc<10:
                continue
            #get WG div
            [wgdist,ngenes]=wgmat[idxi][idxj]
            dist=float(nsnps)/nmc
            if dist<=1:
                mynsites.append(nmc)
                genedivs.append(dist)
                wgdivs.append(wgdist)
                relativedivs.append(float(dist)/wgdist)
                allwg.append(wgdist)
                alldivs.append(dist)

    if len(genedivs)==0:
        return

    
    myncovered=0.5*(float(1.0+math.sqrt(1+8.0*len(genedivs))))
    ncellscovered.append(myncovered)
    meandivs.append(repo.avelist(genedivs))
    mediandivs.append(statistics.median(genedivs))
    maxdivs.append(max(genedivs))
    meanrats.append(repo.avelist(relativedivs))
    medianrats.append(statistics.median(relativedivs))
    maxrats.append(max(relativedivs))
    diamratio.append(float(max(genedivs))/max(wgdivs))
    wgdiam.append(max(wgdivs))
    genediam.append(max(genedivs))

    
    if plotbool and myncovered>12:
        
                
        if 1:
            if min(wgdivs)<0.02:
                #if diamratio[-1]>2:

                if repo.avelist(mynsites)>30:
                    #myannot=dict_annotation[genenum]
                    #print('Genen %i %s' %(genenum,myannot))
                    myannot=''
                    plt.plot(wgdivs,genedivs,'o',alpha=0.3)
                    plt.xlabel('Whole Genome %s Divergence' %meanmedian)
                    plt.ylabel('Flexible Gene Divergence')
                    plt.title('Flexible Gene %i %s Divergence\n(%i cell pairs, %.2f sites average)\n%s' %(genenum,tfna,len(wgdivs),repo.avelist(mynsites),myannot))
                    x=np.linspace(0,max([max(genedivs),max(wgdivs)]),100)
                    plt.plot(x,x,'-')
                    plt.axis('square')
                    plt.show()
                    #t=input('t')
        if 1:#distancemat
            #if diamratio[-1]>2 and repo.avelist(mynsites)>30:
            if min(wgdivs)<0.02:
                
                plotrepo2.plot_single_gene_dendro(celllist,genenum,groupname,tfna,'average',0,infilematrices)
            
   
    
    


###########################

def make_files_for_orthogroups(celllist,groupname,sizecut,qcut,pidcut):
    group_cells,group_ids=load_ortho_group_dict_from_clustering(celllist,groupname,sizecut,qcut,pidcut)
    #print(group_ids)
    filenamelist=[]
    #for each item in group_ids, title the file FlexibleGene_(item+2000)_ORFnames.txt. add this name to list. then write the orfs in the file
    for item in group_ids:
        if groupname!='C1':
            mytitle='FlexibleGene_%i_ORFnames.txt' %(int(item)+2000)
        if groupname=='C1':
            mytitle='FlexibleGene_%i_ORFnames.txt' %(int(item)+3000)
        filenamelist.append(mytitle)
        with open(mytitle,'w') as outf:
            for orf in group_ids[item]:
                outf.write('%s\n' %orf)
    with open('Flexible_Gene_File_Name_List.txt','w') as outf:
        for i in range(len(filenamelist)):
            outf.write('%s\n' %(filenamelist[i]))

def load_ortho_group_dict_from_clustering(celllist,groupname,sizecut,qcut,pidcut):
    f='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    groups_cells,groups_ids={},{}

    lower_cut=4
    lineindex=0
    for line in lines:
        group=ast.literal_eval(line)
        if len(group)<lower_cut:
            continue    
        if len(group)>sizecut:
            continue
        newgroup=[]
        newid=[]
        #now find if each pair is in :(
        for c in group:
            idx=c.find('_')
            if groupname=='C1':
                mynewid=c[:idx]
                cell=c[idx+1:]
            if groupname!='C1':
                cell=c[:idx]
                mynewid=c[idx+1:]
            if cell not in celllist:
                print(cell)
                print(celllist)
                t=input('t')
                continue
            newgroup.append(cell)
            newid.append(mynewid)
        groups_cells[lineindex]=newgroup
        groups_ids[lineindex]=newid
        lineindex+=1
    return groups_cells,groups_ids


    
##########################

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

    if prog==0:
        pidcut,qcut=80,80
        make_files_for_orthogroups(celllist,groupname,sizecut,qcut,pidcut)
      
    if prog==1:
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

    if prog==2:
        get_annotations_hlii()
