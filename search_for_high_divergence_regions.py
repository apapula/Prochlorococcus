

import csv
import matplotlib.pyplot as plt
import math
from os.path import exists
import re
import ast
import statistics
import random
from matplotlib import gridspec, colors
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.cluster.hierarchy as hac
import scipy.spatial.distance as ssd
from Kashtan_lists import CN2, Non_CN2, C1, clique1,cl2,cl3,cl5,cl6,fourcliques,C3, C1_without_ribotype
from matplotlib.backends.backend_pdf import PdfPages

'''
This script uses matrices with gene divergences made with make_distance_matrices_for_each_gene.py

The divergences of cell pairs in celllist are taken at a given gene and clustered using fcluster with a divergence cutoff.  If multiple clusters are found, the gene is flagged and various statistics recorded in outfile (gene number, length, top merge divergence, number of clusters, and cells in each cluster).  Synteny is checked to reduce alignment error influence.

After this the csv file can be checked for co-located clusters along the genome, indicating a transfer event spanning more than 1 gene

This file also contains programs for plotting recombined divergence distributions, and checking local synteny
'''
    
def get_celllist_allhlii():
    myfile='../input_data/HLII_clonelist.txt'
    with open(myfile,'r') as infile:
        celllist=[line.rstrip() for line in infile if len(line)>0]
    for  item in CN2:
        celllist.append(item)
    for item in Non_CN2:
        celllist.append(item)
    return celllist


def Produce_Clusters_File(infile_matrices,celllist,overallcelllist,outfile,divcut,groupname,celllist2,clusteringmethod='single',plotbool=0,pdf='',colorcelllists=[]):
    print('hello')
    cell_indices=[]
    fieldnames=['GeneNum','Top_Merge_Div','NClusters']
    for i in range(len(celllist)):
        fieldnames.append(celllist[i])
        cell_indices.append(overallcelllist.index(celllist[i]))
    
    with open(outfile,'w') as myoutfile:
        outcsv=csv.writer(myoutfile,delimiter='\t')
        outcsv.writerow(fieldnames)
        cluster_genes_and_record(infile_matrices,celllist,cell_indices,overallcelllist,outcsv,divcut,clusteringmethod,plotbool,pdf,groupname,colorcelllists,celllist2)
    if plotbool:
        pdf.close()
        
def cluster_genes_and_record(infile_matrices,celllist,cell_indices,overallcelllist,outcsv,divcut,clusteringmethod,plotbool,pdf,groupname,colorcelllists,celllist2):
    pairwise_divs=[]
    covered_genes=[]
    with open(infile_matrices,'r') as myinfile:
        l=myinfile.readline()       
        mat=[]
        genenum=1
        while len(l.rstrip())!=0:
            l=myinfile.readline()
            if len(l)==0:
                break
            if l[0]=='#':
                ind1=l.find('(')
                print('Gene %i' %genenum)
               
                #now compute on matrix mat
                ###########################
                get_clusters_from_on_matrix(genenum,mat,celllist,cell_indices,overallcelllist,outcsv,divcut,clusteringmethod,plotbool,pdf,groupname,colorcelllists,celllist2)
                ###########################
                
                #genenum=int(l[7:ind1-1])#for NON synonymous
                ind2=l.find(' ',11)
                #genenum=int(l[11:ind2])#for synonymous
                ind=l.find('_')
                genenum=int(l[8:ind])
                mat=[]
                l=myinfile.readline()
            mat.append(list(map(float, l.split())))
    return pairwise_divs,covered_genes


def get_clusters_from_on_matrix(genenum,mat,celllist,cell_indices,overallcelllist,outcsv,divcut,clusteringmethod,plotbool,pdf,groupname,colorcelllists,celllist2):
    #first, check the gene and re-acquire a celllist of acceptable cells
    all_same_contig,uncovered_pair_same_contig,pair_same_contig=Check_Local_Synteny(celllist2,genenum)
    syntenic_celllist=[]
    for item in all_same_contig:
        syntenic_celllist.append(item+'_C1')
    for item in pair_same_contig:
        syntenic_celllist.append(item+'_C1')
    for item in uncovered_pair_same_contig:
        syntenic_celllist.append(item+'_C1')
    
    ##########################
    submat, coveredcells=get_submatrix(mat,cell_indices,celllist,overallcelllist,syntenic_celllist)
    if len(coveredcells)<2:
        return
    clusters,labels,topmerge=cluster_matrix(submat,coveredcells,clusteringmethod,divcut,plotbool,pdf,groupname,genenum,colorcelllists)
    if max(labels)>0:
        write_clusters(genenum,clusters,coveredcells,labels,topmerge,outcsv,celllist)
    
def cluster_matrix(submat,coveredcells,clusteringmethod,divcut,plotbool,pdf,groupname,genenum,colorcelllists):
    myarray=np.array(submat)
    myarray_s=ssd.squareform(myarray)
    clustering=hac.linkage(myarray_s,method='%s' %clusteringmethod,metric='precomputed')
    topmerge=clustering[-1,-2]
    myleaves=hac.leaves_list(clustering).tolist()

    if plotbool:
        plot_dendrogram_of_clustering_in_pdf(pdf,clustering,myleaves,coveredcells,colorcelllists,groupname,genenum,clusteringmethod)
    myclusters=fcluster(clustering,divcut,criterion='distance')
    labels    = fcluster(clustering,divcut,criterion='distance')-1
    clusters  = [[] for i in range(max(labels) + 1)]
    for i in range(len(labels)):
        cellname=coveredcells[i]

        clusters[labels[i]].append(cellname)
    return clusters,labels,topmerge

def plot_dendrogram_of_clustering_in_pdf(pdf,clustering,myleaves,coveredcells,colorcelllists,groupname,genenum,clusteringmethod):
    clustercolors=[]
    colordict={0:'r',1:'b',2:'g',3:'k',4:'m',5:'orange'}
    for i in range(len(myleaves)):
        mycell=coveredcells[myleaves[i]]
        for j in range(len(colorcelllists)):
            if mycell in colorcelllists[j]:
                clustercolors.append(colordict[j])
                continue
        
    dn=dendrogram(clustering,labels=coveredcells,leaf_font_size=3.5,leaf_rotation=90)
    plt.title('%s Gene %s %s Linkage' %(groupname,str(genenum),clusteringmethod))
    ax=plt.gca()
    xlbls=ax.get_xmajorticklabels()
    for j in range(len(xlbls)):
        xlbls[j].set_color(clustercolors[j])
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    
def write_clusters(genenum,clusters,coveredcells,labels,topmerge,outcsv,celllist):
    myrow=[genenum,topmerge,max(labels)+1]
    for i in range(len(celllist)):
        c=celllist[i]
        if c not in coveredcells:
            myrow.append('-')
        if c in coveredcells:
            myidx=coveredcells.index(c)
            myallele=labels[myidx]
            myrow.append(myallele)
    outcsv.writerow(myrow)
    
    
        
def get_submatrix(mat,cell_indices,celllist,overallcelllist,syntenic_celllist):
    mat=np.array(mat)
    mydim=len(overallcelllist)
    m1tot=mat.reshape(mydim,mydim)
    coveredcells=[]
    covered_indices=[]
    for i in range(len(cell_indices)):
        idx=cell_indices[i]
        mycell=overallcelllist[idx]
        if mycell not in syntenic_celllist:
            continue
        mycol=m1tot[idx]
        mycol=list(mycol[:idx])+list(mycol[idx+1:])
        mycol=np.array(mycol)
        if np.all(mycol==2.):
            continue
        covered_indices.append(idx)
        coveredcells.append(overallcelllist[idx])
    ncov=len(coveredcells)
    submat=[[0 for x in range(ncov)] for y in range(ncov)]

    neworder_cells=[]
    covered_indices=sorted(covered_indices)
    for j in range(len(covered_indices)):
        oldj=covered_indices[j]
        neworder_cells.append(overallcelllist[oldj])
        for k in range(j):
            oldk=covered_indices[k]
            submat[j][k]=m1tot[oldj][oldk]
            submat[k][j]=m1tot[oldj][oldk]
    return submat,neworder_cells#neworder_cells should be in a corresponding order to submat



def Plot_Merge_Divs_vs_Background_Official(data_csv,infile_matrices,compare_celllist,overallcelllist,groupname,corebool,divcut,clusteringmethod,celllist,plotsfsbool=0):

    
    background_divs_allpairs,background_divs_averagepergene=get_background_divergences(data_csv,infile_matrices,compare_celllist,overallcelllist,corebool)
    truncated_background_divs_allpairs=[]
    for item in background_divs_allpairs:
        if item>=divcut:
            truncated_background_divs_allpairs.append(item)
    merge_divs=get_column_from_csv(data_csv,'Top_Merge_Div',column_type='float',corebool=corebool)

   
    
    corestring=''
    if corebool:
        corestring='Core '
    
    plt.hist([truncated_background_divs_allpairs,merge_divs],weights=[[len(merge_divs)/float(len(truncated_background_divs_allpairs)) for i in range(len(truncated_background_divs_allpairs))],[1 for i in range(len(merge_divs))]],bins=40,label=['Available Divergences in HLII','Recombined Divergences'])
    # plt.hist([truncated_background_divs_allpairs,merge_divs],bins=40,density=True,label=['Available Divergences','Recombined Divergences'])
    plot_reference_divergences_NT_2()
    plt.xlabel('Nucleotide Divergence',fontsize=25)
    plt.title('Recombined Divergences within Cluster\n vs Available Divergences within HLII',fontsize=25)
    plt.legend()
    ax=plt.gca()
    plt.ylabel('Number of Events',fontsize=25)
    #ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=2))
    ax.set_xticks([0,0.1,0.2,.3,.4])
    ax.set_yticks([0,20,40,60])
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.set_xticklabels(['0%','10%','20%','30%','40%'])
    plt.tight_layout()
    plt.show()


    plt.hist([truncated_background_divs_allpairs,merge_divs],bins=40,density=True,label=['Available Divergences','Recombined Divergences'])
    plt.xlabel('Nucleotide Divergence',fontsize=25)
    plt.title('Recombined Divergences within Ribotypes\n vs Available Divergences within HLII',fontsize=30)
    ax=plt.gca()
    ax.set_xticks([0,0.1,0.2,.3,.4])
    ax.set_xticklabels(['0%','10%','20%','30%','40%'])
    ax.tick_params(axis='both', which='major', labelsize=15)
    plt.legend()
    plt.xlim(0,0.2)
    plt.show()

   


def plot_reference_divergences_NT():#not synonymous
    ax=plt.gca()
    ax.axvspan(0.0, 0.06,hatch='\ ',alpha=0.3,label='CN2',color='k',fill=0)
    ax.axvspan(0.03,0.09, hatch='/',alpha=0.3,label='HLII',color='k',fill=0)
    ax.axvspan(0.10,0.40,alpha=0.3,label='HLI-HLII',color='k',hatch='. ',fill=0)

def plot_reference_divergences_NT_2():#not synonymous
    ax=plt.gca()
    ax.axvspan(0.0, 0.01,alpha=0.2,label='Cluster Backbones',color='k')
    ax.axvspan(0.10,0.40,alpha=0.2,label='HLI-HLII',color='m')
    
    
    
def get_background_divergences(data_csv,infile_matrices,compare_celllist,overallcelllist,corebool):
    #get genenums from data_csv
    #go to infile_matrices. for each gene, load the matrix, get all pairdivs in the compare_celllist. also get the average of these divs.
    genelist=get_column_from_csv(data_csv,'GeneNum',column_type='int',corebool=corebool)
    cell_indices=[]
    for i in range(len(compare_celllist)):
        cell_indices.append(overallcelllist.index(compare_celllist[i]))
    
    background_divs_allpairs,background_divs_averagepergene=get_comparison_divergences_new(infile_matrices,cell_indices,genelist)
    return background_divs_allpairs,background_divs_averagepergene

def get_core_list():
    
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as myinfile:
        mylist=ast.literal_eval(myinfile.read())
    corelist=[]
    for i in range(len(mylist)):
        corelist.append(int(mylist[i]))
    return corelist

def get_column_from_csv(data_csv,column_header,column_type='str',corebool=False):
    if corebool:
        corelist=get_core_list()
    genelist=[]
    with open(data_csv,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            if corebool:
                gnum=int(row['GeneNum'])
                if gnum not in corelist:
                    continue
            if column_type=='int':
                genelist.append(int(row[column_header]))
            if column_type=='str':
                genelist.append(row[column_header])
            if column_type=='float':
                genelist.append(float(row[column_header]))
    return genelist

def get_comparison_divergences(infile_matrices,compare_cell_indices,genelist):
    background_divs_allpairs=[]
    background_divs_averagepergene=[]

    with open(infile_matrices,'r') as myinfile:
         l=myinfile.readline()       
         mat=[]
         genenum=1
         while len(l.rstrip())!=0:
             l=myinfile.readline()
             if len(l)==0:
                 break
             if l[0]=='#':
                 ind1=l.find('(')
                 if genenum in genelist:
                     print('Gene %i' %genenum)
                
                     #now compute on matrix mat
                     ###########################
                     single_gene_divs,single_gene_average=get_divs_from_gene(mat,compare_cell_indices,background_divs_allpairs)
                     update_div_lists(background_divs_allpairs,background_divs_averagepergene,single_gene_divs,single_gene_average)
                     ###########################
                
                 genenum=int(l[7:ind1-1])#for NON synonymous
                 ind2=l.find(' ',11)
                 #genenum=int(l[11:ind2])#for synonymous
                 mat=[]
                 l=myinfile.readline()
             mat.append(list(map(float, l.split())))
    return background_divs_allpairs,background_divs_averagepergene

def get_gene_from_key_to_check(key):
    idx1=key.find('_')
    idx2=key.find('_',idx1+1)
    g=int(key[idx1+1:idx2])
    return g
                
def load_arrays(arrayfile):
    #load into a dict
    genearrays={}
    npfile=np.load(arrayfile,allow_pickle=1)
    mykeys=npfile[()].keys()
    for item in mykeys:
        genearrays[item]=npfile[()][item]
    return genearrays


def get_comparison_divergences_new(infile_matrices,compare_cell_indices,genelist):
    background_divs_allpairs=[]
    background_divs_averagepergene=[]


    genearrays=distmatrepo.load_arrays(infilematrices)
    for item in list(genearrays.keys()):
        mygene=get_gene_from_key_to_check(item)
        if mygene in genelist:
            mat=genearrays[item]
            single_gene_divs,single_gene_average=get_divs_from_gene(mat,compare_cell_indices,background_divs_allpairs)
            update_div_lists(background_divs_allpairs,background_divs_averagepergene,single_gene_divs,single_gene_average)

  
    return background_divs_allpairs,background_divs_averagepergene


def get_divs_from_gene(mat,compare_cell_indices,background_divs_allpairs):
    single_gene_divs=[]
    for i in range(len(compare_cell_indices)):
        idxi=compare_cell_indices[i]
        for j in range(i):
            idxj=compare_cell_indices[j]
            div=mat[idxi][idxj]
            if div<1:
                single_gene_divs.append(div)    
    if len(single_gene_divs)>0:
        single_gene_average=float(sum(single_gene_divs))/len(single_gene_divs)
    if len(single_gene_divs)==0:
        single_gene_average='-'
    return single_gene_divs,single_gene_average

def update_div_lists(background_divs_allpairs,background_divs_averagepergene,single_gene_divs,single_gene_average):
    if single_gene_average=='-':
        return
    background_divs_averagepergene.append(single_gene_average)
    for item in single_gene_divs:
        background_divs_allpairs.append(item)




######################################
#########################################

def Check_for_High_Divergences(data_csv,corebool=1):
    if corebool:
        corelist=get_core_list()
    genelist=[]
    with open(data_csv,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            gnum=int(row['GeneNum'])
            if corebool:
                if gnum not in corelist:
                    continue
            mydiv=float(row['Top_Merge_Div'])
            if mydiv>0.25:
                print('Gene %i Top Merge %f' %(gnum,mydiv))
    return genelist

############
#Synteny checking


def Check_Local_Synteny(celllist,center_gene,gene_range=1):
    print('Gene %i' %center_gene)
    querydict=fill_querydict(celllist,center_gene,gene_range)
    #print(querydict)
    #t=input('t')
    querynum_dict,contig_dict=parse_querydict(querydict)
    #print(querynum_dict)
    #t=input('t')
    #print(contig_dict)
    #t=input('t')
    all_same_contig,uncovered_pair_same_contig,pair_same_contig=check_synteny_of_parsed_querydict(querynum_dict,contig_dict,querydict)
    return all_same_contig,uncovered_pair_same_contig,pair_same_contig

def fill_querydict(celllist,center_gene,gene_range):
    querydict={}
    for cell in celllist:
        querydict[cell]=[]
    for i in list(range(center_gene-gene_range,center_gene+gene_range+1)):
        get_queries_from_genefile(i,querydict)
    return querydict

def get_queries_from_genefile(genenum,querydict):
    celllist=list(querydict.keys())
    #get the queries and put them into querydict (cell:[list of all full query names at genes of interest, or '-' if uncovered]). process the queries later (get contig and genenum)
    genefile='/scratch/groups/dsfisher/Prochlorococcus/Berube_Kashtan_Orthologous_Groups/NEW_SCRIPTS_DOWNSAMPLING_7_12/Aligned_Gene_Groups/Gene%s_Kashtan_and_HLII.fasta' %(str(genenum))
    
    if not exists(genefile):
        return
    get_query_names(genefile,celllist,querydict)

def get_query_names(genefile,celllist,querydict):
    with open(genefile,'r') as myinfile:
        seqs=myinfile.read()
    seqlist=re.split(">",seqs)
    for cell in celllist:
        foundbool=0
        for i in range(len(seqlist)):
            if cell in seqlist[i]:
                queryname=re.split('\n',seqlist[i])[0]
                querydict[cell].append(queryname)
                seqlist.remove(seqlist[i])
                foundbool=1
                break
        if foundbool==0:
            querydict[cell].append('-')

                
    
def parse_querydict(querydict):
    #want a dict with: cell:[[querynum,contig],[querynum2,contig2]...] or ['-','-'] if uncovered
    querynum_dict={}
    contig_dict={}
    for cell in querydict:
        contiglist=[]
        querylist=[]
        oldlist=querydict[cell]
        #print(oldlist)
        for item in oldlist:
            if item=='-':
                querylist.append('-')
                contiglist.append('-')
                continue
            ind1=item.find(' ')
            querynum=int(item[1:ind1])
            #ind2=item.find('_')
            #ind3=item.find('_',ind2+1)
            ind4=item.find('(Ga')
            ind5=item.find('_',ind4)
            ind6=item.find(')',ind5)
            contignum=int(item[ind5+1:ind6])
            #contignum=int(item[ind2+1:ind3])
            querylist.append(querynum)
            contiglist.append(contignum)
        querynum_dict[cell]=querylist
        contig_dict[cell]=contiglist
    return querynum_dict,contig_dict

def check_synteny_of_parsed_querydict(querynum_dict,contig_dict,querydict):
    #check for each cell
    all_same_contig=[]
    pair_same_contig=[]
    uncovered_pair_same_contig=[]#all on same contig but only 2 genes covered
    celllist=list(querynum_dict.keys())
    
    for cell in querynum_dict:
        querylist=querynum_dict[cell]
        contiglist=contig_dict[cell]
        print(cell,contiglist)
        t=input('t')
        if len(contiglist)==0:
            #print('zero length contig list cell %s' %cell)
            #t=input('t')
            continue
        if '-' in querylist:
            #check if at least covered things are syntenic
            #querylist=[querylist[j] for j in range(len(querylist)) if querylist[j]!='-']
            contiglist=[contiglist[j] for j in range(len(contiglist)) if contiglist[j]!='-']
            if len(contiglist)<2:
                continue
        if contiglist.count(contiglist[0])==len(contiglist):
            if len(contiglist)>2:
                all_same_contig.append(cell)
            if len(contiglist)>=2:
                uncovered_pair_same_contig.append(cell)
            
        
        if contiglist.count(contiglist[0])!=len(contiglist):
            if len(contiglist)<len(querylist):
                continue
            if contiglist.count(contiglist[1])>=2:
                pair_same_contig.append(cell)
        #print(contiglist,all_same_contig,uncovered_pair_same_contig,pair_same_contig)
        #t=input('t')
    return all_same_contig,uncovered_pair_same_contig,pair_same_contig


if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB')
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-d','--divcut',type=float)
    parser.add_argument('-cm','--clusteringmethod',type=str,default='average')
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    n=args.n
    divcut=args.divcut
    clusteringmethod=args.clusteringmethod
        
    
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype
        
        
    if n==3:
        celllist=clique1_without_ribotype
        groupname='Clique1'
        

   
    
    
    if prog==0:
        
        overallcelllist=get_celllist_allhlii()
        infile_matrices=''#use the gene matrices file here
        outfile=''
        plotbool=False
        pdf=''
        celllist2=celllist
        Produce_Clusters_File(infile_matrices,celllist,overallcelllist,outfile,divcut,groupname,celllist2,clusteringmethod=clusteringmethod,plotbool=plotbool,pdf=pdf,colorcelllists=colorcelllists)
