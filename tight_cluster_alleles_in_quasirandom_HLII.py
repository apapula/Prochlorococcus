import numpy as np
import csv
import matplotlib.pyplot as plt
import math
import statistics
import re
import ast
import networkx as nx
from networkx.algorithms.approximation import clique
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import dendrogram, linkage,fcluster
import scipy.cluster.hierarchy as hac
from matplotlib.backends.backend_pdf import PdfPages
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,AllKashBATS_without_ribotype, BB_closegrp_1,CN2_without_ribotype,Berube_HOT,Berube_BATS,BB_HOT,BB_BATS,Atlantic,Pacific, C3_without,C2_without,C4_without,C5_without, C1,C2,C3,C4,C5, HLII_NoClose
import matplotlib.pyplot as plt
import Repo_site_catalogues_single_pair as repoi
import Get_R_From_Cell_Pair_Matrices as rrepo
import Plotting_Repo as plotrepo
from mpl_toolkits.axes_grid1 import make_axes_locatable


'''
This program searches for alles from c1-c5 in quasi-random hlii cells, or from one of the 3 7-cell clusters within the quasirandom dataset (depending on the boolean toggle at the bottom of the script)
'''


def check_for_Berube_SAGs_Clustered_with_ribos(infilematrices,tfna,cluster_cells,group_cells,berube_celllist,groupname_list):
    if tfna in ['t','f','a']:
        genearrays=rrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=rrepo.load_two_arrays(infilematrices)
    overalllist=rrepo.overall_order_list()
    
    berube_indexlist,cluster_indexlist,group_indexlist=get_indexlists(berube_celllist,cluster_cells,group_cells,overalllist)
    berube_ribotypes,berube_ribogenes,berube_coveredgenes,clustering_genes,uncovered_gene_list=loop_matrices(genearrays,group_indexlist,cluster_indexlist,berube_indexlist,berube_celllist,overalllist,groupname_list)
    #above:
    #berube_ribotypes: a dictionary with cell: C1-C5
    #berube_ribogenes: dic with berube cell: gene list where clustered
    #berube_covered_genes:dict with berube cell: list of all genes covered that are considered here (not just ribogenes I think?)
    #clustering_genes: a list of the HLII core genes where ribos are clear clusters
    
    print('%i genes uncovered by ribos' %len(uncovered_gene_list))
    print(uncovered_gene_list)
    print('Genes that cluster along Ribotype Boundaries:')
    print(clustering_genes)

    print_berube_ribotypes(berube_ribotypes, clustering_genes,groupname_list)
    print_genes_for_each_cell(berube_ribotypes,clustering_genes,groupname_list,berube_ribogenes)
    found_genes,found_gene_cells=print_genes(clustering_genes,berube_ribotypes,berube_ribogenes,groupname_list)
    #found_genes: genenum of clustering gene if at least 1 cell ribocell: ribotypes of cells
    #found_gene_cells: gennum: cells that are ribocells, corresponding

    
    #want to get: mean number of clustering genes per sample location
    #get_mean_clustering_genes_per_sample(berube_celllist,berube_ribotypes)
    totncells=len(berube_celllist)
    make_plots(clustering_genes,berube_ribotypes,berube_ribogenes,found_genes,found_gene_cells,totncells)


def make_plots(clustering_genes,berube_ribotypes,berube_ribogenes,found_genes,found_gene_cells,totncells):
    ribocells=['AG-355-K20','AG-355-P11','AG-355-K13','AG-355-K23','AG-355-L02','AG-355-L20','AG-355-N23']
    #plot histogram of number of genes per cell, total and striated
    gpctot, gpc1,gpc2,gpc3,gpc4,gpc5=[],[],[],[],[],[]#genes per cell
    cpgtot,cpg1,cpg2,cpg3,cpg4,cpg5=[],[],[],[],[],[]#cells per gene
    for item in berube_ribotypes:
        l=berube_ribotypes[item]
        gpctot.append(len(l))
        gpc1.append(l.count('C1'))
        gpc2.append(l.count('C2'))
        gpc3.append(l.count('C3'))
        gpc4.append(l.count('C4'))
        gpc5.append(l.count('C5'))

    for item in found_genes:
        if item in ribocells:
            continue
        l=found_genes[item]
       
        cpgtot.append(len(l))
        cpg1.append(l.count('C1'))
        cpg2.append(l.count('C2'))
        cpg3.append(l.count('C3'))
        cpg4.append(l.count('C4'))
        cpg5.append(l.count('C5'))


    print('Cells per gene')
    print(cpgtot)
    myl=[item for item in cpgtot]
    for i in range(len(clustering_genes)-len(cpgtot)):
         myl.append(0)
    print('Mean %.4f median %.4f of %i genes when including genes with 0 ribocells' %(repo.avelist(cpgtot)*(float(len(cpgtot))/len(clustering_genes)),statistics.median(myl),len(clustering_genes)))
    print('Mean %.4f median %.4f of %i genes when EXcluding genes with 0 ribocells' %(repo.avelist(cpgtot),statistics.median(cpgtot),len(cpgtot)))

        
    if 1:#sample striation
        sampledict={}#sample numL list of number of ribogenes per cell
        for item in berube_ribotypes:
            if item in ribocells:
                continue
            samp=item[3:6]
            if samp in sampledict:
                sampledict[samp].append(len(berube_ribotypes[item]))
            if samp not in sampledict:
                sampledict[samp]=[len(berube_ribotypes[item])]
        mysamps=list(sampledict.keys())
        mydists=[]
        labs=[]
        for i in range(len(mysamps)):
            labs.append('%s\n(%i Cells)' %(mysamps[i],len(sampledict[mysamps[i]])))
            mydists.append(sampledict[mysamps[i]])
        plotseps=list(range(len(mydists)))
        plt.violinplot(mydists,positions=plotseps,showmedians=True)
        ax=plt.gca()
        ax.set_xticks(plotseps)
        ax.set_xticklabels(labs)
        plt.xlabel('Berube Sample Number')
        plt.ylabel('Number of C1-C5 genes per SAG from the Sample')
        plt.title('Number of C1-C5 genes (of %i examined) per Berube SAG as a function of sample\nExcluding 7 C1-C5 Berube SAGs' %len(clustering_genes))
        plt.show()

    print(berube_ribotypes)
    
    if 1:#ocean striation
        atlanticlist=['355','363','388','412','429','432','436']+['409','418','420','424','442','444']
        pacificlist=['347','402','335','459','469','449','455']+['311','315','316','321','323','331','341','345','450','463','670','673','676','679','683','686']
        oceandictlists=[[],[]]#atl, pac
        batslist=[]
        atlnotbats=[]
        for item in berube_ribotypes:
            if item in ribocells:
                continue
            samp=item[3:6]
            if samp in ['355','363']:
                batslist.append(len(berube_ribotypes[item]))
            if samp in atlanticlist:
                if samp not in ['355','363']:
                    atlnotbats.append(len(berube_ribotypes[item]))
                oceandictlists[0].append(len(berube_ribotypes[item]))
            if samp in pacificlist:
                oceandictlists[1].append(len(berube_ribotypes[item]))
        labs=['Atlantic\n(%i Cells)' %len(oceandictlists[0]),'Pacific\n(%i Cells)' %len(oceandictlists[1])]
       
        plotseps=list(range(len(oceandictlists)))
        plt.violinplot(oceandictlists,positions=plotseps,showmeans=True)
        ax=plt.gca()
        ax.set_xticks(plotseps)
        ax.set_xticklabels(labs)
        plt.xlabel('Sample Location',fontsize=25)
        plt.ylabel('Number of C1-C5 genes per SAG from the Sample')
        plt.title('Number of C1-C5 genes (of %i examined) per Berube SAG as a function of sample\nExcluding 7 C1-C5 Berube SAGs (Mean Shown)' %len(clustering_genes))
        plt.show()

        newlabs=['BATS\n(%i Cells)' %(len(batslist)),'Atlantic, not BATS\n(%i Cells)' %len(atlnotbats),'Pacific\n(%i Cells)' %len(oceandictlists[1])]
        plotseps=list(range(3))
        plt.violinplot([batslist,atlnotbats,oceandictlists[1]],positions=plotseps,showmeans=True)
        ax=plt.gca()
        ax.set_xticks(plotseps)
        ax.set_xticklabels(newlabs)
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.xlabel('Sample Location',fontsize=25)
        plt.ylabel('No. C1-C5 genes per SAG',fontsize=25)
        #plt.title('Number of C1-C5 genes (of %i examined) per Berube SAG as a function of sample\nExcluding 7 C1-C5 Berube SAGs (Mean Shown)' %len(clustering_genes))
        plt.show()

        
    if 1:
        mymax=max(gpctot)
        mybins=[-.5]
        for i in range(mymax+1):
            mybins.append(i+.5)
        ribocells=[item for item in gpctot if item <20]
        plt.hist(gpctot,bins=mybins)
        plt.title('Total Number of Genes per Berube HLII SAG that Cluster\nwith one of C1-C5\n Mean %.2f Median %.2f Cells per gene Excluding 7 ribocells' %(repo.avelist(ribocells),statistics.median(ribocells)))#\n(%i Cells with 0 Clustering Genes Not shown)' %(totncells-len(list(berube_ribotypes.keys()))))
        plt.xlabel('Number of Genes (of %i Used)' %(len(clustering_genes)))
        plt.ylabel('Number of Berube HLII SAGs (167 total)')
        plt.show()

        plt.hist(gpctot,bins=mybins)
        plt.title('Total Number of Genes per Berube HLII SAG that Cluster\nwith one of C1-C5')#\n(%i Cells with 0 Clustering Genes Not shown)' %(totncells-len(list(berube_ribotypes.keys()))))
        plt.xlabel('Number of Genes (of %i Used)' %(len(clustering_genes)))
        plt.ylabel('Number of Berube HLII SAGs (167 total)')
        plt.yscale('log')
        plt.show()

        plt.hist([gpc1,gpc2,gpc3,gpc4,gpc5],label=['C1','C2','C3','C4','C5'],bins=mybins)
        plt.title('Number of Genes per Berube HLII SAG that Cluster\nwith one of C1-C5\n Mean %.2f Median %.2f Cells per gene Excluding 7 ribocells' %(repo.avelist(ribocells),statistics.median(ribocells)))#\n(%i Cells with 0 Clustering Genes Not shown)' %(totncells-len(list(berube_ribotypes.keys()))))
        plt.legend()
        plt.xlabel('Number of Genes (of %i Used)' %(len(clustering_genes)))
        plt.ylabel('Number of Berube HLII SAGs (167 total)')
        plt.show()

        plt.hist([gpc1,gpc2,gpc3,gpc4,gpc5],label=['C1','C2','C3','C4','C5'],bins=mybins)
        plt.title('Number of Genes per Berube HLII SAG that Cluster\nwith one of C1-C5\n Mean %.2f Median %.2f Cells per gene Excluding 7 ribocells' %(repo.avelist(ribocells),statistics.median(ribocells)))#\n(%i Cells with 0 Clustering Genes Not shown)' %(totncells-len(list(berube_ribotypes.keys()))))
        plt.legend()
        plt.yscale('log')
        plt.xlabel('Number of Genes (of %i Used)' %(len(clustering_genes)))
        plt.ylabel('Number of Berube HLII SAGs (167 total)')
        plt.show()

    if 1:#genes per cell 2
        plt.plot(gpctot,gpc1,'o',label='C1',alpha=.5)
        plt.plot(gpctot,gpc2,'o',label='C2',alpha=.5)
        plt.plot(gpctot,gpc3,'o',label='C3',alpha=.5)
        plt.plot(gpctot,gpc4,'o',label='C4',alpha=.5)
        plt.plot(gpctot,gpc5,'o',label='C5',alpha=.5)
        plt.legend()
        plt.xlabel('Total Number of Genes Clustering with C1-C5 per cell')
        plt.ylabel('Number of genes clustering with each of C1-C5 per cell')
        plt.title('Berube HLII SAGs Number of Genes Clustering with C1-C5 per cell')
        plt.show()

    if 1:

       

        
        mymax=max(cpgtot)
        mybins=[-.5]
        for i in range(mymax+1):
            mybins.append(i+.5)
        #plot cells per gene
        plt.hist(cpgtot,bins=mybins)
        plt.title('Total Number of Berube HLII SAGs per Gene that Cluster\nwith one of C1-C5, Excluding 7 C1-C5 Cells\n(%i Genes of %i with 0 Clustering Cells Not shown)' %(len(clustering_genes)-len(list(found_genes.keys())),len(clustering_genes)))
        plt.xlabel('Number of non-C1-C5 HLII Cells (of 160)')       
        plt.ylabel('Number of Genes with Clear C1-C5 Clusters (%i total)' %(len(clustering_genes)))
        plt.show()

       
        plt.hist(cpgtot,bins=mybins)
        plt.title('Total Number of Berube HLII SAGs per Gene that Cluster\nwith one of C1-C5, Excluding 7 C1-C5 Cells\n(%i Genes of %i with 0 Clustering Cells Not shown)' %(len(clustering_genes)-len(list(found_genes.keys())),len(clustering_genes)))
        plt.xlabel('Number of non-C1-C5 HLII Cells (of 160)')       
        plt.ylabel('Number of Genes with Clear C1-C5 Clusters (%i total)' %(len(clustering_genes)))
        plt.yscale('log')
        plt.show()

        plt.hist([cpg1,cpg2,cpg3,cpg4,cpg5],label=['C1','C2','C3','C4','C5'],bins=mybins)
        plt.title('Number of Berube HLII SAGs per Gene that Cluster\nwith one of C1-C5\nwith one of C1-C5, Excluding 7 C1-C5 Cells\n(%i Genes of %i with 0 Clustering Cells Not shown)' %(len(clustering_genes)-len(list(found_genes.keys())),len(clustering_genes)))
        plt.xlabel('Number of non-C1-C5 HLII Cells (of 160)')
        plt.legend()
        plt.ylabel('Number of Genes with Clear C1-C5 Clusters (%s total)' %(len(clustering_genes)))
        plt.show()
    

def loop_matrices(genearrays,group_indexlist,cluster_indexlist,berube_indexlist,berube_celllist,matrixorderlist,group_namelist):#group_indexlist is index of ribotypes of which we want to check membership for each gene, for berube cells. cluster_indexlist is ribotypes that we want to form separated clusters at each gene to consider the gene
    berube_ribotypes,berube_ribogenes,berube_coveredgenes=initiate_dictionaries(berube_celllist)
    clustering_genes=[]
    uncovered_gene_list=[]
    for i in range(1908):
        mykey='Gene_%i_%s_Pairs' %(i,tfna)
        #print(mykey)
        #print(list(genearrays.keys()))
        if mykey not in genearrays:
            continue
        print(i)
        mat=genearrays[mykey]
        compute_on_individual_matrix(mat,i,berube_ribotypes,berube_ribogenes,berube_coveredgenes,matrixorderlist,berube_indexlist,group_indexlist,cluster_indexlist,clustering_genes,group_namelist,uncovered_gene_list)

              
    return berube_ribotypes,berube_ribogenes,berube_coveredgenes,clustering_genes,uncovered_gene_list


    
def print_berube_ribotypes(berube_ribotypes, clustering_genes,groupname_list):
    print('%i Clustering Genes' %len(clustering_genes))
    no_genes=0
    for item in berube_ribotypes:
        l=berube_ribotypes[item]
        if len(l)==0:
            no_genes+=1
            continue
        mys=item+' %i Genes total '%len(l)
        for i in range(len(groupname_list)):
            name=groupname_list[i]
            mys+=' %i %s Genes; ' %(l.count(name),name)
        print(mys)
    print('%i cells with no clustering genes' %no_genes)
    
def print_genes(clustering_genes,berube_ribotypes,berube_ribogenes,groupname_list):
    ribocells=['AG-355-K20','AG-355-P11','AG-355-K13','AG-355-K23','AG-355-L02','AG-355-L20','AG-355-N23']
    #number of genes with any cells (and fraction)
    #list of cells per gene x ribo
    found_genes={}#genenum: ribotypes of genes
    found_gene_cells={}#genenum:cells
    for item in berube_ribotypes:
        if item in ribocells:
            continue
        l=berube_ribotypes[item]
        l2=berube_ribogenes[item]
        if len(l)==0:
            continue
        for g in range(len(l2)):
            ge=l2[g]
            if ge in found_genes:
                found_genes[ge].append(l[g])
                found_gene_cells[ge].append(item)
            if ge not in found_genes:
                found_genes[ge]=[l[g]]
            if ge not in found_gene_cells:
                found_gene_cells[ge]=[item]
    print(sorted(list(found_genes.keys())))
    for item in sorted(list(found_genes.keys())):
        if len(found_genes[item])<3:
            continue
        print('Gene %s' %str(item))
        print(found_gene_cells[item])
        print(found_genes[item])
    print('%i genes used, %i have ribo-cells (%f)' %(len(clustering_genes),len(list(found_genes.keys())),float(len(list(found_genes.keys())))/len(clustering_genes)))
    return found_genes,found_gene_cells

def print_genes_for_each_cell(berube_ribotypes,clustering_genes,groupname_list,berube_ribogenes):
    print('%i Clustering Genes' %len(clustering_genes))
    for item in berube_ribotypes:
        l=berube_ribotypes[item]
        l2=berube_ribogenes[item]
        if len(l)==0:
            continue
        mys=item+' %i Genes total '%len(l)
        print(mys)
        for i in range(len(groupname_list)):
            name=groupname_list[i]
            if l.count(name)>0:
                glist=[]
                for k in range(len(l)):
                    if l[k]==name:
                        glist.append(l2[k])
                print(' %i %s Genes; ' %(l.count(name),name),glist)
        print(mys)

def get_indexlists(berube_celllist,cluster_cells,group_cells,matrixorderlist):
    print(matrixorderlist)
    berube_indexlist=[]
    for i in range(len(berube_celllist)):
        berube_indexlist.append(matrixorderlist.index(berube_celllist[i]))
    group_indexlist=[]
    for i in range(len(group_cells)):
        l=group_cells[i]
        idxl=[]
        for j in range(len(l)):
            if l[j] not in matrixorderlist:
                #if l[j][:-3] not in matrixorderlist:
                print('%s Not in matrix order list' %l[j])
                t=input('t')
                continue
            idxl.append(matrixorderlist.index(l[j]))
        group_indexlist.append(idxl)
    cluster_indexlist=[]
    for i in range(len(cluster_cells)):
        l=cluster_cells[i]
        idxl=[]
        for j in range(len(l)):
            if l[j] not in matrixorderlist:
                print('%s Not in matrix order list' %l[j])
                t=input('t')
                continue
            idxl.append(matrixorderlist.index(l[j]))
        cluster_indexlist.append(idxl)
    return berube_indexlist,cluster_indexlist,group_indexlist


def get_indexlists_old(berube_celllist,cluster_cells,group_cells,matrixorderlist):
    print(matrixorderlist)
    berube_indexlist=[]
    for i in range(len(berube_celllist)):
        berube_indexlist.append(matrixorderlist.index(berube_celllist[i]))
    group_indexlist=[]
    for i in range(len(group_cells)):
        l=group_cells[i]
        idxl=[]
        for j in range(len(l)):
            if l[j][:-3] not in matrixorderlist:
                #if l[j][:-3] not in matrixorderlist:
                print('%s Not in matrix order list' %l[j])
                t=input('t')
                continue
            idxl.append(matrixorderlist.index(l[j][:-3]))
        group_indexlist.append(idxl)
    cluster_indexlist=[]
    for i in range(len(cluster_cells)):
        l=cluster_cells[i]
        idxl=[]
        for j in range(len(l)):
            if l[j][:-3] not in matrixorderlist:
                print('%s Not in matrix order list' %l[j])
                t=input('t')
                continue
            idxl.append(matrixorderlist.index(l[j][:-3]))
        cluster_indexlist.append(idxl)
    return berube_indexlist,cluster_indexlist,group_indexlist


def initiate_dictionaries(berube_celllist):
    #want: dictoinary that for each berube cell, is list of its 'ribotype', and a list of the locations of these genes
    berube_ribotypes={}
    berube_ribogenes={}
    berube_coveredgenes={}
    for item in berube_celllist:
        berube_ribotypes[item]=[]
        berube_ribogenes[item]=[]
        berube_coveredgenes[item]=[]
    return berube_ribotypes,berube_ribogenes,berube_coveredgenes

def compute_on_individual_matrix(mat,genenum,berube_ribotypes,berube_ribogenes,berube_coveredgenes,matrixorderlist,berube_indexlist,group_indexlist,cluster_indexlist,clustering_genes,group_namelist,uncovered_gene_list):
    print('Gene %i' %genenum)
  
    check_coverage(mat,berube_coveredgenes,berube_indexlist,matrixorderlist,genenum)
   
    new_index_list,clusteredbool=check_if_clustered_with_outliers(mat,clustering_genes,cluster_indexlist,genenum,uncovered_gene_list,matrixorderlist)
    if not clusteredbool:
        return
    clustering_genes.append(genenum)
   
    check_if_berube_cells_in_ribotypes(mat,genenum,berube_ribotypes,berube_ribogenes,berube_indexlist,new_index_list,matrixorderlist,group_namelist)

def check_if_berube_cells_in_ribotypes(mat,genenum,berube_ribotypes,berube_ribogenes,berube_indexlist,group_indexlist,matrixorderlist,group_namelist):#grouplist and group_indexlist are in same order. 
    #want the berube cell max distance to be less than the group diam
    for i in range(len(berube_indexlist)):
        idx=berube_indexlist[i]
        for j in range(len(group_indexlist)):
            gindex=group_indexlist[j]
            if len(gindex)==0:
                continue
            mydiam=mymax(get_group_distances(mat,gindex,gindex))
            if mydiam=='-':
                continue
            berubemax=mymax(get_group_distances(mat,gindex,[idx]))
            if berubemax=='-':
                continue#should not change results
            if berubemax<=mydiam:
                groupname=group_namelist[j]
                berubecell=matrixorderlist[idx]
                berube_ribotypes[berubecell].append(groupname)
                berube_ribogenes[berubecell].append(genenum)
                break



                          
def mymax(mylist):
    #returns max unless list is empty, in which case returns '-'
    if len(mylist)==0:
        return '-'
    return max(mylist)

def mymin(mylist):
    #returns max unless list is empty, in which case returns '-'
    if len(mylist)==0:
        return '-'
    return min(mylist)

def check_coverage(mat,berube_coveredgenes,berube_indexlist,matrixorderlist,genenum):
    for i in range(len(berube_indexlist)):
        idx=berube_indexlist[i]
        mycol=mat[idx]
        mycol=list(mycol[:idx])+list(mycol[idx:])
        mycol=np.array(mycol)
        if np.all(mycol==2.):
            continue
        else:
            berubecell=matrixorderlist[idx]
            berube_coveredgenes[berubecell].append(genenum)


def check_if_clustered_with_outliers(matrix,clustering_genes,cluster_indexlist,genenum,uncovered_gene_list,overalllist,divcut=0.02,minority_sizelist=[1,2,3,4,5], majority_minority_ratio=2):
    new_cluster_indexlist=[]
    #for each group, make matrix and cluster it and cut it. if >1 cluster, check if meets minority-size and maj/min ration criteria
    for i in range(len(cluster_indexlist)):
        clist=cluster_indexlist[i]
        main_cluster=return_ribotype_clusters_with_divcut(matrix,clist,divcut,minority_sizelist,majority_minority_ratio,overalllist)
        if main_cluster==['-']:
            new_cluster_indexlist.append([])
            continue
            #continue
        #if len(main_cluster)==0:
            #print('No main cluster')
            #return [], False
        new_cluster_indexlist.append(main_cluster)
    if new_cluster_indexlist==[[] for i in range(len(cluster_indexlist))]:
        print('Uncovered')
        return [], False
    clusteredbool=check_if_clustered(matrix,clustering_genes,new_cluster_indexlist,genenum,uncovered_gene_list)
    return new_cluster_indexlist,clusteredbool

def return_ribotype_clusters_with_divcut(matrix,clusterlist,divcut,minority_sizelist,majority_minority_ratio,overalllist,clusteringmethod='average'):
    #get covered cells
    
    
    wgbool=0
    lengthcut=11
    #cell_indices=[overalllist.index(clusterlist[i]) for i in range(len(clusterlist))]
    cluscells=[overalllist[clusterlist[i]] for i in range(len(clusterlist))]
    submat,coveredcells=plotrepo.get_submatrix(matrix,clusterlist,cluscells,overalllist,lengthcut,wgbool)#only works for discrete
    if len(coveredcells)<2:
        return ['-']
    covered_indices=[overalllist.index(coveredcells[i]) for i in range(len(coveredcells))]
    
    myarray=np.array(submat)
    print(myarray)
    myarrays=ssd.squareform(myarray)
    clustering=hac.linkage(myarrays,method='%s' %clusteringmethod,metric='precomputed')
    myclusters=fcluster(clustering,divcut,criterion='distance')
    labels    = fcluster(clustering,divcut,criterion='distance')-1
    #print('Labels', labels)
    clusters  = [[] for i in range(max(labels) + 1)]
    for i in range(len(labels)):
        cellidx=covered_indices[i]
        clusters[labels[i]].append(cellidx)
    #print('Clusters',clusters)
    #now check if meets criteria
    main_cluster=check_if_ribotype_clusters_with_divcut_meet_criteria(clusters,minority_sizelist,majority_minority_ratio)
    return main_cluster

    
def check_if_ribotype_clusters_with_divcut_meet_criteria(clusters,minority_sizelist,majority_minority_ratio):
    if len(clusters)==1:
        return clusters[0]
    print('%i clusters' %len(clusters))  
    cluster_sizes=[]
    for i in range(len(clusters)):
        cluster_sizes.append(len(clusters[i]))
    sorted_cluster_sizes=sorted(cluster_sizes)
    print(sorted_cluster_sizes)
    #t=input('t')
    if sorted_cluster_sizes[-2]>max(minority_sizelist):
        print('Subcluster too large')
        return []
    if sorted_cluster_sizes[-2]*majority_minority_ratio>sorted_cluster_sizes[-1]:
        print('Main cluster too small')
        return []
    mylen=sorted_cluster_sizes[-1]
    for item in clusters:
        if len(item)==mylen:
            return item


def check_if_clustered(matrix,clustering_genes,cluster_indexlist,genenum,uncovered_gene_list,uncovered_cutoff=.9):#if more than 0.6 of groups are uncovered, return
    #print(cluster_indexlist)
    clusteredbool=1
    diamlist=[]
    mininterlist=[]
    uncovered_groups=0
    print('%i clustering genes' %(len(clustering_genes)))
    for i in range(len(cluster_indexlist)):
        c=cluster_indexlist[i]
        if len(c)==0:
            continue
        mydiam=mymax(get_group_distances(matrix,c,c))
        if mydiam=='-':
            uncovered_groups+=1
            if uncovered_groups>uncovered_cutoff*len(cluster_indexlist):
                uncovered_gene_list.append(genenum)
                return 0
            continue
        diamlist.append(mydiam)
        alter_list=[]
        for j in range(len(cluster_indexlist)):
            if j!=i:
                for item in cluster_indexlist[j]:
                    alter_list.append(item)
                intermin=mymin(get_group_distances(matrix,c,alter_list))
                if intermin=='-':
                    print('No interribo distances')
                    #t=input('t')
                    return 0#unsure what to do here. only 1 ribotype is covered.
                mininterlist.append(intermin)
                if intermin<mydiam:
                    print(intermin,mydiam)
                    #t=input('t')
                    clusteredbool=0
                    return clusteredbool
    return clusteredbool

def get_group_distances(matrix,cellindices1,cellindices2):#cellindices1,2 could be same if doing intra-group. cellindices1 could be a single (berube) cell
    distances=[]
    for i in cellindices1:
        for j in cellindices2:
            if i==j:
                continue
            [n,l]=matrix[i][j]
            dist=float(n)/l
            if dist<1 and l>11:
                distances.append(dist)
    return distances


##################################################

def check_for_Berube_SAGs_Clustered_with_Berube_clusters(infilematrices,tfna,cluster_cells,celllist_berube,cluster_groupname,noclosecells):
    outgroup_cells=[item for item in noclosecells if item not in cluster_cells]#celllist_berube if item not in cluster_cells]
    if tfna in ['t','f','a']:
        genearrays=rrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=rrepo.load_two_arrays(infilematrices)
    overalllist=rrepo.overall_order_list()
    outfile='../input_data/Transferred_Alleles_%s.csv' %cluster_groupname
    fieldnames=['GeneNum','Cluster_Diameter','Clustering_BerubeCells']
    genes,diams,clusteringcells=[],[],[]
    cell_indices=[overalllist.index(celllist_berube[i]) for i in range(len(celllist_berube))]
    #get submatrix
    #get cluster diam
    #all outgroup cells: get their max div among all cluster cells
    for i in range(1908):
        mykey='Gene_%i_%s_Pairs' %(i,tfna)
        if mykey not in genearrays:
            continue
        print(i)
        mat=genearrays[mykey]
        lengthcut=11
        wgbool=False
        submat,neworder_cells=plotrepo.get_submatrix(mat,cell_indices,celllist_berube,overalllist,lengthcut,wgbool)
        mydiam=get_group_diameter(submat,neworder_cells,cluster_cells)
        if mydiam=='-':
            continue
        genes.append(i)
        diams.append(mydiam)
        #now get clusteringcells
        clusteringcells.append(get_clustering_cells(submat,neworder_cells,cluster_cells,mydiam,outgroup_cells))
    #now write
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=fieldnames)
        outwriter.writeheader()
        for i in range(len(genes)):
            outwriter.writerow({'GeneNum':genes[i],'Cluster_Diameter':diams[i],'Clustering_BerubeCells':str(clusteringcells[i])})


        
def get_clustering_cells(submat,neworder_cells,cluster_cells,diam,outcells):
    covccells=list(set(neworder_cells).intersection(set(cluster_cells)))
    covocells=list(set(neworder_cells).intersection(set(outcells)))
    clusteringcells=[]
    for i in range(len(covocells)):
        idxo=neworder_cells.index(covocells[i])
        cello=covocells[i]
        thesedivs=[]
        for j in range(len(covccells)):
            idxc=neworder_cells.index(covccells[j])
            div=submat[idxc][idxo]
            thesedivs.append(div)
        if max(thesedivs)<=diam:
            clusteringcells.append(cello)
    return clusteringcells

def get_group_diameter(submat,neworder_cells,cluster_cells):
    divs=[]
    for i in range(len(cluster_cells)):
        celli=cluster_cells[i]
        if celli not in neworder_cells:
            continue
        myindex=neworder_cells.index(celli)
        for j in range(i):
            cellj=cluster_cells[j]
            if cellj not in neworder_cells:
                continue
            myindexj=neworder_cells.index(cellj)
            div=submat[myindex][myindexj]
            if div<1:
                divs.append(div)
    if len(divs)<2:
        return '-'
    return max(divs)

def get_group_diameter_med(submat,neworder_cells,cluster_cells):
    divs=[]
    for i in range(len(cluster_cells)):
        celli=cluster_cells[i]
        if celli not in neworder_cells:
            continue
        myindex=neworder_cells.index(celli)
        for j in range(i):
            cellj=cluster_cells[j]
            if cellj not in neworder_cells:
                continue
            myindexj=neworder_cells.index(cellj)
            div=submat[myindex][myindexj]
            if div<1:
                divs.append(div)
    if len(divs)<2:
        return '-','-'
    return max(divs),statistics.median(divs)
###################################################


####################################################



   
##############
tfna='a'

suffix=''
if tfna=='a':
    outfilewhole='Matrices_Berube_Kash_All_NT%s' %suffix
outfilewhole='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfilewhole

outfilewhole='../input_data/'+outfilewhole

infilematrices=outfilewhole+'.npy'



if 0:
    myfile='../input_data/HLII_clonelist2.txt'
    with open(myfile,'r') as infile:
        berube_celllist=[line.rstrip() for line in infile if len(line)>0]
    
    cluster_cells=[C1_without_ribotype,C2_without,C3_without,C4_without,C5_without]
        
    group_cells=cluster_cells
    groupname_list=['C1','C2','C3','C4','C5']
    
    check_for_Berube_SAGs_Clustered_with_ribos(infilematrices,tfna,cluster_cells,group_cells,berube_celllist,groupname_list)



if 0:

    myfile='../input_data/HLII_clonelist2.txt'
    with open(myfile,'r') as infile:
        celllist_berube=[line.rstrip() for line in infile if len(line)>0]
    if 0:
        cluster_cells=['AG-355-N23', 'AG-355-P11', 'AG-418-O02', 'AG-449-D22', 'AG-355-O19', 'AG-355-K20', 'AG-355-L20']#contains 355K20, 355P11,355n23, which are ribocells
        cluster_groupname='NonBB1_mostlyAtl'

    if 0:
        cluster_cells=['AG-355-A02', 'AG-424-A14', 'AG-432-K16', 'AG-418-J19', 'AG-424-G03', 'AG-355-L21', 'AG-355-B18']
        cluster_groupname='NonBB2_allAtl'

    if 1:
        cluster_cells=['AG-355-N22', 'AG-347-N19', 'AG-347-M08', 'AG-355-M02', 'AG-355-L22', 'AG-424-J22', 'AG-429-E20']
        cluster_groupname='BB1_mixed'

    cluster_groupname+='_NoClose'
    
    check_for_Berube_SAGs_Clustered_with_Berube_clusters(infilematrices,tfna,cluster_cells,celllist_berube,cluster_groupname,HLII_NoClose)#celllist_berube,cluster_groupname)





