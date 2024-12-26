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
import Repo_site_catalogues_single_pair as repo
import Make_Synonymous_Distance_Matrices_Berube_Kashtan_Using_All_Sites_Per_Pair_and_HalfGene_Record_sites_SNPs as distmatrepo
import Plotting_Repo as plotrepo

from mpl_toolkits.axes_grid1 import make_axes_locatable





def check_for_Berube_SAGs_Clustered_with_ribos(infilematrices,tfna,cluster_cells,group_cells,berube_celllist,groupname_list):
    if tfna in ['t','f','a']:
        genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=distmatrepo.load_two_arrays(infilematrices)
    overalllist=distmatrepo.overall_order_list()
    
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
        genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=distmatrepo.load_two_arrays(infilematrices)
    overalllist=distmatrepo.overall_order_list()
    #here, we can't have the clustering criterion of e.g. all 5 ribotypes cluster. we could ask: the cluster is less than some certain diameter (e.g. its WG divergence?) but this will mess up for conserved genes. at first pass we can just output the cluster diameter and weed out genes with high diameter. then any berube cells within the cluster are considered having the cluster's allele
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

def analyze_other_clusters(closegroup_name,diamcut,totcells,clustercells,corebool):
    #['GeneNum','Cluster_Diameter','Clustering_BerubeCells']
    diameters=[]#plot diameter distribution
    outfile='../input_data/Transferred_Alleles_%s.csv' %closegroup_name
    print(outfile)
    cells_per_gene=[]
    genes_per_cell={}
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]  

    for item in totcells:
        if item not in clustercells:
            genes_per_cell[item]=0
    with open(outfile,'r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            g=int(row['GeneNum'])
            if corebool:
                if g not in mycores:
                    continue
            diameters.append(float(row['Cluster_Diameter']))
            if diameters[-1]>diamcut:
                continue
            mycells=ast.literal_eval(row['Clustering_BerubeCells'])
            cells_per_gene.append(len(mycells))
            for item in mycells:
                genes_per_cell[item]+=1
            
    if corebool:
        closegroup_name+=' Core'
    if not corebool:
        closegroup_name+=' Including NonCore'
    if 0:
        plt.hist(diameters,bins=80)
        plt.title('%s Diameters of MIT9301 Genes' %closegroup_name)
        plt.xlabel('Number of Genes')
        plt.ylabel('Nucleotide Divergence Diameter')
        plt.show()
    #ribocells per gene
    #ribogenes per outcell
    #atl pac: ribogenes per cell atl-nonbats, bats, pac
    if 1:
        mybins=[-.5]
        for i in range(max(cells_per_gene)+1):
            mybins.append(i+.5)
        plt.hist(cells_per_gene,bins=mybins)
        plt.xlabel('Number of Cells Clustering with %s per Gene' %closegroup_name)
        plt.ylabel('Number of Genes (of %i)' %len(cells_per_gene))
        plt.yscale('log')
        plt.title('Alleles from %s in HLII:\nCells per gene (Mean %.1f Median %.1f)' %(closegroup_name,repo.avelist(cells_per_gene),statistics.median(cells_per_gene)))
        plt.show()

        gpclist=[genes_per_cell[item] for item in list(genes_per_cell.keys())]

        mybins=[-.5]
        for i in range(max(gpclist)+1):
            mybins.append(i+.5)
        plt.hist(gpclist,bins=40)#mybins)
        plt.xlabel('Number of Genes Clustering with %s per HLII Cell (of %i)' %(closegroup_name,len(cells_per_gene)))
        plt.ylabel('Number of HLII Cells')
        plt.title('Alleles from %s in HLII:\nGenes per cell (Mean %.1f Median %.1f)' %(closegroup_name,repo.avelist(gpclist),statistics.median(gpclist)))
        plt.show()

    if 1:#geo plot
        #batslist,atlnotbats,pac=[],[],[]
        atlanticlist=['355','363','388','412','429','432','436']+['409','418','420','424','442','444']
        pacificlist=['347','402','335','459','469','449','455']+['311','315','316','321','323','331','341','345','450','463','670','673','676','679','683','686']
        oceandictlists=[[],[]]#atl, pac
        batslist=[]
        atlnotbats=[]
        for item in genes_per_cell:
            samp=item[3:6]
            if samp in ['355','363']:
                batslist.append(genes_per_cell[item])
            if samp in atlanticlist:
                if samp not in ['355','363']:
                    atlnotbats.append(genes_per_cell[item])
                oceandictlists[0].append(genes_per_cell[item])
            if samp in pacificlist:
                oceandictlists[1].append(genes_per_cell[item])
        labs=['Atlantic\n(%i Cells)' %len(oceandictlists[0]),'Pacific\n(%i Cells)' %len(oceandictlists[1])]
       
        plotseps=list(range(len(oceandictlists)))
        plt.violinplot(oceandictlists,positions=plotseps,showmeans=True)
        ax=plt.gca()
        ax.set_xticks(plotseps)
        ax.set_xticklabels(labs)
        plt.xlabel('Sample Ocean')
        plt.ylabel('Number of %s genes per SAG from the Sample' %closegroup_name)
        plt.title('Number of %s genes (of %i examined) per Berube SAG as a function of sample\nExcluding 7 C1-C5 Berube SAGs (Mean Shown)' %(closegroup_name,len(cells_per_gene)))
        plt.show()

        newlabs=['BATS\n(%i Cells)' %(len(batslist)),'Atlantic, not BATS\n(%i Cells)' %len(atlnotbats),'Pacific\n(%i Cells)' %len(oceandictlists[1])]
        plotseps=list(range(3))
        plt.violinplot([batslist,atlnotbats,oceandictlists[1]],positions=plotseps,showmeans=True)
        ax=plt.gca()
        ax.set_xticks(plotseps)
        ax.set_xticklabels(newlabs)
        plt.xlabel('Sample Ocean')
        plt.ylabel('Number of %s genes per SAG from the Sample' %closegroup_name)
        plt.title('Number of %s genes (of %i examined) per Berube SAG as a function of sample\nExcluding 7 C1-C5 Berube SAGs (Mean Shown)' %(closegroup_name,len(cells_per_gene)))
        plt.show()
####################################################

def get_all_diams_and_meds(celllist,tfna,infilematrices,corebool,groupname):
    if tfna in ['t','f','a']:
        genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=distmatrepo.load_two_arrays(infilematrices)
    overalllist=distmatrepo.overall_order_list()
    diams=[]
    meds=[]
    closegenes=[]
    genes,ncells=[],[]
    cell_indices=[overalllist.index(celllist[i]) for i in range(len(celllist))]
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]
    for i in range(1908):
        if corebool:
            if i not in mycores:
                continue
        mykey='Gene_%i_%s_Pairs' %(i,tfna)
        if mykey not in genearrays:
            continue
        print(i)
        mat=genearrays[mykey]
        lengthcut=11
        wgbool=False
        submat,neworder_cells=plotrepo.get_submatrix(mat,cell_indices,celllist,overalllist,lengthcut,wgbool,uncoveredinput=1.5)
        mydiam,mymed=get_group_diameter_med(submat,neworder_cells,celllist)
        if mydiam=='-':
            continue
        meds.append(mymed)
        diams.append(mydiam)
        ncells.append(len(neworder_cells))
        genes.append(i)
        if mydiam<.1:
            print('Gene %i Diam %f' %(i,mydiam))
            closegenes.append(i)
  

    print('Close genes')
    print(closegenes)
    
    plt.plot(meds,diams,'o',alpha=.3)
    plt.xlabel('Median %s Divergence' %tfna)
    plt.ylabel('Maximum %s Divergence' %tfna)
    if not corebool:
        groupname+=' Including NonCore'
    if corebool:
        groupname+=' Core'
    plt.suptitle('%s %s Gene Diameters vs Median' %(groupname,tfna))
    plt.show()

    plt.plot(genes,diams,'.-')
    plt.title('%s %s Diameter along the genome' %(groupname,tfna))
    plt.xlabel('Gene Number')
    plt.ylabel('Maximum %s Divergence' %tfna)
    plt.show()
    
    plt.plot(diams,ncells,'o',alpha=.3)
    plt.xlabel('Maximum %s Divergence' %tfna)
    plt.ylabel('Number of cells covered' )
    plt.suptitle('%s %s Gene Coverage vs Diameters' %(groupname,tfna))
    plt.show()
    
    fig,ax=plt.subplots()
    plt.plot(meds, diams,'o',alpha=.3)
    plt.xlabel('Median %s Divergence' %tfna)
    plt.ylabel('Maximum %s Divergence' %tfna)
    divider=make_axes_locatable(ax)
    ax_histx=divider.append_axes("top",1.2,pad=0.1,sharex=ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)

    # make some labels invisible
    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)

    binsx,binsy=40,40
    
    ax_histx.hist(meds, bins=binsx)
    ax_histy.hist(diams, bins=binsy, orientation='horizontal')

    #ax_histx.set_yscale('log')
    #ax_histy.set_xscale('log')
    plt.suptitle('%s %s Gene Diameters vs Median' %(groupname,tfna))
    #fig.tight_layout()
   
    #plt.tight_layout()
    plt.show()

###############


def Official_Plot():

    berube_ribotypes={'AG-347-E23': [],'AG-347-G18': [],'AG-347-I04': ['C5'], 'AG-347-I06': ['C3'], 'AG-347-I19': [], 'AG-347-I21': [], 'AG-347-J14': [], 'AG-347-J19': ['C1'], 'AG-347-J20': [], 'AG-347-J21': ['C1', 'C1'], 'AG-347-K16': ['C1', 'C3'], 'AG-347-K17': ['C3', 'C5'], 'AG-347-K18': [], 'AG-347-K19': ['C1', 'C1', 'C1', 'C1', 'C3'], 'AG-347-K20': [], 'AG-347-L02': ['C3'], 'AG-347-L17': ['C1', 'C5'], 'AG-347-L20': ['C1', 'C3'], 'AG-347-M08': ['C1', 'C1'], 'AG-355-A09': [], 'AG-355-B18': [], 'AG-355-I04': ['C3', 'C5', 'C1', 'C1', 'C1', 'C4', 'C3'], 'AG-355-J09': ['C1', 'C1', 'C1'], 'AG-355-M02': ['C3'], 'AG-355-N02': ['C5', 'C1', 'C3'], 'AG-355-N16': ['C1', 'C1', 'C1', 'C1', 'C1'], 'AG-355-N22': ['C3'], 'AG-355-O17': ['C1'], 'AG-355-P07': [], 'AG-355-P16': ['C1', 'C5'], 'AG-388-A01': [], 'AG-402-I23': [], 'AG-402-K21': [], 'AG-402-L23': [], 'AG-402-M15': [], 'AG-402-N08': [], 'AG-402-O16': [], 'AG-412-J13': [], 'AG-418-D13': [], 'AG-418-I20': [], 'AG-418-M21': [], 'AG-424-E18': ['C1'], 'AG-442-B03': [], 'AG-442-N07': ['C5', 'C1'], 'AG-335-I15': [], 'AG-335-J02': [], 'AG-347-B08': [], 'AG-347-J05': [], 'AG-347-K21': ['C1'], 'AG-355-J17': ['C1'], 'AG-355-J21': ['C3'], 'AG-355-L21': [], 'AG-355-N21': [], 'AG-355-O19': ['C5', 'C3'], 'AG-355-P23': [], 'AG-363-A05': [], 'AG-402-O23': ['C5', 'C1'], 'AG-418-L19': [], 'AG-424-L22': ['C3'], 'AG-429-A02': [], 'AG-429-C19': [], 'AG-429-E20': [], 'AG-429-P21': [], 'AG-432-G10': [], 'AG-432-K16': [], 'AG-449-D16': [], 'AG-449-D22': ['C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1'], 'AG-449-K21': [], 'AG-455-E04': ['C4'], 'AG-459-P20': [], 'AG-469-M13': [], 'AG-670-J21': [], 'AG-347-B23': ['C1'], 'AG-347-C10': [], 'AG-347-E03': [], 'AG-347-G20': ['C3'], 'AG-347-G22': ['C3'], 'AG-347-I15': ['C3'], 'AG-347-I22': [], 'AG-347-I23': ['C5'], 'AG-347-J06': ['C1', 'C1'], 'AG-347-J23': ['C2'], 'AG-347-K02': [], 'AG-347-K10': [], 'AG-347-K15': [], 'AG-347-K22': [], 'AG-347-L13': ['C1', 'C1'], 'AG-347-L19': ['C1', 'C3'], 'AG-347-L21': [], 'AG-347-M15': ['C1'], 'AG-347-M18': [], 'AG-347-M23': ['C1', 'C1', 'C4'], 'AG-347-N19': ['C3'], 'AG-347-N23': ['C1'], 'AG-347-O22': ['C1', 'C1'], 'AG-355-A02': [], 'AG-355-A18': [], 'AG-355-B23': ['C1', 'C3'], 'AG-355-G23': [], 'AG-355-I20': [], 'AG-355-J04': ['C3'], 'AG-355-J23': ['C5', 'C5', 'C5', 'C1', 'C1', 'C1'],  'AG-355-K15': ['C1', 'C1', 'C1', 'C3'], 'AG-355-L22': ['C3'], 'AG-355-M18': [], 'AG-355-N18': ['C1', 'C1'], 'AG-355-P15': [], 'AG-355-P18': ['C1'], 'AG-388-F11': [], 'AG-402-A04': ['C1'], 'AG-402-C22': [], 'AG-402-F05': ['C1'], 'AG-402-G23': [], 'AG-402-K16': ['C5', 'C1', 'C1', 'C1'], 'AG-402-K22': [], 'AG-402-N17': [], 'AG-402-N23': [], 'AG-412-A14': [], 'AG-412-C21': [], 'AG-418-B17': [], 'AG-418-C17': ['C1'], 'AG-418-F08': ['C1'], 'AG-418-F16': [], 'AG-418-G18': ['C1', 'C1', 'C1', 'C1'], 'AG-418-G23': [], 'AG-418-I21': [], 'AG-418-J17': ['C1', 'C1', 'C1'], 'AG-418-J19': ['C1', 'C3'], 'AG-418-K17': [], 'AG-418-O02': ['C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1'], 'AG-418-O03': ['C3'], 'AG-418-P06': [], 'AG-418-P13': ['C1', 'C1', 'C1'], 'AG-424-A03': ['C4', 'C4'], 'AG-424-A14': [], 'AG-424-G03': [], 'AG-424-J22': ['C1', 'C1', 'C1', 'C3', 'C3', 'C4'], 'AG-424-M03': [], 'AG-424-P16': ['C1', 'C1', 'C1'], 'AG-424-P18': [], 'AG-424-P23': [], 'AG-436-A04': [], 'AG-436-D21': [], 'AG-436-E22': [], 'AG-449-C14': ['C1'], 'AG-449-G23': ['C1', 'C1', 'C1', 'C1'],'AG-449-J16': [], 'AG-449-O05': [], 'AG-455-E15': ['C3'], 'AG-459-A02': ['C3'], 'AG-459-B06': [], 'AG-459-D04': [], 'AG-459-E08': [], 'AG-459-J14': [], 'AG-459-N19': [], 'AG-459-O09': [], 'AG-670-M15': [], 'AG-670-M18': ['C1', 'C1', 'C1', 'C3']}

    ribocells=['AG-355-K20','AG-355-P11','AG-355-K13','AG-355-K23','AG-355-L02','AG-355-L20','AG-355-N23']
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

##############
tfna='a'

if tfna=='a':
    #outfilewhole='Matrices_Berube_Kash_All_NT_Discrete_ACTUALLY_Including_Non_Core'#this actually has noncore
    outfilewhole='Matrices_Berube_Kash_All_NT_Discrete_Including_Non_Core'#this doesn't have noncore
outfilewhole='Cell_Pair_Matrices/Each_Cell_Pair_Sites/%s' %outfilewhole

outfilewhole='../input_data/'+outfilewhole

infilematrices=outfilewhole+'.npy'


if 1:
    Official_Plot()
if 0:
    myfile='../input_data/HLII_clonelist2.txt'
    with open(myfile,'r') as infile:
        berube_celllist=[line.rstrip() for line in infile if len(line)>0]
    
    cluster_cells=[C1_without_ribotype,C2_without,C3_without,C4_without,C5_without]
    #cluster_cells=[C1,C2,C3,C4,C5]
        
    group_cells=cluster_cells#[C1, C2, C3, C4, C5]
    groupname_list=['C1','C2','C3','C4','C5']
    
    check_for_Berube_SAGs_Clustered_with_ribos(infilematrices,tfna,cluster_cells,group_cells,berube_celllist,groupname_list)



if 0:

    myfile='../input_data/HLII_clonelist2.txt'
    with open(myfile,'r') as infile:
        celllist_berube=[line.rstrip() for line in infile if len(line)>0]
    if 0:
        cluster_cells=['AG-355-N23', 'AG-355-P11', 'AG-418-O02', 'AG-449-D22', 'AG-355-O19', 'AG-355-K20', 'AG-355-L20']
        cluster_groupname='NonBB1_mostlyAtl'

    if 0:
        cluster_cells=['AG-355-A02', 'AG-424-A14', 'AG-432-K16', 'AG-418-J19', 'AG-424-G03', 'AG-355-L21', 'AG-355-B18']
        cluster_groupname='NonBB2_allAtl'

    if 1:
        cluster_cells=['AG-355-N22', 'AG-347-N19', 'AG-347-M08', 'AG-355-M02', 'AG-355-L22', 'AG-424-J22', 'AG-429-E20']
        cluster_groupname='BB1_mixed'

    cluster_groupname+='_NoClose'
    
    check_for_Berube_SAGs_Clustered_with_Berube_clusters(infilematrices,tfna,cluster_cells,celllist_berube,cluster_groupname,HLII_NoClose)#celllist_berube,cluster_groupname)


if 0:
    myfile='../input_data/HLII_clonelist2.txt'
    with open(myfile,'r') as infile:
        celllist_berube=[line.rstrip() for line in infile if len(line)>0]
    totcells=HLII_NoClose#celllist_berube
    
    if 1:
        cluster_groupname='NonBB2_allAtl'
        diamcut=.02
        cluster_cells=['AG-355-A02', 'AG-424-A14', 'AG-432-K16', 'AG-418-J19', 'AG-424-G03', 'AG-355-L21', 'AG-355-B18']
    if 0:
        cluster_groupname='NonBB1_mostlyAtl'
        diamcut=.02#arbitrayr
        cluster_cells=['AG-355-N23', 'AG-355-P11', 'AG-418-O02', 'AG-449-D22', 'AG-355-O19', 'AG-355-K20', 'AG-355-L20']
    if 1:
        cluster_cells=['AG-355-N22', 'AG-347-N19', 'AG-347-M08', 'AG-355-M02', 'AG-355-L22', 'AG-424-J22', 'AG-429-E20']
        cluster_groupname='BB1_mixed'
        diamcut=.02#75

    corebool=1
    analyze_other_clusters(cluster_groupname+'_NoClose',diamcut,totcells,cluster_cells,corebool)
