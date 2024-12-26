import numpy as np
import math
import ast
import csv
import statistics
from os.path import exists
import statistics
import argparse
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,AllKashBATS_without_ribotype, BB_closegrp_1,BB_closegrp_2,BB_closegrp_3,BB_closegrp_4,BB_closegrp_5,BB_closegrp_6,BB_closegrp_7,BB_closegrp_8,cl2,cl3,cl4,cl5,cl6
import matplotlib.pyplot as plt
import Repo_site_catalogues_single_pair as repo
import random
import Make_Synonymous_Distance_Matrices_Berube_Kashtan_Using_All_Sites_Per_Pair_and_HalfGene_Record_sites_SNPs as distmatrepo
import Get_R_From_Cell_Pair_Matrices as rrepo
import Partial_Sweeps as sweeprepo
import Make_Synonymous_Distance_Matrices_REFERENCES_Berube_Kashtan as referencerepo
import networkx as nx
#import Kashtan_Alignments_Pairwise_Event_Lengths as kashrepo
#import Comparing_Close_Cliques_to_Poisson as poissonrepo

def wg_v_lambda(celllist,groupname):#0.01 pcut
    if groupname=='C1' or '7-Cell C1 Clique':
        celllist=[item+'_C1' for item in celllist]
    pcut=.01
    chunklength=1000
    lambdas=[]
    lambdadict={}
    wgdict={}
    wgs=[]
    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
    for i in range(len(celllist)):
        print('%i of %i' %(i,len(celllist)))
        for j in range(i):
            ci=celllist[i]
            cj=celllist[j]
            mypair=sorted([celllist[i],celllist[j]])
            mypair=(mypair[0][:-3],mypair[1][:-3])
            mylambda,mywg=poissonrepo.prob_comparer_get_poisson_single_pair_for_flexgenes(i,j,celllist,chunk_seqdict,chunklength,pcut,returnwg=True)
            if mylambda=='-':
                continue
            
            lambdas.append(mylambda/float(chunklength))
            lambdadict[mypair]=lambdas[-1]
            wgs.append(mywg)
    plt.plot(lambdas,wgs,'o',alpha=.3)
    plt.xlabel('Asexual Divergence')
    plt.ylabel('Mean whole genome core divergence')
    plt.title('%s Mean WG Div vs Low Div' %(groupname))
    plt.show()


    
def get_wg_div_continuousfile_multiplegroups(infilematrices,tfna,celllist_list,groupname_list,outfile,refbool):
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]   
    if tfna in ['t','f','a','aa','n']:
        genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=distmatrepo.load_two_arrays(infilematrices)
    tfna_string=repo.return_tfna_string(tfna)
    discretebool=0
    if tfna=='a' and refbool:
        discretebool=1
    if tfna=='n':
        discretebool=1
    if tfna in ['t','f','tf','aa']:
        discretebool=1
    if not refbool:
        overalllist=distmatrepo.overall_order_list()
    if refbool:
        overalllist=referencerepo.overall_order_list(True)
    #make dict of dicts
    meds,means={},{}
    for item in groupname_list:
        meds[item]=[]
        means[item]=[]
    #lazy way out: just do the program again for each gene list.
    for i in range(len(groupname_list)):
        genedivsdict={}
        totaldivdict={}
        gn=groupname_list[i]
        initiate_dicts(genedivsdict,totaldivdict,celllist_list[i])
        for item in list(genearrays.keys()):
            mygene=repo.get_gene_from_key_to_check(item)
            if mygene not in mycores:
                continue
            print(mygene)
            if discretebool:
                get_divs_single_gene(genearrays[item],overalllist,celllist_list[i],totaldivdict,genedivsdict)
            if not discretebool:#tfna=='a' and not refbool:#continuous
                get_divs_single_gene_continuous(genearrays[item],overalllist,celllist_list[i],genedivsdict)
            
        for item in genedivsdict:
            if len(genedivsdict[item])==0:
                continue
            meds[gn].append(statistics.median(genedivsdict[item]))
            means[gn].append(repo.avelist(genedivsdict[item]))
        print('Done %s' %(groupname_list[i]))
    #now plot
    maxlist=[max(means[groupname_list[i]]) for i in range(len(groupname_list))]
                                
    mymax=max(maxlist)
    binwidth=0.001
    mybins=[binwidth*i for i in range(math.ceil(mymax/binwidth))]


    maxlistmed=[max(meds[groupname_list[i]]) for i in range(len(groupname_list))]
                                
    mymaxmed=max(maxlistmed)
    binwidth=0.002
    mybinsmed=[binwidth*i for i in range(math.ceil(mymaxmed/binwidth))]

    '''
    plt.hist([meds[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=mybinsmed,label=[groupname_list[i]+' (Mean %f)' %(repo.avelist(meds[groupname_list[i]])) for i in range(len(groupname_list))])
    plt.legend()
    plt.xlabel('Median Gene Divergence')
    plt.ylabel('Number of Cell Pairs')
    plt.title('%s Median Gene Divergences' %(tfna_string))
    plt.show()


    plt.hist([meds[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=mybinsmed,label=[groupname_list[i]+' (Mean %f)' %(repo.avelist(meds[groupname_list[i]])) for i in range(len(groupname_list))],density=True)
    plt.legend()
    plt.xlabel('Median Gene Divergence')
    plt.ylabel('Fraction of Cell Pairs per group')
    plt.title('%s Median Gene Divergences' %(tfna_string))
    plt.show()

    '''
    
    plt.hist([means[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=mybins,label=[groupname_list[i]+' (Mean %f)' %(repo.avelist(means[groupname_list[i]])) for i in range(len(groupname_list))],linewidth=3.5)
    plt.legend(fontsize=15)
    plt.xlabel('Mean Gene Divergence',fontsize=20)
    plt.ylabel('Fraction of Cell Pairs per group',fontsize=20)
    if tfna!='a':
        plt.title('%s Mean Gene Divergences' %(tfna_string),fontsize=25)
    if tfna=='a':
        plt.title('Mean Gene Nucleotide Divergences' ,fontsize=25)
    ax1=plt.gca()
    ax1.tick_params(axis='both', which='major', labelsize=25)
    plt.tight_layout()
    plt.show()


    plt.hist([means[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=mybins,label=[groupname_list[i]+' (Mean %f)' %(repo.avelist(means[groupname_list[i]])) for i in range(len(groupname_list))],density=True)
    plt.legend(fontsize=15)
    plt.xlabel('Mean Gene Divergence',fontsize=20)
    plt.ylabel('Fraction of Cell Pairs per group',fontsize=20)
    plt.title('%s Mean Gene Divergences' %(tfna_string),fontsize=25)
    ax1=plt.gca()
    ax1.tick_params(axis='both', which='major', labelsize=25)
    plt.show()

    
    if len(celllist_list)==1:
        write_wg_divs(outfile,genedivsdict,tfna,refbool)
            
def get_divs_single_gene_continuous(mat,overalllist,celllist,genedivsdict):
      for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            #d=mat[idxi][idxj]
            #print(d)
            #if d>1:
                #continue
            [d,l]=mat[idxi][idxj]
            if l<100:
                continue
            d=float(d)/l
            mypair=sorted([celllist[i],celllist[j]])
            pairkey=(mypair[0],mypair[1])
            genedivsdict[pairkey].append(d)
           
#######################################################


def get_wg_div_continuousfile_BETWEENgroups(infilewg,tfna,celllist_list,groupnames):
    if len(celllist_list)!=2:
        print('For this program, need exactly 2 grps of cells, put %i' %(len(celllist_list)))
        return
    
    wgarray=distmatrepo.load_arrays(infilewg)
    tfna_string=repo.return_tfna_string(tfna)
    overalllist=distmatrepo.overall_order_list()
   
    #make dict of dicts
    meds,means={},{}
    meds['00']=[]
    meds['11']=[]
    meds['01']=[]
    means['00']=[]
    means['11']=[]
    means['01']=[]

    l0=celllist_list[0]
    l1=celllist_list[1]

    matmean=wgarray['%sWG_%s' %('Mean',tfna)]
    matmed=wgarray['%sWG_%s' %('Median',tfna)]
    for i in range(len(l0)):
        idxi=overalllist.index(l0[i])
        for j in range(i):
            idxj=overalllist.index(l0[j])
            [distmean,nmean]=matmean[idxi][idxj]
            if nmean>10:
                means['00'].append(distmean)
            [distmed,nmed]=matmed[idxi][idxj]
            if nmed>10:
                meds['00'].append(distmed)
        for k in range(len(l1)):
            idxk=overalllist.index(l1[k])
            [distmean,nmean]=matmean[idxi][idxk]
            if nmean>10:
                means['01'].append(distmean)
            [distmed,nmed]=matmed[idxi][idxk]
            if nmed>10:
                meds['01'].append(distmed)
    for i in range(len(l1)):
        idxi=overalllist.index(l1[i])
        for j in range(i):
            idxj=overalllist.index(l1[j])
            [distmean,nmean]=matmean[idxi][idxj]
            if nmean>10:
                means['11'].append(distmean)
            [distmed,nmed]=matmed[idxi][idxj]
            if nmed>10:
                meds['11'].append(distmed)
    
   
    #now plot
    
                                
    mymax=max([max(means['00']),max(means['01']),max(means['11'])])
    binwidth=0.001
    mybins=[binwidth*i for i in range(math.ceil(mymax/binwidth))]


                                
    mymaxmed=max([max(meds['00']),max(meds['01']),max(meds['11'])])
    binwidth=0.002
    mybinsmed=[binwidth*i for i in range(math.ceil(mymaxmed/binwidth))]

    groupname_list=['%s-%s' %(groupnames[0],groupnames[0]),'%s-%s' %(groupnames[1],groupnames[1]),'%s-%s' %(groupnames[0],groupnames[1])]


    meanlabs=['%s (Mean %.3f)' %(groupname_list[0],repo.avelist(means['00'])),'%s (Mean %.3f)' %(groupname_list[1],repo.avelist(means['11'])),'%s (Mean %.3f)' %(groupname_list[2],repo.avelist(means['01'])) ]
    medlabs=['%s (Mean %.3f)' %(groupname_list[0],repo.avelist(meds['00'])),'%s (Mean %.3f)' %(groupname_list[1],repo.avelist(meds['11'])),'%s (Mean %.3f)' %(groupname_list[2],repo.avelist(meds['01'])) ]
    
    plt.hist([means[a] for a in ['00','11','01']],histtype='step',bins=mybins,label=meanlabs,linewidth=3.5)#[groupname_list[i]+' (Mean %f)' %(repo.avelist(means[a])) for a in ['00','11','01']],linewidth=3.5)              
    plt.legend(fontsize=15,loc='upper left')
    plt.xlabel('Mean Gene Divergence',fontsize=20)
    plt.ylabel('Number of Cell Pairs',fontsize=20)
    plt.title('Mean Gene Divergences' ,fontsize=25)
    ax1=plt.gca()
    ax1.tick_params(axis='both', which='major', labelsize=25)
    plt.tight_layout()
    plt.show()


    plt.hist([meds[a] for a in ['00','11','01']],histtype='step',bins=mybins,label=medlabs,linewidth=3.5)#[groupname_list[i]+' (Mean %f)' %(repo.avelist(meds[a])) for a in ['00','11','01']])
    plt.legend(fontsize=15,loc='upper left')
    plt.xlabel('Median Gene Divergence',fontsize=20)
    plt.ylabel('Number of Cell Pairs',fontsize=20)
    plt.title('Median Gene Divergences' ,fontsize=25)
    ax1=plt.gca()
    ax1.tick_params(axis='both', which='major', labelsize=25)
    plt.tight_layout()
    plt.show()
    
             
              


#########################################################
def get_wg_div_discretefile_singlegroup(infilematrices,celllist,tfna,groupname,outfile,corebool=1):
    print(infilematrices)
    #plot mean of all genes, median of all genes, and mean snps/length
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]     
    if tfna!='tf':
        genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=distmatrepo.load_two_arrays(infilematrices)
    overalllist=distmatrepo.overall_order_list()
    genedivsdict={}#list of all gene divs for a cell pair
    totaldivdict={}#for each cell pair: pair:[nsnpstotal,nmctotal]
    initiate_dicts(genedivsdict,totaldivdict,celllist)
    for item in list(genearrays.keys()):
        mygene=repo.get_gene_from_key_to_check(item)
        if corebool:
            if mygene not in mycores:
                continue
        print(mygene)
        get_divs_single_gene(genearrays[item],overalllist,celllist,totaldivdict,genedivsdict)
    #plot_divs(groupname,tfna,totaldivdict,genedivsdict)
    write_wg_divs(outfile,genedivsdict,tfna,False)



    
def write_wg_divs(outfile,genedivsdict,tfna,refbool):
   
    if not refbool:
        overalllist=distmatrepo.overall_order_list()
    if refbool:
        overalllist=referencerepo.overall_order_list(True)
    mydim=len(overalllist)
    medmat=[[[0,0] for i in range(mydim)] for j in range(mydim)]#div, number of genes covered
    meanmat=[[[0,0] for i in range(mydim)] for j in range(mydim)]#div, number of genes covered
    for item in genedivsdict:
        c1=item[0]
        c2=item[1]
        idx1=overalllist.index(c1)
        idx2=overalllist.index(c2)
        l=genedivsdict[item]
        if len(l)<2:
            continue
        medmat[idx1][idx2]=[statistics.median(genedivsdict[item]),len(genedivsdict[item])]
        medmat[idx2][idx1]=[statistics.median(genedivsdict[item]),len(genedivsdict[item])]
        #print([repo.avelist(genedivsdict[item]),len(genedivsdict[item])])
        meanmat[idx1][idx2]=[repo.avelist(genedivsdict[item]),len(genedivsdict[item])]
        meanmat[idx2][idx1]=[repo.avelist(genedivsdict[item]),len(genedivsdict[item])]
    mydict={'MeanWG_%s' %tfna:meanmat,'MedianWG_%s' %tfna :medmat}
    np.save(outfile,mydict)
        
def plot_divs(groupname,tfna,totaldivdict,genedivsdict):
    tfna_string=repo.return_tfna_string(tfna)
    meds,genemeans,totmeans=[],[],[]
    for item in totaldivdict:
        lt=totaldivdict[item]
        #print(item,lt)
        if lt[1]==0:
            print(item, "ZERO COVERAGE")
            continue
        totmeans.append(float(lt[0])/lt[1])
        meds.append(statistics.median(genedivsdict[item]))
        genemeans.append(repo.avelist(genedivsdict[item]))
    mymax=max([max(meds),max(genemeans),max(totmeans)])
    binwidth=0.001
    mybins=[binwidth*i for i in range(math.ceil(mymax/binwidth))]
    plt.hist(meds,bins=mybins)
    plt.title('%s %s Divergences:\nMedian Gene Divergence per cell pair\nMean %f Var %f' %(groupname,tfna_string,repo.avelist(meds),statistics.variance(meds)))
    plt.xlabel('Median Gene Divergence')
    plt.ylabel('Number of Cell Pairs')
    plt.tight_layout()
    plt.show()

    plt.hist(genemeans,bins=mybins)
    plt.title('%s %s Divergences:\nMean Gene Divergence per cell pair\nMean %f Variance %f' %(groupname,tfna_string,repo.avelist(genemeans),statistics.variance(genemeans)))
    plt.xlabel('Mean Gene Divergence')
    plt.ylabel('Number of Cell Pairs')
    plt.tight_layout()
    plt.show()

    plt.hist(totmeans,bins=mybins)
    plt.title('%s %s Divergences:\nMean Core Genome Divergence per cell pair (Not Weighted by Gene)\nMean %f Variance %f' %(groupname,tfna_string,repo.avelist(totmeans),statistics.variance(totmeans)))
    plt.xlabel('Mean Core Genome Divergence')
    plt.ylabel('Number of Cell Pairs')
    plt.tight_layout()
    plt.show()

    
def initiate_dicts(genedivsdict,totaldivdict,celllist):
     for i in range(len(celllist)):
        ci=celllist[i]
        for j in range(i):
            cj=celllist[j]
            mypair=sorted([ci,cj])
            pairkey=(mypair[0],mypair[1])
            genedivsdict[pairkey]=[]
            totaldivdict[pairkey]=[0,0]


def get_divs_single_gene(mat,overalllist,celllist,totaldivdict,genedivsdict):
      for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            [n,l]=mat[idxi][idxj]
            if l<=10:
                continue
            mypair=sorted([celllist[i],celllist[j]])
            pairkey=(mypair[0],mypair[1])
            genedivsdict[pairkey].append(float(n)/l)
            totaldivdict[pairkey][0]+=n
            totaldivdict[pairkey][1]+=l

#########################################################################
#coverage
#for each cell group, a hist.  the hist is: one entry per gene, the fraction of the cell group covered. do also number



def get_coverage_multiplegroups(infilematrices,tfna,celllist_list,groupname_list):
    #genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna in ['t','f','a']:
        genearrays=distmatrepo.load_arrays(infilematrices)
    if tfna=='tf':
        genearrays=distmatrepo.load_two_arrays(infilematrices)
    tfna_string=repo.return_tfna_string(tfna)
    overalllist=distmatrepo.overall_order_list()
    #make dict of dicts
    ncellscov={}
    fraccov={}
    for item in groupname_list:
        ncellscov[item]=[]
        fraccov[item]=[]
    #lazy way out: just do the program again for each gene list.
    for i in range(len(groupname_list)):
        genedivsdict={}
        totaldivdict={}
        gn=groupname_list[i]
        initiate_dicts(genedivsdict,totaldivdict,celllist_list[i])
        for item in list(genearrays.keys()):
            mygene=repo.get_gene_from_key_to_check(item)
            get_coverage_single_gene(genearrays[item],overalllist,celllist_list[i],ncellscov,gn,fraccov)
        print('Done %s' %(groupname_list[i]))
    #now plot
                                
    #binwidth=0.002
    #mybins=[binwidth*i for i in range(math.ceil(1.0/binwidth))]
    
    plt.hist([fraccov[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=50,label=['%s (Mean %f)' %(groupname_list[i],repo.avelist(fraccov[groupname_list[i]])) for i in range(len(groupname_list))])
    plt.legend()
    plt.xlabel('Fraction of Cells Covered')
    plt.ylabel('Number of Genes')
    plt.title('Fraction of Cells Covered per HLII Core Gene')
    plt.show()


def get_coverage_single_gene(mat,overalllist,celllist,ncellscov,groupname,fraccov):
    ncov=0
    for i in range(len(celllist)):
        idx=overalllist.index(celllist[i])
        mycol=mat[idx]
        for j in range(len(mycol)):
            [n,l]=mycol[j]
            if l>10:
                ncov+=1
                break
    ncellscov[groupname].append(ncov)
    fraccov[groupname].append(float(ncov)/len(celllist))

#########################################################
#WG all sites

def all_sites_wg_fromlist(infilecsv,celllist_list,groupname_list):
    meandict,meddict={},{}
    for item in groupname_list:
        meandict[item]=[]
        meddict[item]=[]
    with open(infilecsv,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            c1=row['Cell1']
            c2=row['Cell2']
            med=float(row['Median'])
            mean=float(row['Mean'])
            for i in range(len(celllist_list)):
                cl=celllist_list[i]
                if c1 in cl and c2 in cl:
                    meandict[groupname_list[i]].append(mean)
                    meddict[groupname_list[i]].append(med)
    mybins=40
    plt.hist([meddict[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=mybins,label=[groupname_list[i] for i in range(len(groupname_list))])
    plt.legend()
    plt.xlabel('Median Gene Divergence')
    plt.ylabel('Number of Cell Pairs')
    plt.title('All Sites Median Gene Divergences')
    plt.show()


    plt.hist([meandict[groupname_list[i]] for i in range(len(groupname_list))],histtype='step',bins=mybins,label=[groupname_list[i] for i in range(len(groupname_list))])
    plt.legend()
    plt.xlabel('Mean Gene Divergence')
    plt.ylabel('Fraction of Cell Pairs per group')
    plt.title('All Sites Mean Gene Divergences')
    plt.show()


##########################################################
#where are narB cells from and do they cluster on divergence matrix?
##########################################################################
#get group of outlier cells, div >8% WG
def get_outlier_group(infilewg,orig_celllist,cutoff,meanmedian,tfna):
    wgarray=distmatrepo.load_arrays(infilewg)
    overalllist=distmatrepo.overall_order_list()
    edges=[]
    mat=wgarray['%sWG_%s' %(meanmedian,tfna)]
    for i in range(len(orig_celllist)):
        idxi=overalllist.index(orig_celllist[i])
        for j in range(i):
            idxj=overalllist.index(orig_celllist[j])
            [dist,n]=mat[idxi][idxj]
            # print(dist,n)
            if dist>=cutoff:
                edges.append((orig_celllist[i],orig_celllist[j]))
    ncells,maxclique=sweeprepo.get_cliques(edges)
    print(maxclique)


def get_close_groups(infilewg,orig_celllist,cutoff,meanmedian,tfna):
    wgarray=distmatrepo.load_arrays(infilewg)
    overalllist=distmatrepo.overall_order_list()
    edges=[]
    mat=wgarray['%sWG_%s' %(meanmedian,tfna)]
    for i in range(len(orig_celllist)):
        idxi=overalllist.index(orig_celllist[i])
        for j in range(i):
            idxj=overalllist.index(orig_celllist[j])
            [dist,n]=mat[idxi][idxj]
            # print(dist,n)
            if dist<=cutoff:
                edges.append((orig_celllist[i],orig_celllist[j]))
    get_graph_connected_comps(edges)
   

def get_graph_connected_comps(edges):#rcut is mean_median, minsep is divcut
    G=nx.Graph()
    G.add_edges_from(edges)


    orthogroups=[]
    for clq in nx.connected_components(G):#find_cliques(G):#connected_components(G):
        miniclqlist=[]
        for item in clq:
            miniclqlist.append(item)
        print(miniclqlist)
        miniclqlist=sorted(miniclqlist)
        orthogroups.append(miniclqlist)
        for item in miniclqlist:
            if item in BB_list:
                print('BB')
    groupsizes=[len(orthogroups[i]) for i in range(len(orthogroups))]
    print(sorted(groupsizes))
    if 1:
        mylist=list(range(max(groupsizes)+1))
        mybins=[-0.5]
        for i in range(len(mylist)):
            mybins.append(i+0.5)
        plt.hist(groupsizes,bins=mybins)
        
        plt.show()
    
    
###################
def compare_cell_wg_div(celllist,out_cells):
    totdivs=[]
    outdivs=[]
    for cell in celllist:
        divs=get_single_cell_divs(cell,celllist)
        for item in divs:
            totdivs.append(item)
        if cell in out_cells:
            for item in divs:
                outdivs.append(item)
    plt.hist([totdivs,outdivs],bins=100,label=['C1','Cells with singletons'])
    #plt.title('WG nt divergences for C1 vs cells with many singleton events\nC1 ave wg div %.4f Median %.4f; outlier cells mean %.4f median %.4f' %(repo.avelist(totdivs),statistics.median(totdivs),repo.avelist(outdivs),statistics.median(outdivs)))
    plt.title('WG nt divergences for C1 vs %s\nC1 ave wg div %.4f Median %.4f; outlier cells mean %.4f median %.4f' %(str(out_cells),repo.avelist(totdivs),statistics.median(totdivs),repo.avelist(outdivs),statistics.median(outdivs)))
    plt.show()
            
    
def get_single_cell_divs(cell,celllist):
    arrayfile= '../input_data/Cell_Pair_Matrices/Berube_Kash_WG_Matrices_a.npy'#Berube_Kash_6refs_REFRENCES_WG_Matrices_a.npy'
    genearrays=distmatrepo.load_arrays(arrayfile)
    overalllist=referencerepo.overall_order_list(True)
    divs=[]
    distmat=genearrays['MeanWG_a']
    idxi=overalllist.index(cell)
    for i in range(len(celllist)):
        idxj=overalllist.index(celllist[i])
        if idxi==idxj:
            continue
        [d,l]=distmat[idxi][idxj]
        divs.append(d)
    return divs


#################################

def get_high_div_cells(celllist,cut):
    arrayfile= '../input_data/Cell_Pair_Matrices/Berube_Kash_6refs_REFRENCES_WG_Matrices_a.npy'
    genearrays=distmatrepo.load_arrays(arrayfile)
    overalllist=referencerepo.overall_order_list(True)
    divs=[]
    distmat=genearrays['MeanWG_a']
    celldict={}#cell:number of high div pts
    for i in range(len(celllist)):
        idxi=overalllist.index(celllist[i])
        for j in range(i):
            idxj=overalllist.index(celllist[j])
            [d,l]=distmat[idxi][idxj]
            if d>cut:
                print('Cells %s %s %f' %(celllist[i],celllist[j],d))
                if i in celldict:
                    celldict[i]+=1
                if i not in celldict:
                    celldict[i]=1
                if j in celldict:
                    celldict[j]+=1
                if j not in celldict:
                    celldict[j]=1
                '''
                if celllist[i] in celldict:
                    celldict[celllist[i]]+=1
                if celllist[i] not in celldict:
                    celldict[celllist[i]]=1
                if celllist[j] in celldict:
                    celldict[celllist[j]]+=1
                if celllist[j] not in celldict:
                    celldict[celllist[j]]=1
                '''
    print('Cells with WG divs > %f' %cut)
    print(celldict)
########################

def print_wg_distances(list1,list2):
    #arrayfile= '../input_data/Cell_Pair_Matrices/Berube_Kash_6refs_REFRENCES_WG_Matrices_a.npy'
    arrayfile= '../input_data/Cell_Pair_Matrices/Berube_Kash_WG_Matrices_a.npy'
    genearrays=distmatrepo.load_arrays(arrayfile)
    overalllist=referencerepo.overall_order_list(True)
    mat=[[0 for i in range(1+len(list2))] for j in range(1+len(list1))]
    justdists=[[0 for i in range(len(list2))] for j in range(len(list1))]
    distmat=genearrays['MeanWG_a']
    for i in range(len(list1)):
        idxi=overalllist.index(list1[i])
        mat[i+1][0]=list1[i]
        for j in range(len(list2)):
            mat[0][j+1]=list2[j]
            idxj=overalllist.index(list2[j])
            [d,l]=distmat[idxi][idxj]
            mat[i+1][j+1]=d
            justdists[i][j]=d
    
    #print(list1)
    #print(list2)
    print(mat)
    arr=np.array(justdists)
    np.save('../input_data/HLII_WG_Matrices_allNT.npy',arr)

                                                                                

