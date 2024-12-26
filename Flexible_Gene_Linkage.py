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
import Flexible_Gene_Shared as flexrepo
import ast
import re


#Atlantic=['355','363','388','412','429','432','436']
#N_Pacific=['347','402']
#S_Pacific=['335','459','469','449','455']
#Pacific=N_Pacific+S_Pacific
Atlantic=['355','363','388','412','429','432','436']+['409','418','420','424','442','444']
Pacific=['347','402','335','459','469','449','455']+['311','315','316','321','323','331','341','345','450','463','670','673','676','679','683','686']



def analyze_samps(samps):#have a list of cells, check if they're from same ocean or sample
    #get plurality sample
    plu=''
    count=0
    for item in samps:
        if samps.count(item)>count:
            plu=item
            count=samps.count(item)
    if count>1:
        print('Plurality sample %s: %i of %i' %(plu,count,len(samps)))
    pac,atl=0,0
    for item in samps:
        if item in Atlantic:
            atl+=1
        elif item in Pacific:
            pac+=1
    print('%f from Atl, %f from Pac (of %i)' %(float(atl)/len(samps),float(pac)/len(samps),len(samps)))

def get_core_coverage(groupname,celllist):
    infile='../input_data/%s_Cells_Covered_Per_Orthogroup.txt' %groupname
    with open(infile,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    with open('../input_data/CoreGenes_MIT9301_AS9601_MIT9215_MIT9312_MIT0604.txt','r') as cf:
        mycores=[int(line.rstrip()) for line in cf if len(line)>0]  
    covs=[]
    for item in lines:
        idx=item.find(':')
        filename=item[:idx]
        num=int(item[idx+1:])
        idx2=filename.find('_')
        genenum=int(filename[4:idx2])
        if genenum in mycores:
            covs.append(num)

    print(sorted(covs))
    print('%s (%i cells): %i genes found.  Mean %f Cells Median %f Min %i Max %i' %(groupname,len(celllist),len(covs),float(sum(covs))/len(covs),statistics.median(covs),min(covs),max(covs)))
    
def linkage_matrix_flex_genes(celllist,groupname,sizecut,qcut,pidcut):
    #countdictpair=flexrepo.get_number_flex_genes_shared(groupname,celllist,qcut,pidcut)
    #countdictpair is cellpair:ngenesshared
    #need for each orf, list of cells with and without

    if groupname=='HLII':
        c=0.71856
    if groupname=='Blue_Basin_NoClose':
        c=0.71856#coverage. just guessed from above

    grouplist=load_ortho_groups(celllist,groupname,sizecut,qcut,pidcut)
    #grouplist=load_ortho_groups_from_pfam()
    print('%i groups' %len(grouplist))
    t=input('t')
    r2,filist,fjlist,fijlist=[],[],[],[]
    nilist,njlist=[],[]
    edges=[]
    #could: for each gene pair, get f12, f1, f2 observed, and use coverage to estimate the the real f12, f1 f2, then do r2
    for i in range(len(grouplist)):
      
        listi=grouplist[i]#listi is a list of cells, not orf ids
        fio=float(len(listi))/len(celllist)# o for observed
        for j in range(i):
            listj=grouplist[j]
            fjo=float(len(listj))/len(celllist)
            nilist.append(len(listi))
            njlist.append(len(listj))
            nij=len(list(set(listi).intersection(set(listj))))
            fijo=float(nij)/len(celllist)
            n_i_only=len(list(set(listi).difference(set(listj))))
            n_j_only=len(list(set(listj).difference(set(listi))))
            fij,fi,fj=convert_observed_to_theoretical_freqs(fijo,float(n_i_only)/len(celllist),float(n_j_only)/len(celllist),c)
            #myr2=calculate_r2(fij,fi,fj)
            myr2=calculate_r2(fijo,fio,fjo)
            if myr2>=1:
                print(i,j,'r2 %f ni only %i nj only %i nij %i' %(myr2, n_i_only,n_j_only,nij))
            r2.append(myr2)
            filist.append(fi)
            fjlist.append(fj)
            fijlist.append(fij)
            if myr2>0.3:
                edges.append((i,j))
    if 1:
        highly_linked_groups(edges,grouplist)
        
    if 1:#plot r2 distribution
        plt.hist(r2,bins=100)
        plt.xlabel('r^2 values using estimated theoretical frequencies')
        plt.ylabel('Number of flexible gene group pairs')
        plt.title('%s Flexible Gene Groups (Pid %i Qcut %i)\nLinkage distribution (4-%i cells per grp)' %(groupname,pidcut,qcut,sizecut))
        #plt.title('%s Flexible Gene Groups from Pfam\nLinkage distribution (4-%i cells per grp)' %(groupname,sizecut))
        plt.yscale('log')
        plt.show()

    if 1:
        make_2d_hist(nilist,njlist,r2,sizecut)

def get_gff_annotation(idlist):
    with open('../input_data/All_HLII_Gff.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    for myid in [idlist[0]]:
        for line in lines:
            if re.search(myid,line):
                #print(line)
                idx=line.index('product')
                print(line[idx+8:])
                break

def return_gff_annotation(idlist):
    with open('../input_data/All_HLII_Gff.txt','r') as myinfile:
        lines=[line.rstrip() for line in myinfile if len(line)>0]
    for myid in [idlist[0]]:
        for line in lines:
            if re.search(myid,line):
                #print(line)
                if 'product' not in line:
                    return 'no product'
                idx=line.index('product')
                return line[idx+8:]
               

    
def clustering_orthogroups_get_annotations(celllist,groupname,sizecut,qcut,pidcut):
    #countdictpair=flexrepo.get_number_flex_genes_shared(groupname,celllist,qcut,pidcut)
    #countdictpair is cellpair:ngenesshared
    #need for each orf, list of cells with and without

    if groupname=='HLII' or groupname=='Blue_Basin_NoClose':
        c=0.71856

    groups_cells,groups_ids=load_ortho_group_dict_from_clustering(celllist,groupname,sizecut,qcut,pidcut)#groups_cells and _ids are dicts with : [line index of flex gene file]:list of cells or ORF ids
    
    #grouplist=load_ortho_groups_from_pfam()
    print('%i groups' %len(list(groups_cells.keys())))
    
    r2,filist,fjlist,fijlist=[],[],[],[]
    nilist,njlist=[],[]
    edges=[]
    #could: for each gene pair, get f12, f1, f2 observed, and use coverage to estimate the the real f12, f1 f2, then do r2
    mymin,sizecut=4,34
    mykeys=list(groups_cells.keys())
    mat=[[np.nan for i in range(len(mykeys))] for x in range(len(mykeys))]
    for i in range(len(list(groups_cells.keys()))):
        myi=mykeys[i]
        listi=groups_cells[myi]
        fio=float(len(listi))/len(celllist)# o for observed
        for j in range(i):
            listj=groups_cells[mykeys[j]]
            fjo=float(len(listj))/len(celllist)
            nilist.append(len(listi))
            njlist.append(len(listj))
            nij=len(list(set(listi).intersection(set(listj))))
            fijo=float(nij)/len(celllist)
            n_i_only=len(list(set(listi).difference(set(listj))))
            n_j_only=len(list(set(listj).difference(set(listi))))
            fij,fi,fj=convert_observed_to_theoretical_freqs(fijo,float(n_i_only)/len(celllist),float(n_j_only)/len(celllist),c)
            #myr2=calculate_r2(fij,fi,fj)
            myr2=calculate_r2(fijo,fio,fjo)
            if myr2>=1:
                print(i,j,'r2 %f ni only %i nj only %i nij %i' %(myr2, n_i_only,n_j_only,nij))
                if 0:#nij>5:
                    print('Gene 1')
                    get_gff_annotation(groups_ids[myi])
                    print('\nGene 2')
                    get_gff_annotation(groups_ids[mykeys[j]])
                    t=input('t')
            r2.append(myr2)
            filist.append(fi)
            fjlist.append(fj)
            fijlist.append(fij)
            mat[myi][mykeys[j]]=myr2
            mat[mykeys[j]][myi]=myr2
            if myr2>0.3:
                edges.append((i,j))
                #check_if_linkage_is_presence_absence(listi,listj)
    if 0:
        #highly_linked_groups_from_clustering(edges,groups_ids)
        grouplist=load_ortho_groups(celllist,groupname,sizecut,qcut,pidcut)
        highly_linked_groups(edges,grouplist)
    if 1:
        xticklabs=list(range(mymin,sizecut+1))
        xtickpos=list(range(len(mat)))
        plt.imshow(mat)
        plt.colorbar(spacing='uniform')
        plt.xlabel('Flexible Gene Index (Arbitrary)')
        ax=plt.gca()
        ax.set_xticks(xtickpos)
        ax.set_yticks(xtickpos)
        ax.set_xticklabels(xticklabs)
        ax.set_yticklabels(xticklabs)
        plt.title('%s r^2 of Flexible Gene Presence/Absence\nfor each Gene Pair\nUsing Connected Component Orthogroups' %groupname)
        plt.show()

        refine_matrix(mat,0.2,groupname)
        refine_matrix(mat,0.3,groupname)

def check_if_linkage_is_presence_absence(listi,listj):
    ni=len(listi)
    nj=len(listj)
    nij=len(list(set(listi).intersection(set(listj))))
    if nij<0.7*min([ni,nj]):
        print('absence linkage: ni %i nj %i nij %i' %(ni,nj,nij))
        t=input('t')
        
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
                mynewid=cell[:idx]
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


def load_ortho_groups_from_pfam():
    groups=[]
    with open('../input_data/HLII_Flexible_Gene_Groups.csv','r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            l=ast.literal_eval(row['Cells'])
            if len(l)<4:
                continue
            newl=[]
            for item in l:
                if 'Screened' in item:
                    idx=item.find(' ')
                    item=item[:idx]
                newl.append(item)
            groups.append(newl)
    return groups

def make_2d_hist(nilist,njlist,r2,sizecut):
    mymin=min(nilist)
    mat=[[[] for i  in range(mymin,sizecut+1)] for j in range(mymin,sizecut+1)]
    for i in range(len(nilist)):
        myni=nilist[i]-mymin
        mynj=njlist[i]-mymin
        mat[myni][mynj].append(r2[i])
        mat[mynj][myni].append(r2[i])
    #now average
    avemat=[[0 for i in range(mymin,sizecut+1)] for j in range(mymin,sizecut+1)]
    nummat=[[0 for i in range(mymin,sizecut+1)] for j in range(mymin,sizecut+1)]
    for i in range(len(mat)):
        for j in range(len(mat)):
            mylist=mat[i][j]
            if len(mylist)==0:
                #print(i,j,'zero entries')
                #t=input('t')
                continue
            nummat[i][j]=len(mylist)
            #nummat[j][i]=len(mylist)
            avemat[i][j]=float(sum(mylist))/len(mylist)

    xticklabs=list(range(mymin,sizecut+1))
    xtickpos=list(range(len(mat)))
    plt.imshow(avemat)
    plt.colorbar(spacing='uniform')
    plt.xlabel('Number of Cells with Flexible Gene')
    ax=plt.gca()
    ax.set_xticks(xtickpos)
    ax.set_yticks(xtickpos)
    ax.set_xticklabels(xticklabs)
    ax.set_yticklabels(xticklabs)
    plt.title('<r^2> of Flexible Gene Presence/Absence\nBinned by frequency')
    plt.show()

    plt.imshow(nummat,norm=colors.LogNorm())
    plt.colorbar()#spacing='uniform')
    plt.xlabel('Number of Cells with Flexible Gene')
    ax=plt.gca()
    ax.set_xticks(xtickpos)
    ax.set_yticks(xtickpos)
    ax.set_xticklabels(xticklabs)
    ax.set_yticklabels(xticklabs)
    plt.title('<r^2> of Flexible Gene Presence/Absence\nBinned by frequency')
    plt.show()
    
def convert_observed_to_theoretical_freqs(fboth,fionly,fjonly,c):#c is coverage?
    #p_ij=f_ij c^2
    #p_i_only=fi*c-f_ij *C^2
    #p_neither=1-f_i-f_j+f_ij +f_i *(1-c)^2 + (1-c)*(fi-f_ij+fj-fij)
    fij=float(fboth)/(c**2)
    fi=(fionly+fij*c*c)/float(c)
    fj=(fjonly+fij*c*c)/float(c)
    return fij,fi,fj

def calculate_r2(fij,fi,fj):
    r2=((float(fi*fj-fij)**2))/(fi*(1.0-fi)*(1.0-fj)*fj)
    #if r2>=1:
        #print('r2 %f fi %f fj %f fij %f' %(r2,fi,fj,fij))
        #t=input('t')
    return r2

            

def load_ortho_groups(celllist,groupname,sizecut,qcut,pidcut):#returns cells. want ORF ids
    f='../input_data/%s_Orthologous_Gene_Group_Qcut%s_pid%s.txt' %(groupname,str(qcut),str(pidcut))
    with open(f,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    groups=[]

    lower_cut=4
    for line in lines:
        group=ast.literal_eval(line)
        if len(group)<lower_cut:
            continue    
        if len(group)>sizecut:
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
        groups.append(newgroup)
    return groups



###########################################################


def make_HLII_flex_orf_list():
    f1='../input_data/All_Berube_HLII_ORF_Names.txt'
    f2='../input_data/All_HLII_Pfam.txt'
    f3='../input_data/All_HLII_Gff.txt'
    #make a network of pfams. get networks size 34 or fewer, nonsingleon
    #then get the cells
    #then write csv with: pfram number, annotation, cells

    sizecut=34
    pfamdict={}#pfam:[list of orfs]
    pfamname={}#pfam number: name
    with open(f2,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            pfam=row['pfam_id']
            if pfam=='pfam_id':
                continue
            if pfam not in pfamname:
                name=row['pfam_name']
                pfamname[pfam]=name
                pfamdict[pfam]=[]
            orf=int(row['gene_oid'])
            pfamdict[pfam].append(orf)
    #now that's done
    #prune to size

    with open(f1,'r') as myinfile:
        lines=myinfile.readlines()
        
    pfamcells={}
    count=0
    for item in pfamdict:
        #if len(list(pfamcells.keys()))==2:
            #break
        count+=1
        if count%10==0:
            print('count %i of %i' %(count,len(list(pfamdict.keys()))))
        l=pfamdict[item]
        if len(l)==1:
            continue
        if len(l)>sizecut:
            continue
        pfamcells[item]=[]
        for orf in l:
            cell=search_for_cell(lines,orf)
            pfamcells[item].append(cell)


    with open(f3,'r') as myinfile:
        lines=myinfile.readlines()
        
    pfamgff={}
    count=0
    for item in pfamcells:
        
        count+=1
        if count%10==0:
            print('count %i of %i' %(count,len(list(pfamdict.keys()))))
        l=pfamdict[item]
        for orf in [l[0]]:
            prod=search_for_product(lines,orf)
            pfamgff[item]=prod
    outfile='../input_data/HLII_Flexible_Gene_Groups.csv'
    outfields=['Pfam_ID','Pfam_Name','Annotation','Cells']
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,fieldnames=outfields,delimiter='\t')
        outwriter.writeheader()
        for item in pfamcells:
            outwriter.writerow({'Pfam_ID':item,'Pfam_Name':pfamname[item],'Annotation':pfamgff[item],'Cells':str(pfamcells[item])})

def search_for_cell(lines,orf):
    for i, line in enumerate(lines):
        if re.search(str(orf),line):
            idx=line.find('AG-')
            cell=line[idx:-2]
            return cell
        
    
def search_for_product(lines,orf):
    for i, line in enumerate(lines):
        if re.search(str(orf),line):
            idx=line.find('product')
            prod=line[idx+8:]
            return prod
#######################

def get_linkage_from_pfam_file(celllist,groupname,sizecut):
    f1='../input_data/HLII_Flexible_Gene_Groups.csv'
    #celllist=hliilist2
    mymin=4
    #sizecut=34

    r2=[]
    edges=[]
    cutoff=0.2
    #get cell linkage
    #can bin by freq
    #can compare to other orthogroups
    #can annotations of linked things
    groupdict=dict_ortho_groups_from_pfam()
    groups=list(groupdict.keys())
    print('%i groups' %(len(groups)))
    t=input('t')
    mat=[[np.nan for i in range(len(groups))] for x in range(len(groups))]
    for i in range(len(groups)):
        groupi=groupdict[groups[i]]
        for j in range(i):
            groupj=groupdict[groups[j]]
            ni,nj=len(groupi),len(groupj)
            nij=len(list(set(groupi).intersection(set(groupj))))
            myr2=calculate_r2(float(nij)/len(celllist),float(ni)/len(celllist),float(nj)/len(celllist))
            mat[i][j]=myr2
            mat[j][i]=myr2
            r2.append(myr2)
            if myr2>.3:
                print('ni %i nj %i nij %i r2 %f %s %s' %(ni,nj,nij,myr2,groups[i],groups[j]))
            if myr2>.3:#cutoff:
                edges.append((groups[i],groups[j]))

    if 1:
        grouplist=load_ortho_groups_from_pfam()##new nov 2024
        highly_linked_groups(edges,groupdict)#list)
    if 1:
        plt.hist(r2,bins=40)
        plt.title('HLII flexible orthogroups presence/absence linkage\nUsing Pfram groups 4-34 cells')
        plt.xlabel('r**2')
        plt.ylabel('Number of ORF pairs')
        plt.yscale('log')
        plt.show()
        
        xticklabs=list(range(mymin,sizecut+1))
        xtickpos=list(range(len(mat)))
        plt.imshow(mat)
        plt.colorbar(spacing='uniform')
        plt.xlabel('Flexible Gene Index (Arbitrary)')
        ax=plt.gca()
        ax.set_xticks(xtickpos)
        ax.set_yticks(xtickpos)
        ax.set_xticklabels(xticklabs)
        ax.set_yticklabels(xticklabs)
        plt.title('r^2 of Flexible Gene Presence/Absence\nfor each Gene Pair')
        plt.show()

        refine_matrix(mat,0.3,groupname)
        
def refine_matrix(mat,cutoff,groupname):
    #keep only rows and columns with at least 1 entry>cutoff
    indices_to_keep=[]
    for i in range(len(mat)):
        col=mat[i]
        col=col[:i]+col[i+1:]
        if max(col)>cutoff:
            indices_to_keep.append(i)
    newmat=[[np.nan for i in range(len(indices_to_keep))] for j in range(len(indices_to_keep))]
    for i in range(len(indices_to_keep)):
        for j in range(len(indices_to_keep)):
            newmat[i][j]=mat[indices_to_keep[i]][indices_to_keep[j]]
    
   
    plt.imshow(newmat)
    plt.colorbar(spacing='uniform')
    plt.xlabel('Flexible Gene Index (Arbitrary)')
    plt.title('%s r^2 of Flexible Gene Presence/Absence\nfor Gene Pairs with r2>%s' %(groupname,cutoff))
    plt.show()


def highly_linked_groups(edges,grouplist):#grouplist is list of cells with the gene?
    G=nx.Graph()
    G.add_edges_from(edges)

    ngenes,avencells=[],[]
    orthogroups=[]

    mylinkedcells=[]
    for clq in nx.connected_components(G):
        miniclqlist=[]
        for item in clq:
            miniclqlist.append(item)
        miniclqlist=sorted(miniclqlist)
        orthogroups.append(miniclqlist)
        print(len(miniclqlist),miniclqlist)
        #get ave number of cells per gene
        if 1:
            ncells=[]
            #print(grouplist)
            mydict={}#cell:ngenes
            for item in miniclqlist:
                ncells.append(len(grouplist[item]))
                for e in grouplist[item]:
                    if e in mydict:
                        mydict[e]+=1
                    if e not in mydict:
                        mydict[e]=1
            print('average %f cells' %(float(sum(ncells))/len(ncells)))
            ngenes.append(len(miniclqlist))
            avencells.append(float(sum(ncells))/len(ncells))
            if 1:
                #if avencells[-1]>=10:#ngenes[-1]>=10:#avencells[-1]>=10:
                if 1:#ngenes[-1]>10:
                    print('%i genes' %ngenes[-1])
                    #print(mydict)
                    #t=input('t')
                    littlelist=[cell for cell in mydict if mydict[cell]>=3]
                    if len(littlelist)>1:
                        mylinkedcells.append(littlelist)
                        #check if littlelist is interesting. from same sample or ocean?
                        samps=[cc[3:6] for cc in littlelist]
                        analyze_samps(samps)
                        t=input('t')
    if 0:#analysize littlelist
        #get overlap of entries
        for i in range(len(mylinkedcells)):
            for j in range(i):
                li=mylinkedcells[i]
                lj=mylinkedcells[j]
                lij=list(set(li).intersection(set(lj)))
                print('%i cells, %i cells, %i mutual cells' %(len(li),len(lj),len(lij)))
                print(lij)
                t=input('t')

        
    mylens=[len(item) for item in orthogroups]
    print('Group lengths:')
    print(mylens)
    plt.plot(ngenes,avencells,'o')
    plt.xlabel('Ngenes linked in grp')
    plt.ylabel('Ave n cells covered per gene')
    plt.show()
    #print(orthogroups)

                


def highly_linked_groups_from_clustering(edges,groups_ids,grouplist):
    G=nx.Graph()
    G.add_edges_from(edges)
    pfamdict=load_pfam_file()
    orthogroups=[]
    for clq in nx.connected_components(G):
        miniclqlist=[]
        for item in clq:
            miniclqlist.append(item)
        miniclqlist=sorted(miniclqlist)
        orthogroups.append(miniclqlist)
        #print(miniclqlist)
        #now get the annotations of all members of the edges
        cluster_annotation(pfamdict,groups_ids,miniclqlist)

        
    mylens=[len(item) for item in orthogroups]
    print('Group lengths:')
    print(mylens)
    #print(orthogroups)


def cluster_annotation(pfamdict,groups_ids,miniclq):
    annotations=[]
    for item in miniclq:
        myids=groups_ids[item]
        i=0
        for i in range(len(myids)):
            if myids[i] in pfamdict:
                annotations.append(pfamdict[myids[i]])
                break
            if myids[i] not in pfamdict:
                if i==len(myids)-1:
                    #annotations.append('Unannotated')
                    annotations.append(return_gff_annotation(myids))
    print(len(annotations),annotations)
    
def load_pfam_file():
    pfamdict={}
    with open('../input_data/All_HLII_Pfam.txt','r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            pfamdict[row['gene_oid']]=row['pfam_name']
    return pfamdict

def dict_ortho_groups_from_pfam():
    groups={}
    with open('../input_data/HLII_Flexible_Gene_Groups.csv','r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            l=ast.literal_eval(row['Cells'])
            if len(l)<4:
                continue
            newl=[]
            for item in l:
                if 'Screened' in item:
                    idx=item.find(' ')
                    item=item[:idx]
                newl.append(item)
            if row['Annotation'] in groups:
                index=2
                while '%s_%i' %(row['Annotation'],index) in groups:
                    index+=1
                groups['%s_%i' %(row['Annotation'],index)]=newl
            if row['Annotation'] not in groups:
                groups[row['Annotation']]=newl
    return groups
    
#######################################################

def compare_flexible_groups_graph_pfam():#inconclusive
    sizecut=34
    pidcut,qcut,groupname,celllist=80,80,'HLII',hliilist2
    groupgraph=load_ortho_groups(celllist,groupname,sizecut,qcut,pidcut)
    grouppfam=load_ortho_groups_from_pfam()
    sizesgraph=[len(item) for item in groupgraph]
    sizespfam=[len(item) for item in grouppfam]
    mybins=[-.5]
    for i in range(35):
        mybins.append(i+.5)
    plt.hist([sizesgraph,sizespfam],label=['Graph','Pfam'],histtype='step',bins=mybins)
    plt.legend()
    plt.xlabel('Number of cells with flexible orf\n%i groups graph, %i groups pfam' %(len(sizesgraph),len(sizespfam)))
    plt.show()

    #grouppfam is list of cells per group
    #so is groupgraph
    #can we find matches?
    matches,nomatch=0,0
    for item in grouppfam:
        print('%i matches %i nomatch' %(matches,nomatch))
        matchbool=0
        for g in groupgraph:
            print(item)
            print(g)
            print(len(list(set(item).intersection(set(g)))),len(item),len(g))
            #t=input('t')
            if len(list(set(item).intersection(set(g))))>len(item)-2 and len(list(set(item).intersection(set(g))))>len(g)-2 :
                matches+=1
                matchbool=1
                t=input('t')
                break
            if matchbool:
                continue
        if not matchbool:
            nomatch+=1
    print('%i pfams, %i matches, %i no match' %(len(grouppfam),matches,nomatch))
###################################################


def r2_null(celllist,groupname,sizecut,qcut,pidcut):
    #put in coverage expectation
    #we have f_observed
    #perfect linkage: freqs are the same
    cov=0.71856

    #need: number of cells per orthogroup

    grouplist=load_ortho_groups(celllist,groupname,sizecut,qcut,pidcut)
    #grouplist=load_ortho_groups_from_pfam()
   
    r2dict={}#fi: list of r2 vals
   
    #could: for each gene pair, get f12, f1, f2 observed, and use coverage to estimate the the real f12, f1 f2, then do r2
    for i in range(len(grouplist)):
      
        listi=grouplist[i]#listi is a list of cells, not orf ids
        fio=float(len(listi))/len(celllist)# o for observed
        for j in range(i):
            listj=grouplist[j]
            fjo=float(len(listj))/len(celllist)
            nij=len(list(set(listi).intersection(set(listj))))
            fijo=float(nij)/len(celllist)
            myr2=calculate_r2(fijo,fio,fjo)
            if myr2=='-':
                continue
            ni=len(listi)
            if ni in r2dict:
                r2dict[ni].append(myr2)
            if ni not in r2dict:
                r2dict[ni]=[myr2]
    count=0
    for item in r2dict:
        l=r2dict[item]
        if count==0:
            plt.plot([item for i in range(len(l))],l,'o',alpha=0.03,color='blue',label='Full Data')
            #plt.violinplot([l],positions=[item])#,label='Full Data')
            plt.plot(item,avelist(l),'o',alpha=1,color='k',label='Data Average')
            #now get max r2
            maxr2=get_max_r2(cov,len(celllist),item)
            plt.plot(item,maxr2,'o',color='red',label='Theoretical Maximum')
            count+=1
        if count>0:
            plt.plot([item for i in range(len(l))],l,'o',alpha=0.03,color='blue')
            #plt.violinplot([l],positions=[item])#,color='blue')
            plt.plot(item,avelist(l),'o',alpha=1,color='k')
            #now get max r2
            maxr2=get_max_r2(cov,len(celllist),item)
            plt.plot(item,maxr2,'o',color='red')
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.xlabel('Number of SAGs with flexible gene', fontsize=25)
    plt.ylabel('r**2 for flexible gene\npresence/absence',fontsize=25)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.show()

    count=0
    for item in r2dict:
        l=r2dict[item]
        if count==0:
            plt.plot(item,avelist(l),'o',alpha=1,color='k',label='Data Average')
            #now get max r2
            maxr2=get_max_r2(cov,len(celllist),item)
            plt.plot(item,maxr2,'o',color='red',label='Theoretical Maximum')
            count+=1
        if count>0:
            plt.plot(item,avelist(l),'o',alpha=1,color='k')
            #now get max r2
            maxr2=get_max_r2(cov,len(celllist),item)
            plt.plot(item,maxr2,'o',color='red')
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.xlabel('Number of SAGs with flexible gene', fontsize=25)
    plt.ylabel('r**2 for flexible gene\npresence/absence',fontsize=25)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.show()


    count=0
    mypos=[item for item in r2dict]
    myls=[r2dict[mypos[i]] for i in range(len(mypos))]
    plt.violinplot(myls,positions=mypos,showmeans=True,showextrema=False)
    for item in r2dict:
        l=r2dict[item]
        if count==0:
            #plt.plot(item,avelist(l),'o',alpha=1,color='k',label='Data Average')
            #now get max r2
            maxr2=get_max_r2(cov,len(celllist),item)
            plt.plot(item,maxr2,'o',color='red',label='Theoretical Maximum')
            count+=1
        if count>0:
            #plt.plot(item,avelist(l),'o',alpha=1,color='k')
            #now get max r2
            maxr2=get_max_r2(cov,len(celllist),item)
            plt.plot(item,maxr2,'o',color='red')
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.xlabel('Number of SAGs with flexible gene', fontsize=25)
    plt.ylabel('r**2 for flexible gene\npresence/absence',fontsize=25)
    plt.legend(fontsize=15)
    plt.tight_layout()
    #plt.yscale('log')
    plt.show()


def avelist(ll):
    return float(sum(ll))/len(ll)


def get_max_r2_freal(coverage,ncellstot,freq):#freq is actually a number of cells
    #for a single gene, <n_obs>=n_real*cov
    #for a pair, <n_obs>=? n_real*cov**2 for f_ab
    f1o=freq*coverage/float(ncellstot)
    f12o=(freq*coverage**2)/float(ncellstot)
    r2=float((f12o-f1o**2)**2)/((f1o*(1.0-f1o)**2))
    return r2
           
def get_max_r2(coverage,ncellstot,freq):#freq is actually a number of cells
    #here, freq is observed, so already has the factor of *coverave built in
    #for a single gene, <n_obs>=n_real*cov
    #for a pair, <n_obs>=? n_real*cov**2 for f_ab
    f1o=freq/float(ncellstot)
    f12o=(freq*coverage)/float(ncellstot)
    r2=float((f12o-f1o**2)**2)/((f1o*(1.0-f1o)**2))
    return r2
           
    
##########################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-n','--n',type=int)
    parser.add_argument('-q','--qcut',type=int,default=80)
    parser.add_argument('-pid','--pidcut',type=int,default=80)
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    n=args.n
    qcut=args.qcut
    pidcut=args.pidcut

    if n==0:
        groupname='HLII'
        celllist=hliilist2
        sizecut=75#34
    if n==2:
        groupname='Blue_Basin_NoClose'
        celllist=BB_list
        sizecut=14
   


    if prog==0:
       
        
        #get_core_coverage(groupname,celllist)


        #HLII (167 cells): 1253 genes found.  Mean 117.572227 Cells Median 120.000000 Min 0 Max 138
        #use median: 120./167=72%
        

        linkage_matrix_flex_genes(celllist,groupname,sizecut,qcut,pidcut)

    if prog==1:
        make_HLII_flex_orf_list()

    if prog==2:
        #compare_flexible_groups_graph_pfam()
        get_linkage_from_pfam_file(celllist,groupname,sizecut)


    if prog==3:
        clustering_orthogroups_get_annotations(celllist,groupname,sizecut,qcut,pidcut)
    if prog==4:
        r2_null(celllist,groupname,sizecut,qcut,pidcut)
'''
using 0.3 cutoff,pfam

['phage prohead protease, HK97 family/phage major capsid protein, HK97 family,TIGR01554\n', 'phage prohead protease, HK97 family/phage major capsid protein, HK97 family,TIGR01554\n_2']
['alcohol dehydrogenase/S-(hydroxymethyl)glutathione dehydrogenase / alcohol dehydrogenase\n', 'alcohol dehydrogenase/S-(hydroxymethyl)glutathione dehydrogenase / alcohol dehydrogenase\n_2']
['VTC domain-containing protein\n', 'protein of unknown function (DUF4956)\n']
['CHAT domain-containing protein\n', 'Protein of unknown function (DUF3307)\n', 'adenylate cyclase\n', 'adenylate cyclase\n_2']
['EF hand\n', 'RNA polymerase sigma factor, sigma-70 family\n']
['Protein of unknown function (DUF1800)\n', 'Uncharacterized conserved protein, DUF1501 family\n']
['leader peptidase (prepilin peptidase) / N-methyltransferase\n', 'leader peptidase (prepilin peptidase) / N-methyltransferase\n_2', 'type IV pilus assembly protein PilM\n', 'type IV pilus assembly protein PilQ\n']
['cyanate lyase\n', 'hypothetical protein\n']
['SNF2 family N-terminal domain-containing protein\n', 'protein of unknown function (DUF1998)\n']
['aminopeptidase-like domain-containing protein\n', 'aminopeptidase-like domain-containing protein\n_2', 'aminopeptidase-like domain-containing protein\n_3']
['Superfamily I DNA and/or RNA helicase\n', 'Superfamily I DNA and/or RNA helicase\n_2', 'Superfamily I DNA and/or RNA helicase\n_3']
['Methyltransferase domain-containing protein\n_2', 'Methyltransferase domain-containing protein\n_3']
['Predicted kinase\n', 'uridine kinase\n']
[' T/G mismatch-specific endonuclease\n', 'AIPR protein\n', 'Putative PD-(D/E)XK family member, (DUF4420)\n', 'Z1 domain-containing protein\n']
['Coenzyme PQQ synthesis protein D (PqqD)\n', 'Transglutaminase-like superfamily protein\n', 'Uncharacterised nucleotidyltransferase\n']
['Dual specificity phosphatase, catalytic domain\n', 'capsular exopolysaccharide family\n']
['Heparinase II/III N-terminus\n', 'Heparinase II/III N-terminus\n_2']
['Immunity protein 17\n', 'type IV pilus assembly protein PilQ\n_2']
[' outer membrane transport energization protein ExbB\n', ' outer membrane transport energization protein ExbD\n']
['Choline kinase\n', 'Choline kinase\n_2', 'putative glutamine amidotransferase\n']
['Ribbon-helix-helix protein, copG family\n', 'Uncharacterized membrane protein\n']
['arsenite transporter, ACR3 family\n', 'iron(III) transport system ATP-binding protein\n']
[' methylthioribose-1-phosphate isomerase\n', 'Rhamnose utilisation protein RhaD, predicted bifunctional aldolase and dehydrogenase\n']
['nitrile hydratase\n', 'nitrile hydratase\n_2']
['N-6 DNA Methylase\n', 'type I restriction enzyme M protein\n', 'type I restriction enzyme, R subunit\n_2', 'type I restriction enzyme, R subunit\n_3', 'type I restriction enzyme, R subunit\n_4', 'type I restriction enzyme, S subunit\n']
['4Fe-4S dicluster domain-containing protein\n', 'HEAT repeat-containing protein\n', 'conserved hypothetical integral membrane protein\n', 'taurine transport system substrate-binding protein\n']
['Metallo-peptidase family M12B Reprolysin-like\n_2', 'Tetratricopeptide repeat-containing protein\n_2']
[' amino acid ABC transporter substrate-binding protein, PAAT family\n', 'Na+/H+-dicarboxylate symporter\n']
[' L-glutaminase\n', ' L-glutaminase\n_2']
['Concanavalin A-like lectin/glucanases superfamily protein\n', 'ribonucleoside-diphosphate reductase alpha chain\n']

Group lengths:
[2, 2, 2, 4, 2, 2, 4, 2, 2, 3, 3, 2, 2, 4, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 6, 4, 2, 2, 2, 2]


using 80 80 constructed orthogroups
['Unannotated', 'Unannotated', 'Unannotated', 'Unannotated']
['Unannotated', 'Unannotated']
['Pyr_redox_3', 'Unannotated', 'CN_hydrolase', 'Radical_SAM', 'Acetyltransf_1', 'AIRS_C', 'Unannotated']
['Unannotated', 'Unannotated', 'PsbN', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated']
['OprB', '2OG-FeII_Oxy_5', 'Sulfotransfer_3', 'Aconitase_C', 'DUF3769']
['Unannotated', 'Unannotated', 'LysM', 'Unannotated', 'Sigma70_r4_2', 'Unannotated', 'Inositol_P', 'Unannotated', 'TPR_8', 'Unannotated', 'Carbam_trans_C', 'Unannotated', 'DUF3769', 'Unannotated', 'Unannotated', 'Unannotated', 'EF-hand_5', 'Unannotated', 'N_methyl', 'Unannotated', 'OprB', 'Unannotated', 'Unannotated', 'Wzy_C', 'Unannotated', 'Methyltransf_25', 'DegT_DnrJ_EryC1', 'Bac_transf', 'Transketolase_C', 'E1_dh', 'Glycos_transf_1', 'Polysacc_synt_2', 'Unannotated', 'Unannotated']
['Unannotated', 'Unannotated']
['YjgP_YjgQ', 'Unannotated', 'GtrA']
['HlyD_3', 'ABC_tran', 'Unannotated', 'Unannotated']
['Unannotated', 'Unannotated']
['Unannotated', 'Unannotated', 'Unannotated', 'Unannotated']
['cobW', 'CobW_C']
['Porin_4', 'Sulfotransfer_3', 'Unannotated']
['Unannotated', 'Unannotated']
['Unannotated', 'Unannotated', 'Unannotated']
['ResIII', 'Unannotated', 'Unannotated', 'T5orf172']
['MFS_1', 'HTH_5']
['GCV_T_C', 'Unannotated']
['Unannotated', 'Unannotated', 'YjgP_YjgQ']
['Unannotated', 'cNMP_binding', 'Guanylate_cyc', 'POTRA_2', 'DUF3307', 'Unannotated', 'Unannotated', 'OTCace', 'GDP_Man_Dehyd', 'DUF3769', 'NTP_transferase', 'Glycos_transf_2', 'Bac_transf', 'Wzy_C', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Guanylate_cyc']
['Pribosyltran', 'Sulfotransfer_3', 'Porin_4', 'Sulfotransfer_3', 'Unannotated']
['DUF1028', 'M20_dimer', 'Asn_synthase', 'OTCace']
['Sulfotransfer_1', 'Unannotated', 'TPR_11']
['Unannotated', 'Unannotated', 'DUF1330']
['Unannotated', 'Autotransporter']
['Fer4_7', 'Succ_DH_flav_C', 'Unannotated']
['Unannotated', 'Unannotated']
['Unannotated', 'Unannotated']
['Unannotated', 'Unannotated', 'Unannotated', 'DNA_methylase']
['Glycos_transf_2', 'Glyco_transf_11', 'Abi']
['ZnuA', 'NMT1_2', 'BPD_transp_1', 'ABC_tran', 'Cyanate_lyase', 'Unannotated', 'DUF1348', 'Unannotated']
['Unannotated', 'Unannotated']
['ABC_tran', 'Glycos_transf_1', 'Wzy_C', 'PMT_2', 'Glycos_transf_2', 'Unannotated', 'Chloroa_b-bind', 'Unannotated']
['Unannotated', 'OprB', 'Unannotated', 'Gly_transf_sug', 'Unannotated', 'Glycos_transf_1', 'HlyD_3']
['Glycos_transf_1', 'Bac_transf']
['Unannotated', 'Unannotated', 'Unannotated']
['Unannotated', 'Unannotated', 'Unannotated', 'PDDEXK_1']
['Unannotated', 'Thiamine_BP']
['Unannotated', 'OprB']
['HAD_2', 'Unannotated', 'Polysacc_lyase', 'HTH_Crp_2', 'Phosphonate-bd', 'ABC_tran', 'BPD_transp_1', 'BPD_transp_1', 'PDZ', 'Unannotated', 'Unannotated', 'Unannotated', 'SBF', 'Unannotated', 'Unannotated', 'SBP_bac_6', 'TOBE_2', 'BPD_transp_1', 'HAD_2', 'DAO', 'HD', 'NUDIX', 'Unannotated', 'Amidohydro_2', 'Unannotated']
['Unannotated', 'OprB', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'TP_methylase', 'Inositol_P', 'ParA', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'DUF1049', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'UDPG_MGDP_dh_C', 'Wzy_C', 'Bac_transf', 'GDP_Man_Dehyd', 'M20_dimer', 'BPD_transp_1', 'NMT1', 'Amidase', 'ABC_tran', 'NMT1_2', 'BPD_transp_1', 'ABC_tran', 'HTH_1', 'Succ_DH_flav_C', 'Fer4_9', 'HEAT_2', 'Unannotated', 'Unannotated', 'LysR_substrate', 'Unannotated', 'Unannotated', 'Unannotated', 'OTCace']
['Unannotated', 'Unannotated']
['Unannotated', 'DUF1330']
['SdiA-regulated', 'Pentapeptide']
['Unannotated', 'Unannotated', 'LptC']
['YjgP_YjgQ', 'GtrA', 'DUF2029']
['UDPG_MGDP_dh_C', 'Wzy_C', 'Methyltransf_21']
['Mg_chelatase', 'VWA_2']
['Unannotated', 'Unannotated']
['Unannotated', 'Unannotated']
['Unannotated', 'OprB']
['Vsr', 'Unannotated', 'DNA_methylase']
['ABC_tran', 'Phosphonate-bd', 'BPD_transp_1', '2-Hacid_dh_C', 'Acetyltransf_10', 'PhyH', 'HD', 'OsmC', 'Unannotated']
['Hydrolase_3', 'DAHP_synth_1', 'CBS', 'Acetyltransf_5', 'Capsule_synth', 'NTP_transferase', 'Acyltransferase']
['Unannotated', 'Glyco_transf_41', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'Sulfotransfer_1', 'Unannotated', 'Unannotated', 'Unannotated']
['Unannotated', 'Autotransporter']
['Unannotated', 'Ysc84']
['Resolvase', 'YCII']
['DUF1651', 'Unannotated', 'Unannotated', 'DUF1330', 'Unannotated', 'Unannotated', 'Unannotated', 'Unannotated', 'PMT_2', 'GtrA', 'Unannotated', 'Unannotated', 'OprB', 'PsbY', 'Unannotated', 'NmrA', 'Unannotated', 'Unannotated', 'Unannotated', 'Ycf9', 'Unannotated']
['Unannotated', 'BPD_transp_1']
['Unannotated', 'VirE']
['Glycos_transf_2', 'POTRA']
['Peptidase_S13', 'Peptidase_S13']
['Abi', 'Glyco_transf_11', 'Glycos_transf_2']
['DUF2834', 'Unannotated']
['Met_synt_B12', 'Pterin_bind']
['Carbam_trans_C', 'Unannotated', 'Unannotated']
['Unannotated', 'DUF3721']
['Unannotated', 'DUF839']
['Unannotated', 'Unannotated']
['MannoseP_isomer', 'UDPG_MGDP_dh_C']
['Unannotated', 'DNA_ligase_aden', 'Unannotated', 'VirE']
Group lengths:
[4, 2, 7, 7, 5, 34, 2, 3, 4, 2, 4, 2, 3, 2, 3, 4, 2, 2, 3, 23, 5, 4, 3, 3, 2, 3, 2, 2, 4, 3, 8, 2, 8, 7, 2, 3, 4, 2, 2, 25, 42, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 3, 9, 7, 11, 2, 2, 2, 21, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 2, 2, 4]




some groups of cells identified using our graph orthogroups:

7 genes linked
{'AG-347-M18': 6, 'AG-355-K10': 7, 'AG-347-I06': 7, 'AG-347-J14': 7, 'AG-418-O03': 7, 'AG-355-O19': 7, 'AG-355-I20': 7, 'AG-347-L21': 7, 'AG-459-B06': 7, 'AG-355-J17': 7, 'AG-347-N23': 7, 'AG-347-M15': 7, 'AG-347-E23': 7, 'AG-418-B17': 6, 'AG-355-A18': 4, 'AG-355-I04': 7, 'AG-355-B23': 7, 'AG-347-K19': 6, 'AG-347-I23': 7, 'AG-347-G22': 3, 'AG-355-L22': 7, 'AG-347-B08': 7, 'AG-402-N08': 7, 'AG-347-L02': 7, 'AG-355-J09': 7, 'AG-402-N17': 6, 'AG-355-P15': 7, 'AG-429-C19': 7, 'AG-449-J16': 6, 'AG-347-L13': 5, 'AG-459-J14': 4}


7 genes
{'AG-347-B08': 5, 'AG-347-I04': 6, 'AG-347-E03': 3, 'AG-418-B17': 4, 'AG-418-O02': 7, 'AG-347-L20': 7, 'AG-347-L13': 7, 'AG-442-B03': 7, 'AG-347-K17': 4, 'AG-347-I15': 3, 'AG-418-G18': 6, 'AG-347-G20': 7, 'AG-402-N23': 6, 'AG-347-L02': 4, 'AG-402-K16': 6, 'AG-402-G23': 7, 'AG-418-J17': 7, 'AG-347-E23': 4, 'AG-347-I22': 4, 'AG-347-N23': 4, 'AG-347-G18': 6, 'AG-347-M08': 6, 'AG-347-I19': 6, 'AG-449-G23': 7, 'AG-418-F08': 6, 'AG-347-I23': 3, 'AG-347-C10': 2, 'AG-402-L23': 2, 'AG-402-M15': 2, 'AG-412-J13': 2, 'AG-418-L19': 1, 'AG-355-K10': 3, 'AG-418-M21': 3, 'AG-347-K15': 1, 'AG-402-I23': 1, 'AG-355-M18': 1, 'AG-429-E20': 1, 'AG-355-N22': 1, 'AG-355-N02': 1}



BELOW: ocean-specific gene grps

25 genes
average 24.360000 cells
Plurality sample 355: 29 of 41
1.000000 from Atl, 0.000000 from Pac (of 41)


3 genes
Plurality sample 355: 3 of 5
1.000000 from Atl, 0.000000 from Pac (of 5)

'''



'''
blue basin , 0.3 cut

all hypothetical grps deleted, and all but 1 hypo deleted



3 ['OprB', '2OG-FeII_Oxy_5', 'Sulfotransfer_3']
5 [ 'Fer2', 'Epimerase']
2 ['MFS_1', 'HTH_5']
2 ['GCV_T_C']
2 ['OprB', '2OG-FeII_Oxy_3']
6 ['POTRA_2', 'DUF3307',  'ZnuA',]
4 ['DUF1028', 'M20_dimer', 'Asn_synthase', 'OTCace']
4 ['hypothetical protein', 'hypothetical protein', 'DUF1330'']
29 ['DUF3303',, 'Autotransporter', 'HAD_2', 'Phosphonate-bd', 'ABC_tran', 'Amidohydro_2', 'BPD_transp_1', 'HAD_2', 'DAO', 'SBP_bac_6', 'ABC_tran', 'HD', 'NUDIX',  'SBF']
6 ['SLBB',  'OprB', 'UDPG_MGDP_dh_C', 'Wzy_C', 'Methyltransf_21']
2 ['Glycos_transf_1', 'Bac_transf']
2 [ 'Thiamine_BP']
4 ['GtrA',  'YjgP_YjgQ']
7 ['OTCace', 'NTP_transferase', 'Glycos_transf_2', 'Bac_transf', 'Wzy_C']
4 ['YjgP_YjgQ', 'GtrA', 'DUF2029', 'Polysacc_lyase']
2 ['Mg_chelatase', 'VWA_2']
2 ['SIMPL', 'Phage_integrase']
2 ['Hydrolase_3', 'CBS']
2 ['CobW_C', 'cobW']
7 ['ZnuA', 'NMT1_2', 'BPD_transp_1', 'ABC_tran', 'Cyanate_lyase', 'hypothetical protein', 'DUF1348']
'''


'''
hlii r2>.3


7 ['Pyr_redox_3', 'MSMEG_0572 family protein', 'CN_hydrolase', 'Radical_SAM', 'Acetyltransf_1', 'AIRS_C', 'MSMEG_0570 family protein']
5 ['OprB', '2OG-FeII_Oxy_5', 'Sulfotransfer_3', 'Aconitase_C', 'DUF3769']
34 ['LysM', 'hypothetical protein', 'Sigma70_r4_2',  'Inositol_P', 'hypothetical protein', 'TPR_8', 'hypothetical protein', 'Carbam_trans_C', 'hypothetical protein', 'DUF3769', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'EF-hand_5', 'hypothetical protein', 'N_methyl', 'hypothetical protein', 'OprB', 'hypothetical protein', 'hypothetical protein', 'Wzy_C', 'hypothetical protein', 'Methyltransf_25', 'DegT_DnrJ_EryC1', 'Bac_transf', 'Transketolase_C', 'E1_dh', 'Glycos_transf_1', 'Polysacc_synt_2', 'hypothetical protein', 'hypothetical protein']
3 ['YjgP_YjgQ', 'hypothetical protein', 'GtrA']
4 ['HlyD_3', 'ABC_tran', 'hypothetical protein', 'hypothetical protein']
2 ['cobW', 'CobW_C']
3 ['Porin_4', 'Sulfotransfer_3', 'hypothetical protein']
4 ['ResIII',  'T5orf172']
2 ['MFS_1', 'HTH_5']
2 ['GCV_T_C', 'no product']
23 ['hypothetical protein', 'cNMP_binding', 'Guanylate_cyc', 'POTRA_2', 'DUF3307', 'hypothetical protein', 'hypothetical protein', 'OTCace', 'GDP_Man_Dehyd', 'DUF3769', 'NTP_transferase', 'Glycos_transf_2', 'Bac_transf', 'Wzy_C', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'Guanylate_cyc']
5 ['Pribosyltran', 'Sulfotransfer_3', 'Porin_4', 'Sulfotransfer_3', 'hypothetical protein']
4 ['DUF1028', 'M20_dimer', 'Asn_synthase', 'OTCace']
3 ['Sulfotransfer_1', 'hypothetical protein', 'TPR_11']
3 ['Fer4_7', 'Succ_DH_flav_C', ' succinate dehydrogenase subunit C']
3 ['Glycos_transf_2', 'Glyco_transf_11', 'Abi']
8 ['ZnuA', 'NMT1_2', 'BPD_transp_1', 'ABC_tran', 'Cyanate_lyase', 'hypothetical protein', 'DUF1348', 'hypothetical protein']
2 ['hypothetical protein', 'hypothetical protein']
8 ['ABC_tran', 'Glycos_transf_1', 'Wzy_C', 'PMT_2', 'Glycos_transf_2', 'hypothetical protein', 'Chloroa_b-bind', 'hypothetical protein']
7 ['hypothetical protein', 'OprB', 'hypothetical protein', 'Gly_transf_sug', 'hypothetical protein', 'Glycos_transf_1', 'HlyD_3']
2 ['Glycos_transf_1', 'Bac_transf']
25 ['HAD_2', 'hypothetical protein', 'Polysacc_lyase', 'HTH_Crp_2', 'Phosphonate-bd', 'ABC_tran', 'BPD_transp_1', 'BPD_transp_1', 'PDZ', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'SBF', 'hypothetical protein', 'hypothetical protein', 'SBP_bac_6', 'TOBE_2', 'BPD_transp_1', 'HAD_2', 'DAO', 'HD', 'NUDIX', 'hypothetical protein', 'Amidohydro_2', 'hypothetical protein']
42 ['hypothetical protein', 'OprB', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'TP_methylase', 'Inositol_P', 'ParA', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'DUF1049', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'UDPG_MGDP_dh_C', 'Wzy_C', 'Bac_transf', 'GDP_Man_Dehyd', 'M20_dimer', 'BPD_transp_1', 'NMT1', 'Amidase', 'ABC_tran', 'NMT1_2', 'BPD_transp_1', 'ABC_tran', 'HTH_1', 'Succ_DH_flav_C', 'Fer4_9', 'HEAT_2', 'hypothetical protein', 'hypothetical protein', 'LysR_substrate', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'OTCace']
2 ['hypothetical protein', 'hypothetical protein']
2 ['hypothetical protein', 'DUF1330']
2 ['SdiA-regulated', 'Pentapeptide']
3 ['hypothetical protein', 'hypothetical protein', 'LptC']
3 ['YjgP_YjgQ', 'GtrA', 'DUF2029']
3 ['UDPG_MGDP_dh_C', 'Wzy_C', 'Methyltransf_21']
2 ['Mg_chelatase', 'VWA_2']
3 ['Vsr', 'hypothetical protein', 'DNA_methylase']
9 ['ABC_tran', 'Phosphonate-bd', 'BPD_transp_1', '2-Hacid_dh_C', 'Acetyltransf_10', 'PhyH', 'HD', 'OsmC', 'hypothetical protein']
7 ['Hydrolase_3', 'DAHP_synth_1', 'CBS', 'Acetyltransf_5', 'Capsule_synth', 'NTP_transferase', 'Acyltransferase']
11 ['hypothetical protein', 'Glyco_transf_41', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'Sulfotransfer_1', 'DNA-binding transcriptional regulator, MarR family', 'hypothetical protein', 'hypothetical protein']
2 ['hypothetical protein', 'Autotransporter']
2 ['hypothetical protein', 'Ysc84']
2 ['Resolvase', 'YCII']
21 ['DUF1651', 'hypothetical protein', 'hypothetical protein', 'DUF1330', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'PMT_2', 'GtrA', 'hypothetical protein', 'hypothetical protein', 'OprB', 'PsbY', 'hypothetical protein', 'NmrA', 'hypothetical protein', 'hypothetical protein', 'hypothetical protein', 'Ycf9', 'hypothetical protein']
2 ['phosphonate transport system permease protein', 'BPD_transp_1']
2 ['Glycos_transf_2', 'POTRA']
2 ['Peptidase_S13', 'Peptidase_S13']
3 ['Abi', 'Glyco_transf_11', 'Glycos_transf_2']
2 ['DUF2834', 'hypothetical protein']
2 ['Met_synt_B12', 'Pterin_bind']
3 ['Carbam_trans_C', 'hypothetical protein', 'hypothetical protein']
2 ['MannoseP_isomer', 'UDPG_MGDP_dh_C']
4 ['DNA_ligase_aden',  'VirE']
'''
