import numpy as np
import ast
import csv
import math
from scipy import stats
import networkx as nx
import random
from os.path import exists
import statistics
import argparse
import matplotlib.pyplot as plt
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list
import Repo_site_catalogues_single_pair as repo
import Make_Synonymous_Distance_Matrices_Berube_Kashtan_Using_All_Sites_Per_Pair_and_HalfGene_Record_sites_SNPs as distmatrepo
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mtick
from matplotlib import colors


def make_hists_INTERgene(outfile,infilepair,downsamplesize,bins, agtc=''):#bins, intra or inter...?
    print(infilepair)
    gene1,site1,gene2,site2,myns=repo.load_INTERgene_pair_file(infilepair)
    #may want to separated AG vs TC
    #don't forget monomorphic
    #for each site, get into AA, AG,GA,GG plus TC if agtc==''
    newns,newtrash,newgene1,newgene2=make_aa_ag_ga_gg(myns,agtc,site1,gene1,gene2)
    print('DOne')
    INTERgene_now_make_hists(bins,newns,newgene1,newgene2,downsamplesize,outfile)
    
def INTERgene_now_make_hists(bins,newns,gene1,gene2,downsamplesize,outfile):
    mydict={}
    for i in range(len(bins)-1):
        print(i)
        mydict['%s_%s' %(str(bins[i]),str(bins[i+1]))]=[[[0 for i in range(downsamplesize+1)] for j in range(downsamplesize+1)] for k in range(downsamplesize+1)]
    for i in range(len(gene1)):
        if i%1000==0:
            print('%f (%i of %i)' %(float(i)/len(gene1),i,len(gene1)))
        g1,g2=gene1[i],gene2[i]
        mykey=find_bin_key(bins,g1,g2)
        #print('key',bins,s1,s2,mykey)
        
        if mykey=='-':
            continue
        mynums=[newns[j][i] for j in range(4)]
        #print('mynums',mynums)
        mydict[mykey][mynums[0]][mynums[1]][mynums[2]]+=1
        #print(mydict[mykey])
        #print( mydict[mykey][mynums[0]][mynums[1]][mynums[2]])
        #mydict[mykey][mynums[0]][mynums[1]][mynums[2]][mynums[3]]+=1
    #now save
    np.save(outfile,mydict)

########################################################################
def make_hists_intragene(outfile,infilepair,downsamplesize,bins, agtc='',mediancut=0):#bins, intra or inter...?
    if mediancut:
        infilepair=infilepair[:-4]+'_MEDIAN_LENGTH_CUT.npz'
        outfile+='_MEDIAN_LENGTH_CUT'
    gene,site1,site2,myns=repo.load_pair_file(infilepair)#intragene
    #may want to separated AG vs TC
    #don't forget monomorphic
    #for each site, get into AA, AG,GA,GG plus TC if agtc==''
    newns,newgene,newsite1,newsite2=make_aa_ag_ga_gg(myns,agtc,gene,site1,site2)
    print(len(gene),len(newgene))
    #print(newns)
    print('New sites 1 sample: ',newsite1[:30])
    #t=input('t')
    print('DOne')
    make_hists(bins,newns,newgene,newsite1,newsite2,downsamplesize,outfile)
    
def make_hists(bins,newns,gene,site1,site2,downsamplesize,outfile):
    mydict={}
    for i in range(len(bins)-1):
        print(i)
        mydict['%s_%s' %(str(bins[i]),str(bins[i+1]))]=[[[0 for i in range(downsamplesize+1)] for j in range(downsamplesize+1)] for k in range(downsamplesize+1)]
        #mydict['%s_%s' %(str(bins[i]),str(bins[i+1]))]=[[[[0 for i in range(downsamplesize+1)] for j in range(downsamplesize+1)] for k in range(downsamplesize+1)] for m in range(downsamplesize+1)]
    #indicies: mat[aa][ag][ga][gg]
    for i in range(len(gene)):
        if i%1000==0:
            print('%f (%i of %i)' %(float(i)/len(gene),i,len(gene)))
        s1,s2=site1[i],site2[i]
        mykey=find_bin_key(bins,s1,s2)
        #print('key',bins,s1,s2,mykey)
        
        if mykey=='-':
            continue
        mynums=[newns[j][i] for j in range(4)]
        #print('mynums',mynums)
        mydict[mykey][mynums[0]][mynums[1]][mynums[2]]+=1
        #print(mydict[mykey])
        #print( mydict[mykey][mynums[0]][mynums[1]][mynums[2]])
        #mydict[mykey][mynums[0]][mynums[1]][mynums[2]][mynums[3]]+=1
    #now save
    np.save(outfile,mydict)

        
def find_bin_key(bins,site1,site2):
    sep=abs(site1-site2)
    for i in range(len(bins)-1):
        if sep>=bins[i] and sep<bins[i+1]:
            mykey='%s_%s' %(bins[i],bins[i+1])
            return mykey
    return '-'
    
def make_aa_ag_ga_gg(myns,agtc,gene,site1,site2):
    arr=[[],[],[],[]]
    ngene,nsite1,nsite2=[],[],[]
    for i in range(len(gene)):
        
        if agtc=='f':#fourfold. add gc at
            l=[myns[0][i]+myns[1][i]+myns[4][i]+myns[5][i],myns[2][i]+myns[3][i]+myns[6][i]+myns[7][i],myns[9][i]+myns[8][i]+myns[12][i]+myns[13][i],myns[10][i]+myns[11][i]+myns[14][i]+myns[15][i]]
            if sum(l)==0:
                continue
            ngene.append(gene[i])
            nsite1.append(site1[i])
            nsite2.append(site2[i])
            for j in range(4):
                arr[j].append(l[j])
        if agtc=='ag':
            l=[myns[0][i],myns[3][i],myns[12][i],myns[15][i]]
            if sum(l)==0:
                continue
            ngene.append(gene[i])
            nsite1.append(site1[i])
            nsite2.append(site2[i])
            for j in range(4):
                arr[j].append(l[j])
        if agtc=='tc':
            l=[myns[5][i],myns[6][i],myns[9][i],myns[10][i]]
            if sum(l)==0:
                continue
            ngene.append(gene[i])
            nsite1.append(site1[i])
            nsite2.append(site2[i])
            for j in range(4):
                arr[j].append(l[j])
        if agtc=='':
            '''
            print(gene[i],site1[i],site2[i],'trying')
            l=[myns[0][i]+myns[5][i],myns[6][i]+myns[3][i],myns[12][i]+myns[9][i],myns[15][i]+myns[10][i]]
            if sum(l)==0:#this is from comparing a pair of sites that has one site AG and one site TC. this gives AT entries, AC entries, etc. but we don't want to exclude these, so change
                print('problem',site1[i],site2[i],l)
                t=input('t')
                continue
            '''
            #print(gene[i],site1[i],site2[i],'trying')
            l=[myns[0][i]+myns[5][i]+myns[1][i]+myns[4][i],myns[2][i]+myns[7][i]+myns[6][i]+myns[3][i],myns[12][i]+myns[9][i]+myns[8][i]+myns[13][i],myns[11][i]+myns[14][i]+myns[15][i]+myns[10][i]]
            if sum(l)==0:#this is from comparing a pair of sites that has one site AG and one site TC. this gives AT entries, AC entries, etc. but we don't want to exclude these, so change
                print('problem',site1[i],site2[i],l)
                #t=input('t')
                continue
            ngene.append(gene[i])
            nsite1.append(site1[i])
            nsite2.append(site2[i])
            for j in range(4):
                arr[j].append(l[j])
    return arr,ngene,nsite1,nsite2


##############################
def make_shuffled_matrix_from_matrix(matrixfile,downsamplesize):
    #rng = np.random.default_rng()

    outfile=matrixfile[:-4]+'_SHUFFLED'
    matarrays=distmatrepo.load_arrays(matrixfile)
    mykeys=list(matarrays.keys())
    shufmats={}
    for item in mykeys:
        print(item)
        shufmat=get_shuffle_single_mat(matarrays[item],downsamplesize)
        shufmats[item]=shufmat
    np.save(outfile,shufmats)
        
def get_shuffle_single_mat(inmat,downsamplesize):

    newmat=[[[0 for i in range(downsamplesize+1)] for j in range(downsamplesize+1)] for k in range(downsamplesize+1)]
    for  i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                weight=inmat[i][j][k]
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                if weight>0 and downsamplesize-i-j-k<0:
                    print('Error: negative freqs: ',i,j,k,downsamplesize-i-j-k,' weight ',myweight)
                if weight<=0:
                    continue
                fg1=ga+gg
                fg2=ag+gg
                #print('fg1 %i fg2 %i' %(fg1,fg2))
                list1=['G' for a in range(fg1)]+['A' for a in range(downsamplesize-fg1)]
                list2=['G' for a in range(fg2)]+['A' for a in range(downsamplesize-fg2)]
                for w in range(weight):
                    l=[0,0,0,0]#new aa, ag,ga gg
                    alleles1=list(np.random.choice(np.array(list1),replace=False,size=downsamplesize))
                    #if fg1<downsamplesize-1 and fg1>0:
                        #print(alleles1)
                        #print(alleles1.index('A'))
                        #t=input('t')
                    alleles2=list(np.random.choice(np.array(list2),replace=False,size=downsamplesize))
                    for c in range(len(alleles1)):
                        if alleles1[c]=='A':
                            if alleles2[c]=='A':
                                l[0]+=1
                                continue
                            if alleles2[c]=='G':
                                l[1]+=1
                                continue
                        if alleles1[c]=='G':
                            if alleles2[c]=='A':
                                l[2]+=1
                                continue
                            if alleles2[c]=='G':
                                l[3]+=1
                                continue
                    #now add this single point to the new matrix
                    newmat[l[0]][l[1]][l[2]]+=1
    #now return the matrix, save the dictionary of matrices
    return newmat
                    
    

    
#############################################################################
def test_matrices(matrixfile,downsamplesize,outfilesingle):
    matarrays=distmatrepo.load_arrays(matrixfile)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        test_individual_matrix(mat,downsamplesize,item,outfilesingle,agtc)
        t=input('t')

def sum_3d_mat(mat,downsamplesize):
    mysum=0
    for i in range(downsamplesize):
        for j in range(downsamplesize):
            myarr=mat[i][j]
            mysum+=sum(myarr)
    print('sum',mysum)
    
    return mysum

def test_individual_matrix(mat,downsamplesize,mykey,outfilesingle,agtc):
    fa=[]
    fg=[]
    weights=[]
    for i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                weight=mat[i][j][k]
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                if weight>0 and downsamplesize-i-j-k<0:
                    print('Error: negative freqs: ',i,j,k,downsamplesize-i-j-k,' weight ',myweight)
                #print(weight,aa,ag,ga,gg)
                fa.append(aa+ag)
                weights.append(weight)
                fa.append(ga+aa)
                weights.append(weight)
                fg.append(ga+gg)
                fg.append(ag+gg)


   
        
    mybins=[]
    for i in range(downsamplesize+2):
        mybins.append(i-0.5)
    plt.hist(fg,weights=weights,bins=mybins,label='Site Pairs',density=1)
    plot_single_site_sfs(outfilesingle,downsamplesize,agtc)
    plt.yscale('log')
    plt.legend()
    plt.title('freq G from matrices %s' %mykey)
    plt.show()

def plot_single_site_sfs(outfilesingle,downsamplesize,agtc):
    covcut=downsamplesize
         
    mybins=[]
    for i in range(downsamplesize+2):
        mybins.append(i-0.5)
        
    na,nt,nc,ng,gene,site=repo.load_single_sitefile(outfilesingle)
    na,nt,nc,ng=repo.check_if_twofold_degen(na,nt,nc,ng,gene,site)
    print('%i AG sites, %i TC Sites' %(len(na),len(nt)))
    nmcag=[na[i]+ng[i] for i in range(len(na))]
    nmctc=[nc[i]+nt[i] for i in range(len(nt))]
    alist=[float(na[i]) for i in range(len(na)) if nmcag[i]>=covcut]
    glist=[float(ng[i]) for i in range(len(ng)) if nmcag[i]>=covcut]
    tlist=[float(nt[i]) for i in range(len(nt)) if nmctc[i]>=covcut]
    clist=[float(nc[i]) for i in range(len(nc)) if nmctc[i]>=covcut]
    if agtc=='ag':
        plt.hist(glist,bins=mybins,density=1,label='Single Site')
    if agtc=='tc':
        plt.hist(clist,bins=mybins,density=1,label='Single Site')
    if agtc=='':
        plt.hist(clist+glist,bins=mybins,density=1,label='Single Site')

#########################################################
def get_r2_from_matricesfile(matrixfileintra,matrixfileinter,downsamplesize,agtcstring,groupname,mono):
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfileintra=matrixfileintra[:-4]+'_SHUFFLED.npy'
    shuffledfileinter=matrixfileinter[:-4]+'_SHUFFLED.npy'
    r2,sep,r2shuf=[],[],[]
    r2intra,sepintra,r2shufintra=[],[],[]
    matarrays=distmatrepo.load_arrays(matrixfileintra)
    shufarrays=distmatrepo.load_arrays(shuffledfileintra)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sepintra.append(get_sep_from_key(item))
        r2intra.append(get_r2_single_matrix(mat,downsamplesize,mono))
        r2shufintra.append(get_r2_single_matrix(shufmat,downsamplesize,mono))
        #sep.append(get_sep_from_key(item))
        #r2.append(get_r2_single_matrix(mat,downsamplesize,mono))
        #r2shuf.append(get_r2_single_matrix(shufmat,downsamplesize,mono))

    matarrays=distmatrepo.load_arrays(matrixfileinter)
    shufarrays=distmatrepo.load_arrays(shuffledfileinter)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        print(item)
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item)*600)
        r2.append(get_r2_single_matrix(mat,downsamplesize,mono))
        r2shuf.append(get_r2_single_matrix(shufmat,downsamplesize,mono))
    plt.plot(sepintra,r2intra,'*',color='blue',label='Intragene Data',markersize=9)
    plt.plot(sepintra,r2shufintra,'*',color='orange',markersize=9)#label='Random Shuffle',markersize=9)
    print(sepintra)
    print(sep)
    plt.plot(sep,r2,'o',label='Inter-gene Data',color='blue',markersize=9)
    plt.plot(sep,r2shuf,'o',label='Random Shuffle',color='orange',markersize=9)
    plt.xscale('log')
    plt.legend(prop={'size': 15})
    #plt.title('%s %s Twofold-Site <r**2>\nDownsampled to %i Cells' %(groupname,agtcstring,downsamplesize))

    if groupname=='HLII':#tick labels
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)
    if groupname=='BB_NoClose':
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)
        
    plt.xlabel('Separation (bp)',fontsize=25,labelpad=10)
    plt.ylabel('<r**2>',fontsize=25,labelpad=10)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.gca().set_ylim(top=max(r2intra)+.01)
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    plt.show()


def get_r2_from_matricesfile_justintra(matrixfileintra,downsamplesize,agtcstring,groupname,mono,matrixfile2='',agtcstring2=''):
    print(matrixfileintra)
    print(matrixfile2)
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfileintra=matrixfileintra[:-4]+'_SHUFFLED.npy'
    
   
    r2intra,sepintra,r2shufintra=[],[],[]
    matarrays=distmatrepo.load_arrays(matrixfileintra)
    shufarrays=distmatrepo.load_arrays(shuffledfileintra)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sepintra.append(get_sep_from_key(item))
        r2intra.append(get_r2_single_matrix(mat,downsamplesize,mono))
        r2shufintra.append(get_r2_single_matrix(shufmat,downsamplesize,mono))
        #sep.append(get_sep_from_key(item))
        #r2.append(get_r2_single_matrix(mat,downsamplesize,mono))
        #r2shuf.append(get_r2_single_matrix(shufmat,downsamplesize,mono))

    if matrixfile2!='':
        shuffledfile2=matrixfile2[:-4]+'_SHUFFLED.npy'
        r2,sep,r2shuf=[],[],[]

        matarrays=distmatrepo.load_arrays(matrixfile2)
        shufarrays=distmatrepo.load_arrays(shuffledfile2)
        mykeys=list(matarrays.keys())
        for item in mykeys:
            mat=matarrays[item]
            shufmat=shufarrays[item]
            if sum_3d_mat(mat,downsamplesize)==0:
                continue
           
            sep.append(get_sep_from_key(item))
            r2.append(get_r2_single_matrix(mat,downsamplesize,mono))
            r2shuf.append(get_r2_single_matrix(shufmat,downsamplesize,mono))
        plt.plot(sep,r2,'o',color='red',label='Intragene %s' %agtcstring2,markersize=9)
        plt.plot(sep,r2shuf,'o',color='green',label='Shuffle %s' %agtcstring2,markersize=9)
    plt.plot(sepintra,r2intra,'*',color='blue',label='Intragene %s' %agtcstring,markersize=9)
    plt.plot(sepintra,r2shufintra,'*',color='orange',markersize=9,label='Shuffle %s' %agtcstring)#label='Random Shuffle',markersize=9)
    #plt.plot(sep,r2,'o',label='Inter-gene Data',color='blue',markersize=9)
    #plt.plot(sep,r2shuf,'o',label='Random Shuffle',color='orange',markersize=9)
    plt.xscale('log')
    plt.legend(prop={'size': 15})
    if 'MEDIAN' in matrixfileintra:
        groupname+=' Genes >= Median Length\n'
    plt.title('%s %s Sites <r**2>\nDownsampled to %i Cells' %(groupname,agtcstring,downsamplesize),fontsize=25)

    if groupname=='HLII':#tick labels
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)
    if groupname=='BB_NoClose':
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)
        
    plt.xlabel('Separation (bp)',fontsize=25,labelpad=10)
    plt.ylabel('<r**2>',fontsize=25,labelpad=10)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.gca().set_ylim(top=max([max(r2intra),max(r2)])+.01)
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    plt.show()


def get_r2_from_matricesfile_intragene(matrixfile,downsamplesize,agtc,groupname):
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfile=matrixfile[:-4]+'_SHUFFLED.npy'
    r2,sep,r2shuf=[],[],[]
    matarrays=distmatrepo.load_arrays(matrixfile)
    shufarrays=distmatrepo.load_arrays(shuffledfile)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item))
        r2.append(get_r2_single_matrix(mat,downsamplesize))
        r2shuf.append(get_r2_single_matrix(shufmat,downsamplesize))
    plt.plot(sep,r2,'o',label='Data')
    plt.plot(sep,r2shuf,'o',label='Random Shuffle')
    plt.xscale('log')
    plt.legend()
    plt.gca().set_ylim(bottom=0)
    plt.show()



        
def get_sep_from_key(mykey):
    idx=mykey.find('_')
    sep1=int(float(mykey[:idx]))
    sep2=int(float(mykey[idx+1:]))
    midpoint=0.5*(sep1+sep2-1)
    return midpoint

def get_r2_single_matrix(mat,downsamplesize,mono):
    r2_list=[]
    weights=[]
    for i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                weight=mat[i][j][k]
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                if aa+ag==0 or aa+ga==0 or gg+ag==0 or gg+ga==0:
                    if not mono:
                        continue#monomorphic
                    if mono:
                        r2_list.append(0)
                        weights.append(weight)
                        continue
                if weight>0 and downsamplesize-i-j-k<0:
                    print('Error: negative freqs: ',i,j,k,downsamplesize-i-j-k,' weight ',weight)
                weights.append(weight)
                r2=float(((ga+gg)*(ag+gg)-downsamplesize*(gg))**2)/((ga+gg)*(downsamplesize-ga-gg)*(ag+gg)*(downsamplesize-ag-gg))
                r2_list.append(r2)
    #now, get average r2
    aver2=float(sum([r2_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    return aver2


def get_sigmaD2_single_matrix(mat,downsamplesize,mono):
    num_list=[]
    denom_list=[]
    weights=[]
    for i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                weight=mat[i][j][k]
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                if aa+ag==0 or aa+ga==0 or gg+ag==0 or gg+ga==0:
                   continue#monomorphic
                if weight>0 and downsamplesize-i-j-k<0:
                    print('Error: negative freqs: ',i,j,k,downsamplesize-i-j-k,' weight ',myweight)
                weights.append(weight)
                num=float(((ga+gg)*(ag+gg)-downsamplesize*(gg))**2)
                denom=((ga+gg)*(downsamplesize-ga-gg)*(ag+gg)*(downsamplesize-ag-gg))
                num_list.append(num)
                denom_list.append(denom)
    #multiply numerator and denom by weight july 1 2024
    #aver2=float(sum([r2_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    avenum=float(sum([num_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    avedenom=float(sum([denom_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    avesD=float(avenum)/avedenom
    return avesD
                

#################################
def freq_freq_corr(matrixfile,downsamplesize,groupname,agtc,freqbins,pdffile,intragenebool):
    #mykeys will be used in the plot titles
   
    freqbins_withmono=[0]
    for item in freqbins:
        freqbins_withmono.append(item)
    freqbins_withmono.append(downsamplesize)
    r2dict,enrichdict,nsnpsdict,incodict,fulllinkdict={},{},{},{},{}#dicts of matrices for various separations
    r2avedict,incoavedict,flavedict={},{},{}
    matarrays=distmatrepo.load_arrays(matrixfile)
    mykeys=list(matarrays.keys())
    for i in range(len(mykeys)):
        item=mykeys[i]
        mat=matarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        #want: r2 mat, freqfreq corr(SNP enrichment), nsnsps, incompatibilities
        #put average of each (e.g. r2, inco) in plot title
        r2mat,enrichmat,nsnpsmat,incomat,fullylinkmat,aver2,aveinco,avefl=get_r2_inco_SNPenrich_from_mat(mat,downsamplesize,freqbins,freqbins_withmono)
        r2dict[item]=r2mat
        enrichdict[item]=enrichmat
        incodict[item]=incomat
        nsnpsdict[item]=nsnpsmat
        fulllinkdict[item]=fullylinkmat
        flavedict[item]=avefl
        r2avedict[item]=aver2
        incoavedict[item]=aveinco
    #now plot
    make_matrices_pdf(pdffile,groupname,agtc,freqbins,freqbins_withmono,r2dict,enrichdict,incodict,nsnpsdict,fulllinkdict,r2avedict,incoavedict,flavedict,downsamplesize,intragenebool)
    
def make_matrices_pdf(pdffile,groupname,agtc,freqbins,freqbins_withmono,r2dict,enrichdict,incodict,nsnpsdict,fulllinkdict,r2avedict,incoavedict,flavedict,downsamplesize,intragenebool):
    agtcstring='AG and TC'
    agdict={'ag':'AG','tc':'TC'}
    if agtc in agdict:
        agtcstring=agdict[agtc]
    pdf=PdfPages(pdffile)
    myseps=list(r2dict.keys())
    xtickpos=[]
    for ff in range(len(freqbins)-1):
        if freqbins[ff]==freqbins[ff+1]-1:
            xtickpos.append(ff)
        else:
            xtickpos.append(ff-0.5)
    xtickpos.append(-1.5+len(freqbins))
    xticklabs=[]
    for f in range(len(freqbins)):
        xticklabs.append('%s' %(str(math.floor(freqbins[f]))))
    xtickposmono=[]
    for ff in range(len(freqbins_withmono)-1):
        if freqbins_withmono[ff]==freqbins_withmono[ff+1]-1:
            xtickposmono.append(ff)
        else:
            xtickposmono.append(ff-0.5)
    xtickposmono.append(-1.5+len(freqbins_withmono))
    xticklabsmono=[]
    for f in range(len(freqbins_withmono)):
        xticklabsmono.append('%s' %(str(math.floor(freqbins_withmono[f]))))
        
    for i in range(len(myseps)):
        key=myseps[i]
        fig,axes=plt.subplots(nrows=3,ncols=2)
        fig.set_size_inches([7,10])
       
        for j in range(3):
            for k in range(2):
                axes[j,k].tick_params(axis='both', which='major', labelsize=5)
                ax=axes[j,k]
                if [j,k]==[0,0]:
                    ax.title.set_text('<r**2> (Ave %f)' %(r2avedict[key]))
                    im=ax.imshow(r2dict[key])
                    ax.set_xticks(xtickpos)
                    ax.set_yticks(xtickpos)
                    ax.set_xticklabels(xticklabs)
                    ax.set_yticklabels(xticklabs)
                if [j,k]==[0,1]:
                    ax.title.set_text('Incompatibilities (Ave %f)' %(incoavedict[key]))
                    im=ax.imshow(incodict[key])
                    ax.set_xticks(xtickpos)
                    ax.set_yticks(xtickpos)
                    ax.set_xticklabels(xticklabs)
                    ax.set_yticklabels(xticklabs)
                if [j,k]==[1,0]:
                    ax.title.set_text('SNP Enrichment')
                    im=ax.imshow(enrichdict[key])
                    ax.set_xticks(xtickposmono)
                    ax.set_yticks(xtickposmono)
                    ax.set_xticklabels(xticklabsmono)
                    ax.set_yticklabels(xticklabsmono)
                if [j,k]==[2,0]:
                    ax.title.set_text('Fully Linked')
                    im=ax.imshow(fulllinkdict[key])
                    ax.set_xticks(xtickpos)
                    ax.set_yticks(xtickpos)
                    ax.set_xticklabels(xticklabs)
                    ax.set_yticklabels(xticklabs)
                if [j,k]==[1,1]:
                    ax.title.set_text('Number of Sites')
                    im=ax.imshow(nsnpsdict[key],norm=colors.LogNorm())
                    ax.set_xticks(xtickposmono)
                    ax.set_yticks(xtickposmono)
                    ax.set_xticklabels(xticklabsmono)
                    ax.set_yticklabels(xticklabsmono)
                if [j,k]==[2,1]:
                    continue
                ax.set_xlabel('Frequency C/G (of %i)' %downsamplesize)
                ax.set_ylabel('Frequency C/G (of %i)' %downsamplesize)
        
                divider=make_axes_locatable(ax)
                cax=divider.append_axes('right',size='5%',pad=0.05)
                fig.colorbar(im,cax=cax,orientation='vertical')
        if intragenebool:
            fig.suptitle('%s %s Twofold Sites\nSNP Pairs Separated by %s Sites' %(groupname,agtcstring,key))
        if not intragenebool:
            fig.suptitle('%s %s Twofold Sites\nSNP Pairs Separated by %s Genes' %(groupname,agtcstring,key),y=.98)
        #plt.tight_layout()
        #fig.subplots_adjust(top=0.8)
        #fig.tight_layout()
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        pdf.savefig()
        plt.close()
    pdf.close()
     
def get_r2_inco_SNPenrich_from_mat(mat,downsamplesize,freqbins,freqbins_withmono):
    #freqbins=[1,2,3,.....,n-3,n-3,n-1]
    r2mat=[[[] for i in range(len(freqbins))] for y in range(len(freqbins))]
    incomat=[[[] for i in range(len(freqbins))] for y in range(len(freqbins))]
    fullylinkmat=[[[] for i in range(len(freqbins))] for y in range(len(freqbins))]
    #bin by Freq G, C
   
    
    nsnpsmat=[[0 for i in range(len(freqbins_withmono))] for y in range(len(freqbins_withmono))]
    enrichmat=[[np.nan for i in range(len(freqbins_withmono))] for y in range(len(freqbins_withmono))]

    
    for i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                #print(np.shape(np.array(mat)),i,j,k)
                weight=mat[i][j][k]
                if weight==0:
                    continue
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                fg1=ga+gg
                fg2=ag+gg
                nbin1,nbin2=get_mat_index(fg1,freqbins_withmono),get_mat_index(fg2,freqbins_withmono)
                if downsamplesize-i-j-k<0:
                    print('Negative freqs')
                    print(i,j,k,downsamplesize-i-j-k,weight)
                #print(i,j,k,downsamplesize-i-j-k,fg1,fg2,nbin1,nbin2)
                nsnpsmat[nbin1][nbin2]+=weight
                if aa+ag==0 or aa+ga==0 or gg+ag==0 or gg+ga==0:
                    continue#monomorphic
                l=[aa,ag,ga,gg]
                #get r2. add this, times the weight, to the appropriate r2 bins
                r2=float(((ga+gg)*(ag+gg)-downsamplesize*(gg))**2)/((ga+gg)*(downsamplesize-ga-gg)*(ag+gg)*(downsamplesize-ag-gg))
               
                bin1=get_mat_index(fg1,freqbins)
                bin2=get_mat_index(fg2,freqbins)
                for z in range(weight):
                    r2mat[bin1][bin2].append(r2)
                    if l.count(0)==0 and fg1>1 and fg2>1 and fg1<downsamplesize-1 and fg2<downsamplesize-1:
                        incomat[bin1][bin2].append(1)
                    if l.count(0)>=1 and fg1>1 and fg2>1 and fg1<downsamplesize-1 and fg2<downsamplesize-1:
                        incomat[bin1][bin2].append(0)
                    if l.count(0)==2:
                        fullylinkmat[bin1][bin2].append(1)
                    if l.count(0)>2:
                        fullylinkmat[bin1][bin2].append(0)

    #now get snp number and snp enrichment mats. make sure there's monomorphic included
    get_snp_enrichment_from_number_mat(nsnpsmat,enrichmat)
    r2mat,aver2=average_matrix(r2mat)
    incomat,aveinco=average_matrix(incomat,incobool=1)
    fullylinkmat,avefl=average_matrix(fullylinkmat)
    return r2mat,enrichmat,nsnpsmat,incomat,fullylinkmat,aver2,aveinco,avefl

def average_matrix(mat,incobool=0):#if incobool, don't include singletons
    mydim=len(mat)
    avemat=[[np.nan for i in range(mydim)] for j in range(mydim)]
    totcount,totquant=0,0
    for i in range(mydim):
        for j in range(mydim):
            l=mat[i][j]
            if len(l)==0:
                print('No entries ',i,j)
                #print(mat)
                continue
            
            totcount+=len(l)
            totquant+=sum(l)
            avemat[i][j]=repo.avelist(l)
    avequant=float(totquant)/totcount
    return avemat,avequant
   

def get_snp_enrichment_from_number_mat(nsnpsmat,enrichmat):
    totsnps=0
    mydim=len(nsnpsmat)
    for i in range(mydim):
        totsnps+=sum(nsnpsmat[i])
    for i in range(mydim):
        for j in range(i+1):
            sumi=sum(nsnpsmat[i])
            sumj=sum(nsnpsmat[j])
            freq=nsnpsmat[i][j]
            if sumi>0 and sumj>0:
                newfreq=freq*totsnps/(sumi*sumj)
                enrichmat[i][j]=newfreq
                enrichmat[j][i]=newfreq
                
def get_mat_index(freq,freqbins):
    for i in range(len(freqbins)-1):
        if freq>=freqbins[i] and freq<freqbins[i+1]:
            return i
    if freq==freqbins[-1]:
        return len(freqbins)-1


def get_sigmaD2_from_matricesfile(matrixfileintra,matrixfileinter,downsamplesize,agtcstring,groupname,mono):
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfileintra=matrixfileintra[:-4]+'_SHUFFLED.npy'
    shuffledfileinter=matrixfileinter[:-4]+'_SHUFFLED.npy'
    if agtcstring=='Fourfold':
        shuffledfileinter='../input_data/Intragene_SNP_Pair_Catalogs/Matrices/%s_POSSIBLYFOURFOLD_INTERgene_Matrices_AG_and_TC_SHUFFLED.npy' %groupname
    r2,sep,r2shuf=[],[],[]
    r2intra,sepintra,r2shufintra=[],[],[]
    matarrays=distmatrepo.load_arrays(matrixfileintra)
    shufarrays=distmatrepo.load_arrays(shuffledfileintra)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sepintra.append(get_sep_from_key(item))
        r2intra.append(get_sigmaD2_single_matrix(mat,downsamplesize,mono))
        r2shufintra.append(get_sigmaD2_single_matrix(shufmat,downsamplesize,mono))

    matarrays=distmatrepo.load_arrays(matrixfileinter)
    shufarrays=distmatrepo.load_arrays(shuffledfileinter)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        print(item)
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item)*600)
        r2.append(get_sigmaD2_single_matrix(mat,downsamplesize,mono))
        r2shuf.append(get_sigmaD2_single_matrix(shufmat,downsamplesize,mono))
    plt.plot(sepintra,r2intra,'*',color='blue',label='Intragene Data',markersize=9)
    print(r2intra)
    plt.plot(sepintra,r2shufintra,'*',color='orange',markersize=9)
    plt.plot(sep,r2,'o',label='Data',markersize=9)
    plt.plot(sep,r2shuf,'o',label='Random Shuffle',markersize=9)
    plt.xscale('log')
    plt.legend(prop={'size': 15})

    if 0:
        monostr=''
        if agtcstring!='Fourfold':
            plt.title('%s %s Twofold-Site <sigma_D**2>\nDownsampled to %i Cells%s' %(groupname,agtcstring,downsamplesize,monostr))
        if agtcstring=='Fourfold':
            plt.title('%s %s AT vs CG <sigma_D**2>\nDownsampled to %i Cells%s' %(groupname,agtcstring,downsamplesize,monostr))
    if 0:#tick labels
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)

    plt.xlabel('Separation (bp)',fontsize=25,labelpad=10)
    #plt.ylabel('<sigma_D**2>',fontsize=25,labelpad=10)
    plt.ylabel(r'$\sigma_D^2$',fontsize=25,labelpad=10)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.gca().set_ylim(top=max(r2intra)+.01)
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    plt.show()

def get_expression_intragene_compare_medianlengthcut(matrixfileintra_cut,matrixfileintra,downsamplesize,agtcstring,groupname,mono,expressionname):
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfileintracut=matrixfileintra_cut[:-4]+'_SHUFFLED.npy'
    shuffledfileintra=matrixfileintra[:-4]+'_SHUFFLED.npy'

    print(matrixfileintra)
    print('cut', matrixfileintra_cut)
    
    r2,sep,r2shuf=[],[],[]
    cr2,csep,cr2shuf=[],[],[]#Cut median length
    matarrays=distmatrepo.load_arrays(matrixfileintra)
    shufarrays=distmatrepo.load_arrays(shuffledfileintra)

    print('uncut')
    matarrayscut=distmatrepo.load_arrays(matrixfileintra_cut)
    shufarrayscut=distmatrepo.load_arrays(shuffledfileintracut)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item))
        r2.append(get_expression_single_matrix(mat,downsamplesize,mono,expressionname))
        r2shuf.append(get_expression_single_matrix(shufmat,downsamplesize,mono,expressionname))

   

    ##########################
    print('cut')
    for item in list(matarrayscut.keys()):
        mat=matarrayscut[item]
        shufmat=shufarrayscut[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        csep.append(get_sep_from_key(item))
        cr2.append(get_expression_single_matrix(mat,downsamplesize,mono,expressionname))
        cr2shuf.append(get_expression_single_matrix(shufmat,downsamplesize,mono,expressionname))

   

        
    plt.plot(sep,r2,'o',label='Data',markersize=9)
    plt.plot(sep,r2shuf,'o',label='Random Shuffle',markersize=9)
    plt.plot(sep,cr2,'o',label='Median Cut Data',markersize=9)
    plt.plot(sep,cr2shuf,'o',label='Median Cut Random Shuffle',markersize=9)
    plt.xscale('log')
    plt.legend()#prop={'size': 15})
   
    monostr=''
    if agtcstring!='Fourfold':
        plt.title('%s %s Twofold-Site %s\nDownsampled to %i Cells%s' %(groupname,agtcstring,expressionname,downsamplesize,monostr))
    if agtcstring=='Fourfold':
        plt.title('%s %s AT vs CG %s\nDownsampled to %i Cells%s' %(groupname,agtcstring,expressionname,downsamplesize,monostr))
    if 0:#tick labels
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)

    plt.xlabel('Separation (bp)',fontsize=25,labelpad=10)
    plt.ylabel('%s' %expressionname,fontsize=25,labelpad=10)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.gca().set_ylim(top=max(r2)+.01)
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    plt.show()

def get_expression_single_matrix(mat,downsamplesize,mono,expressionname):
    if expressionname=='sigmaD2':
        l=get_sigmaD2_single_matrix(mat,downsamplesize,mono)
        return l
    if expressionname=='2<f_AB f_ab + f_Ab f_aB>':
        l=get_exp1_single_matrix(mat,downsamplesize,mono)
        return l
    if expressionname=='<(1-f_a)*(1-f_b)*(f_ab-f_a*f_b)>':
        l=get_exp2_single_matrix(mat,downsamplesize,mono)
        return l
########expression 1

def get_exp1_from_matricesfile(matrixfileintra,matrixfileinter,downsamplesize,agtcstring,groupname,mono):
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfileintra=matrixfileintra[:-4]+'_SHUFFLED.npy'
    shuffledfileinter=matrixfileinter[:-4]+'_SHUFFLED.npy'
    if agtcstring=='Fourfold':
        shuffledfileinter='../input_data/Intragene_SNP_Pair_Catalogs/Matrices/%s_POSSIBLYFOURFOLD_INTERgene_Matrices_AG_and_TC_SHUFFLED.npy' %groupname
    r2,sep,r2shuf=[],[],[]
    matarrays=distmatrepo.load_arrays(matrixfileintra)
    shufarrays=distmatrepo.load_arrays(shuffledfileintra)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item))
        r2.append(get_exp1_single_matrix(mat,downsamplesize,mono))
        r2shuf.append()

    matarrays=distmatrepo.load_arrays(matrixfileinter)
    shufarrays=distmatrepo.load_arrays(shuffledfileinter)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        print(item)
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item)*600)
        r2.append(get_exp1_single_matrix(mat,downsamplesize,mono))
        r2shuf.append(get_exp1_single_matrix(shufmat,downsamplesize,mono))
    plt.plot(sep,r2,'o',label='Data',markersize=9)
    plt.plot(sep,r2shuf,'o',label='Random Shuffle',markersize=9)
    plt.xscale('log')
    plt.legend(prop={'size': 15})
    monostr='\nDimorphic Only'
    if mono:
        monostr='\nMonomorphic Included, Subtracting Pi**2'
    if agtcstring!='Fourfold':
        plt.title('%s %s Twofold-Site 2<f_AB f_ab + f_Ab f_aB>\nDownsampled to %i Cells%s' %(groupname,agtcstring,downsamplesize,monostr))
    if agtcstring=='Fourfold':
        plt.title('%s %s AT vs CG 2<f_AB f_ab + f_Ab f_aB>\nDownsampled to %i Cells%s' %(groupname,agtcstring,downsamplesize,monostr))
   

    if 0:#tick labels
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)

    plt.xlabel('Separation (bp)',fontsize=25,labelpad=10)
    plt.ylabel('2<f_AB f_ab + f_Ab f_aB>',fontsize=25,labelpad=10)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #plt.gca().set_ylim(top=max(r2)+.01)
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    plt.show()


def get_exp1_single_matrix(mat,downsamplesize,mono):#probability different at both sites
    val_list=[]
    weights=[]
    pisqtot,pisqweight=0,0
    for i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                weight=mat[i][j][k]
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                #how to get pi squared average...can we take, for each pair of sites, the pi squared for each site and average?
                pisqweight+=weight+weight
                fg1=float(ga+gg)/downsamplesize
                fg2=float(gg+ag)/downsamplesize
                pisqtot+=weight*(fg1*(1.0-fg1)+fg2*(1.0-fg2))*2
                #square before or after? AFTER, ASSUMES INDEPENDENCE. average first, then square the average
                if aa+ag==0 or aa+ga==0 or gg+ag==0 or gg+ga==0:
                    if mono:
                        weights.append(weight)
                        val_list.append(0)
                    continue#monomorphic
                if weight>0 and downsamplesize-i-j-k<0:
                    print('Error: negative freqs: ',i,j,k,downsamplesize-i-j-k,' weight ',myweight)
                weights.append(weight)
                val=2*(gg*aa+ag*ga)/(downsamplesize**2)
                #subtract off pi^2, the unlinked probabitliy diff at both sites
                val_list.append(val)
    #multiply numerator and denom by weight july 1 2024
    #aver2=float(sum([r2_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    aveval=float(sum([val_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    if mono:
        avepi=float(pisqtot)/pisqweight
        aveval-=avepi**2
    return aveval



########expression 2

def get_exp2_from_matricesfile(matrixfileintra,matrixfileinter,downsamplesize,agtcstring,groupname,mono):
    #will want: r2, separation
    #what about shuffle? fa1, fa2 fixed.  fa1a2, fa1g2 vary e.g.
    shuffledfileintra=matrixfileintra[:-4]+'_SHUFFLED.npy'
    shuffledfileinter=matrixfileinter[:-4]+'_SHUFFLED.npy'
    if agtcstring=='Fourfold':
        shuffledfileinter='../input_data/Intragene_SNP_Pair_Catalogs/Matrices/%s_POSSIBLYFOURFOLD_INTERgene_Matrices_AG_and_TC_SHUFFLED.npy' %groupname
    r2,sep,r2shuf=[],[],[]
    matarrays=distmatrepo.load_arrays(matrixfileintra)
    shufarrays=distmatrepo.load_arrays(shuffledfileintra)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item))
        r2.append(get_exp2_single_matrix(mat,downsamplesize,mono))
        r2shuf.append(get_exp2_single_matrix(shufmat,downsamplesize,mono))

    matarrays=distmatrepo.load_arrays(matrixfileinter)
    shufarrays=distmatrepo.load_arrays(shuffledfileinter)
    mykeys=list(matarrays.keys())
    for item in mykeys:
        print(item)
        mat=matarrays[item]
        shufmat=shufarrays[item]
        if sum_3d_mat(mat,downsamplesize)==0:
            continue
        sep.append(get_sep_from_key(item)*600)
        r2.append(get_exp2_single_matrix(mat,downsamplesize,mono))
        r2shuf.append(get_exp2_single_matrix(shufmat,downsamplesize,mono))
    plt.plot(sep,r2,'o',label='Data',markersize=9)
    plt.plot(sep,r2shuf,'o',label='Random Shuffle',markersize=9)
    plt.xscale('log')
    plt.legend(prop={'size': 15})
    monostr='\nDimorphic Only'
    if mono:
        monostr='\nMonomorphic Included'
    if agtcstring!='Fourfold':
        plt.title('%s %s Twofold-Site <(1-f_a)*(1-f_b)*(f_ab-f_a*f_b)>\nDownsampled to %i Cells%s' %(groupname,agtcstring,downsamplesize,monostr))
    if agtcstring=='Fourfold':
        plt.title('%s %s AT vs CG <(1-f_a)*(1-f_b)*(f_ab-f_a*f_b)>\nDownsampled to %i Cells%s' %(groupname,agtcstring,downsamplesize,monostr))

    if 0:#tick labels
        yticklabs,ytickpos=['0.00','0.02','0.04','0.06'],[0,.02,.04,.06]
        plt.yticks(ticks=ytickpos,labels=yticklabs)

    plt.xlabel('Separation (bp)',fontsize=25,labelpad=10)
    plt.ylabel('<(1-f_a)*(1-f_b)*(f_ab-f_a*f_b)>',fontsize=25,labelpad=10)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #plt.gca().set_ylim(top=max(r2)+.01)
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    plt.show()


def get_exp2_single_matrix(mat,downsamplesize,mono):
    val_list=[]
    weights=[]
    for i in range(downsamplesize+1):
        for j in range(downsamplesize+1):
            for k in range(downsamplesize+1):
                weight=mat[i][j][k]
                [aa,ag,ga,gg]=[i,j,k,downsamplesize-i-j-k]
                if aa+ag==0 or aa+ga==0 or gg+ag==0 or gg+ga==0:
                    if mono:
                        val_list.append(0)
                        weights.append(weight)
                    continue#monomorphic
                if weight>0 and downsamplesize-i-j-k<0:
                    print('Error: negative freqs: ',i,j,k,downsamplesize-i-j-k,' weight ',myweight)
                weights.append(weight)
                f1=(ga+gg)/float(downsamplesize)
                f2=(ag+gg)/float(downsamplesize)
                fgg=gg/float(downsamplesize)
                val=(1.0-2*f1)*(1.0-2*f2)*(fgg-f1*f2)
                val_list.append(val)
    #multiply numerator and denom by weight july 1 2024
    #aver2=float(sum([r2_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    aveval=float(sum([val_list[i]*weights[i] for i in range(len(weights))]))/sum(weights)
    return aveval

#############################################################
