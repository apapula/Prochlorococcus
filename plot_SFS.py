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
from Cell_Lists import C1_without_ribotype,clique1_without_ribotype,hliilist2,BB_list,HLII_NoClose
import Repo_site_catalogues_single_pair as repo
from scipy.special import gamma,binom

'''
This is a script for plotting the data SFS and various model site frequency spectra.
'''


    
def make_sfs_plot(infile,covcut,tfna,greenbool=0):
    na,nt,nc,ng,gene,site=repo.load_single_sitefile('../input_data/%s' %infile)
    get_monomorphic_singletons(na,nt,nc,ng)
    ncg=[nc[i]+ng[i] for i in range(len(ng))]
    arraylist=[ncg]
    freqnamelist=['G+C']
    plot_sfs2(arraylist[0],covcut,groupname,freqnamelist[0],tfna,greenbool)

def get_monomorphic_singletons(na,nt,nc,ng):
    mono=[0,0,0,0]
    single=[0,0,0,0]#single[0] is the number of sites with singletons and the rest A
    for i in range(len(na)):
        ncells=na[i]+nt[i]+nc[i]+ng[i]
        if na[i]==ncells:
            mono[0]+=1
        if na[i]==ncells-1:
            single[0]+=1
        if nt[i]==ncells:
            mono[1]+=1
        if nt[i]==ncells-1:
            single[1]+=1
        if nc[i]==ncells:
            mono[2]+=1
        if nc[i]==ncells-1:
            single[2]+=1
        if ng[i]==ncells:
            mono[3]+=1
        if ng[i]==ncells-1:
            single[3]+=1
    print('A: %i monomorphic, %i singleton (%f)' %(mono[0],single[0],float(single[0])/mono[0]))
    print('T: %i monomorphic, %i singleton (%f)' %(mono[1],single[1],float(single[1])/mono[1]))
    print('C: %i monomorphic, %i singleton (%f)' %(mono[2],single[2],float(single[2])/mono[2]))
    print('G: %i monomorphic, %i singleton (%f)' %(mono[3],single[3],float(single[3])/mono[3]))
    print('Overall: %i monomorphic, %i singleton (%f)' %(sum(mono),sum(single),float(sum(single))/(sum(mono))))
    print('Overall: %i sites, %f mono, %f single' %(len(na),float(sum(mono))/len(na),float(sum(single))/len(na)))
def twofold_partner(allele):
    print('allele',allele)
    if allele=='A':
        return 'G'
    if allele=='G':
        return 'A'
    if allele=='T':
        return 'C'
    if allele=='C':
        return 'T'
    
def plot_sfs2(myarray,covcut,groupname,freqname,tfna,greenbool):
    
    intmybins=list(range(covcut+1))
    mybins=[intmybins[i]-0.5 for i in range(len(intmybins))]
    mybins.append(covcut+0.5)
        
    xticklabs=[]
    xtickpos=[]
    #want about 10 labels
    nticks=5
    perc=round(100/nticks,2)
    binwidth=float(covcut)/nticks
    #for i in range(nticks+1):
        #xtickpos.append(i*binwidth)
        #xticklabs.append('%s%%' %(str(perc*i)))
    xtickpos=[0,float(covcut)/5,2*float(covcut)/5,3*float(covcut)/5,4*float(covcut)/5,covcut]
    xticklabs=['0%','20%','40%','60%','80%','100%']

    sitetype_string=''
    if tfna=='t':
        sitetype_string='%s Twofold Degenerate Sites' %(freqname+twofold_partner(freqname))
    if tfna=='t':
        sitetype_string='Fourfold Degenerate Sites'
    if tfna=='n':
        sitetype_string='Nonsynonymous Sites (Second Bp)'
    if tfna=='a':
        sitetype_string='All Sites'

   
        
    mymids=[0.5*(mybins[i]+mybins[i+1]) for i in range(len(mybins)-1)]

    myweight=120./len(myarray)#want to sum to 1. the x axis is labeled as percentage, or 1/120 bins. so, myweight*(1/120)*len(myarray)=1?
    histcounts,histbins,_=plt.hist(myarray,bins=mybins,label='Data',weights=[myweight for i in range(len(myarray))])
    print(list(histcounts))

    
    
   
    if 0:#plot for mut rate
        mymaxx=45
        print(len(mymids),len(histcounts))
        plt.plot(mymids[:mymaxx],histcounts[:mymaxx],'o-',color='orange',label='Majority A/T')
        plt.plot(mymids[len(histcounts)-mymaxx:],histcounts[len(histcounts)-mymaxx:],'o-',color='black',label='Majority G/C')
        if greenbool:
            weights=[float(sum(histcounts[len(histcounts)-mymaxx:]))/sum(gsims) for i in range(len(gsims))]
            gsims=[gsims[i]*weights[0] for i in range(len(gsims))]
            plt.plot(mymids[len(histcounts)-mymaxx:],gsims[::-1],'o-',color='green',label='Convolution Simulation',markeredgecolor='white')
    plt.xlabel('Frequency %s' %freqname,fontsize=25)


    u,v=0,0
    if groupname=='HLII' and tfna=='f':
        u=.17953#Nmu_AT_to_GC
        v=.808373
        ncells=covcut
        #plot_neutral_discrete(u,v,ncells,weight=120)
        #####plot_neutral(u,v,ncells,weight=120)
        plot_DSF_hitchhiking(weight=120)
        #plot_expansion()
        plt.legend(fontsize=15)#prop={'size': 20})
    if tfna=='n':
        plt.title('Polarized Nonsynonymous\nSite Frequency Spectrum',fontsize=25)
    plt.xticks(ticks=xtickpos,labels=xticklabs,fontsize=25)
    plt.yticks(fontsize=25)
    plt.yscale('log')
    plt.tight_layout()
    ax=plt.gca()
    ax.set_ylim(bottom=.003)
    #plt.ylabel('Fraction (of %i Sites)' %(len(myarray)))
    #plt.title('%s %s Core Site %s Frequency Spectrum\nDownsampled to %i Cells' %(groupname,sitetype_string,freqname,covcut))
    plt.tight_layout()
    plt.show()

   
def plot_expansion():
    weight=120
    with open('../input_data/SFS4fold120Expansion.txt','r') as myinfile:
        SFS=[float(line.rstrip())*weight for line in myinfile if len(line)>0]
    plt.plot(list(range(121)),SFS,'-',label='Exponential Expansion',color='orange',linewidth=3)

        


def plot_neutral(u,v,ncells,weight=1):
    #want to plot a curve
    x=np.linspace(0,ncells,200)
    norm=gamma(u+v)/(gamma(u)*gamma(v))*(1./ncells)*weight
    rhof=((x/ncells)**(u-1))*((1-x/ncells)**(v-1))*norm
    plt.plot(x,rhof,label='Neutral Drift Theory',color='red',linewidth=3)#color was k, changed at DSF's request 10-29-24

def plot_neutral_discrete2(u,v,ncells,weight=1):
    #want to plot a curve
    x=np.linspace(0,ncells,ncells+1)#200)
    norm=gamma(u+v)/(gamma(u)*gamma(v))*(1./ncells)*weight
    rhof=((x/ncells)**(u-1))*((1-x/ncells)**(v-1))*norm
    discreterhof=norm*binom(n,x)*(gamma(x+u)*gamma(n-x+v))/(gamma(n+u+v))
    for x in range(ncells+1):
        print(x,rhof,discreterhof)
    plt.plot(x,discreterhof,label='Neutral Drift Theory',color='k',linewidth=3)

def plot_neutral_discrete(u,v,ncells,weight=1,mycolor='k'):
    #want to plot a curve
    #x=np.linspace(0,ncells,ncells+1)#200)
    
    norm=gamma(u+v)/(gamma(u)*gamma(v))*ncells#*(1./ncells)*weight
    myx,mydiscrete=[],[]
    for x in range(ncells+1):
        myx.append(x)
        #rhof=((x/ncells)**(u-1))*((1-x/ncells)**(v-1))*norm
        discreterhof=norm*binom(ncells,x)*(gamma(x+u)*gamma(ncells-x+v))/(gamma(ncells+u+v))
        print(binom(ncells,x)*(gamma(x+u)*gamma(ncells-x+v))/(gamma(ncells+u+v)))
        mydiscrete.append(discreterhof)
        #print(x,discreterhof)
    plt.plot(myx[1:-1],mydiscrete[1:-1],'-',label='Neutral Drift Theory',color=mycolor,linewidth=3)#was 'x'
    #plt.show()
    print(mydiscrete[0],mydiscrete[-1])
    plt.plot([0],[mydiscrete[0]],'_',color=mycolor)#'red')#,markersize=10)
    plt.plot([ncells],[mydiscrete[-1]],'_',color=mycolor)#'red')#,markersize=10)
    ax=plt.gca()
    ax.axvline(0.5 , ls='--', lw=0.5, c=mycolor)#'red')
    ax.axvline((ncells - 0.5) , ls='--', lw=0.5, c=mycolor)#'red')#color was k, changed at dsf's request 10-29-24
    #plt.show()
    
def plot_DSF_hitchhiking(weight=1):
    with open('../input_data/SFS4fold120HitchHLII.txt','r') as myinfile:
        SFS=[float(line.rstrip())*weight for line in myinfile if len(line)>0]
    plt.plot(list(range(121)),SFS,'-',label='Hitchhiking Model',color='r',linewidth=3)#green',linewidth=3), color was red, changed at DSF's request 10-29-24


    

def get_convolution(glist,ncells,nconvolve,downsamplesize):
    covcut=downsamplesize
    alist=[downsamplesize-glist[i] for i in range(len(glist))]
    intmybins=list(range(covcut+1))
    mybins=[intmybins[i]-0.5 for i in range(len(intmybins))]
    mybins.append(downsamplesize+0.5)
        
    mycounts,mbins,_=plt.hist(alist,bins=mybins)#must be downsampled
    plt.close()
    #print(mycounts)
    #t=input('t')
    counts_amaj=mycounts[::-1][:ncells+1]
    counts_gmaj=mycounts[:ncells+1]
    my_pdf=[float(counts_amaj[i])/sum(counts_amaj) for i in range(len(counts_amaj))]
    freqs=list(range(ncells+1))#len(mycounts)))
    number_of_draws=math.ceil((nconvolve+2)*sum(counts_gmaj))               
    draws=np.random.choice(freqs,number_of_draws,p=my_pdf)
    print(my_pdf)
    gsims=[0 for i in range(len(counts_gmaj))]#counts of freqs, not freqs themselves
    i=0
   
    while sum(gsims)<sum(counts_gmaj):
        mysim=sum(draws[i*nconvolve:(i+1)*nconvolve])
        if mysim<len(gsims):
            gsims[mysim]+=1
        i+=1
    scaled_countsa=[float(counts_amaj[i])*sum(counts_gmaj)/sum(counts_amaj) for i in range(len(counts_amaj))]
    freqs=list(range(ncells+1))
    return gsims



    
        
if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--n',type=int,help='0: HLII; 1: C1; 2: BB')
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    parser.add_argument('-ds','--downsamplebool',type=int,default=1)
    parser.add_argument('-tfna','--tfna',type=str,help='Twofold Fourfold Nonsyn (Second nt) or All SNPs')
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
    n=args.n
    ds=args.downsamplebool
    
    if n==0:
        groupname='HLII'
        celllist=hliilist2
        downsamplesize=120
    if n==1:
        groupname='C1'
        celllist=C1_without_ribotype
        downsamplesize=37#math.floor(0.7*len(celllist))
    if n==2:
        groupname='BB_NoClose'
        celllist=BB_list#HLII_BB_NoCloseCells_Median04
        downsamplesize=48#math.floor(0.7*len(celllist))

    if n==3:
        celllist=clique1_without_ribotype
        downsamplesize=7#math.floor(0.7*len(celllist))
        groupname='Clique1'

    if n==5:
        celllist=HLII_NoClose
        groupname='HLII_NoClose'
        downsamplesize=math.floor(0.7*len(celllist))


  
    suffix2='All_Same_AA_without_6fold'
    if n==5:
        suffix1=suffix2
    tfna=args.tfna
    if tfna=='t':       
        outfile2='%s_TWOFOLDDEGENERATE_Single_Sites_NEWDS_%s.npz' %(groupname,suffix2)
    if tfna=='f':
        outfile2='%s_FOURFOLDDEGENERATE_Single_Sites_NEWDS_%s.npz' %(groupname,suffix2)
        outfile=outfile2
    if tfna=='n':
        outfile='%s_NONSYN_SECONDNT_Single_Sites_NEWDS.npz' %(groupname)
    if tfna=='a':
        outfile='%s_All_Single_Sites_NEWDS.npz' %(groupname)
        
    if prog==0:
        covcut=downsamplesize

        if tfna in ['n','a']:
            outfile='Single_Site_Catalogs/%s' %outfile
            #outfile='../input_data/Single_Site_Catalogs/%s' %outfile
            outfileds=outfile[:-4]+'_DownSample_%iCells.npz' %downsamplesize
            print('outfile %s' %outfile)
            
            if not ds:
                make_sfs_plot(outfile,covcut,tfna)
            if ds:
                make_sfs_plot(outfileds,covcut,tfna)

            
        if tfna in ['t','f']:

           

            outfile2='../input_data/Single_Site_Catalogs/%s' %outfile2
            outfileds2=outfile2[:-4]+'_DownSample_%iCells.npz' %downsamplesize

            
            print('outfile downsampled: %s' %outfileds2)
            if ds:
                make_sfs_plot(outfileds2,covcut,tfna,greenbool=1)
            if not ds:
                make_sfs_plot(outfile2,covcut,tfna,greenbool=1)
   
