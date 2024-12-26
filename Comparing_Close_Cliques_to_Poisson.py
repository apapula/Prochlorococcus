

import csv
import numpy as np
import ast
import csv
import math
from scipy import stats
from os.path import exists
import statistics
import argparse
import matplotlib.pyplot as plt
import re
import Load_Sequences_Repo as load_sequences_repo
from difflib import SequenceMatcher
from Bio import Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.cluster.hierarchy as hac
import Kashtan_Alignments_Pairwise_Event_Lengths as kashrepo
import Kashtan_Alignments_Record_Gaps_Between_Highdiv_Events as kashclusterrepo
from Cell_Lists import clique1_without_ribotype,C1_without_ribotype,clique1,cl2,cl3,cl4,cl5,cl6
from scipy.stats import linregress, pearsonr


def go_through_clique_pairs_and_compare_to_poisson(chunklength=1000,gabriel_method=True):
    fieldnames=['Cell1','Cell2','Manual_Reco_Length','Poisson_Excess_Length','Lambda','Freco_Fitting','F_reco_excess','Clique']
    totmc,totmc_highdiv=[],[]
    labels=[]
    recodict=read_in_close_clique_events()
    mystr=''#_NoCovCut'
    outfile='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp_Gabriel_methods%s.csv' %(chunklength,mystr)
    if not gabriel_method:
         outfile='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp%s.csv' %(chunklength,mystr)
    cliques=['1','2','3','4','5','6']
    celllist_list=[clique1,cl2,cl3,cl4,cl5,cl6]
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,fieldnames=fieldnames,delimiter='\t')
        outwriter.writeheader()
        for i in range(len(cliques)):
            mc,mc_highdiv=do_one_clique(outwriter,celllist_list[i],cliques[i],chunklength,recodict,gabriel_method)
            #totmc.append(mc)
            #totmc_highdiv.append(mc_highdiv)
            #totmc.append(mc_highdiv)
            #labels.append('Cl %s All Chunks' %cliques[i])
            #labels.append('Cl %s High Div Chunks' %cliques[i])
            
            for item in mc:
                totmc.append(item)
            for item in mc_highdiv:
                totmc_highdiv.append(item)
    #plt.hist(totmc,label=labels,bins=100)
    
    plt.hist(np.array([totmc,totmc_highdiv],dtype=object),label=['All chunks','High-Div Chunks'],bins=30,density=True)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Number of bp mutually covered, of %i' %chunklength)
    plt.ylabel('Number of %i-bp Chunks' %chunklength)
    plt.title('Comparing Mutual Coverage for cell pairs in all regions vs high divergence regions\nWhat mutual coverage cutoff should we set?')
    plt.show()

def do_one_clique(outwriter,celllist,clique,chunklength,recodict,gabriel_method):
    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
    totmc,totmc_highdiv=[],[]
    for i in range(len(celllist)):
        for j in range(i):
            recoamount=get_recombined_length_for_cellpair_from_recodict(recodict,celllist,i,j,clique)
            mc,mc_highdiv=get_poisson_single_pair(i,j,celllist,chunk_seqdict,chunklength,recoamount,clique,outwriter,gabriel_method,show_or_write='show')
            for item in mc:
                totmc.append(item)
            for item in mc_highdiv:
                totmc_highdiv.append(item)
    return totmc,totmc_highdiv

def get_recombined_length_for_cellpair_from_recodict(recodict,celllist,i,j,clique):
    mydict=recodict[clique]
    cellpair=sorted([celllist[i],celllist[j]])
    cellpair=(cellpair[0],cellpair[1])
    if cellpair not in mydict:
        return 0
    if cellpair in mydict:
        return mydict[cellpair]
    
            
def get_poisson_single_pair(i,j,celllist,chunk_seqdict,chunksize,recoamount,clique,outwriter,gabriel_method,show_or_write='write'):
    covfrac=80
    covcut=0.01*covfrac*chunksize
    celli=celllist[i]
    cellj=celllist[j]
    seqi=chunk_seqdict[celli]
    seqj=chunk_seqdict[cellj]
    nsnps,mc=[],[]
    uncovered=0
    mc_highdiv=[]
    #do for all kb with e.g. at least 80% of sites, and also only for completely covered kb
    #how many kb have fewer than 80% covered?
    nchunks=len(seqi)
    for c in range(nchunks):
        si=seqi[c]
        sj=seqj[c]
        mydiff,mymc=get_div_chunk(si,sj)
        if mymc<covcut:
            uncovered+=1
            continue
        nsnps.append(mydiff)
        if mydiff>8:
            #print(celli,cellj,c*chunksize,(c+1)*chunksize,'%i snps %i mc' %(mydiff,mymc))
            #t=input('t')
            mc_highdiv.append(mymc)
        if mymc>7:
            mc.append(mymc)
    mymax=max(nsnps)
    mybins=[-.5]
    for i in range(mymax+1):
        mybins.append(i+.5)
    
    poisvals,fasex,mylambda,excess=get_poisson_model(nsnps,chunksize,gabriel_method)
    if show_or_write=='show':
        plt.hist(np.array([nsnps,list(range(mymax+1))],dtype=object),weights=np.array([[1 for i in range(len(nsnps))],[len(nsnps)*poisvals[j] for j in range(len(poisvals))]],dtype=object),label=['Data','Poisson\nF_asex %f; Lambda %f' %(fasex,mylambda)],bins=mybins)
        plt.legend()
        plt.xlabel('SNPs per %i bp' %chunksize)
        plt.ylabel('Number of %i-bp chunks at least %s %% Covered' %(chunksize,str(covfrac)))
        plt.title('%s vs %s SNPs\nper %i-bp along C1 Composite Reference\n%i bp found recombined; %i bp excess of Poisson' %(celli,cellj,chunksize,recoamount,excess))
        plt.yscale('log')
        plt.show()
        #t=input('t')

    if show_or_write=='write':
        outwriter.writerow({'Cell1':celli,'Cell2':cellj,'Manual_Reco_Length':recoamount,'Poisson_Excess_Length':excess,'Lambda':mylambda,'Freco_Fitting':1.0-fasex,'F_reco_excess':float(excess/chunksize)/len(nsnps),'Clique':clique})
    return mc,mc_highdiv

def get_poisson_model(nsnps,chunksize,gabriel_method):
    if not gabriel_method:
        #mypois,my_fasex,my_lambda,excess=get_poisson_model_old_method(nsnps,chunksize)
        mypois,my_fasex,my_lambda,excess=get_poisson_model_mean_method(nsnps,chunksize)
    if gabriel_method:
        mypois,my_fasex,my_lambda,excess=get_poisson_model_gabriel_method(nsnps,chunksize)
    return mypois,my_fasex,my_lambda,excess

def get_poisson_model_gabriel_method(nsnps,chunksize):
    f0=float(nsnps.count(0))/len(nsnps)
    f1=float(nsnps.count(1))/len(nsnps)
    my_fasex=1
    my_lambda=-math.log(f0)
    mymax=max(nsnps)
    mypois=[]
    probcut=1.0/len(nsnps)
    ncut=0
    for n in range(mymax+1):
        poisval=my_fasex*np.exp(-my_lambda)*(my_lambda**n)/math.factorial(n)
        mypois.append(poisval)
        if ncut==0:
            if poisval<=probcut:
                ncut=n
    #for excess here, we're going to take a prob cutoff =1/(len(nsnps)) and find the minimum number of snps per chunksize whose probability is less than this.  any chunk with this number of snps or more is counted as excess
    excess=0
    for item in nsnps:
        if item>=ncut:
            excess+=chunksize
    ###dec 27 BUT NOW ALSO WANT TO SUBTRACT HOW MANY OF THIS LENGTH WE'D EXPECT FROM THE POISSON
    expected_excess=get_poisson_expectation_above_cutoff_to_subtract(my_lambda,ncut,nsnps,chunksize)

    excess-=expected_excess
    ###
    
    return mypois,my_fasex,my_lambda,excess

def get_poisson_expectation_above_cutoff_to_subtract(my_lambda,ncut,nsnps,chunklength):
    #sum up the poisson probability from 0 to ncut-1
    #then, subtract this from 1
    #then, multiply by len(nsnps)*chunklength
    prob=0
    for i in range(ncut):
        poisval=np.exp(-my_lambda)*(my_lambda**i)/math.factorial(i)
        prob+=poisval
    expected_excess=len(nsnps)*chunklength*(1.0-prob)
    return expected_excess

def get_poisson_model_mean_method(nsnps,chunksize):#to show how bad the fit is if we include the high div regions
    my_lambda=float(sum(nsnps))/(len(nsnps))
    
    my_fasex=1
    
    mymax=max(nsnps)
    mypois=[]
    for n in range(mymax+1):
        poisval=my_fasex*np.exp(-my_lambda)*(my_lambda**n)/math.factorial(n)
        mypois.append(poisval)
    print('F asex %f; lambda %f' %(my_fasex,my_lambda))
    excess=get_excess_over_poisson(nsnps,mypois,chunksize)
    print('%i bp excess of poisson (%f asex)' %(excess, 1.0-float(excess/chunksize)/len(nsnps)))

    if 0:
        poisnorm=sum(mypois)*len(nsnps)
        print([len(nsnps)*mypois[i] for i in range(len(mypois))])
        print('excess %f; nsnps-poisnorm %f' %(excess,chunksize*(len(nsnps)-poisnorm)))
        t=input('t')
    return mypois,my_fasex,my_lambda,excess
    
def get_poisson_model_old_method(nsnps,chunksize):
    f0=float(nsnps.count(0))/len(nsnps)
    f1=float(nsnps.count(1))/len(nsnps)
    my_lambda=f1/f0
    my_fasex=f0*np.exp(my_lambda)
    #my_fasex=1
    #my_lambda=-math.log(f0)
    mymax=max(nsnps)
    mypois=[]
    for n in range(mymax+1):
        poisval=my_fasex*np.exp(-my_lambda)*(my_lambda**n)/math.factorial(n)
        mypois.append(poisval)
    print('F asex %f; lambda %f' %(my_fasex,my_lambda))
    excess=get_excess_over_poisson(nsnps,mypois,chunksize)
    print('%i bp excess of poisson (%f asex)' %(excess, 1.0-float(excess/chunksize)/len(nsnps)))

    if 0:
        poisnorm=sum(mypois)*len(nsnps)
        print([len(nsnps)*mypois[i] for i in range(len(mypois))])
        print('excess %f; nsnps-poisnorm %f' %(excess,chunksize*(len(nsnps)-poisnorm)))
        t=input('t')
    return mypois,my_fasex,my_lambda,excess

def get_excess_over_poisson(nsnps,mypois,chunksize):
    excess=0
    for i in range(3,max(nsnps)+1):
        poisval=mypois[i]*len(nsnps)
        datval=nsnps.count(i)
        if datval<int(poisval):
            print('Data value is less than poisson value: data %i , pois %f, snps %i' %(datval,poisval,i))
        if datval>poisval:
            excess+=chunksize*(datval-poisval)
    return excess

def get_div_chunk(seq1,seq2):
    ndiff=sum(1 for a, b in zip(seq1,seq2) if a != b and a!='-' and b!='-' and a!='N' and b!='N')
    nmc=sum(1 for a, b in zip(seq1,seq2) if a!='-' and b!='-' and a!='N' and b!='N')
    return ndiff,nmc

def read_in_close_clique_events():
    #want dict with: clique num string: {dict of cell pair: number of bp where they are relatively recombined}
    recodict={'4':{}}
    
    infile='../input_data/Close_Cliques_Manual_Event_Full_Data.csv'
    with open(infile,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            clique=row['Clique_Num']
            clusters=ast.literal_eval(row['Clusters'])
            length=float(row['Length'])
            add_event_to_recodict(recodict,clique,clusters,length)
    return recodict

def add_event_to_recodict(recodict,clique,clusters,length):
    if clique not in recodict:
        recodict[clique]={}
    #there are always 2 clusters here
    c1=clusters[0]
    c2=clusters[1]
    for i in range(len(c1)):
        celli=c1[i]
        for j in range(len(c2)):
            cellj=c2[j]
            cellpair=sorted([celli,cellj])
            cellpair=(cellpair[0],cellpair[1])
            if cellpair in recodict[clique]:
                recodict[clique][cellpair]+=length
            if cellpair not in recodict[clique]:
                recodict[clique][cellpair]=length


##############################################

def plot_results(chunklength=1000,gabriel_method=True):
    manual,excess=[],[]#lengths of reco for a pair from manual, vs poisson excess
    f_reco_fit,f_reco_excess=[],[]#1-f_asex from poisson fit, vs from poisson excess
    lambdas=[]#lambdas, to see if there is a correlation with fasex or excess
    mystr='_NoCovCut'
    outfile='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp_Gabriel_methods%s.csv' %(chunklength,mystr)
    if not gabriel_method:
         outfile='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp%s.csv' %(chunklength,mystr)
    with open(outfile,'r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            manual.append(float(row['Manual_Reco_Length']))
            lambdas.append(float(row['Lambda'])/chunklength)
            excess.append(float(row['Poisson_Excess_Length']))
            f_reco_fit.append(float(row['Freco_Fitting']))
            f_reco_excess.append(float(row['F_reco_excess']))
    plt.plot(manual,excess,'o',alpha=0.3)
    plot_R(manual,excess)
    plt.xlabel('Length Recombined from Manually Identified Events')
    plt.ylabel('Length in excess divergence of Poisson')
    plt.title('High Divergence Length from Poisson vs Manually-Found\nRecombination Events, for Close Cliques%s' %mystr)
    plt.show()

    plt.plot(lambdas,f_reco_excess,'o',alpha=0.3)
    plot_R(lambdas,f_reco_excess)
    plt.xlabel('Divergence from Fitting Poisson')
    plt.ylabel('Fraction of Genome in excess divergence of Poisson')
    plt.title('High Divergence Fraction of Genome from Poisson vs Lambda\nfor Close Cliques%s' %mystr)
    plt.show()

    plt.plot(f_reco_fit,f_reco_excess,'o',alpha=0.3)
    plot_R(f_reco_fit,f_reco_excess)
    plt.xlabel('Fraction Recombined from Fitting Poisson')
    plt.ylabel('Fraction of Genome in excess divergence of Poisson')
    plt.title('High Divergence Fraction of Genome from Poisson Fit vs Excess\nfor Close Cliques%s' %mystr)
    plt.show()
    
def plot_R(x,y):
    slope, intercept, r_value, p_value, std_err = linregress(x,y)
    myx=np.linspace(0,max(x),100)#(max(x),max(y)),100)
    plt.plot(myx,slope*myx+intercept,label='y=%.3fx+%.3f\nR=%f' %(slope,intercept,r_value))
    plt.legend()

################################

def compare_results(chunklength=1000):
    mystr=''
    lamgab,lamold,fgab,fold,exgab,exold=[],[],[],[],[],[]#lambda, freco using gabriels method an my old method
    outfilegab='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp_Gabriel_methods.csv' %(chunklength)#,mystr)
    outfileold='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp_Gabriel_methods.csv' %(300)#,mystr)
    #outfileold='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp.csv' %(chunklength)#,mystr)
    #outfileold='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp_2SNPs_and_up_excess.csv' %(chunklength)#,mystr)
    #outfileold='../input_data/Conservative_Event_Calling_Estimation_for_Close_Cliques_%ibp_3SNPs_and_up_excess.csv' %(chunklength)#,mystr)
    with open(outfilegab,'r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            lamgab.append(float(row['Lambda'])/chunklength)
            exgab.append(float(row['Poisson_Excess_Length']))
            fgab.append(float(row['F_reco_excess']))
    with open(outfileold,'r') as myoutfile:
        csvreader=csv.DictReader(myoutfile,delimiter='\t')
        for row in csvreader:
            lamold.append(float(row['Lambda'])/300)#chunklength)
            exold.append(float(row['Poisson_Excess_Length']))
            fold.append(float(row['F_reco_excess']))
    plt.plot(lamgab,lamold,'o',alpha=0.3)
    myx=np.linspace(0,max(lamgab),100)
    plt.plot(myx,myx)
    plt.xlabel('Divergence Using Method 1')
    plt.ylabel('Divergence Using Method 2')
    plt.title('Asexual Divergences Gauged using Method 1 vs 2 for Poisson\nfor 6 Close Cliques')
    plt.show()

    plt.plot(fgab,fold,'o',alpha=0.3)
    myx=np.linspace(0,max(fgab),100)
    plt.plot(myx,myx)
    plt.xlabel('F_reco Using Method 1')
    plt.ylabel('F_reco Using Method 2')
    plt.title('Fraction Recombined Gauged using Method 1 vs 2 for Poisson\nfor 6 Close Cliques')
    plt.show()

    plt.plot(exgab,exold,'o',alpha=0.3)
    myx=np.linspace(0,max(exgab),100)
    plt.plot(myx,myx)
    plt.xlabel('Length Recombined Using Method 1')
    plt.ylabel('Length Recombined Using Method 2')
    plt.title('Length Recombined Gauged using Method 1 vs 2 for Poisson\nfor 6 Close Cliques')
    plt.show()
###################################

def coverage_investigator(chunklength=1000):
    '''
    plot the coverage for chunks with 0,1,2,3,4-8,9-15,15+ snps
    '''
    mc=[[],[],[],[],[],[],[]]#0 snps, 1 snp,2,3,4-8,9-15,15+
    cliques=['1','2','3','4','5','6']
    celllist_list=[clique1,cl2,cl3,cl4,cl5,cl6]
    for i in range(len(cliques)):
        mc_one_clique(celllist_list[i],cliques[i],chunklength,mc)
    plt.hist(mc,label=['0 SNPs','1 SNP','2 SNPs','3 SNPs','4-8 SNPs','9-15 SNPs','>=15 SNPs'],bins=30,density=True,histtype='step')
    plt.legend()
    plt.xlabel('Number of bases (of %i) mutually covered for pair' %(chunklength))
    plt.ylabel('Number of %i-bp chunks (over all pairs in all cliques)' %chunklength)
    plt.title('Mutually Covered Length per %i bp binned by chunk divergence\nfor cell pairs in close cliques;\nIs there a natural scale to set the coverage cutoff?' %chunklength)
    plt.yscale('log')
    plt.show()

    plt.hist([mc[0]+mc[1]+mc[2]+mc[3],mc[4],mc[5]+mc[6]],label=['0-3 SNPs','4-8 SNPs','>=9 SNPs'],bins=50,density=True,histtype='step')
    plt.legend(loc=2)
    plt.xlabel('Number of bases (of %i) mutually covered for pair' %(chunklength))
    plt.ylabel('Number of %i-bp chunks (over all pairs in all cliques)' %chunklength)
    plt.title('Mutually Covered Length per %i bp binned by chunk divergence\nfor cell pairs in close cliques;\nIs there a natural scale to set the coverage cutoff?' %chunklength)
    plt.yscale('log')
    plt.show()

def mc_one_clique(celllist,clique,chunklength,mc):
    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
   
    for i in range(len(celllist)):
        for j in range(i):
            celli=celllist[i]
            cellj=celllist[j]
            seqi=chunk_seqdict[celli]
            seqj=chunk_seqdict[cellj]
            nchunks=len(seqi)
            for c in range(nchunks):
                si=seqi[c]
                sj=seqj[c]
                mydiff,mymc=get_div_chunk(si,sj)
                if mymc<15:# or mymc==chunklength:
                    continue
                if mydiff==0:
                    mc[0].append(mymc)
                if mydiff==1:
                    mc[1].append(mymc)
                if mydiff==2:
                    mc[2].append(mymc)
                if mydiff==3:
                    mc[3].append(mymc)
                if mydiff>=4 and mydiff<9:
                    mc[4].append(mymc)
                if mydiff>=9 and mydiff<15:
                    mc[5].append(mymc)
                if mydiff>=15:
                    mc[6].append(mymc)
#################################################


def compare_results_for_various_probability_cuts_and_violinplot(chunklength=1000):
    #want 3 plots:
    #1: on x axis, the probability cut. use f0 only for lambda.  on y axis, violin plots of the fraction recombined for all pairs (or box whisker) for that prob cut
    #2: lambda on x axis (same for all), and perc reco on y axis, with a different color for each prob cut (3-4 at a time)
    #3: x axis is length reco from events (same for all). y axis is length reco estimated using that prob cut, different color for all prob cuts
   
    
    pcuts=[.001,.01,.1]
    #pcuts=[0.01]

    
    cliques=['1','2','3','4','5','6']
    celllist_list=[clique1,cl2,cl3,cl4,cl5,cl6]
    
    
    #lambdas, eventlen are one entry per cell pair.  lenreco,freco are lists of length pcuts, with each list being in the same order as lambdas and eventlen, one entry per pair
    lambdas=[]#same for all
    freco=[[] for i in range(len(pcuts))]
    lenreco=[[] for i in range(len(pcuts))]
    eventlen=[]#from found events in recodict
    
    recodict=read_in_close_clique_events()
    for i in range(len(cliques)):
        prob_comparer_do_one_clique(celllist_list[i],cliques[i],chunklength,recodict,pcuts,lambdas,freco,lenreco,eventlen)

    #now plot
    prob_comparer_violin_plot(pcuts,freco)
    prob_comparer_plot_scatters(lambdas,freco,pcuts,'Lambda from Poisson','Fraction Recombined','Fraction Recombined for Close Cliques vs Asex Div\nEstimated for Various Probability Cutoffs')
    prob_comparer_plot_scatters(eventlen,lenreco,pcuts,'Recombined Length Estimated from Events','Recombined Length Estimated from Poisson','Recombined Lengths from Poisson vs Identified Events\nfor Close Clique Pairs')
    
def prob_comparer_violin_plot(pcuts,distributions):
    
    plotseps=list(range(1,len(pcuts)+1))
    
    plt.violinplot(distributions,positions=plotseps,showmedians=True)
    ax=plt.gca()
    ax.set_xticks(plotseps)
    ax.set_xticklabels(pcuts)
    plt.xlabel('Probability Cutoff Used for Poisson')
    plt.ylabel('Fraction of Genome Labeled Recombined')
    plt.title('Fraction of Genome Labeled Recombined for Close Clique Pairs\nvs the Probability Cutoff Used for Poisson')
    plt.show()

def prob_comparer_plot_scatters(xdist,ydists,pcuts,xlabel,ylabel,title):
    for i in range(len(pcuts)):
        plt.plot(xdist,ydists[i],'o',alpha=0.3,label='%s Cutoff' %str(pcuts[i]))
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def get_probs_per_kb(probs,nsnps,my_lambda):
    #for item in nsnps, want to get its probability under the poisson with lambda, then append to probs (which is for all cell pairs)
    for ns in nsnps:
        myprob=np.exp(-my_lambda)*(my_lambda**ns)/math.factorial(ns)
        probs.append(myprob)

def prob_comparer_do_one_clique(celllist,clique,chunklength,recodict,pcuts,lambdas,freco,lenreco,eventlen,returnprobs=False):
    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
    totmc,totmc_highdiv=[],[]
    probs=[]#foreach kb for each pair, the probability of that kb under the poisson
    for i in range(len(celllist)):
        for j in range(i):
            if clique!='C1':
                recoamount=get_recombined_length_for_cellpair_from_recodict(recodict,celllist,i,j,clique)#unchanged
                eventlen.append(recoamount)
            mylambda=prob_comparer_get_poisson_single_pair(i,j,celllist,chunk_seqdict,chunklength,clique,pcuts,lambdas,freco,lenreco,probs)
            if mylambda=='-':
                continue
            lambdas.append(mylambda/float(chunklength))
            if lambdas[-1]>0.004:
                print(celllist[i],celllist[j])
    if returnprobs:
        return probs
    
def prob_comparer_get_poisson_single_pair(i,j,celllist,chunk_seqdict,chunksize,clique,pcuts,lambdas,freco,lenreco,probs,returnwg=False):
    covfrac=80
    covcut=0.01*covfrac*chunksize
    celli=celllist[i]
    cellj=celllist[j]
    seqi=chunk_seqdict[celli]
    seqj=chunk_seqdict[cellj]
    nsnps,mc=[],[]
    uncovered=0
    mc_highdiv=[]
    #do for all kb with e.g. at least 80% of sites, and also only for completely covered kb
    #how many kb have fewer than 80% covered?
    nchunks=len(seqi)
    for c in range(nchunks):
        si=seqi[c]
        sj=seqj[c]
        mydiff,mymc=get_div_chunk(si,sj)
        if mymc<covcut:
            uncovered+=1
            continue
        nsnps.append(mydiff)
    #now get poisson model, for various pcuts
    for k in range(len(pcuts)):
        prob_cut=pcuts[k]
        thislambda,excess,thisfreco=prob_comparer_get_poisson_model_gabriel_method(nsnps,chunksize,prob_cut)
        if thislambda=='-':
            continue
        freco[k].append(thisfreco)
        lenreco[k].append(excess)
    get_probs_per_kb(probs,nsnps,thislambda)
    return thislambda

def prob_comparer_get_poisson_single_pair_for_flexgenes(i,j,celllist,chunk_seqdict,chunksize,pcut,returnwg=False):
    covfrac=80
    covcut=0.01*covfrac*chunksize
    celli=celllist[i]
    cellj=celllist[j]
    seqi=chunk_seqdict[celli]
    seqj=chunk_seqdict[cellj]
    nsnps,mc=[],[]
    uncovered=0
    mc_highdiv=[]
    #do for all kb with e.g. at least 80% of sites, and also only for completely covered kb
    #how many kb have fewer than 80% covered?
    nchunks=len(seqi)
    for c in range(nchunks):
        si=seqi[c]
        sj=seqj[c]
        mydiff,mymc=get_div_chunk(si,sj)
        if mymc<covcut:
            uncovered+=1
            continue
        nsnps.append(mydiff)
    #now get poisson model, for various pcuts
    if 1:
        prob_cut=pcut
        thislambda,excess,thisfreco=prob_comparer_get_poisson_model_gabriel_method(nsnps,chunksize,prob_cut)
        
        
       
    if not returnwg:
        return thislambda
    if returnwg:
        return thislambda,float(sum(nsnps))/(len(nsnps)*chunksize)


def prob_comparer_get_poisson_model_gabriel_method(nsnps,chunksize,probcut):
    f0=float(nsnps.count(0))/len(nsnps)
    if f0==0:
        return '-','-','-'
    my_lambda=-math.log(f0)
    mymax=max(nsnps)
    mypois=[]
    ncut=0
    for n in range(100):#mymax+1):
        poisval=np.exp(-my_lambda)*(my_lambda**n)/math.factorial(n)
        #print(n,poisval,my_lambda)
        mypois.append(poisval)
        if ncut==0:
            if poisval<=probcut:
                ncut=n
                break
    
    excess=0
    for item in nsnps:
        if item>=ncut:
            excess+=chunksize

    #12-27-23 edit
    expected_excess=get_poisson_expectation_above_cutoff_to_subtract(my_lambda,ncut,nsnps,chunksize)
    excess-=expected_excess
    ###
    
    freco=float(excess)/(chunksize*len(nsnps))

    if freco==1:
        print('%i chunks , lambda %f, prob %f, ncut %i' %(len(nsnps),my_lambda,probcut,ncut))
        t=input('t')
        print(nsnps)
        t=input('t')
    return my_lambda,excess,freco


########################################################3
def C1_official():
    freco=[]
    lambdas=[]
    myoutfile='../input_data/C1_freco_vs_lambda.csv'
    with open(myoutfile,'r') as outfile:
        csvreader=csv.DictReader(outfile,delimiter='\t')
        for row in csvreader:
            freco.append(float(row['Freco']))
            lambdas.append(float(row['Lambda']))

    plt.plot(lambdas,freco,'o',alpha=0.3) 
    ax=plt.gca()
    ax.set_xticks([0,0.001,0.002,0.003])
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.set_xticklabels(['0%','0.1%','0.2%','0.3%'])
    ax.set_yticks([0,0.05,0.1,0.15,0.2])
    ax.set_yticklabels(['0%','5%','10%','15%','20%'],fontsize=25)
    plt.xlabel('Asexual Divergence',fontsize=20)#,pad=40)
    plt.ylabel('Fraction Recombined',fontsize=20)#,pad=40)
    #plt.title('Fraction of Genome Recombined vs\nCell Pair Asexual Divergence',fontsize=25,pad=20)
   
    plt.tight_layout()
    plt.show()

    
def do_C1():
    pcuts=[0.01]
    chunklength=500
    groupname='C1'
    celllist=[item+'_C1' for item in C1_without_ribotype][:5]
    lambdas=[]#same for all
    freco=[[] for i in range(len(pcuts))]
    lenreco=[[] for i in range(len(pcuts))]
    eventlen=[]#from found events in recodict
    recodict={}
    probs=prob_comparer_do_one_clique(celllist,groupname,chunklength,recodict,pcuts,lambdas,freco,lenreco,eventlen,returnprobs=True)

    myoutfile='../input_data/C1_freco_vs_lambda____.csv'
    outfields=['Freco','Lambda']
    with open(myoutfile,'w') as outfile:
        csvreader=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
        csvreader.writeheader()
        for i in range(len(lambdas)):
            csvreader.writerow({'Freco':freco[0][i],'Lambda':lambdas[i]})

    print(lambdas)
    print(freco[0])
    t=input('t')
    for i in range(len(pcuts)):
        plt.plot(lambdas,freco[i],'o',alpha=0.3)
   
   
    ax=plt.gca()
    ax.set_xticks([0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035])
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xticklabels(['0%','0.05%','0.1%','0.15%','0.2%','0.25%','0.3%','0.35%'])
    ax.set_yticks([0,0.05,0.1,0.15,0.2])
    ax.set_yticklabels(['0%','5%','10%','15%','20%'])
    plt.xlabel('Asexual Divergence',fontsize=25)#,pad=40)
    plt.ylabel('Fraction Recombined',fontsize=25)#,pad=40)
    #plt.title('Fraction of Genome Recombined vs\nCell Pair Asexual Divergence',fontsize=25,pad=20)
   
    plt.tight_layout()
    plt.show()

    return

    plt.hist(probs,bins=100)
    plt.xlabel('Probability of SNPs for each kb under Poisson model')
    plt.ylabel('Number of kb x cell pairs')
    plt.title('Probability of each kb for each cell pair\nUnder Poisson model for each cell pair asexual div')
    plt.show()

    logbins=np.logspace(np.log10(min(probs)),np.log10(.1),100)
    plt.hist(probs,bins=logbins)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Probability of SNPs for each kb under Poisson model')
    plt.ylabel('Number of kb x cell pairs')
    plt.title('Probability of each kb for each cell pair\nUnder Poisson model for each cell pair asexual div')
    plt.show()

    logbins2=np.logspace(-15,np.log10(.1),100)
    plt.hist(probs,bins=logbins2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Probability of SNPs for each kb under Poisson model (down to 10^(-15))')
    plt.ylabel('Number of kb x cell pairs')
    plt.title('Probability of each kb for each cell pair\nUnder Poisson model for each cell pair asexual div')
    plt.show()


############################################
def C1_SNPs_official():
    freco=[]
    lambdas=[]
    myoutfile='../input_data/C1_freco_vs_lambda_with_SNPs.csv'
    with open(myoutfile,'r') as outfile:
        csvreader=csv.DictReader(outfile,delimiter='\t')
        for row in csvreader:
            freco.append(float(row['SNPsasex']))
            lambdas.append(float(row['SNPsreco']))

    plt.plot(freco,lambdas,'o',alpha=0.3) 
    ax=plt.gca()
    #ax.set_xticks([0,0.001,0.002,0.003])
    ax.tick_params(axis='both', which='major', labelsize=25)
    #ax.set_xticklabels(['0%','0.1%','0.2%','0.3%'])
    #ax.set_yticks([0,0.05,0.1,0.15,0.2])
    #ax.set_yticklabels(['0%','5%','10%','15%','20%'],fontsize=25)
    plt.xlabel('Number of SNPs in Asexual Backbone',fontsize=20)
    plt.ylabel('Number of SNPs\n from Recombination',fontsize=20)
   
    plt.tight_layout()
    plt.show()
    
def do_C1_frac_snps_reco():
    pcut=0.01
    chunklength=500
    groupname='C1'
    celllist=[item+'_C1' for item in C1_without_ribotype]
    lambdas=[]#same for all
    freco=[]
    lenreco=[]
    eventlen=[]#from found events in recodict
    recodict={}
    snpsasex,snpsreco=[],[]


    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
    totmc,totmc_highdiv=[],[]
    for i in range(len(celllist)):
        for j in range(i):
            if groupname!='C1':
                recoamount=get_recombined_length_for_cellpair_from_recodict(recodict,celllist,i,j,clique)#unchanged
                eventlen.append(recoamount)
            mylambda=single_pair_snps(i,j,celllist,chunk_seqdict,chunklength,clique,pcut,lambdas,freco,lenreco,snpsasex,snpsreco)
            if mylambda=='-':
                continue
            lambdas.append(mylambda/float(chunklength))
            if lambdas[-1]>0.004:
                print(celllist[i],celllist[j])
   
    
    #single_pair_snps(celllist,groupname,chunklength,recodict,pcut,lambdas,freco,lenreco,eventlen,snpsasex,snpsreco)
    
    myoutfile='../input_data/C1_freco_vs_lambda_with_SNPs.csv'
    outfields=['Freco','Lambda','SNPsasex','SNPsreco']
    with open(myoutfile,'w') as outfile:
        csvreader=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
        csvreader.writeheader()
        for i in range(len(lambdas)):
            csvreader.writerow({'Freco':freco[i],'Lambda':lambdas[i],'SNPsasex':snpsasex[i],'SNPsreco':snpsreco[i]})
            
    
      
    
    print(lambdas)
    print(freco)
    t=input('t')
    plt.plot(lambdas,freco,'o',alpha=0.3)
    ax=plt.gca()
    ax.set_xticks([0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035])
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xticklabels(['0%','0.05%','0.1%','0.15%','0.2%','0.25%','0.3%','0.35%'])
    ax.set_yticks([0,0.05,0.1,0.15,0.2])
    ax.set_yticklabels(['0%','5%','10%','15%','20%'])
    plt.xlabel('Asexual Divergence',fontsize=25)#,pad=40)
    plt.ylabel('Fraction Recombined',fontsize=25)#,pad=40)
    #plt.title('Fraction of Genome Recombined vs\nCell Pair Asexual Divergence',fontsize=25,pad=20)
   
    plt.tight_layout()
    plt.show()


    plt.plot(snpsasex,snpsreco,'o',alpha=.3)
    plt.xlabel('Number of SNPs in Asexual Backbone',fontsize=20)
    plt.ylabel('Number of SNPs\n from Recombination',fontsize=20)
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.tight_layout()
    plt.show()
    


def single_pair_snps(i,j,celllist,chunk_seqdict,chunksize,clique,pcut,lambdas,freco,lenreco,snpsasex,snpsreco):
    covfrac=80
    covcut=0.01*covfrac*chunksize
    celli=celllist[i]
    cellj=celllist[j]
    seqi=chunk_seqdict[celli]
    seqj=chunk_seqdict[cellj]
    nsnps,mc=[],[]
    uncovered=0
    mc_highdiv=[]
    #do for all kb with e.g. at least 80% of sites, and also only for completely covered kb
    #how many kb have fewer than 80% covered?
    nchunks=len(seqi)
    for c in range(nchunks):
        si=seqi[c]
        sj=seqj[c]
        mydiff,mymc=get_div_chunk(si,sj)
        if mymc<covcut:
            uncovered+=1
            continue
        nsnps.append(mydiff)
    #now get poisson model, for various pcuts
    
    prob_cut=pcut
    thislambda,excess,thisfreco,asex,reco=prob_comparer_get_poisson_snps(nsnps,chunksize,prob_cut)
    if thislambda=='-':
        return '-'
    freco.append(thisfreco)
    lenreco.append(excess)
    snpsasex.append(asex)
    snpsreco.append(reco)
    return thislambda

def prob_comparer_get_poisson_snps(nsnps,chunksize,prob_cut):
    asex,reco=0,0#number of snps from each
    f0=float(nsnps.count(0))/len(nsnps)
    if f0==0:
        return '-','-','-'
    my_lambda=-math.log(f0)
    mymax=max(nsnps)
    mypois=[]
    ncut=0
    for m in range(100):#mymax+1):
        poisval=np.exp(-my_lambda)*(my_lambda**m)/math.factorial(m)
        #print(n,poisval,my_lambda)
        mypois.append(poisval)
        if ncut==0:
            if poisval<=prob_cut:
                ncut=m
                break
    
    excess=0
    for item in nsnps:
        if item>=ncut:
            excess+=chunksize
            reco+=item
        if item<ncut:
            asex+=item
    #now change by amount poisson would expect
    for m in range(ncut,100):
        poisval=np.exp(-my_lambda)*(my_lambda**m)/math.factorial(m)
        if round(poisval*len(nsnps))>=1:
            asex+=m*round(poisval*len(nsnps))
            reco-=m*round(poisval*len(nsnps))

    #12-27-23 edit
    expected_excess=get_poisson_expectation_above_cutoff_to_subtract(my_lambda,ncut,nsnps,chunksize)
    excess-=expected_excess
    ###
    
    freco=float(excess)/(chunksize*len(nsnps))

    if freco==1:
        print('%i chunks , lambda %f, prob %f, ncut %i' %(len(nsnps),my_lambda,probcut,ncut))
        t=input('t')
        print(nsnps)
        t=input('t')
    return my_lambda,excess,freco,asex, reco
###################################################

def plot_single_pair_poisson_Official(cell1,cell2,chunklength):
    cell1=cell1+'_C1'
    cell2+='_C1'
    celllist=[cell1,cell2]
    chunk_seqdict=kashrepo.get_region_lists(celllist,False,chunklength)
    totmc,totmc_highdiv=[],[]
    clique=''
    freco,lambdas,lenreco,probs=[],[],[],[]
  

    covfrac=80
    covcut=0.01*covfrac*chunklength

    gabriel_method=True
    
    seqi=chunk_seqdict[cell1]
    seqj=chunk_seqdict[cell2]
    nsnps,mc=[],[]
    uncovered=0
    mc_highdiv=[]
    #do for all kb with e.g. at least 80% of sites, and also only for completely covered kb
    #how many kb have fewer than 80% covered?
    nchunks=len(seqi)
    for c in range(nchunks):
        si=seqi[c]
        sj=seqj[c]
        mydiff,mymc=get_div_chunk(si,sj)
        if mymc<covcut:
            uncovered+=1
            continue
        nsnps.append(mydiff)
        if mydiff>8:
            #print(celli,cellj,c*chunksize,(c+1)*chunksize,'%i snps %i mc' %(mydiff,mymc))
            #t=input('t')
            mc_highdiv.append(mymc)
        if mymc>7:
            mc.append(mymc)
    mymax=max(nsnps)
    mybins=[-.5]
    for i in range(mymax+1):
        mybins.append(i+.5)
    
    poisvals,fasex,mylambda,excess=get_poisson_model(nsnps,chunklength,gabriel_method)
   
    plt.hist(np.array([nsnps,list(range(mymax+1))],dtype=object),weights=np.array([[1 for i in range(len(nsnps))],[len(nsnps)*poisvals[j] for j in range(len(poisvals))]],dtype=object),label=['Data','Poisson (Lambda=%f)' %(mylambda/float(chunklength))],bins=mybins)
    plt.legend(fontsize=18)
    #if chunklength==1000:
        #plt.xlabel('SNPs per kb',fontsize=20)
        #plt.ylabel('Number of Mutually-Covered kb',fontsize=20)
        
    if chunklength!=1000:
        plt.xlabel('SNPs per %i bp' %chunklength)
        plt.ylabel('Number of %i-bp chunks at least %s %% Covered' %(chunklength,str(covfrac)))
    #plt.title('%s vs %s SNPs\nper %i-bp along C1 Composite Reference\n%i bp found recombined; %i bp excess of Poisson' %(cell1,cell2,chunklength,recoamount,excess))
    ax=plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    





