import re
import statistics
import numpy as np
from os.path import exists

def load_single_gene_sequences_total(infile,celllist,tfna,gapsbool,verbosebool,returnAAlistbool=0):
    seqdict,lengthdict=load_single_gene_sequences_step1(infile,celllist,tfna)
   
    if len(seqdict.keys())==0:
        if returnAAlistbool:
            return seqdict,[],[]
        if not returnAAlistbool:
            return seqdict,[]
    seqdict_length=filter_seqdict_for_length(seqdict,lengthdict)
    #print('%i cells from load sequences ' %(len(list(seqdict_length.keys()))))##edit
    if len(list(seqdict_length.keys()))==0:
        if returnAAlistbool:
            return {},[],[]
        if not returnAAlistbool:
            return seqdict_length,[]
    seqdict_tfna,positionlist,AAlist=get_tfna_from_seqdict(seqdict_length,tfna,gapsbool,verbosebool)
    #print('%i cells from load sequences 2' %(len(list(seqdict_tfna.keys()))))##edit
    if returnAAlistbool:
        return seqdict_tfna,positionlist,AAlist
    if not returnAAlistbool:
        return seqdict_tfna,positionlist
    
def load_single_gene_sequences_step1(infile,celllist,tfna):
    #output: a dictionary with cell:sequence for all covered cells in celllist within a certain length cutoff.
    totseqdict={}#all nucleotides
    lengthdict={}
    if not exists(infile):
        print('File %s does not exist' %(infile))
        return {},{}
    with open(infile,'r') as myinfile:
        lines=myinfile.read()
    cellseqs=re.split('>',lines)
    seqlist=[]
    ngaps=[]
    for i in range(len(cellseqs)):
        sublist=re.split('\n',cellseqs[i])
        if not any(cell in sublist[0] for cell in celllist):
            continue
        #now find which cell is this
        mycell=''
        for cell in celllist:
            if cell in sublist[0]:
                mycell=cell
                break
        if mycell=='':
            print('Error: cell not found')
            t=input('t')
        seq=''.join(sublist[1:])
        if tfna in ['t','f','tf','n','a','atvscg']:
            seq=seq.replace("N","-")
        lengthdict[mycell]=len(seq)-seq.count('-')
        totseqdict[mycell]=seq
    if len(totseqdict)==0:
        return {},{}
    return totseqdict,lengthdict





def filter_seqdict_for_length(seqdict,lengthdict):#returns new seq dict with all seqs within .9 and 1.1 of median length
    lengthlist=[lengthdict[item] for item in list(lengthdict.keys())]
    if len(lengthlist)<2:
        return {}
    length_lower=0.9*statistics.median(lengthlist)
    length_upper=1.1*statistics.median(lengthlist)
    newseqdict={}
    for item in seqdict:
        if lengthdict[item]>=length_lower and lengthdict[item]<=length_upper:
            newseqdict[item]=seqdict[item]
    return newseqdict
    
    
def get_tfna_from_seqdict(seqdict,tfna,gapsbool,verbosebool):
    #take in seqdict. tfna is twofold, fourfold, nonsyn (second bp), or all. return proper thing. gapsbool is whether to include intermediate gaps for e.g. twofold sites (preserving the original site separation). it is set to True by default
    if tfna=='a':
        return seqdict,[],[]
    if tfna=='n':#just taking second sites
        newseqdict,positionlist=simple_nonsyns_seqs(seqdict,gapsbool)
        AAlist=[]
    if tfna=='t':
        newseqdict,positionlist,AAlist=twofold_seqs(seqdict,gapsbool,verbosebool)
    if tfna=='f':
        newseqdict,positionlist,AAlist=fourfold_seqs(seqdict,gapsbool,verbosebool)
    if tfna=='tf':
        newseqdict,positionlist,AAlist=tf_seqs(seqdict,gapsbool,verbosebool)
        #newseqdict4,positionlist4,AAlist4=fourfold_seqs(seqdict,gapsbool,verbosebool)
        #newseqdict={}
        #for item in newseqdict2:
            #newseqdict[item]=newseqdict2[item]+newseqdict4[item]
        #positionlist=positionlist2+positionlist4
        #AAlist=AAlist2+AAlist4
    if tfna=='aminoacid':
        return seqdict,[],[]
    if tfna=='atvscg':
        newseqdict,positionlist,AAlist=at_vs_gc_f_seqs(seqdict)
    if verbosebool:
        print('Twofold (t), Fourfold (f), Nonsyn second codons (n): %s' %tfna)
        print('Position List')
        print(positionlist)
    return newseqdict,positionlist,AAlist

def at_vs_gc_f_seqs(seqdict):
    gapsbool=1
    genelength=len(seqdict[list(seqdict.keys())[0]])
    if genelength%3!=0:
        if verbosebool:
            print('Length of Gene is not a multiple of 3 : %i '%genelength)
        #t=input('t')
        return {},[],[]
    ncodons=int(genelength/3)
    positionlist=[]
    AAlist=[]
    for c in range(ncodons):
        myAA=get_AA_if_twofold_fourfold(seqdict,c,'f')#changed from tf 8-22-24
        if myAA=='-':
            continue
        else:
            positionlist.append(c*3+2)
            AAlist.append(myAA)
    newseqdict={}
    for item in seqdict:
        myseq=seqdict[item]
        synseq=''
        for c in range(len(positionlist)):
            mycodon=myseq[positionlist[c]-2:positionlist[c]+1]
            if '-' in mycodon or 'N' in mycodon:
                synseq+='-'*(positionlist[c]-len(synseq)+1)
            if '-' not in mycodon and 'N' not in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq))
                synseq+=myseq[positionlist[c]]
        converted_seq=convert_myseq_to_at_vs_gc(synseq)
        newseqdict[item]=converted_seq
    return newseqdict,positionlist,AAlist


def convert_myseq_to_at_vs_gc(myseq):
    newseq=''
    for i in range(len(myseq)):
        b=myseq[i]
        if b not in ['A','T','C','G','a','t','c','g']:
            newseq+=b
            continue
        if b in ['A','a','T','t']:
            newb='A'
        if b in ['g','G','c','C']:
            newb='G'
        newseq+=newb
    return newseq
            
def tf_seqs(seqdict,gapsbool,verbosebool):
    #need to get length of alignment
    genelength=len(seqdict[list(seqdict.keys())[0]])
    if genelength%3!=0:
        if verbosebool:
            print('Length of Gene is not a multiple of 3 : %i '%genelength)
        #t=input('t')
        return {},[],[]
    ncodons=int(genelength/3)
    positionlist=[]
    AAlist=[]
    for c in range(ncodons):
        myAA=get_AA_if_twofold_fourfold(seqdict,c,'tf')
        if myAA=='-':
            continue
        else:
            positionlist.append(c*3+2)
            AAlist.append(myAA)
    ###edi
    #print('genelength ',genelength)
    #print('positions',positionlist)
    #t=input('t')

    
    newseqdict={}
    for item in seqdict:
        myseq=seqdict[item]
        synseq=''
        for c in range(len(positionlist)):
            mycodon=myseq[positionlist[c]-2:positionlist[c]+1]
            if '-' in mycodon or 'N' in mycodon:
                synseq+='-'*(positionlist[c]-len(synseq)+1)
            if '-' not in mycodon and 'N' not in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq))
                synseq+=myseq[positionlist[c]]
        newseqdict[item]=synseq
    return newseqdict,positionlist,AAlist

def fourfold_seqs(seqdict,gapsbool,verbosebool):
    #need to get length of alignment
    genelength=len(seqdict[list(seqdict.keys())[0]])
    if genelength%3!=0:
        if verbosebool:
            print('Length of Gene is not a multiple of 3 : %i '%genelength)
        #t=input('t')
        return {},[],[]
    ncodons=int(genelength/3)
    positionlist=[]
    AAlist=[]
    for c in range(ncodons):
        myAA=get_AA_if_twofold_fourfold(seqdict,c,'f')
        if myAA=='-':
            continue
        else:
            positionlist.append(c*3+2)
            AAlist.append(myAA)
    ###edi
    #print('genelength ',genelength)
    #print('positions',positionlist)
    #t=input('t')

    
    newseqdict={}
    for item in seqdict:
        myseq=seqdict[item]
        synseq=''
        for c in range(len(positionlist)):
            mycodon=myseq[positionlist[c]-2:positionlist[c]+1]
            if '-' in mycodon or 'N' in mycodon:
                synseq+='-'*(positionlist[c]-len(synseq)+1)
            if '-' not in mycodon and 'N' not in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq))
                synseq+=myseq[positionlist[c]]
        newseqdict[item]=synseq
    return newseqdict,positionlist,AAlist






def twofold_seqs(seqdict,gapsbool,verbosebool):
    #need to get length of alignment
    genelength=len(seqdict[list(seqdict.keys())[0]])
    if genelength%3!=0:
        if verbosebool:
            print('Length of Gene is not a multiple of 3 : %i '%genelength)
        #t=input('t')
        return {},[],[]
    ncodons=int(genelength/3)
    positionlist=[]
    AAlist=[]
    for c in range(ncodons):
        myAA=get_AA_if_twofold_fourfold(seqdict,c,'t')
        if myAA=='-':
            continue
        else:
            AAlist.append(myAA)
            positionlist.append(c*3+2)
    newseqdict={}
    for item in seqdict:
        myseq=seqdict[item]
        #print('original %i' %len(myseq))
        synseq=''
        for c in range(len(positionlist)):
            mycodon=myseq[positionlist[c]-2:positionlist[c]+1]
            if '-' in mycodon or 'N' in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq)+1)
                if not gapsbool:
                    synseq+='-'
                #print('bad')
            if '-' not in mycodon and 'N' not in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq))
                synseq+=myseq[positionlist[c]]
        newseqdict[item]=synseq
        #print(len(synseq))
        #t=input('t')
    return newseqdict,positionlist,AAlist

def fourfold_seqs_old_with_6fold(seqdict,gapsbool):
    #need to get length of alignment
    genelength=len(seqdict[list(seqdict.keys())[0]])
    if genelength%3!=0:
        print('Length of Gene is not a multiple of 3 : %i '%genelength)
        t=input('t')
        return {},[],[]
    ncodons=int(genelength/3)
    positionlist=[]
    AAlist=[]
    for c in range(ncodons):
        #mycod=get_consensus_codon(seqdict,c)
        #if mycod[:2] in fourfoldlist:
            #positionlist.append(c*3+2)
        mycod=get_if_codons_are_all_synonymous_and_return_consensus_codon(seqdict,c,'f')
        if mycod=='-':
            continue
        if mycod[:2] in fourfoldlist:
            positionlist.append(c*3+2)
            AAlist.append(codon_dict[mycod])
            
    newseqdict={}
    for item in seqdict:
        myseq=seqdict[item]
        synseq=''
        for c in range(len(positionlist)):
            mycodon=myseq[positionlist[c]-2:positionlist[c]+1]
            if '-' in mycodon or 'N' in mycodon:
                synseq+='-'*(positionlist[c]-len(synseq)+1)
            if '-' not in mycodon and 'N' not in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq))
                synseq+=myseq[positionlist[c]]
        newseqdict[item]=synseq
    return newseqdict,positionlist,AAlist


def twofold_seqs_old_with_6fold(seqdict,gapsbool):
    #need to get length of alignment
    genelength=len(seqdict[list(seqdict.keys())[0]])
    if genelength%3!=0:
        print('Length of Gene is not a multiple of 3 : %i '%genelength)
        t=input('t')
        return {},[],[]
    ncodons=int(genelength/3)
    positionlist=[]
    AAlist=[]
    for c in range(ncodons):
        #mycod=get_consensus_codon(seqdict,c)
        mycod=get_if_codons_are_all_synonymous_and_return_consensus_codon(seqdict,c,'t')
        if mycod=='-':
            continue
        if mycod in twofoldlist:
            AAlist.append(codon_dict[mycod])
            positionlist.append(c*3+2)
    newseqdict={}
    for item in seqdict:
        myseq=seqdict[item]
        #print('original %i' %len(myseq))
        synseq=''
        for c in range(len(positionlist)):
            mycodon=myseq[positionlist[c]-2:positionlist[c]+1]
            if '-' in mycodon or 'N' in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq)+1)
                if not gapsbool:
                    synseq+='-'
                #print('bad')
            if '-' not in mycodon and 'N' not in mycodon:
                if gapsbool:
                    synseq+='-'*(positionlist[c]-len(synseq))
                synseq+=myseq[positionlist[c]]
        newseqdict[item]=synseq
        #print(len(synseq))
        #t=input('t')
    return newseqdict,positionlist,AAlist

def simple_nonsyns_seqs(seqdict,gapsbool):
    nonsyndict={}
    positionlist=[]
    for item in seqdict:
        nseq=''
        oldseq=seqdict[item]
        for i in range(int(len(oldseq)/3)):
            #this does each codon, not each bp
            positionlist.append(i*3+1)
            if gapsbool:
                nseq+='-'+oldseq[i*3+1]+'-'
            if not gapsbool:
                nseq+=oldseq[i*3+1]
        nonsyndict[item]=nseq
    return nonsyndict,positionlist

def get_AA_if_twofold_fourfold(seqdict,codon_num,tfna):
    codons=[]
    AAs=[]
    for item in seqdict:
        mycodon=seqdict[item][codon_num*3:(codon_num+1)*3]
        if '-' in mycodon or 'N' in mycodon:
            continue
        codons.append(seqdict[item][codon_num*3:(codon_num+1)*3])
        AAs.append(codon_dict[mycodon])

    
    if len(codons)==0:
        return '-'
    if len(AAs)!=AAs.count(AAs[0]):
        return '-'
    AA=AAs[0]
    if tfna=='t':
        if AA in twofoldAA:
            return AA
        else:
            return '-'
    if tfna=='f':
        if AA in fourfoldAA:
            return AA
        else:
            return '-'
    if tfna=='tf':
        if AA in twofoldAA or AA in fourfoldAA:
            return AA
        else:
            return '-'

def get_if_codons_are_all_synonymous_and_return_consensus_codon(seqdict,codon_num,tfna):
    codons=[]
    for item in seqdict:
        mycodon=seqdict[item][codon_num*3:(codon_num+1)*3]
        if '-' in mycodon or 'N' in mycodon:
            continue
        codons.append(seqdict[item][codon_num*3:(codon_num+1)*3])
    if len(codons)==0:
        return '-'
    #now check if all synonymous
    if tfna=='f':
        first_two_letters=[cod[:2] for cod in codons]
        if first_two_letters.count(first_two_letters[0])==len(first_two_letters):
            #if they're fourfold degen, they're all the same AA. we'll check if they're fourfold degen in the fourfold_seqs program
            return most_frequent(codons)
        else:
            return '-'
    if tfna=='t':
        first_two_letters=[cod[:2] for cod in codons]
        if first_two_letters.count(first_two_letters[0])!=len(first_two_letters):
            return '-'
        #first, check if first two letters are the same. note that things can happen where most frequent is leucine TTA but things like CTC are also leucine, and then have AC sites
        
        AAs=[]
        for item in codons:
            AAs.append(codon_dict[item])
        if len(AAs)==AAs.count(AAs[0]):
            return most_frequent(codons)
        else:
            return '-'
        

    
def get_consensus_codon(seqdict,codon_num):
    codons=[]
    for item in seqdict:
        codons.append(seqdict[item][codon_num*3:(codon_num+1)*3])
    #get majority codon
    #then will want to check if twofold or fourfold
    if 0:
        print(codons)
        print('most common of above: ',most_frequent(codons))
        t=input('t')
    return most_frequent(codons)

def most_frequent(List):
    unique, counts = np.unique(List, return_counts=True)
    index = np.argmax(counts)
    return unique[index]






fourfoldlist=['CT','GT','TC','CC','GC','AC','GG','CG']
twofold_list_ct=['TTT','TTC','TAC','TAT','CAT','CAC','AAT','AAC','GAT','GAC','AGT','AGC','TGT','TGC']
twofold_list_ag=['TTA','TTG','CAA','CAG','AAA','AAG','GAA','GAG','AGA','AGG']
twofoldlist=twofold_list_ct+twofold_list_ag
#all fourfold degenerate codons have same first two np
#for twofold degenerate codons that are twofold at the third site (that is, excluding L TTG->CTG, TTA->CTA and R CGA->AGA, CGG->AGG, : they have the same first two sites.

AA_dict={'F':['TTT','TTC'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],'I':['ATT','ATC','ATA'],'M':['ATG'],'V':['GTT','GTA','GTG','GTC'],'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'P':['CCT','CCC','CCA','CCG'],'T':['ACT','ACC','ACA','ACG'],'A':['GCT','GCC','GCA','GCG'],'Y':['TAT','TAC'],'STOP':['TAA','TAG','TGA'],'H':['CAT','CAC'],'Q':['CAA','CAG'],'N':['AAT','AAC'],'K':['AAA','AAG'],'D':['GAT','GAC'],'E':['GAA','GAG'],'C':['TGT','TGC'],'W':['TGG'],'R':['CGT','CGC','CGA','CGG','AGG','AGA'],'G':['GGT','GGC','GGA','GGG']}
codon_dict={}
for item in AA_dict:
    mylist=AA_dict[item]
    for c in mylist:
        codon_dict[c]=item
        
twofoldAA=['F','Y','H','Q','N','K','D','E','C']#excluding 6fold R, L, A
fourfoldAA=['V','P','T','A','G']#excluding 6fold R, L, S
