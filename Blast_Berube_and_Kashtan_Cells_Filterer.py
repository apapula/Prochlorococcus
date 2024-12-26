import csv
from collections import Counter


#objectives: MIT9301 gene: cell gene.  each 9301 gene can only have best match.  each cell gene can only have best match. require mutual best

def execute():
    pidcutoff=40.
    qcovcutoff=55.
    #with open('OutfileList.txt','r') as myinfile:
    with open('Berube_Blast_Output_List.txt','r') as myinfile:
        outfiles=[line.rstrip() for line in myinfile if len(line)>0]

    
    for i in range(len(outfiles)):
        myf=outfiles[i]
        if myf[-4:]!='.out':
            continue
        newf='New_New_Filters/%s' %myf#Filtered
        #myf='Blast_Output/'+myf
        myf='Old_Filtered_with_HLII/'+myf
        filter_cell_blastfile(myf,newf,qcovcutoff,pidcutoff)


def filter_cell_blastfile(infile,outfile,qcovcutoff,pidcutoff):
    
   
    cellhitdict={}#cell hit:[ref gene,score]
    refhitdict={}#ref gene: [celgene,score]# score is pid*qcov

    mybest=[]#for a given reference/query gene: cell gene, pid, qcov, length
    refgene=''#when it changes, have done all the hits for that reference/query gene
    print(infile)##
    
    with open(infile,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            
            q=row['qseqid']
            s=row['sseqid']
            pid=float(row['pident'])
            qcov=float(row['qcovhsp'])
            score=pid*qcov
            if q not in refhitdict and s not in cellhitdict:
                refhitdict[q]=[s,score,pid,qcov]
                cellhitdict[s]=[q,score]
                continue
            if q in refhitdict and s not in cellhitdict:
                oldscore=refhitdict[q][1]
                olds=refhitdict[q][0]
                if score>oldscore:
                    refhitdict[q]=[s,score,pid,qcov]
                    del cellhitdict[olds]
                    cellhitdict[s]=[q,score]
                continue
            if q not in refhitdict and s in cellhitdict:
                oldscore=cellhitdict[s][1]
                oldq=cellhitdict[s][0]
                if score>oldscore:
                    refhitdict[q]=[s,score,pid,qcov]
                    cellhitdict[s]=[q,score]
                    del refhitdict[oldq]
                continue
            if q in refhitdict and s in cellhitdict:
                oldscoreq=refhitdict[q][1]
                olds=refhitdict[q][0]
                oldscores=cellhitdict[s][1]
                oldq=cellhitdict[s][0]
                if score>oldscores and score>oldscoreq:
                    refhitdict[q]=[s,score,pid,qcov]
                    del cellhitdict[olds]
                    cellhitdict[s]=[q,score]
                    del refhitdict[oldq]
                continue
                    

   
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=['qseqid','sseqid','pident','qcovhsp'])
        outwriter.writeheader()
        reflist=sorted(list(refhitdict.keys()))
        for i in range(len(reflist)):
            r=reflist[i]
            rl=refhitdict[r]
            outwriter.writerow({'qseqid':r,'sseqid':rl[0],'pident':rl[2],'qcovhsp':rl[3]})
            
                    
    
def old_filter_cell_blastfile(infile,outfile,qcovcutoff,pidcutoff):
    #here, mit9301 (reference) is query and cell is subject
    #gather all the hits to a query (ref) gene. get best. keep a dictionary: cell (subject) hit: reference (query) hit.  AND reference hit: [subject hit, pid, qcov, length]
    cellhitdict={}#cell hit:ref gene
    refhitdict={}#ref gene: [cell gene, pid, qcov, length]

    mybest=[]#for a given reference/query gene: cell gene, pid, qcov, length
    refgene=''#when it changes, have done all the hits for that reference/query gene
    print(infile)##
    
    with open(infile,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            
            q=row['qseqid']
           
            if q!=refgene and len(mybest)>0:
                #pick best hit
                #one could do product of qcov and pid
                hitindex=-1
                score=0
                for i in range(len(mybest)):
                    myl=mybest[i]
                    myscore=myl[1]*myl[2]
                    if myscore>score:
                        hitindex=i
                bestl=mybest[hitindex]
                if bestl[0] in cellhitdict:
                    #have a duplicate hit in cell. need to find the better of the hits.
                    othercellhitlist=refhitdict[cellhitdict[bestl[0]]]
                    otherscore=othercellhitlist[1]*othercellhitlist[2]
                    ##if otherscore>score:
                        #old hit wins
                        ##continue
                    if otherscore<=score:
                        #new hit
                        cellhitdict[bestl[0]]=refgene
                        refhitdict[refgene]=bestl
                if bestl[0] not in cellhitdict:
                    cellhitdict[bestl[0]]=refgene
                    refhitdict[refgene]=bestl
                
                
                refgene=q
                mybest=[]
            if q!=refgene and len(mybest)==0:
                refgene=q
            if q==refgene:
                #add entry to besthit if qcov and pid high enough
                pid=float(row['pident'])
                if pid<pidcutoff:
                    continue
                qcov=float(row['qcovhsp'])
                if qcov<qcovcutoff:
                    continue
                s=row['sseqid']
                leng=float(row['length'])
                mybest.append([s,pid,qcov,leng])

    #do for final gene
    if len(mybest)>0:
        hitindex=-1
        score=0
        for i in range(len(mybest)):
            myl=mybest[i]
            myscore=myl[1]*myl[2]
            if myscore>score:
                hitindex=i
        bestl=mybest[hitindex]
        if bestl[0] in cellhitdict:
            #have a duplicate hit in cell. need to find the better of the hits.
            othercellhitlist=refhitdict[cellhitdict[bestl[0]]]
            otherscore=othercellhitlist[1]*othercellhitlist[2]
            ##if otherscore>score:
                #old hit wins
                ##continue
            if otherscore<=score:
                #new hit
                cellhitdict[bestl[0]]=refgene
                refhitdict[refgene]=bestl
        if bestl[0] not in cellhitdict:
            cellhitdict[bestl[0]]=refgene
            refhitdict[refgene]=bestl
                
    #now, write results
   
    with open(outfile,'w') as myoutfile:
        outwriter=csv.DictWriter(myoutfile,delimiter='\t',fieldnames=['qseqid','sseqid','pident','qcovhsp','length'])
        outwriter.writeheader()
        reflist=sorted(list(refhitdict.keys()))
        for i in range(len(reflist)):
            r=reflist[i]
            rl=refhitdict[r]
            outwriter.writerow({'qseqid':r,'sseqid':rl[0],'pident':rl[1],'qcovhsp':rl[2],'length':rl[3]})
            
                                 

                
###########
execute()
