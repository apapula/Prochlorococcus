import csv
from collections import Counter



'''
THis file filters through the blast hits of each cell to the database (mit9301 ref genome here) and finds mutual best hits---not always the first hit listed by blast

Requires OutfileList.txt, a list of all the output blast files from step 1
'''


    
        
def filter_cell_blastfile(infile,outfile,qcovcutoff,pidcutoff):
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
            
                                 

                

##########################

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--prog',type=int,help='Program to run. 0: write SNPs on cluster; 1: plot SNPs locally',default=0)
    args = parser.parse_args()
    print(args)
    
    prog=args.prog
  

   
        
    if prog==0:
        with open('OutfileList.txt','r') as myinfile:
            outfiles=[line.rstrip() for line in myinfile if len(line)>0]

    
        for i in range(len(outfiles)):
            myf=outfiles[i]
            if myf[-4:]!='.out':
                continue
            newf='Filtered/%s' %myf
            myf='Blast_Output/'+myf
            filter_cell_blastfile(myf,newf,qcovcutoff,pidcutoff)
