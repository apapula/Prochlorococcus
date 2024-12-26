import matplotlib.pyplot as plt
import csv
import math
import numpy as np
import statistics

import pandas as pd
from shapely.geometry import Point

import geopandas as gpd
from geopandas import GeoDataFrame



latlongdict={'311':(-20.08,-70.8),'315':(-20.08,-71),'316':(-20.08,-71),'321':(-23.46,-88.77),'323':(-23.46,-88.77),'331':(-23.46,-89),'335':(-26.25,-103.96),'341':(-26.24,-104),'345':(-26.24,-104),'347':(23.75,-158),'402':(23.75,-158),'355':(31,-64),'363':(31,-65),'388':(36,-53.3),'412':(24.7,-67),'409':(9.55,-50.5),'418':(24,-22),'420':(24,-22),'424':(38.3,-69),'429':(26.14,-44.8),'432':(22.33,-36),'436':(17.4,-24.5),'442':(-35,16),'444':(-35,16),'449':(-30,156),'450':(-30,156),'455':(-30,174),'459':(-32.5,-170),'463':(-32.5,-170),'469':(-32.5,-154),'670':(28.14,-158),'673':(29,-158),'676':(32.7,-158),'679':(32,-158),'683':(36.57,-158),'686':(36,-158)}

depthdict={'311':20,'315':55,'321':14,'331':112,'335':14,'341':180,'347':5,'402':100,'355':10,'363':100,'388':8,'412':119,'409':100,'418':58,'424':91,'429':90,'432':99.7,'436':71.6,'442':21.2,'449':50.6,'455':50.4,'459':76,'463':203,'469':100,'670':5,'673':90,'676':5,'679':60,'683':5,'686':65}#meters

ocean_ranges_lat={'North Atlantic':[-1,68.6],'South Atlantic':[-60,0],'Southern':[-86,-60],'Mediterranean':[30,48],'North Pacific':[0,67],'South Pacific':[-60,4],'Indian':[-60,32],'Red Sea':[12,28]}
ocean_ranges_long={'North Atlantic':[-98,12],'South Atlantic':[-70,20],'Southern':[-180,180],'Mediterranean':[-6,43],'North Pacific':[117,-77],'South Pacific':[130,-67],'Indian':[20,147],'Red Sea':[33,43]}


def return_ocean(samplenumber):
    if samplenumber not in latlongdict:
        return '-'
    ll=latlongdict[samplenumber]
    lat=ll[0]
    longitude=ll[1]
    for ocean in ocean_ranges_lat:
        if cyclic_range(ocean_ranges_lat[ocean],lat):
            if cyclic_range(ocean_ranges_long[ocean],longitude):
                return ocean

 
def plot_cell_group_on_map(cellgrouplist):#plots the samply locations of each cell in all the groups in grouplist as different color scatter plots on map
    #make dataframe for each group

    world=gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    for cellgroup in cellgrouplist:
        lats=[]
        longs=[]
        for item in cellgroup:
            sample=item[3:6]
            lats.append(latlongdict[sample][0])
            longs.append(latlongdict[sample][1])
            df=pd.DataFrame(list(zip(lats,longs)),columns=['Latitude','Longitude'])
            geometry=[Point(xy) for xy in zip(df['Longitude'],df['Latitude'])]
            gdf=GeoDataFrame(df,geometry=geometry)
            #print(gdf)
            #gdf.plot(ax=world,marker='o',markersize=15)
            gdf.plot(ax=world.boundary.plot(figsize=(10,6),color='k'),marker='o',markersize=15)
    plt.show()

def plot_cell_group_on_map2(cellgrouplist):#plots the samply locations of each cell in all the groups in grouplist as different color scatter plots on map
    #make dataframe for each group
    #fig,ax=plt.subplots()
    #ax.set_aspect('equal')
    #world=gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    #world.plot(ax=ax,color='k')
    
    for cellgroup in cellgrouplist:
        lats=[]
        fig,ax=plt.subplots()
        ax.set_aspect('equal')
        world=gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        world.plot(ax=ax,color='k')
        longs=[]
        for item in cellgroup:
            sample=item[3:6]
            lats.append(latlongdict[sample][0])
            longs.append(latlongdict[sample][1])
        df=pd.DataFrame(list(zip(lats,longs)),columns=['Latitude','Longitude'])
        geometry=[Point(xy) for xy in zip(df['Longitude'],df['Latitude'])]
        gdf=GeoDataFrame(df,geometry=geometry)
        gdf.plot(ax=ax,marker='o',label='%i cells' %len(lats))
            
        plt.title('Berube HLII Close Cliques\n%s' %str(cellgroup))
        plt.legend()
        plt.show()
    
berube_close_cliques=[['AG-412-A14', 'AG-449-D16'],['AG-449-D22', 'AG-418-O02', 'AG-442-N07', 'AG-469-M13', 'AG-449-C14'],['AG-449-D22', 'AG-418-O02', 'AG-442-N07', 'AG-418-P13', 'AG-455-E04'], ['AG-449-D22', 'AG-418-J17', 'AG-469-M13', 'AG-449-C14'], ['AG-418-L19', 'AG-670-M18'], ['AG-459-J14', 'AG-469-M13', 'AG-418-J17'], ['AG-459-J14', 'AG-429-E20'], ['AG-424-G03', 'AG-424-A14', 'AG-418-J19', 'AG-432-K16'], ['AG-429-A02', 'AG-442-B03'], ['AG-412-C21', 'AG-418-F16'], ['AG-418-P06', 'AG-449-C14'], ['AG-424-J22', 'AG-429-E20'], ['AG-429-C19', 'AG-455-E15'], ['AG-429-C19', 'AG-388-F11'], ['AG-429-C19', 'AG-418-F16']]

mylist=['AG-424-E18','AG-355-J21']#BATS and GFST westerlies gulf stream
plot_cell_group_on_map([mylist])#berube_close_cliques[0])
plot_cell_group_on_map2([mylist])#berube_close_cliques)
    
def cyclic_range(mytuple, samplevalue):#for latitude and longitudes, e.g. North Pac longitude range
    lower=mytuple[0]
    upper=mytuple[1]
    if lower<upper:
        if samplevalue>=lower and samplevalue<=upper:
            return True
    if upper<lower:
        if samplevalue>lower and samplevalue<180:
            return True
        if samplevalue<upper and samplevalue>-180:
            return True
    else:
        return False
    
#https://www.marineregions.org/gazetteer.php?p=details&id=1912&from=rss
def plotlatlong(l1,l2,l3,l4,gnum):
    lat1=[]
    long1=[]
    lat2=[]
    long2=[]
    lat3=[]
    long3=[]
    lat4=[]
    long4=[]
    for item in l1:
        num=item[3:6]
        ll=latlongdict[num]
        lat1.append(ll[0])
        long1.append(ll[1])
    for item in l2:
        num=item[3:6]
        ll=latlongdict[num]
        lat2.append(ll[0])
        long2.append(ll[1])
    for item in l3:
        num=item[3:6]
        ll=latlongdict[num]
        lat3.append(ll[0])
        long3.append(ll[1])
    for item in l4:
        num=item[3:6]
        ll=latlongdict[num]
        lat4.append(ll[0])
        long4.append(ll[1])
    jitteredlat1=lat1+5*np.random.rand(len(lat1))-2.5
    jitteredlong1=long1+5*np.random.rand(len(long1))-2.5
    jitteredlat2=lat2+5*np.random.rand(len(lat2))-2.5
    jitteredlong2=long2+5*np.random.rand(len(long2))-2.5

    jitteredlat3=lat3+5*np.random.rand(len(lat3))-2.5
    jitteredlong3=long3+5*np.random.rand(len(long3))-2.5
    jitteredlat4=lat4+5*np.random.rand(len(lat4))-2.5
    jitteredlong4=long4+5*np.random.rand(len(long4))-2.5
   
    plt.scatter(jitteredlat1,jitteredlong1,alpha=0.3,color='r',label='Cluster 1')
    plt.scatter(jitteredlat2,jitteredlong2,alpha=0.3,color='b',label='Cluster 2')
    plt.scatter(jitteredlat3,jitteredlong3,alpha=0.3,color='g',label='Cluster 3')
    plt.scatter(jitteredlat4,jitteredlong4,alpha=0.3,color='orange',label='Cluster 4')
    plt.xlabel('Latitude')
    plt.ylabel('Longitude')
    plt.legend()
    plt.title('Lat/Long Positions of Cells in the Clusters of Gene %s' %str(gnum))
    plt.show()




def plotdepth(l1,l2,l3,l4,lab1,lab2,lab3,lab4):
    depth1=[]
    depth2=[]
    depth3=[]
    depth4=[]
    for item in l1:
        num=item[3:6]
        d=depthdict[num]
        depth1.append(d)
    for item in l2:
        num=item[3:6]
        d=depthdict[num]
        depth2.append(d)
    for item in l3:
        num=item[3:6]
        d=depthdict[num]
        depth3.append(d)
    for item in l4:
        num=item[3:6]
        d=depthdict[num]
        depth4.append(d)
    jitteredd1=depth1+5*np.random.rand(len(depth1))-2.5
    jitteredd2=depth2+5*np.random.rand(len(depth2))-2.5
    jitteredd3=depth3+5*np.random.rand(len(depth3))-2.5
    jitteredd4=depth4+5*np.random.rand(len(depth4))-2.5

    
    if len(depth3)>0 and len(depth4)==0:
        x=[depth1,depth2,depth3]
        labs=['%s (%i cells)'  %(lab1,len(l1)), '%s (%i cells)' %(lab2, len(l2)),'%s (%i cells)' %(lab3,len(l3))]
    if len(depth4)==0 and len(depth3)==0:
        x=[depth1,depth2]
        labs=['%s (%i cells)' %(lab1,len(l1)), '%s (%i cells)' %(lab2,len(l2))]
   
    
  
    plt.hist(x,label=labs,bins=20,density=True)
    #plt.hist(depth1,bins=40,label='%s' %lab1)
    #plt.hist(depth2,bins=40,label='%s' %lab2)
    #if len(depth3)>0:
        #plt.hist(depth3,bins=40,label='%s' %lab3)
    #if len(depth4)>0:
        #plt.hist(depth4,bins=40,label='%s' %lab4)
    
    plt.legend()
    plt.title('Depth (m) of Cells (Normalized)')
    plt.show()

    plt.hist(x,label=labs,bins=20)
    #plt.hist(depth1,bins=40,label='%s' %lab1)
    #plt.hist(depth2,bins=40,label='%s' %lab2)
    #if len(depth3)>0:
        #plt.hist(depth3,bins=40,label='%s' %lab3)
    #if len(depth4)>0:
        #plt.hist(depth4,bins=40,label='%s' %lab4)
    
    plt.legend()
    plt.title('Depth (m) of Cells')
    plt.show()

def wgdiv_vs_location(infilewgdiv,celllist):
    #open the file.  get all WG distances and categorize by same or different ocean.  can bin by oceans (w.g. north vs south atlantic, then atlantic vs pacific...get closest pair in each ocean, and between each pair of oceans.  plot somehow...
    celloceans={}#cell:ocean
    oceancount={}#ocean:count of cells
    for cell in celllist:
        num=cell[3:6]
        ll=latlongdict[num]
        lat=ll[0]
        longitude=ll[1]
        foundbool=False
        for ocean in ocean_ranges_lat:
            if cyclic_range(ocean_ranges_lat[ocean],lat):
                if cyclic_range(ocean_ranges_long[ocean],longitude):
                    foundbool=True
                    celloceans[cell]=ocean
                    if ocean in oceancount:
                        oceancount[ocean]+=1
                    if ocean not in oceancount:
                        oceancount[ocean]=1
                    break
    print(oceancount)
    #print(celloceans)
    
    ocean_means={}#[ocean,ocean]:[means]
    ocean_medians={}
    ocean_means_pairs={}
    ocean_medians_pairs={}
    with open(infilewgdiv,'r') as myinfile:
        csvreader=csv.DictReader(myinfile,delimiter='\t')
        for row in csvreader:
            c1=row['Cell1']
            c2=row['Cell2']
            if c1 not in celllist or c2 not in celllist:
                continue
            mean=float(row['Mean'])
            median=float(row['Median'])
            oc1=celloceans[c1]
            oc2=celloceans[c2]
            ocpair=sorted([oc1,oc2])
            ocpair=(ocpair[0],ocpair[1])
            if ocpair in ocean_means:
                ocean_means[ocpair].append(mean)
                ocean_medians[ocpair].append(median)
                ocean_means_pairs[ocpair].append((c1,c2))
                ocean_medians_pairs[ocpair].append((c1,c2))
            if ocpair not in ocean_means:
                ocean_means[ocpair]=[mean]
                ocean_medians[ocpair]=[median]
                ocean_means_pairs[ocpair]=[(c1,c2)]
                ocean_medians_pairs[ocpair]=[(c1,c2)]
    for item in ocean_means:
        #print('%s Closest Pair Mean: %f ; Median %f ; Furthest Pair Mean %f; Median %f' %(item,min(ocean_means[item]),min(ocean_medians[item]),max(ocean_means[item]),max(ocean_medians[item])))
        print('%s Closest Pair Mean: %f %s ; Median %f %s' %(item,min(ocean_means[item]),ocean_means_pairs[item][ocean_means[item].index(min(ocean_means[item]))],min(ocean_medians[item]),ocean_medians_pairs[item][ocean_medians[item].index(min(ocean_medians[item]))]))




def gene_vs_location(infile,celllist):
    celloceans={}#cell:ocean
    oceancount={}#ocean:count of cells
    for cell in celllist:
        num=cell[3:6]
        ll=latlongdict[num]
        lat=ll[0]
        longitude=ll[1]
        foundbool=False
        for ocean in ocean_ranges_lat:
            if cyclic_range(ocean_ranges_lat[ocean],lat):
                if cyclic_range(ocean_ranges_long[ocean],longitude):
                    foundbool=True
                    celloceans[cell]=ocean
                    if ocean in oceancount:
                        oceancount[ocean]+=1
                    if ocean not in oceancount:
                        oceancount[ocean]=1
                    break
    print(oceancount)
    #print(celloceans)
    
    ocean_divs={}#[ocean,ocean]:[means]
    ocean_divs_pairs={}
    seqdict=extract_seqs(infile,celllist)


    coveredcells=list(seqdict.keys())
    for i in range(len(coveredcells)):
        c1=coveredcells[i]
        seqi=seqdict[c1]
        for j in range(i):
            c2=coveredcells[j]
            seqj=seqdict[c2]
            ndiff=sum(1 for a, b in zip(seqi, seqj) if a != b and a!='N' and b!='N')
            nmc=sum(1 for a, b in zip(seqi,seqj) if a!='N' and b!='N')
            mean=float(ndiff)/nmc
            oc1=celloceans[c1]
            oc2=celloceans[c2]
            ocpair=sorted([oc1,oc2])
            ocpair=(ocpair[0],ocpair[1])
            if ocpair in ocean_divs:
                ocean_divs[ocpair].append(mean)
                ocean_divs_pairs[ocpair].append((c1,c2))
            if ocpair not in ocean_divs:
                ocean_divs[ocpair]=[mean]
                ocean_divs_pairs[ocpair]=[(c1,c2)]
    for item in ocean_divs:
        #print('%s Closest Pair Mean: %f ; Median %f ; Furthest Pair Mean %f; Median %f' %(item,min(ocean_means[item]),min(ocean_medians[item]),max(ocean_means[item]),max(ocean_medians[item])))
        print('%s Closest Pair Mean: %f %s ' %(item,min(ocean_divs[item]),ocean_divs_pairs[item][ocean_divs[item].index(min(ocean_divs[item]))]))

def extract_seqs(infile,celllist):
    seqdict={}
    lengthdict={}
    with open(infile,'r') as myinfile:
        lines=[line.rstrip() for line in myinfile]
    myseq=''
    mycell=''
    for line in lines:
        if line[0]=='>':
            if mycell!='':
                #if mycell not in celllist:
                    #print(mycell,'not in HLI?',mynewhandlist.index(mycell))
                if mycell in celllist:
                    seqdict[mycell]=myseq
                    lengthdict[mycell]=len(myseq)-myseq.count('-')
                myseq=''
            ind1=line.find('AG-')
            mycell=line[ind1:ind1+10]
            if mycell=='AG-676-P15':
                mycell='AG-676-P15-1'
        if line[0]!='>':
            myseq+=line
    mylengths=[]
    for item in lengthdict:
        mylengths.append(lengthdict[item])
    mycutoff=0.9*statistics.median(mylengths)
    print('At least %i AA/NT' %mycutoff)
    print(len(seqdict))
    for item in lengthdict:
        if lengthdict[item]<mycutoff:
            seqdict.pop(item)
    print('After applying length cutoff:', len(seqdict))
    return seqdict




