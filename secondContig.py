#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

#### new second contigs set constructing
#read genomes data
genomes_file = open(r'./secondContig/genomes.txt', 'r')
genomes =  genomes_file.readlines()
genomes = genomes[2:]

idxs = []
for idx, line in enumerate(genomes): 
    if 'ChrNumber' in line:
        idxs.append(idx)
        print(idx)

######################
######################
##get generalized adjacencies of weight more than 2 but not included
##in the ancestral contigs

#adj before mwm:

adjs = []
adjs_copy = []
for t in open(r'./secondContig/1stOld/W7InputTreeNode3_50_10'
              '.txt').read().split("),"):
    a, b, c = t.strip('()').split(',')
    if int(c) > 0:
        adjs.append((int(a), int(b), int(c)))
        adjs_copy.append((int(a), int(b)))


#adj after mwm
mwm = open(r"./secondContig/1stOld/W7TreeNode3_50_10.txt", 'r')
mwm = mwm.readlines()
for l in range(0, len(mwm)):
    mwm[l] = tuple(int(item) for item in mwm[l].split())

#all mwm input adjs whixh are not in the mwm output
adj2 = list(set(adjs).difference(mwm))

######################
######################



##get all duplicated adjs in each genome:
dupsAdj = []
     
for i in range(0, len(idxs)):
    # i=0
    if i == len(idxs)-1:
        data = genomes[idxs[i]+1:][1::2]
    else:
        data = genomes[idxs[i]+1:idxs[i+1]][1::2]
    print(i)
    ##get int contigs number
    for z in range(0, len(data)):     
        data[z] = data[z].split()
        data[z] = list(map(int, data[z]))
        
    allAdjs = []
    for p in range(0, len(data)):
        alist=data[p]
        #print(alist)
        z = 0
        while z < len(alist) - 1:
            #testAdj = []
            t1 = alist[z]
            if t1 > 0:
                left1 = 2*t1 - 2
                right1 = 2*t1 -1
            else:
                t1 = abs(t1)
                left1 = 2*t1 - 1
                right1 = 2*t1 - 2
            
            t2 = alist[z+1]
            if t2 > 0:
                left2 = 2*t2 - 2
                right2 = 2*t2 -1
            else:
                t2 = abs(t2)
                left2 = 2*t2 - 1
                right2 = 2*t2 - 2
            testAdj = sorted([right1, left2])
            # if right1==left2==0:
            #     print(p, z)
            allAdjs.append(testAdj)
            z += 1
            
            
            
    ##get duplicate adjacencies in one genome
    dups = {tuple(x) for x in allAdjs if allAdjs.count(x)>1}
    dups = sorted(list(dups))
    #dups1 = sorted([list(x) for x in dups])
    #find the weight for each duplicated adj
    
    for dup in dups:
        try: 
            index = adjs_copy.index(dup)
            # print(dup, index)
            weight = adjs[index][2]
            dupsAdj.append(tuple([dup[0], dup[1], weight]))
        except ValueError:
            continue
   
dupsAdj = set(dupsAdj) #remove copies
dupsAdj2 = list(dupsAdj.intersection(mwm)) ##get all duplicated adjs in mwm

dupsAdj2 = list(dupsAdj.intersection(mwm)) ##get all duplicated adjs in mwm

geneDup = []
for num in dupsAdj2:
    if num[0] %2 == 0:
        geneDup.append(int(num[0]/2 +1))
    else:
        geneDup.append(int((num[0]+1)/2))
        
    if num[1] %2 == 0:
        geneDup.append(int(num[1]/2 +1))
    else:
        geneDup.append(int((num[1]+1)/2))
            
geneDup = sorted(set(geneDup))       
#pd.DataFrame(sorted(geneDup)).to_csv('geneDup.txt', header=None, index=False)
gfs = pd.read_csv(r"./secondContig/1stOld/CleanedGF.csv", header=None)   
for g in geneDup:
    # g = 10
    # gfs2 = copy.deepcopy(gfs)
    # gfs.loc[(gfs[1] == g) & (gfs[8] == 25734), 2] = 26
    # gfs.loc[(gfs[1] == g) & (gfs[8] == 33018), 2] = 17
    # gfs.loc[(gfs[1] == g) & (gfs[8] == 33908), 2] = 11
    # gfs.loc[(gfs[1] == g) & (gfs[8] == 51051), 2] = 22
    # gfs.loc[(gfs[1] == g) & (gfs[8] == 51364), 2] = 21
    # gfs.loc[(gfs[1] == g) & (gfs[8] == 54711), 2] = 13
    gfs = gfs[gfs[1] != g]
    # gfdata = gfs[gfs[1] == g].reset_index(drop=True)
    # gfdata.loc[gfdata[8]==25734, 2] = 26
    
    
gfs.to_csv(r"./secondContig/contigs2nd/CleanedGFanc3.csv", header=None, index=False)           

##get all candidate adjacencies for the second run
adjsAll = set(sorted(adj2+dupsAdj2))
# pd.DataFrame(sorted(dupsAdj2)).to_csv('adjsDup.txt', header=None, index=False, sep=',')







with open (r'./secondContig/W7InputTreeNode3_50_10.txt', 'w') as f:
    print("maxWeightMatching([", file=f)
    for j in adjsAll:
        print("("+str(j[0])+", "+str(j[1])+", "+str(j[2])+"),", end='', file=f)
        # if j == len(adjsAll)-1:
    print("])", file=f)
    #print(adjsAll, sep=',', file =f)



path = os.path.abspath(os.getcwd())

fin = open (r'./secondContig/Sample_MWM_Output.py',"rt")
fout = open (r'./secondContig/mwmInputW7_50_10_TreeNode3.py',"wt")
for line in fin:
    fout.write(line.replace('/Users/qiaojixu/Desktop/6Genomes_Project/TreeNode/MWMOutput/S4.txt',path+'/secondContig/2ndMWM/W7TreeNode3_50_10.txt'))

fin.close()
fout.close()

with open (r'./secondContig//W7InputTreeNode3_50_10.txt') as f:
    with open (r'./secondContig/mwmInputW7_50_10_TreeNode3.py','a') as f1:
        for line in f:
            f1.write(line)
        print('end = time.perf_counter()'+'\n'+'print (end - start)',file =f1)
   

exec(open(r'./secondContig/mwmInputW7_50_10_TreeNode3.py').read())

path2 = os.path.join(path, "secondContig")
##get second contig file
os.system('java -jar GetContig.jar ' + path2 + ' W7TreeNode3_50_10.txt')


        
