#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os


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
for t in open(r'./secondContig/old/W7InputTreeNode2_50_10'
              '.txt').read().split("),"):
    a, b, c = t.strip('()').split(',')
    if int(c) > 0:
        adjs.append((int(a), int(b), int(c)))
        adjs_copy.append((int(a), int(b)))


#adj after mwm
mwm = open(r"./secondContig/old/W7TreeNode2_50_10.txt", 'r')
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

   

##get all candidate adjacencies for the second run
adjsAll = set(sorted(adj2+dupsAdj2))
# pd.DataFrame(sorted(dupsAdj2)).to_csv('adjsDup.txt', header=None, index=False, sep=',')







with open (r'./secondContig/W7InputTreeNode2_50_10.txt', 'w') as f:
    print("maxWeightMatching([", file=f)
    for j in adjsAll:
        print("("+str(j[0])+", "+str(j[1])+", "+str(j[2])+"),", end='', file=f)
        # if j == len(adjsAll)-1:
    print("])", file=f)
    #print(adjsAll, sep=',', file =f)



path = os.path.abspath(os.getcwd())

fin = open (r'./secondContig/Sample_MWM_Output.py',"rt")
fout = open (r'./secondContig/mwmInputW7_50_10_TreeNode2.py',"wt")
for line in fin:
    fout.write(line.replace('/Users/qiaojixu/Desktop/6Genomes_Project/TreeNode/MWMOutput/S4.txt',path+'/secondContig/MWMOutput/W7TreeNode2_50_10.txt'))

fin.close()
fout.close()

with open (r'./secondContig//W7InputTreeNode2_50_10.txt') as f:
    with open (r'./secondContig/mwmInputW7_50_10_TreeNode2.py','a') as f1:
        for line in f:
            f1.write(line)
        print('end = time.perf_counter()'+'\n'+'print (end - start)',file =f1)
   

exec(open(r'./secondContig/mwmInputW7_50_10_TreeNode2.py').read())   

path2 = os.path.join(path, "secondContig")
##get second contig file
os.system('java -jar GetContig.jar ' + path2 + ' W7TreeNode2_50_10.txt')


        