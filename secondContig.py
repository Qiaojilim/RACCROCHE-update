#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from itertools import groupby
from operator import itemgetter
import numpy as np    
import copy
def find_gf(c,alist):
    for i, sub in enumerate(alist):
        try:
            j = [index for index, value in enumerate(sub) if value == c]
            #j = sub.index(c)
        except ValueError:
            continue
        for t in j:
            yield i, t

def triangular(M, N):
    for j in range(0, len(M)):
        i = j + 1
        while i < N:
            
            M.loc[i,j] = M.loc[j,i]
            i += 1
    return M

# path1 ='/Users/qiaojixu/Desktop/RACCROCHE/project-monocots'


#### new second contigs set constructing
#read genomes data
genomes_file = open(r'./project-monocots/data/GeneFamily/genomes.txt', 'r')
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
nodes=[1,2,3,4]
for tr in nodes:
    treenode= str(tr)
    # treenode= str(4)
    N = 500   
    
    WS=7
    gf1=50
    gf2=10
    TreeNode ='TreeNode'+treenode          
    
    #adj before mwm:
    
    # adjs = []
    # adjs_copy = []
    # for t in open(r'./secondContig/1stOld/W7InputTreeNode4_50_10'
    #               '.txt').read().split("),"):
    #     a, b, c = t.strip('()').split(',')
    #     if int(c) > 0:
    #         adjs.append((int(a), int(b), int(c)))
    #         adjs_copy.append((int(a), int(b)))
    
    adjs = []
    adjs_copy = []
    afile = open(r'./project-monocots/results/InputPyfile/W7InputTreeNode'+treenode+
                          '_50_10.txt', 'r')
    afile = afile.readlines()
    afile = afile[1].replace("])", "")
    # afile = afile.replace("\n", "")
    
    afile = afile.split("),")
    # for t in open(os.path.join(path1,'results/old/InputPyfile/simulated3/W7InputTreeNode'+treenode+
    #                       '_50_10.txt')).read().split("),"):
    for t in afile:    
        # print(t)
        if len(t.strip('()').split(',')) == 3:
            a, b, c = t.strip('()').split(',')
            if int(c) > 0:
                adjs.append((int(a), int(b), int(c)))
                adjs_copy.append((int(a), int(b)))
    
    
    # #adj after mwm
    # mwm = open(r"./secondContig/1stOld/W7TreeNode4_50_10.txt", 'r')
    # mwm = mwm.readlines()
    # for l in range(0, len(mwm)):
    #     mwm[l] = tuple(int(item) for item in mwm[l].split())
    
    # #all mwm input adjs whixh are not in the mwm output
    # adj2 = list(set(adjs).difference(mwm))
    
    
    mwm = open(r'./project-monocots/results/InputPyfile/mwmOutput/W7TreeNode'+
                            treenode+'_50_10.txt', 'r')
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
    dupsAdj = sorted(dupsAdj)
    dupsAdj3 = []
    for num, item in enumerate(mwm): ## num is index and item id the element
        #num[0]
        for num2, item2 in enumerate(dupsAdj):
            if item2[0] % 2 == 0:
                if item[0] - item2[0] == 0 or item[0] - item2[0] == 1:
                    if item2[1] %2 == 0:
                        if item[1] - item2[1] == 0 or item[1] - item2[1] == 1:
                            dupsAdj3.append(item2)
                    else:
                        if item[1] - item2[1] == 0 or item[1] - item2[1] == -1:
                            dupsAdj3.append(item2)
                        
            else:
                if item[0] - item2[0] == 0 or item[0] - item2[0] == -1 :
                    if item2[1] %2 == 0:
                        if item[1] - item2[1] == 0 or  item[1] - item2[1] ==1:
                            dupsAdj3.append(item2)
                    else:
                        if item[1] - item2[1] == 0 or item[1] - item2[1] == -1:
                            dupsAdj3.append(item2) 
                    
                
                
                
            # if abs(item2[0] - item[0]) <= 1 and abs(item2[1] - item[1]) <= 1:
            #     dupsAdj3.append(item2)
    
    geneDup = []
    for num in dupsAdj3:
        if num[0] %2 == 0:
            geneDup.append(int(num[0]/2 +1))
        else:
            geneDup.append(int((num[0]+1)/2))
            
        if num[1] %2 == 0:
            geneDup.append(int(num[1]/2 +1))
        else:
            geneDup.append(int((num[1]+1)/2))
                
                
    geneDup = sorted(set(geneDup))           
    pd.DataFrame(sorted(geneDup)).to_csv(r'./secondContig/contigs2nd/dupGenes/W7TreeNode'+treenode+'geneDup.txt', header=None, index=False)
    # gfs = pd.read_csv("CleanedGF.csv", header=None)   
        
    
    
    ##get all candidate adjacencies for the second run
    adjsAll = set(sorted(adj2+dupsAdj3))
    # pd.DataFrame(sorted(dupsAdj2)).to_csv('adjsDup.txt', header=None, index=False, sep=',')
    
    
    
    
    
    
    
    with open (r'./secondContig/W' + str(WS) + 'Input'+TreeNode + '_'+ str( gf1)+'_'+str(gf2) +'.txt', 'w') as f:
        print("maxWeightMatching([", file=f)
        for j in adjsAll:
            print("("+str(j[0])+", "+str(j[1])+", "+str(j[2])+"),", end='', file=f)
            # if j == len(adjsAll)-1:
        print("])", file=f)
        #print(adjsAll, sep=',', file =f)
    
    
    
    path = os.path.abspath(os.getcwd())
    
    fin = open (r'./secondContig/Sample_MWM_Output.py',"rt")
    fout = open (r'./secondContig/mwmInputW'+str(WS)+'_'+str( gf1)+'_'+str(gf2)+'_' + TreeNode+'.py',"wt")
    for line in fin:
        fout.write(line.replace('/Users/qiaojixu/Desktop/6Genomes_Project/TreeNode/MWMOutput/S4.txt',path+'/secondContig/2ndMWM/W'+str(WS)+TreeNode + '_'+ str( gf1)+'_'+str(gf2) +'.txt'))
    
    fin.close()
    fout.close()
    
    with open (r'./secondContig/W' + str(WS) + 'Input'+TreeNode + '_'+ str( gf1)+'_'+str(gf2) +'.txt') as f:
        with open (r'./secondContig/mwmInputW' +str(WS)+'_'+ str( gf1)+'_'+str(gf2)+'_'+ TreeNode+'.py','a') as f1:
            for line in f:
                f1.write(line)
            print('end = time.perf_counter()'+'\n'+'print (end - start)',file =f1)
       
    
    exec(open(r'./secondContig/mwmInputW'+str(WS)+'_50_10_'+TreeNode+'.py').read())
    
    path2 = os.path.join(path, "secondContig")
    ##get second contig file
    os.system('java -jar GetContig.jar ' + path2 + ' W'+str(WS)+TreeNode+'_50_10.txt')
    
    
    geneDup = pd.read_csv(r'./secondContig/contigs2nd/dupGenes/W7TreeNode'+treenode+'geneDup.txt', header=None)
    geneDup=list(geneDup[0])
    ###########
    #combine new and old contigs as one
    #read old ocntig file
    # contig1 = open('/Users/qiaojixu/Desktop/RACCROCHE/project-monocots/data/'
                   # '/Contig/old/ContigW7TreeNode'+treenode+'_50_10.txt', 'r')
    # contig1 = open(os.path.join(path1, 'data/Contig/old/ContigW7TreeNode'+
    #                             treenode+'_50_10.txt'), 'r')
    contig1 = open(r'./project-monocots/data/Contig/ContigW7TreeNode'+
                                treenode+'_50_10.txt', 'r')
    
    
    contig1 = contig1.readlines()
    contig1 = contig1[1::2]
    size1 = []
    for c in range(0,len(contig1)):
        contig = contig1[c]
        size = len(contig.split())
        size1.append(tuple([1, size, c]))
    
    
    
    contig2 = open(r'./secondContig/contigs2nd/ContigW7TreeNode'+treenode+'_50_10.txt', 'r')
    contig2 = contig2.readlines()
    contig2 = contig2[1::2]
    size2 = []
    for c in range(0,len(contig2)):
        contig = contig2[c]
        size = len(contig.split())
        size2.append(tuple([2, size, c]))
    
    # size12 = size1 + size2
    # # size12.sort(key=lambda x:-x[1])
    
    # #ContigW7TreeNode4_50_10old
    
    # with open (r'./Contig/ContigW7TreeNode'+treenode+'_50_10.txt', 'w') as f:
    #     for z in range(0, len(size12)):
    #         #z=0
    #         ss = size12[z]
    #         if ss[0] == 1:
    #             t = ss[2]
    #             out = contig1[t].strip()
    #         else:
    #             t = ss[2]
    #             out = contig2[t].strip()
    #         if ss[0] ==1 and ss[1]==1:
    #             break
    #         else:
    #             print("contig "+str(z), file = f)
    #             print("  "+out, file = f)
    "combine 250 old and 250 new contigs together"        
    size12 = size1[0:N//2] + size2[0:N//2]
    # size12.sort(key=lambda x:-x[1])
    # pd.DataFrame(size12).to_csv('size.csv', header=None, index=False)
    #ContigW7TreeNode4_50_10old
    
    with open (r'./secondContig/contigs2nd/combineContig/ContigW7TreeNode'+treenode+'_50_10_'+str(N)+'.txt', 'w') as f:
        for z in range(0, len(size12)):
            #z=0
            ss = size12[z]
            if ss[0] == 1:
                t = ss[2]
                out = contig1[t].strip()
            else:
                t = ss[2]
                out = contig2[t].strip()
            if ss[0] ==1 and ss[1]==1:
                break
            else:
                print("contig "+str(z), file = f)
                print("  "+out, file = f)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    r'./Contig/ContigW7TreeNode2_50_10.txt'
    contig12 = open(r'./secondContig/contigs2nd/combineContig/ContigW7TreeNode'+treenode+'_50_10_'+str(N)+'.txt', 'r')
    
    contig12 = contig12.readlines()
    contig12 = contig12[1::2]
    for l in range(0, len(contig12)):
        contig12[l] = [int(item) for item in contig12[l].split()]
    
    # count = 0
    # contigPos = []
    # for g in geneDup:
    #     # g=1
    #     match1 = [match for match in find_gf(g, contig12)] ##contig index
    #     match2 = [match for match in find_gf(-g, contig12)] ## gene index
    #     m1 = []
    #     m2 = []
    #     if len(match1) > 0:
    #         for mm in match1:  
    #             # print(mm)
    #             m1.append(mm[0])
        
    #     if len(match2) > 0:
    #         for mm2 in match1:   
    #             m2.append(mm2[0])
        
    #     collect = set(m1 + m2)
    #     if len(collect) > 0:
    #         count += 1
    #         for item in collect:
    #             contigPos.append(item)
    contig12c = copy.deepcopy(contig12)
    for ind, item in enumerate(contig12c):
        contig12c[ind] = [abs(t) for t in item]
    
    
    # contigPos = list(set(contigPos))
    
    # N = 500
    matrixDup = pd.DataFrame(index=range(0, N), columns=range(0, N))
    
    for i in range(0, N//2):
        # i = 0
        contig_index1 = contig12c[i]
        j = N//2
        while j < N :
            # j = 150
            contig_index2 = contig12c[j]
            
            common= set(contig_index1).intersection(set(contig_index2))
            # common=set([1,2,10])
            if len(common) > 0:
                check =  common.intersection(set(geneDup))
                if len(check) > 0:    
                    matrixDup.loc[i, j] = len(check)
            j += 1
    
    
    
    
    
    
    
    
    
    
    # N=250
    # contigPos.sort()
    # size = []
    # for num in set(contigPos):
    #     if num < N:
    #         l = contigPos.count(num)
    #         size.append((num,l))
    ##create a 250*250 matrix to contain dup genes size
    # matrixDup = pd.DataFrame(index=range(0, N), columns=range(0, N))
    # matrixDup = matrixDup.fillna(0)
    
    # for i in range(0, len(size)):
    #     pair1 = size[i]
    #     print(pair1)
    #     contig_index1 = pair1[0]
    #     l1 = pair1[1]
    #     if i < len(size) - 1:
    #         for j in range(i+1, len(size)):
    #             pair2 = size[j]
    #             contig_index2 = pair2[0]
    #             l2 = pair2[1]
    #             # matrixDup[contig_index1][contig_index2] = l1 + l2
    #             matrixDup[contig_index2][contig_index1] = l1 + l2
    
    
    # matrixDup['Sum'] = matrixDup.sum(axis=1)
    
    matrixDup['mean'] = matrixDup.mean(axis=1)
    matrixDup = matrixDup.fillna(0)
    for i in range(0, len(matrixDup)):
        ave = matrixDup['mean'][i]
        if ave > 0:        
            j = i + 1
            while j < N:
                if 0 < matrixDup.loc[i,j]:
                    matrixDup.loc[i,j] = 0.25
                # if 0 < matrixDup.loc[i,j] <= ave:
                #     matrixDup.loc[i,j] = 0.5
                # if matrixDup.loc[i,j] > ave:
                #     matrixDup.loc[i,j] = 0.25
                j += 1    
    
    del matrixDup['mean']
    matrixDup = triangular(matrixDup,N)
    
    matrixDup =  matrixDup.replace(to_replace=0, value=1)
    matrixDup.to_csv(r'./project-monocots/data/Contig/W7TreeNode'+treenode+
                                  '_50_10dupMatrix'+#str(N)
                                  '.txt', sep=',', index=False )
    
