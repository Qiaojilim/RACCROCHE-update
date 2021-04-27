#!/usr/bin/env Rscript

###################################################
### This program matches each ancestral genome to extant genomes
### and counts the contig co-occurrence on extant chromosomes for every pair of ancestral contigs
### It also summarizes measures for contig matching
###################################################
### input:  1. genome IDs and ancestor tree nodes defined in Genomes.txt 
###         2. extant genome karyotypes defined in "karyotype" folder: karyotype_genomeID_genomeName.txt, 
###             where genomeID and genomeName match the info in Genomes.txt 
###         3. contig gene feature files for each descendent genome in ./data/contigGFF/ContigGFF_gid_W*TreeNode*_*_*.txt
### output: 1. heatmap for each ancestor: results/clustering/AncestorNode_*_heat.pdf
###         2. heatmap reordered contigs for each ancestor: results/clustering/ancestor*_clusters.csv
###         3. summary of measures: results/ancestorStats/block_measures.csv


source("./module3/config.R")
source("./module3/helper.R")
check.packages("psych")
##install package for cor2dist()


## initialize a data frame for stats for each parameter
stats.DF <- data.frame( ws=numeric(), treeNode=numeric(), gf1=numeric(), gf2=numeric(), 
                        gID=numeric(), gName=character(), avgBLKLen=numeric(), BLKLen.N50=numeric(), coverage=numeric(),
                        avgNoChr=numeric(), stringsAsFactors = FALSE)

for (trn in trn.vector){
  
  ## count cooccurrence for each ancestor tree node
  ctgPairs<-c() ## initialize pairs of cooccured contigs
  
  for (gID in gid.vector){ ## loop through each genome to be analyzed
    
    gName <- as.character(genomeCoGeID[genomeCoGeID$genomeID == gID,]$genomeName)
    
    cat("\n---Analyzing", gName, "for ancestor", trn,"\n")
    
    ### read in extant genome karyotype: chr, size
    karyotype <- readIn.karyotype(karyotype.path, gID)
    
    ### read in gene features in contig
    contigGFF <- readIn.contigGFF(gID, ws, trn, gf1, gf2, nctg, contigGFF.path)  
    if(is.null(contigGFF)) {
      warningMessage <- paste0("*",gID,"_W",ws,"TreeNode",trn,"_",gf1,"_",gf2,".txt file doesn't exist under ", contigGFF.path)
      warning(warningMessage); 
      next
      } ## if file doesn't exist, go to next iteration
   
    
    ### merge genes in contigGFF that are within a DIS.threshold distance in extant genome into blocks
    blockDF <- generate.blockDF.2(contigGFF, DIS.threshold, ws) ## blocks with window size ws
    
    
    # add a column lenBLK: length of block
    blockDF <- setDT(blockDF)[, lenBLK := end - start]
    
    ###########
    ###########
    
    ###########
    ###########
    
    ############################################
    ### calculate block measures: 
    ### average block length in bp
    ### block length N50
    ### average number of chromosomes on which a contig produces a "significant" size block  
    ### extant genome coverage
    ############################################
    
    ## average length and N50 of blocks in bp
    avgBlockLen <- mean(blockDF$lenBLK)
    blockLen.N50 <- N50(blockDF$lenBLK)
    
    ## How many blocks total per chromosome  
    ## How many different contigs/colors per chromosome
    chr.st <- blockDF[blockDF$lenBLK>=blockLEN.threshold,] %>% 
      group_by(chr) %>%
      summarise(num_blocks = length(chr), num_diff_blk = length(unique(contig)))
    
    ## For each contig/color, how many chrs is it on
    ## excludes small blocks shorter than blockLEN.threshold
        
    blk.st <- blockDF[blockDF$lenBLK>=blockLEN.threshold,] %>% 
      group_by(contig) %>%
      summarise(no_chr = length(unique(chr))) 
    
    ctg.st <- contigGFF %>% group_by(contig) %>% summarize(no_chr = length(unique(chr)))
    ## average number of chromosomes on which a contig produces a "significant" size block (excludes small blocks shorter than blockLEN.threshold)
    avgNoChr <- mean(ctg.st$no_chr)
    
    ## block coverage over all chromosomes
    coverage <- setDT(merge(karyotype, aggregate(blockDF$lenBLK, by=list(chr=blockDF$chr), FUN=sum), by.x="chr", by.y="chr"))
    coverage[,nCovered := size - x]  
    
    ##QX: read dup genes contig matrix
    matrixDup <- readIn.dupContig(ws, trn, gf1, gf2, nctg, contig.path)
    
    ## calculate extant genome coverage without counting overlapped blocks
    total.coverage <- c()
    for (c in unique(blockDF$chr)){
      coverage.2 <- setDT(data.frame(with(blockDF[chr==c,], sets(start, end, chr, 1)))) ## take union of the intervals
      colnames(coverage.2) <- c("start","end")
      coverage.2 <- coverage.2[,LenCovered:=end-start]
      total.coverage <- c(total.coverage, sum(coverage.2$LenCovered))
    }
    pCoverage.2 <- sum(total.coverage) / sum(karyotype[,2])#
    
    stats.DF <- rbind(stats.DF, data.frame(ws=ws, treeNode=trn, gf1=gf1, gf2=gf2, gID=gID, gName=gName, 
                                           avgBLKLen=avgBlockLen, BLKLen.N50=blockLen.N50, 
                                           coverage=pCoverage.2, avgNoChr=avgNoChr) )
    
    
    ################################################################
    ### counting co-occurrence of contigs from all extant genomes
    ################################################################
    
    ## get uniq rows of (chr, contig), sort by chr then by contig
    ## c(1,4) means combine the first and forth col of blockDF as a new two cols
    uniq <- unique(blockDF[blockDF$lenBLK>=lenBLK.threshold,c(1,4)]) 
    #print(uniq)
    uniq <- uniq[order(uniq$chr,uniq$contig)]
    #print(uniq)
    ## get all combinations of cooccurrence of contigs on the same chromosome
    for (c in unique(uniq$chr)) {
      c <- as.numeric(as.character(c))
      if(length(uniq[uniq$chr==c,]$contig)>=2) {
        combinations <- combn(as.numeric(as.character(uniq[uniq$chr==c,]$contig)),m=2) ## m=2: cooccurrence of two contigs on the same chromosome
        ctgPairs <- cbind(ctgPairs, combinations) ## combine all combinations 
        ## note: by applying cbind, values in pairs are added by 1, starting from 1
        ## contig0 becomes 1
      }
    }
    
  } ## done looping through each genome and counting contig co-occurrence
  
  
  
  ############################################
  ### construct co-occurrence matrix for the current ancestor
  ### gather pairs of contigs that appear on the same chromosome, 
  ### then construct the cooccurrence matrix
  ############################################
  
  if(length(ctgPairs)==0) {
    warningMessage <- paste0("No contig pairs detected (co-occurred contigs) in genome ", gName)
    warning(warningMessage); 
    next
    }
  
  ## count frequencies of uniq ctgPairs
  pairs.freq <- data.frame(t(ctgPairs))
  colnames(pairs.freq) <- c("ctgID1", "ctgID2")
  pairs.freq <- count(pairs.freq, vars = c("ctgID1", "ctgID2")) ## add column "freq" that counts frequency of uniq rows
  pairs.freq <- pairs.freq[order(pairs.freq$ctgID1, pairs.freq$ctgID2),]
  
  mat <- as.matrix(list2dist.1(pairs.freq))
  mat[is.na(mat)] <- 0
  #print(mat)
  correlation <- cor(mat)
  ####QX: new correlation cut down to 1/4
  if (NROW(correlation) != nctg) {
    correlation <- correlation[, order(as.integer(colnames(correlation)))]
    correlation <- correlation[order(as.integer(rownames(correlation))), ]
    rowele <-as.array( as.integer(rownames(correlation)))
    miss <- seq.int(0, nctg-1) %in% rowele
    miss <- as.vector(which(FALSE == miss))
    miss <- miss - 1
    matrixDup <- matrixDup[, -miss]
    matrixDup <- matrixDup[-miss, ]
   # for (q in miss){
      #q=239
      #delq <- seq.int(0, nctg-1)[q]+1
      #matrixDup <- matrixDup[, -delq]
      ##}
    
    mat4 <- correlation*matrixDup
    }
  else {
    correlation <- correlation[, order(as.integer(colnames(correlation)))]
    correlation <- correlation[order(as.integer(rownames(correlation))), ]
    mat4 <- correlation*matrixDup
  }

  #colnames(correlation)
  #NROW(correlation)
  #NROW(matrixDup)
  
  ####correlation <- mat
  #print(correlation)
  ##correlation matrix from co-occurence matrix
  ## construct contig distance matrix from cooccurrence matrix
  max_freq <- sum(genomeCoGeID$numChr)
  mat2 <- (max_freq - mat)/max_freq# distance is between 0 and 1
  diag(mat2) <- 0
  
  mat3 <- cov(mat2,mat2, method = c("pearson"))
  
  #vec <- as.vector(correlation)
  vec <- as.vector(mat4)
  vec2 <- sort(vec, decreasing = FALSE)
  ###print(vec2)
  #get nth value in 9 groups as breaks
  #l <- (sum(correlation)-length(correlation))%/%9
  l <- length(vec2)
  print(l)
  sm <- min(vec2)
  lr <- max(vec2)
  colors=c(sm,
    nth(vec2,as.integer(l*0.5)),
    nth(vec2,as.integer(l*0.65)),
    nth(vec2,as.integer(l*0.75)),
    nth(vec2,as.integer(l*0.81)),
    nth(vec2,as.integer(l*0.85)),
    nth(vec2,as.integer(l*0.89)),
    nth(vec2,as.integer(l*0.93)),
    nth(vec2,as.integer(l*0.995)), lr)#0.995
  
  my_colors <- c(y1 = "#F7F7D0",
      y2 = "#FCFC3A",
      y3 = "#D4D40D",
      b1 = "#40EDEA",
      b2 = "#18B3F0",
      b3 = "#186BF0",
      r1 = "#FA8E8E",
      r2 = "#F26666",
      r1 = "#C70404")
  ##correlation <- correlation + 1
  print(colors)
  ## write distance matrix into file under the directory of results.path
  # distMat.fname <- paste0(results.path, "distanceMat_trn",trn,"_W",ws,"(",gf1,",",gf2,")_",lenBLK.threshold/1000,"kbBLK.csv")
  # write.csv(mat2, file=distMat.fname, row.names=FALSE)
  
  
  ############################################
  ########### generate heatmaps to group contigs into ancestral chromosomes
  ############################################
  
  ## convert distance to -log d
  ## -log d, the larger, the darker blue color
  d3 <- -log(mat2); diag(d3) <- 0
  d3[which(!is.finite(d3))] <- 0
  #print(d3)
  
  
    
  pdf(mat4, file=file.path(results.path, "clustering", paste0("AncestorNode_",trn,"_heat.pdf")), width=40, height=40)
  par(cex.main=4)
  p <- heatmap.2(mat4,
                 main = paste0("Ancestor ",trn),
                 dendrogram="row",     # only draw a row dendrogram/ none
                 hclustfun=hclust,
                 #hclustfun=function(x)
                      #{hclust(x,method="average")}, #hclust, ## defaut method for hclust: "complete"
                 srtCol=0,   adjCol = c(0.5,1),
                 lhei=c(.1,1), lwid=c(.2,1), key.title="Color Key", keysize=0.75,
                 col = brewer.pal(9, "Greys"), breaks = colors,#Greys
                 trace = "none")
  ##print(d3)
  print(p)
  dev.off()
  
  write.csv(rev(rownames(p$carpet)), file=file.path(results.path, "clustering", paste0("cluster_trn",trn,".csv")), row.names=TRUE)
  
  #pdf(file=file.path(results.path, "clustering", paste0("KmeansAncestorNode_",trn,"_heat.pdf")), width=40, height=40)
  #par(cex.main=4)
  #k2 <- kmeans(p, centers = 7, nstart = 25)
  #kcluster <- fviz_cluster(k2, data = p)
  #dev.off()
  
} ## loop through each ancestor

### export analysis stats data 
write.csv(stats.DF, file=file.path(results.path, "ancestorStats", "block_measures.csv"), row.names=FALSE)

message("\n~~~~~Rscript finished analyzing co-occurrence and clustering")
