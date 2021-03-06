#!/usr/bin/env Rscript

###################################################
### This program matches and paints ancestral chromosomes to extant genomes
### It also summarizes measures for choppiness
###################################################
### input:  1. genome IDs and ancestor tree nodes defined in Genomes.txt 
###         2. extant genome karyotypes defined in "karyotype" folder: karyotype_genomeID_genomeName.txt, 
###             where genomeID and genomeName match the info in Genomes.txt 
###         3. contig gene feature files for each descendent genome in ./data/contigGFF/ContigGFF_gid_W*TreeNode*_*_*.txt
###         4. clustering results in ./data/clustering/
### output: 1. painted chromosomes: results/paintedChrs/ancestor[trn]-[genomeName].pdf
###         2. choppinness : results/ancestorStats/choppiness_sum.csv

source("./module3/config.R")
source("./module3/helper.R")

## 7 ancestral chromosomes
myClr <- c("green", "red", "blue", "purple", "yellow", "cyan", "orange" )  


#####################################################################
#####################################################################
### paint extant genomes with ancestral contigs at each ancestor tree node
### then perform choppiness analysis

# methods <- c("CompleteLink", "SingleLink", "AverageLink", "Ward", "k-means")
# for (method in methods) {
  
choppiness <- data.frame(genome=character(), ancestor=character(), chr=character(), t=numeric(), r=numeric(), x=numeric())

## initialize a data frame for stats for each parameter
stats.DF <- data.frame( ws=numeric(), treeNode=numeric(), gf1=numeric(), gf2=numeric(),
                        gID=numeric(), gName=character(), avgBLKLen=numeric(), BLKLen.N50=numeric(), coverage=numeric(),
                        avgNoChr=numeric(), stringsAsFactors = FALSE)


for(trn in trn.vector){
  # trn=2
  # gID <- 19990  # only plot Vitis for testing
  for(gID in gid.vector) {
    
    gName <- genomeCoGeID[genomeCoGeID$genomeID == gID,]$genomeName
    myTrn <- genomeCoGeID[genomeCoGeID$genomeID == gID,]$ancestor
    

    ### read in extant genome karyotype: chr, size
    karyotype <- readIn.karyotype(karyotype.path, gID)
    chrMaxLen <- max(karyotype$size)+100
    
    ## read in contigGFF. Each row is one gene from ancestral contigg, sorted by extant chr #, then positions in extant chr
    # Colnames: chr of extant genome, geneFamilyID, pos on extant chr, ancestral contig#, 
    # start position on extant chr, end position on extant chr, distance between two genes
    contigGFF <- readIn.contigGFF(gID, ws, trn, gf1, gf2, nctg, contigGFF.path)
    
    
    ## generate plotting data from contigGFF !!!
    ## 7 cloumns in contigGFF: "chr",	"geneFamilyID",	"pos",	"contig", "start", "end", "distance"
    ## blockDF columns: chr, start, end, contig --> make sure they are all numeric
    ### merge genes in contigGFF that are within a DIS.threshold distance in extant genome into blocks
    blockDF <- generate.blockDF.2(contigGFF, DIS.threshold, ws)
    
    
    ## read in ancestral chromosomes from the clustering results
    clusterVector <- scan(file.path(results.path, "clustering", paste0("cluster_trn",trn,".txt")))
    # clusterVector <- scan(file.path(results.path, "clustering", paste0("trn",trn, "_", method, ".txt")))
    clusters <- cbind(0:(nctg-1), as.data.frame(clusterVector))
    colnames(clusters)<-c("contig","ancestralChr")
    
    
    
    ## add one column: ancestralChr from corresponding cl$cluster
    blockDF[,"contig"] <- as.numeric(as.character(blockDF$contig)) ## unfactorize column contig
    blockDF <- merge(blockDF, clusters)
    blockDF <- blockDF[order(blockDF$chr,blockDF$start),]
    blockDF <- blockDF[ , c("chr", "start", "end", "contig", "ancestralChr")]
    blockDF[, len := end - start ]
    
    # ### export synteny blocks before merging
    # ### this will be used to order contigs within each ancestral chromosome using LOP
    # write.csv(blockDF, file=file.path(results.path, "clustering", "AncestralSyntenyBlocks.csv"), row.names=FALSE)
    
    ##############################################
    ### generate plots of painted extant genomes
    ##############################################
    
    ## merge adjancent blocks if their distance is within DIS.threshold (i.e. 1 MB)
    ## only choose blocks longer than blockLEN threshold to merge
    mergedDF <- mergeBlockDF(blockDF[(end-start) > blockLEN.threshold,], DIS.threshold)
    mergedDF <- as.data.table(mergedDF)
    mergedDF[, len := end - start ]
    
    
    ############################################
    ### calculate block measures: 
    ### average block length in bp
    ### block length N50
    ### average number of chromosomes on which a contig produces a "significant" size block  
    ### extant genome coverage
    ## !! based on blockDF or mergedDF??
    ############################################
    
    ## average length and N50 of blocks in bp
    avgBlockLen <- mean(blockDF$len)
    blockLen.N50 <- N50(blockDF$len)
    
    ## How many blocks total per chromosome  
    ## How many different contigs/colors per chromosome
    chr.st <- blockDF[blockDF$len>=blockLEN.threshold,] %>% 
      group_by(chr) %>%
      summarise(num_blocks = length(chr), num_diff_blk = length(unique(ancestralChr)))
    
    ## For each contig/color, how many chrs is it on
    ## excludes small blocks shorter than blockLEN.threshold
    blk.st <- blockDF[blockDF$len>=blockLEN.threshold,] %>% 
      group_by(ancestralChr) %>%
      summarise(no_chr = length(unique(chr))) 
    
    ctg.st <- contigGFF %>% group_by(contig) %>% summarize(no_chr = length(unique(chr)))
    ## average number of chromosomes on which an ancestralChr produces a "significant" size block (excludes small blocks shorter than blockLEN.threshold)
    avgNoChr <- mean(ctg.st$no_chr)
    
    ## block coverage over all chromosomes
    coverage <- setDT(merge(karyotype, aggregate(blockDF$len, by=list(chr=blockDF$chr), FUN=sum), by.x="chr", by.y="chr"))
    coverage[,nCovered := size - x]  
    
    ## calculate extant genome coverage without counting overlapped blocks
    total.coverage <- c()
    for (c in unique(blockDF$chr)){
      coverage.2 <- setDT(data.frame(with(blockDF[chr==c,], sets(start, end, chr, 1)))) ## take union of the intervals
      colnames(coverage.2) <- c("start","end")
      coverage.2 <- coverage.2[,LenCovered:=end-start]
      total.coverage <- c(total.coverage, sum(coverage.2$LenCovered))
    }
    pCoverage.2 <- sum(total.coverage) / sum(karyotype[,2])#
    
    stats.DF <- rbind(stats.DF,
                      data.frame(ws=ws, treeNode=trn, gf1=gf1, gf2=gf2, gID=gID, gName=gName,
                                 avgBLKLen=avgBlockLen, BLKLen.N50=blockLen.N50,
                                 coverage=pCoverage.2, avgNoChr=avgNoChr) )
    
    #######################
    ## choppiness analysis
    #######################
    for(c in karyotype$chr){
      df.chr <- mergedDF[mergedDF$chr==c,] ## !! choppiness based on blockDF or mergedDF??
      
      ## T_i = the number of different colours on each chromosome - 1
      t <- length(unique(df.chr$ancestralChr)) - 1
      
      ## R_i = the number of single-colour regions on each chromosome
      r <- length(df.chr$ancestralChr) 
      
      ## the number of stripes less than a certain threshold size 
      x <- nrow(df.chr[df.chr$len <= 2*blockLEN.threshold,])
      
      choppiness <- rbind(choppiness, data.frame(genome=gID, ancestor=trn, 
                                                 chr=as.numeric(as.character(c)), 
                                                 t=as.numeric(as.character(t)), 
                                                 r=as.numeric(as.character(r)), 
                                                 x=as.numeric(as.character(x))))
      
    } ## end of looping through chromosomes
    
    
    ############################################
    ### generate genome plots
    ############################################
    
     # if(trn == myTrn){ # uncomment to only map ancestor to its immediate descendant
      myTrn=trn
      ## output plot file
      pdf(file = file.path(results.path, "paintedChrs", paste0("ancestor", myTrn, "-", gName, ".pdf")),paper = "a4r", width = 0, height = 0)

      ## magic number to draw rectangle
      offset = ifelse (max(karyotype$chr)<=7, 0.05, ifelse (max(karyotype$chr)<=11, 0.07, ifelse (max(karyotype$chr)<=16, 0.11, ifelse (max(karyotype$chr)<=25,0.18,0.3))))

      #########################
      p <- ggplot() +
        geom_segment(data = karyotype,
                     aes(y = chr, yend = chr, x = 0, xend = size),
                     lineend = "round", color = "lightgrey", size = 6) +
        scale_x_continuous("Length (Mbp)", breaks = seq(1,chrMaxLen,by=5000000), labels = Ms2, limits=c(0,chrMaxLen))+
        scale_y_continuous("Chromosome", breaks = karyotype$chr, labels = karyotype$chr ) +

        geom_rect(data=mergedDF, mapping=aes(xmin=start, xmax=end, ymin=as.integer(chr)-offset,
                                             ymax=as.integer(chr)+offset, alpha=0.5),
                  fill=myClr[mergedDF$ancestralChr], show.legend = FALSE,  alpha=0.5)+
        ##### Show ancestral chromosome number
        #geom_text(data=mergedDF, aes(x=start+(end-start)/2, y=chr, label=ancestralChr), size=2) +
        coord_flip()  +
        theme(plot.title = element_text(size=15, face="bold"),
              text = element_text(size=15),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black", size = .5, linetype = "solid")) # enable axis lines #, axis.ticks = element_blank() )

      print(p)

      dev.off()
     # }


   } ## end of looping through all genome ids
} ## end of looping through tree nodes

choppiness.sum <- as.data.table(aggregate(choppiness[,4:6], by=list(genomeID = choppiness$genome, Ancestor.map = choppiness$ancestor), FUN = "sum"))
choppiness.sum[, "R-T" := r-t]
choppiness.sum[, "R-X" := r-x]

choppiness.sum <- merge(genomeCoGeID[,1:2], choppiness.sum, by = "genomeID")
setorder(choppiness.sum, genomeID, Ancestor.map)

###################################
## output choppiness stats to file
## details for debuging
# write.csv(choppiness, file=file.path(results.path, "ancestorStats", "choppiness_raw.csv"), row.names=TRUE) 
## summarized results
write.csv(choppiness.sum, file=file.path(results.path, "ancestorStats", paste0("choppiness_sum_ancestor",trn, method, ".csv")), row.names=TRUE)
write.csv(choppiness.sum, file=file.path(results.path, "ancestorStats", "choppiness_sum.csv"), row.names=TRUE)

### export analysis stats data 
# write.csv(stats.DF, file=file.path(results.path, "ancestorStats", paste0("block_measures_ancestor",trn, method, ".csv")), row.names=FALSE)
write.csv(stats.DF, file=file.path(results.path, "ancestorStats", "block_measures.csv"), row.names=FALSE)

# } ## end of methods

message("\n~~~~~Rscript finished paiting extant chromosomes \n")










#####################################################################
#####################################################################
### match ancestral contigs to each other ancestor
### then perform choppiness analysis

choppiness <- data.frame(genome=character(), ancestor=character(), chr=character(), t=numeric(), r=numeric(), x=numeric())

## initialize a data frame for stats for each parameter
# stats.DF <- data.frame( ws=numeric(), treeNode=numeric(), gf1=numeric(), gf2=numeric(),
#                         gID=numeric(), gName=character(), avgBLKLen=numeric(), BLKLen.N50=numeric(), coverage=numeric(),
#                         avgNoChr=numeric(), stringsAsFactors = FALSE)
# 


for(trn in trn.vector){
    trn2 <- trn + 1
    if (trn2 > max(trn.vector)) {break}
    
    ## read in ancestral chromosomes from the clustering results
    clusterVector <- scan(file.path(results.path, "clustering", paste0("cluster_trn",trn,".txt")))
    
    # ## get gene family features for ancestor
    #genefamilyLength.trn <- getGeneFamilyLength(trn, ws, gf1, gf2, nctg, contigGFF.path)
    
    ancestorGF.trn <- readIn.ancestorGF (contig.path, ws, trn, gf1, gf2, clusterVector)
    #merge(ancestorGF.trn[1:10,], genefamilyLength.trn, by="gfID", all.y=FALSE, sort=FALSE, suffixes = c("", ".trn1"))
    
    
    ancestorGF.trn2 <- readIn.ancestorGF (contig.path, ws, trn2, gf1, gf2, clusterVector)
    
    ## match trn1 to trn2 (trn1 is the ancestor, trn2 is the descendent)
    ancestorGF.trn2 <- merge(ancestorGF.trn2, ancestorGF.trn, by="gfID", all.y=FALSE, sort=FALSE, suffixes = c("", ".trn1"))
    
    ancestorGF.trn2 <- ancestorGF.trn2[order(ancestorGF.trn2$chr,ancestorGF.trn2$contig),]
    
   
    blockDF <- generate.blockDF.3 (ancestorGF.trn2) 
      
    # 
    # blockDF <- data.frame() 
    # for (current.chr in unique(as.numeric(as.character(ancestorGF.trn2$chr)))){
    #   ## loop through every chr
    #   current.ctg.gff <- ancestorGF.trn2[which(ancestorGF.trn2$chr==current.chr),]
    # 
    #   ## start the first block from the first gene
    #   blockDF <- rbind( blockDF, current.ctg.gff[1,c("chr", "contig", "gfID", "chr.trn1")] )
    #   if(nrow(current.ctg.gff)<2) {next; }
    #   for (r in 2:nrow(current.ctg.gff)) {
    #     ## the start of a new chromosome, start a new block
    #     if (current.ctg.gff[r,]$chr != current.ctg.gff[r-1,]$chr) { 
    #       # block.no <- 1;  
    #       # block.vector <- c(block.vector, block.no); 
    #       blockDF <- rbind( blockDF, current.ctg.gff[r, c("chr", "contig", "gfID", "chr.trn1")] )
    #       next;
    #     }
    #   
    #     if (current.ctg.gff[r,]$chr.trn1 ==  current.ctg.gff[r-1,]$chr.trn1) { ## only consider window size
    #       ## when this gene and the previous gene blong to the same block 
    #       ## their distance is less than the threshold,
    #       # block.vector <- c(block.vector, block.no)
    #       blockDF[nrow(blockDF),]$gfID <- current.ctg.gff[r,]$gfID
    #     }  
    #     else {
    #       # block.no <- block.no + 1
    #       # block.vector <- c(block.vector, block.no)
    #       blockDF <- rbind( blockDF, current.ctg.gff[r, c("chr", "contig", "gfID", "chr.trn1")] )
    #     }
    #   }
    #   
    # } ## end of loop through contigs
    # 
    # nGF <- head( c( blockDF$gfID[1], (shift(blockDF$gfID, 1L, type="lag") - blockDF$gfID)), -1)
    # 
    # blockDF <- cbind(blockDF[,-3], nGF)
    
    #######################
    ## choppiness analysis
    #######################
    for(c in sort(unique(ancestorGF.trn2$chr))) {
      df.chr <- blockDF[blockDF$chr==c,] 
      
      ## T_i = the number of different colours on each chromosome - 1
      t <- length(unique(df.chr$chr.trn1)) - 1
      
      ## R_i = the number of single-colour regions on each chromosome
      r <- length(df.chr$chr.trn1) 
      
      ## the number of stripes less than a certain threshold size 
      x <- 0 #nrow(df.chr[df.chr$len <= 2*blockLEN.threshold,])
      
      choppiness <- rbind(choppiness, data.frame(genome=paste0("ancestor",trn2), ancestor=paste0("ancestor",trn), 
                                                 chr=as.numeric(as.character(c)), 
                                                 t=as.numeric(as.character(t)), 
                                                 r=as.numeric(as.character(r)), 
                                                 x=as.numeric(as.character(x))))
      
    } ## end of looping through chromosomes
    
    
    
    # ## build ancestral genome karyotype
    # karyotype.trn <- buildAncestralKaryotype(genefamilyGFF.trn, clusterVector) #results.path, trn)
    
    
} ## end of looping through tree nodes

choppiness.sum <- as.data.table(aggregate(choppiness[,4:6], by=list(genomeID = choppiness$genome, Ancestor.map = choppiness$ancestor), FUN = "sum"))
choppiness.sum[, "R-T" := r-t]
choppiness.sum[, "R-X" := r-x]

choppiness.sum <- merge(genomeCoGeID[,1:2], choppiness.sum, by = "genomeID")
setorder(choppiness.sum, genomeID, Ancestor.map)



