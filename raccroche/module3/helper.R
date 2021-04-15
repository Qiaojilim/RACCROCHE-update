##### read in ancestral genome
readIn.ancestorGF <- function(contig.path, ws, trn, gf1, gf2, clusterVector)
{
  ### read in contigs from file for ancestor trn
  ## initialize contigLen list for the current parameters
  CTG.list <- c()
  ancestor <- data.frame()
  
  # input contigs
  contig.fname <- list.files(contig.path, pattern=paste0("ContigW",ws,"TreeNode",trn,"_",gf1,"_",gf2,".txt"))
  if (length(contig.fname)==0) {next}
  
  ## read in contig file line by line
  con <- file(file.path(contig.path, contig.fname), "r")
  
  while(TRUE){
    
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    ## read in and parse the current contig: contigN
    contigN <- as.numeric(strsplit(line, " ")[[1]][2])
    line = trimws(readLines(con, n = 1))
    #print(line)
    
    # the original contig read in from file
    contig_origin = abs(as.numeric(strsplit(line, " +")[[1]]))
    
    gf.ctg <- cbind(contig_origin, contigN, clusterVector[contigN+1]); colnames(gf.ctg) <- c("gfID", "contig", "chr")
    ancestor <- rbind(if(exists("ancestor")) ancestor, gf.ctg)
    
    CTG.list <- c(CTG.list, length(contig_origin))
    
    if (contigN==(nctg-1)) {break}
  } # end of reading file line by line
  # ancestor <- cbind.all(ancestor, CTG.list[1:nctg])
  
  close(con)
  
  return(ancestor)
}


##### get gene family lengths for ancestor
getGeneFamilyLength <- function (trn, ws, gf1, gf2, nctg, contigGFF.path)
{
  
  
  ### build ancestral genome karyotype for the first ancestor trn: chr, size
  ### using genes in the direct descendents to represent gene family
  gid.vector <- genomeCoGeID[genomeCoGeID$ancestor == trn,]$genomeID
  contigGFF=data.frame()
  for(gID in gid.vector) {
    
    ## read in contigGFF. Each row is one gene from ancestral contig, sorted by extant chr #, then positions in extant chr
    # Colnames: chr of extant genome, geneFamilyID, pos on extant chr, ancestral contig#, 
    # start position on extant chr, end position on extant chr, distance between two genes
    contigGFF.temp <- as.data.table(readIn.contigGFF(gID, ws, trn, gf1, gf2, nctg, contigGFF.path))
    contigGFF <- rbind(if(exists("contigGFF")) contigGFF, contigGFF.temp)
  }
  
  contigGFF[, "bp" := end-start]
  
  ## use the LONGEST gene in the gene family among ALL direct descendents to represent the length of the gene family
  genefamilyLen <- as.data.table(aggregate(contigGFF[,"bp"], by=list(gfID = contigGFF$geneFamilyID, contig = contigGFF$contig), FUN = "max"))
  
  return(genefamilyLen)
}


##### build ancestral karyotype
buildAncestralKaryotype <- function(genefamilyGFF, clusterVector) #results.path, trn)
{
  
  ## summarize contigs by sum of gene families
  karyotype.ctg <- as.data.table(aggregate(genefamilyGFF[,"bp"], by=list(contig = genefamilyGFF$contig), FUN = "sum"))
  missing.ctg <- as.data.table(cbind(which(!0:(nctg-1) %in% karyotype.ctg$contig) - 1, 0))
  colnames(missing.ctg)<-c("contig","bp")
  karyotype.ctg <- rbind(karyotype.ctg, missing.ctg)
  
  
  # ## read in ancestral chromosomes from the clustering results
  # clusterVector <- scan(file.path(results.path, "clustering", paste0("cluster_trn",trn,".txt")))
  karyotype.ctg <- cbind(karyotype.ctg, as.data.frame(clusterVector))
  
  ## build ancestral genome karyotype
  karyotype <- aggregate(karyotype.ctg[,"bp"], by=list(chr = karyotype.ctg$clusterVector), FUN = "sum")
  colnames(karyotype)<-c("chr","size")
  
  return (karyotype)
}


##### read in genomes to be analyzed
readIn.genomes <- function(data.path, fname){
  
  genomes.fname <- file.exists(file.path(data.path, fname))
  print(genomes.fname)
  if (length(genomes.fname)==0) {
    errorMessage <- paste0("Can't find input genomes file ",
                           "\nPlease check if the file exists under ", data.path)
    stop(errorMessage)
  }
  
  genomeCoGeID <- read.delim(file.path(data.path, fname), header=TRUE)
  return(genomeCoGeID)
}


##### read in karyotype chromosome
readIn.karyotype <- function(karyotype.path, gID){
  ## read in karyotype chromosome
  karyotype.fname <- list.files(karyotype.path, pattern=paste0("*_",gID,"_*"))
  print(karyotype.fname)
  if (length(karyotype.fname)==0) {
    errorMessage <- paste0("Can't find karyotype for ", genomeCoGeID[genomeCoGeID$genomeID == gID,]$genomeName,
                           "\nPlease check if the file(s) exist under ", karyotype.path)
    stop(errorMessage)
    }
  karyotype <- read.delim(file.path(karyotype.path, karyotype.fname), header=TRUE)
  
  return(karyotype)
}

##### read in gene features in contig
readIn.contigGFF <- function(gID, ws, trn, gf1, gf2, Nctg, contigGFF.path){
  
    contigGFF.fname <- list.files(contigGFF.path, pattern=paste0("*",gID,"_W",ws,"TreeNode",trn,"_",gf1,"_",gf2,".txt"))
    print(contigGFF.fname)
    if (length(contigGFF.fname)==0) {
      errorMessage <- paste0("Can't find contig GFF file for ancestor", genomeCoGeID[genomeCoGeID$genomeID == gID,]$ancestor,
                             " genome ", genomeCoGeID[genomeCoGeID$genomeID == gID,]$genomeName, 
                             "\nPlease check if the file(s) exist under ", contigGFF.path)
      stop(errorMessage)
      }
    
    contigGFF <- read.delim(file.path(contigGFF.path, contigGFF.fname), header=FALSE, row.names=NULL, stringsAsFactors = FALSE)
    colnames(contigGFF) <- c("chr",	"geneFamilyID",	"pos",	"contig", "start", "end", "geneName")
    
    ## clean up data: remove invalid chromosome numbers
    if(length(which(contigGFF$chr>max(karyotype$chr))) >0) {
      contigGFF <- contigGFF[ - which(contigGFF$chr>max(karyotype$chr)) , ]
    }
    # if(length(which(contigGFF$chr==0 )) >0) {  ## 33827_Chenopodium chromosomes start from 0
    #   contigGFF <- contigGFF[ - which(contigGFF$chr==0) , ]
    # }
    ## remove invalid contig number
    if(length(which(contigGFF$contig>=Nctg)) > 0) { contigGFF <- contigGFF [ - which(contigGFF$contig>=Nctg) , ]}

    
    ## process data and calculate distance between genes
    contigGFF$contig <- as.factor(contigGFF$contig)
    ## sort genes by chromosome number, then by their positions within chromosome
    contigGFF <- as.data.table(contigGFF)[order(contigGFF$chr,contigGFF$pos),]
    
    # Create column distance by subtracting column start by the value from the previous row of column end
    # this will result in the first gene in a chromosome having negative distance
    contigGFF[, distance := as.numeric(start - shift(end, 1L, type="lag"))]
    # change negative distance to NA, then replace all NAs with 0
    # a gene with distance 0 is the first gene in a chromosome.
    contigGFF<-contigGFF[distance < 0, distance := NA]
    contigGFF[is.na(contigGFF)] <- 0
 
    for(c in unique(contigGFF$chr)){
      contigGFF[which(contigGFF$chr == c), geneOrder := 1:nrow(contigGFF[which(contigGFF$chr == c),])]
    }
      
    return(contigGFF)   
    ## 9 cloumns in contig GFF: "chr",	"geneFamilyID",	"contig", "start", "end", "geneName", "distance", "geneOrder"
}


## Generate block data frame by merging only adjacent genes
## The first step initializes a syntenic block by merging two adjacent genes given a distance threshold DIS: merge two genes, g_1 and g_2, forming one ancestral syntenic block on G_i if g_1 and g_2 satisfy the following conditions:
## 1. g_1 and g_2 locate the same chromosome of G_i;
## 2. g_1 and g_2 are adjacent to each other; in other words, there could be a non-coding region but no other gene(s) between g_1 and g_2;
## 3. The distance between g_1 and g_2 must be less than or equal to the distance threshold DIS (i.e. DIS=1 MB)
## The second step extends the above identified ancestral syntenic block by merging flanking gene(s) into the block if the gene(s) satisfies the above three conditions. It stops extending the block if no flanking gene could be merged into the block.

generate.blockDF <- function(contigGFF, DIS.threshold) {
  ## merge genes in contigGFF that are within a DIS.threshold distance in extant genome into blocks
  ## columns: chr, start, end, contig --> make sure they are all numeric
  blockDF <- data.frame( chr=numeric(),  start=numeric(), end=numeric(), contig=numeric(), stringsAsFactors = FALSE) 
  
  ## asign block numbers
  # block.no <- 1
  # block.vector <- c(block.no)
  ## 8 cloumns in contig GFF: "chr",	"geneFamilyID",	"pos",	"contig", "start", "end", "geneName", "distance", "geneOrder", "geneOrder"
  blockDF <- rbind( blockDF, contigGFF[1, c("chr","start","end","contig")] )
  
  for (r in 2:nrow(contigGFF)) {
    if (contigGFF[r]$chr != contigGFF[r-1]$chr) { 
      ## the start of a new chromosome
      # block.no <- 1;  
      # block.vector <- c(block.vector, block.no); 
      blockDF <- rbind( blockDF, contigGFF[r, c("chr","start","end","contig")] )
      next;
    }
    
    if (contigGFF[r]$distance <= DIS.threshold && contigGFF[r]$contig == contigGFF[r-1]$contig) { ## ???? remove DIS.threshold?
      ## when this gene and the previous gene blong to the same block 
      ## their distance is less than the threshold,
      # block.vector <- c(block.vector, block.no)
      blockDF[nrow(blockDF),]$end <- contigGFF[r]$end
    }  
    else {
      # block.no <- block.no + 1
      # block.vector <- c(block.vector, block.no)
      blockDF <- rbind( blockDF, contigGFF[r, c("chr","start","end","contig")] )
    }
  }
  
  return(blockDF)
}

## Generate block data frame by allowing a window size to merge adjacent genes
## The first step initializes a syntenic block by merging two adjacent genes given a distance threshold DIS: merge two genes, g_1 and g_2, forming one ancestral syntenic block on G_i if g_1 and g_2 satisfy the following conditions:
## 1. g_1 and g_2 locate the same chromosome of G_i;
## 2. g_1 and g_2 are adjacent to each other with in WS; in other words, there could be less than WS gene(s) between g_1 and g_2;
## 3. The distance between every pair of adjacnet genes between g_1 and g_2 must be less than or equal to the distance threshold DIS (i.e. DIS=1 MB)
## The second step extends the above identified ancestral syntenic block by merging flanking gene(s) into the block if the gene(s) satisfies the above three conditions. It stops extending the block if no flanking gene could be merged into the block.

generate.blockDF.2 <- function(contigGFF, DIS.threshold, WS) {
  
  ## 9 cloumns in contig GFF now: "chr",	"geneFamilyID",	"pos",	"contig", "start", "end", "geneName", "distance", "geneOrder"
  
  ## merge genes in contigGFF that are within a DIS.threshold distance in extant genome into blocks
  ## columns: chr, start, end, contig --> make sure they are all numeric
  blockDF <- data.frame( chr=numeric(),  start=numeric(), end=numeric(), contig=numeric(), stringsAsFactors = FALSE) 
  
    
    for (current.ctg in unique(as.numeric(as.character(contigGFF$contig)))){
      ## loop through every contig
      # current.ctg <- 238
      # cat("\ncurrent.ctg=",current.ctg)
      current.ctg.gff <- contigGFF[which(contig==current.ctg),]
      current.ctg.gff[, distanceInCtg := as.numeric(start - shift(end, 1L, type="lag"))]
      # change alternating chr to NA, then replace all NAs with 0
      # a gene with distance 0 is the first gene in a chromosome.
      current.ctg.gff<-current.ctg.gff[as.numeric(chr - shift(chr, 1L, type="lag")) > 0, distanceInCtg := NA]
      current.ctg.gff[is.na(current.ctg.gff)] <- 0
      ## 10 cloumns in contig GFF: "chr",	"geneFamilyID",	"pos",	"contig", "start", "end", "geneName", "distance", "geneOrder, "distanceInCtg"
      
      ## start the first block from the first gene
      blockDF <- rbind( blockDF, current.ctg.gff[1, c("chr","start","end","contig")] )
      if(nrow(current.ctg.gff)<2) {next; }
      for (r in 2:nrow(current.ctg.gff)) {
        ## the start of a new chromosome, start a new block
        if (current.ctg.gff[r]$chr != current.ctg.gff[r-1]$chr) { 
          # block.no <- 1;  
          # block.vector <- c(block.vector, block.no); 
          blockDF <- rbind( blockDF, current.ctg.gff[r, c("chr","start","end","contig")] )
          next;
        }
  
        # if (current.ctg.gff[r]$distanceInCtg <= DIS.threshold && current.ctg.gff[r]$geneOrder < WS) { ## consider both distance and window size
        
        if ((current.ctg.gff[r]$geneOrder - current.ctg.gff[r-1]$geneOrder)  < WS) { ## only consider window size
          ## when this gene and the previous gene blong to the same block 
          ## their distance is less than the threshold,
          # block.vector <- c(block.vector, block.no)
          blockDF[nrow(blockDF),]$end <- current.ctg.gff[r]$end
        }  
        else {
          # block.no <- block.no + 1
          # block.vector <- c(block.vector, block.no)
          blockDF <- rbind( blockDF, current.ctg.gff[r, c("chr","start","end","contig")] )
        }
      }
      
    } ## end of loop through contigs
  
  blockDF<-blockDF[order(blockDF$chr,blockDF$start),]
  
  return(blockDF)
}


## Generate block data frame by allowing a window size to merge adjacent genes
## The first step initializes a syntenic block by merging two adjacent genes given a distance threshold DIS: merge two genes, g_1 and g_2, forming one ancestral syntenic block on G_i if g_1 and g_2 satisfy the following conditions:
## 1. g_1 and g_2 locate the same chromosome of G_i;
## The second step extends the above identified ancestral syntenic block by merging flanking gene(s) into the block if the gene(s) satisfies the above three conditions. It stops extending the block if no flanking gene could be merged into the block.

generate.blockDF.3 <- function(ancestorGF) 
{
  ancestorGF$gfID <- 1:length(ancestorGF$gfID)
  
  blockDF <- data.frame() 
  for (current.chr in unique(as.numeric(as.character(ancestorGF$chr)))){
    ## loop through every chr
    current.ctg.gff <- ancestorGF[which(ancestorGF$chr==current.chr),]
    
    ## start the first block from the first gene
    blockDF <- rbind( blockDF, current.ctg.gff[1,c("chr", "contig", "gfID", "chr.trn1")] )
    if(nrow(current.ctg.gff)<2) {next; }
    for (r in 2:nrow(current.ctg.gff)) {
      ## the start of a new chromosome, start a new block
      if (current.ctg.gff[r,]$chr != current.ctg.gff[r-1,]$chr) { 
        # block.no <- 1;  
        # block.vector <- c(block.vector, block.no); 
        blockDF <- rbind( blockDF, current.ctg.gff[r, c("chr", "contig", "gfID", "chr.trn1")] )
        next;
      }
      
      if (current.ctg.gff[r,]$chr.trn1 ==  current.ctg.gff[r-1,]$chr.trn1) { ## only consider window size
        ## when this gene and the previous gene blong to the same block 
        ## their distance is less than the threshold,
        # block.vector <- c(block.vector, block.no)
        blockDF[nrow(blockDF),]$gfID <- current.ctg.gff[r,]$gfID
      }  
      else {
        # block.no <- block.no + 1
        # block.vector <- c(block.vector, block.no)
        blockDF <- rbind( blockDF, current.ctg.gff[r, c("chr", "contig", "gfID", "chr.trn1")] )
      }
    }
    
  } ## end of loop through contigs
  
  nGF <- head( c( blockDF$gfID[1], (shift(blockDF$gfID, 1L, type="lag") - blockDF$gfID)), -1)
  
  blockDF <- cbind(blockDF[,-3], nGF)
  
  return (blockDF)
}

#### merge together blocks within distance threshold that are from the same ancestral chromosome cluster into the same block
#### blockDF ("chr", "start", "end", "contig", "ancestralChr")
mergeBlockDF <- function (blockDF, DIS.threshold) {
  blockDF<-blockDF[order(blockDF$chr,blockDF$start),]  # sort 
  
  ## columns: chr, start, end, ancestralChr --> make sure they are all numeric
  mergedDF <- data.frame( chr=numeric(),  start=numeric(), end=numeric(), ancestralChr=numeric(), stringsAsFactors = FALSE) 
  
  
  mergedDF <- rbind( mergedDF, blockDF[1, c("chr","start","end","ancestralChr")] )
  for(r in 2:nrow(blockDF)) {
    
    if (blockDF[r]$chr != blockDF[r-1]$chr) { 
      ## the start of a new chromosome
      mergedDF <- rbind( mergedDF, blockDF[r, c("chr","start","end","ancestralChr")] )
      next;
    }
    
    if ((blockDF[r]$start-blockDF[r-1]$end) <= DIS.threshold && (blockDF[r]$ancestralChr == blockDF[r-1]$ancestralChr)) { 
      ## when this block and the previous block blong to the same ancestral chromosome
      ## their distance is less than the threshold, merge them together
      mergedDF[nrow(mergedDF),]$end <- blockDF[r]$end
    }  
    else {
      # block.no <- block.no + 1
      # block.vector <- c(block.vector, block.no)
      mergedDF <- rbind( mergedDF, blockDF[r, c("chr","start","end","ancestralChr")] )
    } 
    
  }
  
  mergedDF <- mergedDF[order(mergedDF$chr,mergedDF$start),]
  
  return(mergedDF)
}



plotMyHeatMap <- function(data, title) {
  
  heatmap.2 (data,
             main = title, # heat map title
             lhei=c(3, 10), lwid=c(3,10), 
             cexRow=1,cexCol=1,margins=c(12,8),
             dendrogram="none",     # only draw a row dendrogram
             Rowv=FALSE, Colv=FALSE, 
             srtCol=0,   adjCol = c(0.5,1),
             col = brewer.pal(9, "Reds"), trace = "none",
             cellnote=round(data, 2),
             notecex=1.0,
             notecol="gray",
             na.color=par("bg"))
  
}
