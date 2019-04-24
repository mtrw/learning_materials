library(magrittr)
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(parallel)
library(bnstruct)


fill_missing_genotypes_with_NAs <- function( data ){
  vcf <- copy(data)
  sample_list <- sort(unique(vcf$sample))
  pos_list <- sort(unique(vcf$pos))
  chr_list  <- sort(unique(vcf$chr))
  filler <- as.data.table(  expand.grid(sample_list,pos_list,chr_list)  ) %>% setnames( (c("sample","pos","chr") ))
  vcf <- vcf[ filler , on=.(sample,pos,chr) ]
  vcf
}


generate_pca <- function(data,select_snps=unique(data$snp_id)){
  #High-performance processing/analyses with SNPRelate
  #Read in to SNPRelate's native gds format
  #make a SNP matrix
  
  vcf_as_snp_matrix <- copy(data) #pass in a vcf here. just the cols
  setorder(vcf_as_snp_matrix,snp_id,sample)
  vcf_as_snp_matrix <- dcast(vcf_as_snp_matrix , sample ~ snp_id , value.var = "genotype" ) #it's currently still a data table, just a bit "matrixy"
  
  sample_list <- vcf_as_snp_matrix$sample
  snp_list <- colnames(vcf_as_snp_matrix)[2:ncol(vcf_as_snp_matrix)]
  chr_list <- vcf[,.N,by=.(snp_id,chr)][snp_list , on=.(snp_id)]$chr
  pos_list <- vcf[,.N,by=.(snp_id,pos)][snp_list , on=.(snp_id)]$pos
  
  vcf_as_snp_matrix <- as.matrix( vcf_as_snp_matrix[,2:ncol(vcf_as_snp_matrix)] )
  
  # Create a gds file
  fname <- paste0(".genofile_temp_",sample.int(1e10,1),".gds")
  #unlink(fname)
  #showfile.gds(closeall=TRUE)
  snpgdsCreateGeno (
    gds.fn=fname,
    genmat = vcf_as_snp_matrix,
    sample.id = sample_list,
    snp.id = snp_list,
    snp.chromosome = chr_list,
    snp.position = pos_list,
    snpfirstdim=FALSE
  )
  #Open a connection to the GDS file
  genofile <- snpgdsOpen(fname)
  
  #Run a PCA to assess population structure
  pca <- snpgdsPCA( genofile ,  snp.id=select_snps,  autosome.only=FALSE  )
  
  pca_dt <- data.table(
    pca$eigenvect
  )
  pca_dt %>% setnames(paste0("PC_",1:ncol(pca$eigenvect)))
  pca_dt[ , sample:=pca$sample.id ]
  
  
  pca_varprop <- data.table(
    PC =1:ncol(pca$eigenvect),
    extracted_variance = pca$varprop
  )
  
  list( PCs = pca_dt , extracted_variance = pca_varprop )
}

generate_LD_matrix <- function(data,slide.max.bp=250000,ld.threshold=0.5){
  #High-performance processing/analyses with SNPRelate
  #Read in to SNPRelate's native gds format
  #make a SNP matrix
  
  vcf_as_snp_matrix <- copy(data) #pass in a vcf here. just the cols
  setorder(vcf_as_snp_matrix,snp_id,sample)
  vcf_as_snp_matrix <- dcast(vcf_as_snp_matrix , sample ~ snp_id , value.var = "genotype" ) #it's currently still a data table, just a bit "matrixy"
  
  sample_list <- vcf_as_snp_matrix$sample
  snp_list <- colnames(vcf_as_snp_matrix)[2:ncol(vcf_as_snp_matrix)]
  chr_list <- vcf[,.N,by=.(snp_id,chr)][snp_list , on=.(snp_id)]$chr
  pos_list <- vcf[,.N,by=.(snp_id,pos)][snp_list , on=.(snp_id)]$pos
  
  vcf_as_snp_matrix <- as.matrix( vcf_as_snp_matrix[,2:ncol(vcf_as_snp_matrix)] )
  # vcf_as_snp_matrix[1:10,1:10]
  # 
  # dim(vcf_as_snp_matrix)
  # length(sample_list)
  # length(snp_list)
  # length(chr_list)
  # length(pos_list)
  
  # Create a gds file
  fname <- paste0(".genofile_temp_",sample.int(1e10,1),".gds")
  #unlink(fname)
  #showfile.gds(closeall=TRUE)
  snpgdsCreateGeno (
    gds.fn=fname,
    genmat = vcf_as_snp_matrix,
    sample.id = sample_list,
    snp.id = snp_list,
    snp.chromosome = chr_list,
    snp.position = pos_list,
    snpfirstdim=FALSE
  )
  #Open a connection to the GDS file
  genofile <- snpgdsOpen(fname)
  #demo of how to access
  #read.gdsn(index.gdsn(genofile,"sample.id"))
  #to change would be add.gdsn(genofile,"sample.id", vector_of_new_ids ,replace = TRUE)
  
  #make an LD matrix
  LD_matrix <- snpgdsLDMat(genofile,method="corr",slide=0)
  #LD_matrix %>% str
  #LD_matrix$LD[1:5,1:5]
  LD_matrix$pos <- sub(".*_","",LD_matrix$snp.id) %>% as.numeric
  
  #convert to data.able for plotting WITH positions
  LD_dt <- data.table(
    pos1 = rep(LD_matrix$pos       ,times=length(LD_matrix$pos)),
    snp_id1 = rep(LD_matrix$snp.id ,times=length(LD_matrix$snp.id)),
    pos2 = rep(LD_matrix$pos       ,each=length(LD_matrix$pos)),
    snp_id2 = rep(LD_matrix$snp.id ,each=length(LD_matrix$snp.id)),
    LD_corr = LD_matrix$LD %>% as.numeric
  )
  LD_dt[ , snp_order1 := frank(pos1,ties.method="dense") ]
  LD_dt[ , snp_order2 := frank(pos2,ties.method="dense") ]
  LD_dt[ pos1==pos2 , LD_corr := NA ]
  LD_dt[ , LD_corr_abs := abs(LD_corr) ]
  setorder(LD_dt,LD_corr_abs)
  snps_to_keep <- snpgdsLDpruning( genofile   , method="dprime" , autosome.only=FALSE , slide.max.bp=slide.max.bp , ld.threshold=ld.threshold )
  snpgdsClose(genofile)
  unlink(fname)
  list( LD_data=LD_dt , good_SNPs=snps_to_keep )
}

test_imputation <- function(data,n=round( data[!is.na(genotype),.N]/100) , digits=0 , ... ){
  test_vcf <- copy(data)
  if(!is.null(test_vcf$imputed_genotype)) { test_vcf$imputed_genotype <- NULL }
  testset_indices <- which(!is.na(test_vcf$genotype))[sample.int( length( which(!is.na(test_vcf$genotype))) ,n) ]
  testset_dt <- test_vcf[   testset_indices  , .(sample,snp_id,genotype) ]
  test_vcf <- test_vcf[   testset_indices  , genotype := NA ][]
  test_vcf <- impute_genotypes_using_K_nearest_neighbour_method( test_vcf )
  cmp <- test_vcf[ , .(sample,snp_id,imputed_genotype) ][ testset_dt , on=.(sample,snp_id) ]
  cmp[ , paste0("Tested ",n," random SNPs. Imputation accuracy estimated: ",((sum(imputed_genotype==genotype)/.N) * 100) %>% round(digits=2),"%") ] %>% print
  return() %>% invisible
}


LD_heatmap_pos <- function(data,midpoint=0){
  ggplot(  data  ,  aes(x=pos1,y=pos2,colour=LD_corr_abs)  ) + geom_point(size=.01) + theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + scale_color_gradient2(
    midpoint=midpoint,
    low="blue",
    mid="white",
    high="red",
    space ="Lab"
  )
}


LD_heatmap_order <- function(data,midpoint=0){
  ggplot(  data  ,  aes(x=frank(snp_order1),y=frank(snp_order2),colour=LD_corr_abs)  ) + geom_point(size=.01) + theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + scale_color_gradient2(
    midpoint=midpoint,
    low="blue",
    mid="white",
    high="red",
    space ="Lab"
  )
}

#requires the VCF to have columns including sample, snp_id, genotype.
impute_genotypes_using_K_nearest_neighbour_method <- function(data,k=3){
  vcf <- copy(data) #pass in a vcf here.
  vcf_as_snp_matrix <- copy(data)[ , .( sample , snp_id , genotype )]
  vcf_as_snp_matrix <- dcast(vcf_as_snp_matrix , sample ~ snp_id , value.var = "genotype" )
  .invert_marker <- function(mvec){
    swap( vec = mvec , matches = c(0,2) , names = c(2,0) )
  }
  sample_names <- vcf_as_snp_matrix$sample
  snp_names <- colnames(vcf_as_snp_matrix)[2:ncol(vcf_as_snp_matrix)]
  
  vcf_as_snp_matrix <- vcf_as_snp_matrix[ , 2:ncol(vcf_as_snp_matrix) ] %>% as.matrix
  
  filler_inverted_haplotypes <- apply( vcf_as_snp_matrix , 2 , function(col) .invert_marker(col) )
  colnames( filler_inverted_haplotypes ) <- paste0( colnames(filler_inverted_haplotypes),"_i" )
  
  #just check ...
  dim(vcf_as_snp_matrix)
  dim(filler_inverted_haplotypes)
  vcf_as_snp_matrix_filled <- t(  cbind(vcf_as_snp_matrix,filler_inverted_haplotypes) )
  dim(vcf_as_snp_matrix_filled)
  
  imputed_mat <- knn.impute( vcf_as_snp_matrix_filled , k = k, cat.var = NULL , to.impute = 1:ncol(vcf_as_snp_matrix), using = 1:nrow(vcf_as_snp_matrix_filled))
  imputed_mat %<>% t
  
  dt <- as.data.table(imputed_mat)
  dt[ , sample := sample_names ]
  dt <- melt.data.table(dt , measure.vars = snp_names , id.vars="sample" , value.name="imputed_genotype" , variable.name="snp_id" )
  dt[ , imputed_genotype := as.integer(imputed_genotype) ] #kludge alert
  dt
  #join the imputed snps to the vcf data table
  dt[ vcf , on=.( sample , snp_id ) ]
}



read_annotation <- function(file){
  transcript_info <- fread( file , col.names = c("gene_id","transcript_id","type","gene_type","chr","strand","start","end","gene_description") )
  transcript_info[ , pos := pmean(start,end) ]
  transcript_info[ , length := end-start ]
  transcript_info <- transcript_info[ , .SD[length==max(length)][1] , by=.(gene_id) ]
  transcript_info[ chr=="chr1H" ][ , .(gene_id,chr,gene_description,pos,length) ]
}
