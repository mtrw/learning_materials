library(magrittr)
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(parallel)
library(devtools)

#apply a fun of two vars in every combo, and give the results as a matrix
#grid_ply(1:2,3:1,sum)
grid_ply <- function(rows,cols,FUN) {
  lapply( rows , function(row) lapply(cols , function(col) {
    FUN(row,col)
  }) %>% unlist ) %>% unlist %>%
    matrix(nrow=length(rows),dimnames=list(rows,cols))
}

#apply a fun of two vars in every combo, and give the results as a matrix
#mc_grid_ply(1:2,3:1,sum)
mc_grid_ply <- function(rows,cols,FUN,cores=25,...) {
  mclapply( mc.cores=cores , rows , function(row) lapply(cols , function(col) {
    FUN(row,col)
  }) %>% unlist ) %>% unlist %>%
    matrix(nrow=length(rows),dimnames=list(rows,cols))
}

#scale a list of values to between two points, proportionally spaced as they were originally
#rnorm(100) %>% scale_between(20,29) %>% pd
scale_between <- function(x,lower,upper){
  ( x - min(x) ) / (max(x)-min(x)) * (upper-lower) + lower
}

#easy way to see the spread of your values. plot the density.
#pd(c(NA,rnorm(500),NA))
pd <- function(x,...){
  x %>% density(na.rm=TRUE,...) %>% plot
}

#difference between two values but vectorises more intuitively than diff()
#difff(10,-286)
difff <- function(a,b){
  abs(a-b)
}

#z-transform values in a vec
#z_transform(x = ((runif(100))**2)+20  ) %>% pd
z_transform <- function(x){
  if ( sd(x)==0 ){
    return(rep(0,times=length(x)))
  }
  (mean(x)-x) / sd(x)
}

#number of unique entries in vector
#nu(c(1,1,2,3,4,4))
nu <-function(x){
  unique(x) %>% length
}

#fetch the upper triangular matrix in a vector without all that typing
#upper_tri( matrix(1:9,ncol=3) , diag=T )
upper_tri <- function(mat,diag=F){
  mat[upper.tri(mat,diag=diag)]
}

#paralle mean of two numbers. Great for BLAST search (midpoint of hits = pmean(start_vec,end_vec))
#pmean(1:5,2:6,5:1)
pmean <- function(...){
  invecs <- list(...)
  out <- rep(0, times=length(invecs[[1]]) )
  lapply(invecs,function(x) out <<- out+x )
  out/length(invecs)
}
pmean2 <- pmean #legacy reasons


#read a gff3 format file
read_gff <- function(fname){
  fread( paste("cat",fname,"| grep -v '^#'") ,col.names=c("seqname","source","feature","start","end","score","strand","frame","attribute"))
}

#as above but
#not the parts (eg CDSs ... just the whole genes)
read_gff_genes_pgsb <- function(fname){
  g <- fread( paste("cat",fname,"| grep -v '^#' | awk '$3==", '"gene"' ,"'") ,col.names=c("seqname","source","feature","start","end","score","strand","frame","attribute"))
  g[ , id := sub("^ID=(.*);.*","\\1",attribute) ][]
}
#seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
# source - name of the program that generated this feature, or the data source (database or project name)
# feature - feature type name, e.g. Gene, Variation, Similarity
# start - Start position of the feature, with sequence numbering starting at 1.
# end - End position of the feature, with sequence numbering starting at 1.
# score - A floating point value.
# strand - defined as + (forward) or - (reverse).
# frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

#for reading in my fave blast output
#as produiced by ~/bin/blastn_basic.zsh
read_blastn_basic <- function(fname){
  fread(fname,col.names=c("qseqid", "sseqid" , "slength" , "qlength" , "match_len" , "qstart" , "qend" , "sstart" , "send" , "pct_id" , "evalue" , "bitscore" ))
}

#length, because why spend your life typing "ength" all the time?
l <- function(x){
  length(x)
}

#simultaneously change the names of things to a particular thing if they match a particular string.
#name_by_match(vec=iris$Species,matches = c("set","sicol"),names = c("SETOSA","VERSICOLOR"))
name_by_match <- function(vec,matches,names){
  if(l(matches) != l(names)) { stop("Lengths of `matches` and `names` vectors don't match, you massive knob!")  }
  vec %<>% as.character()
  l_ply(1:l(matches) , function(n){
    vec[grep(matches[n],vec)] <<- names[n]
  })
  vec
}

#simultaneously change the names of things to a particular thing if they match (EXACTLY!) a particular string.
#swap(vec=iris$Species,matches = c("virginica"),names = c("XXXXX"))

swap <- function(vec,matches,names,na.replacement=NA){
  orig_vec <- vec
  #if(sum(! matches %in% names ) > 0 ) { stop("Couldn't find all matches in names") }
  if(l(matches) != l(names)) { stop("Lengths of `matches` and `names` vectors don't match, you old bison!") }
  if(is.factor(vec)) { levels(vec) <- c(levels(vec),names,na.replacement) }
  vec[is.na(orig_vec)] <- na.replacement
  l_ply( 1:l(matches) , function(n){
    vec[orig_vec==matches[n]] <<- names[n]
  })
  vec
}

#length of a vector
#euclidean_dist(c(3,4))
euclidean_dist <-function(x){
  ( sum( x**2 ) )**(1/2)
}

#as the title says. relative probs should be positive (duh) and same order as events (duh)
#t <- random_draws_from_any_discreet_distribution(n=50000,events=LETTERS[1:5],relative_probs=1:5) %>% table; t / max(t) * 5 ; rm(t)
random_draws_from_any_discreet_distribution <- function(n=1,events,relative_probs){
  lapply( 1:n,   function(x){
    events[ which( runif(1)*sum(relative_probs) < (cumsum(relative_probs)) )[1] ]    
  }) %>% unlist
}

#n random (non-repeated!) rows of a data frame
#sample_df(iris,20)
sample_df <- function(df,n=10){
  df[ sample.int(nrow(df),n) , ]
}

p <- function(x){
  print(x)
}

#random_matrix( 5 , 5 , symmetric = T, all_pos=F )
random_matrix <- function(m,n=m,symmetric=F,all_pos=T){
  if(all_pos){
    M <- matrix(runif(m*n),nrow=m)
  } else {
    M <- matrix(runif(m*n)-0.5,nrow=m)
  }
  if(symmetric){
    M <- M+t(M)
  }
  M
}

#the middle of the range of a dataset
#midpoint(c(0,4,4,4,4,4,4,4,4,4,10,NA,10))
midpoint <- function(x){
  mean( c( min(x, na.rm=T) , max(x, na.rm=T) ) )
}

#the expected value this many SDs from the mean of distribution (assumed normal, taken from a vector)
#SDs_from_mean(rnorm(10000),2) #around 2
SDs_from_mean <- function(x,SDs){
  mean(x) + SDs*sd(x,na.rm=T)
}

#echo to the global environment. good warning messager. still doesn't work in mclappy, hashtag doh
#ce("beans ",list("hello "," goodbye")," whatever")
ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }

#distance to the nearest element of each element
nearest_in_set <- function(x) { 
  if(length(x)==1) {return(0)} else  { sapply(1:length(x), function(y) { min (  abs( x[y] - x[setdiff(1:length(x),y)] ) ) } ) } 
}
#nearest_in_set(1:5); nearest_in_set((1:5)**2)



#sigmoid curves
sigmoid <- function(x,L=1,k=1,xo=0){
  L/(1+exp(-k*(x-xo)))
}
#sigmoid(-10:10)

interpolate <- function( x , y , points_x=NULL , points_y=NULL , plot=TRUE ){
  if((is.null(points_x) & is.null(points_y)) | !(is.null(points_y) & !is.null(points_x))){
    stop("Must enter one of points_x and points_y, not both.")
  }
  if(!is.null(points_y)){
    x -> t
    x <- y
    y <- t
    rm(t)
    points_x <- points_y
    rm(points_y)
  }
  
  len <- l(x)
  if(len!=l(y) | len<2){
    stop("x and y must be equal lengths of at least 2.")
  }
  
  o_x <- order(x)
  o_y <- order(y)
  if((l(o_x != o_y))<0){
    stop("x and y must be monotonic.")
  }
  x <- x[o_x]
  y <- y[o_x]
  
  out <- lapply(points_x,function(pt_x){
    #it hits an exact point
    if(l(xeq <- which(x==pt_x))==1){
      y[xeq]
    } else if(l(xeq)>1) {
      warning("WARNING: Falls within run of identical x-values. Estimating by mean(y).")
      mean(y[xeq])
    } else if(pt_x>max(x)){
        warning("WARNING: Extrapolating at high end.")
      if((x[len]-x[len-1])==0){
        warning("WARNING: Falls beyond run of identical x-values. Estimating by mean(y).")
        mean(y[x==x[len]])
      } else {
        y[len] + (pt_x-x[len])*( (y[len]-y[len-1])/(x[len]-x[len-1]) )
      }
    } else if(pt_x<min(x)){
      warning("WARNING: Extrapolating based on last two points at low end.")
      if((x[2]-x[1])==0){
        warning("WARNING: Falls before run of identical x-values. Estimating by mean(y).")
        mean(y[x==x[1]])
      } else {
          y[1] + (pt_x-x[1])*( (y[2]-y[1])/(x[2]-x[1]) )
        }
    } else {
      nearest_under_idx <- last(which(x<pt_x))
      y[nearest_under_idx] + (pt_x-x[nearest_under_idx])*( (y[nearest_under_idx+1]-y[nearest_under_idx])/(x[nearest_under_idx+1]-x[nearest_under_idx]) )
    }
  }) %>% unlist
  if(plot){
    plot(x,y,xlim=c(min(points_x,x),max(points_x,x)))
    points(points_x,out,col="red")
  }
  out
}

# X <- -10:100
# Y <- sigmoid(X,L=100,xo=20)-10
# pts_x <- (-20:200)+.4
# out_x <- pts_x

# out_y <- interpolate(x = X , y = Y, points_x=pts_x , plot = TRUE )
