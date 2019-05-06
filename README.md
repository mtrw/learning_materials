# learning_materials
Variously lifted from my workshops and personal learning materials

This is material designed for a course I ran this year. Out of the course context, it may take a little fiddling to get it running, but I'll happily provide help and the required materials on request.

It is designed to be run in a directory called ~/workspace

This directory contains (at least) the following files/paths:

> ~/workspace/bin/course_functions.R

> > Included in the repo

> ~/workspace/bin/tim_functions.R

> > Included in the repo

> ~/workspace/raw_data/barley_annotation.tsv

> > Available on request. A simplified gff-type file. Some example columns:

> > HORVU1Hr1G000010        HORVU1Hr1G000010.1      transcript      canonical_gene  chr1H   +       41811   45049   RING-H2_finger_protein_2B
> > HORVU1Hr1G000010        HORVU1Hr1G000010.1      exon    canonical_gene  chr1H   +       41811   42213   RING-H2_finger_protein_2B
> > HORVU1Hr1G000010        HORVU1Hr1G000010.1      exon    canonical_gene  chr1H   +       42300   42338   RING-H2_finger_protein_2B

> ~/workspace/raw_data/barley_panel_sample_metadata.csv

> > Meta-data for the above. Also included in the repo, but you'll need to provide your own to use your own VCF file. Some example columns (incl. header):

> > cultivation_status,domestication_status,sample,annual_growth_habit,row_type,bridge_region
> > cultivar,domesticated,Sample_7264,spring,6-rowed,Northern_Europe
> > cultivar,domesticated,Sample_8700,spring,2-rowed,Central_Asia

> ~/workspace/calls/barley_panel_snpcalls.tsv

> > A "tidied" VCF file. It's big but I can find a way to get it to you on request. To make your own though, you can adapt this AWK script which will convert a VCF from STDIN:

\#!/bin/awk -f

BEGIN{
 OFS="\t"
 h["1/1"]=2
 h["0/1"]=1
 h["0/0"]=0
}

/^#CHROM/ && !header {
 for(i = 1; i <= NF; i++)
  s[i]=$i
 header=1
}

/INDEL/ || $4 == "N" || $5 ~ /,/ || /^#/ {
 next
}


{
 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  if(a[3] > 0)
   print $1, $2, $4, $5, $6, s[i], $8, h[a[1]], a[5], a[3], a[4]
 }
}
