# this are sourced by the -makefiles in all sub directories
PATH := $(PATH):../bin:../../bin
PYTHONPATH := ../bin/pythonlib:../../bin/pythonlib
export PATH
export PYTHONPATH

# the dataset name, will be used for all directories
#DATASET=dm3dp3Example
#DATASET=dm3dp3Rnd1
#DATASET=dm3dp3MaxTrain2

# target genome identifier in UCSC-speak
TARGET=dm3

# query genome identifier in UCSC-speak
DATASET=dm3dp3Example
QUERY=dp3
QUERY_UP=Dp3

# UNCOMMENT THE FOLLOWING THREE LINES IF YOU WANT TO RUN
# ON D VIRILIS
#DATASET=dm3droVir3Example
#QUERY=droVir3
#QUERY_UP=DroVir3

# to make downstream scripts easier
QUERY_CAMELCASE=$(QUERY_UP)/$(QUERY)
# only process these target chromosomes 
CHROMS=chr2L,chr2R,chr3L,chr3R,chr4,chrX
# the maximum size of a fragment
MAXFRAGSIZE=15000
# minimum overlap of two fragments
MINOVERLAP=5000
# the maximum hole before we stop extending a fragment
MAXHOLE=2000

#username and hostname of compute cluster
#CLUSTERHOST=mhaeussl@templar.ls.manchester.ac.uk
CLUSTERHOST=danzek
# directory on compute cluster 
#CLUSTERDIR=/fs/nas8/scratch/mhaeussl/contemplateJobs/$(DATASET)
CLUSTERDIR=scratch/contemplateJobs/$(DATASET)

# for the spacer plots and the nora benchmark,
# we need two datasets to compare
DATASET_DP=dm3dp3Example
DATASET_DV=dm3droVir3Example

#DATASET=mm9Chr1
#TARGET=mm9
#QUERY=hg18
#QUERY_UP=Hg18
#QUERY_CAMELCASE=$(QUERY_UP)/$(QUERY)
#CHROMS=chr1
#MINOVERLAP=5000
#MAXHOLE=2000
#MAXFRAGSIZE=15000

#DATASET=mm9otherChrom
#TARGET=mm9
#QUERY=hg18
#QUERY_UP=Hg18
#QUERY_CAMELCASE=$(QUERY_UP)/$(QUERY)
#CHROMS=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY	
