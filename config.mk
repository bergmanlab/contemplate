# these are sourced by the makefiles in all sub directories
PATH := $(PATH):../bin:../../bin
PYTHONPATH := ../bin/pythonlib:../../bin/pythonlib
export PATH
export PYTHONPATH

# target genome identifier in UCSC-speak
TARGET=dm3

# query genome identifier in UCSC-speak
# the dataset name, will be used for all directories
DATASET=dm3dp3Example
QUERY=dp3
QUERY_UP=Dp3

# uncomment the following three lines if you want to run on virilis
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
CLUSTERHOST=danzek
# directory on compute cluster 
CLUSTERDIR=scratch/contemplateJobs/$(DATASET)

# for the spacer plots and the nora benchmark,
# we need two datasets to compare
DATASET_DP=dm3dp3Example
DATASET_DV=dm3droVir3Example