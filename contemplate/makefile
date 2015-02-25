ifneq "$(wildcard ../config.mk*)" ""
	include ../config.mk
endif

ifneq "$(wildcard config.mk*)" ""
	include config.mk
endif

.PHONY: test
export CLASSPATH:=javaSource/lib/nmica.jar:javaSource/lib/biojava.jar:javaSource/lib/bytecode.jar:javaSource/bin
JAVA=java -Xmx8000M 
CONTEMPLATE=$(JAVA) ModelRunner javaSource/run/enhSearch.xml 

# how many fragments should we treat per submitted jobs?
JOBCHUNKSIZE=100

# prefix for UCSC tracks when loading them into mirror, to have alpha/beta versions of tracks
#TRACKPREFIX=MaxTrain
#TRACKPREFIX=maxRnd1
#TRACKPREFIX=MaxTrain

# -- MAIN STEPS 
copyOver: toCluster

runJobs:
	@echo remove first line from makefile
	@echo then run make submitJobs

convert: fromCluster convertBlocks

# -- DETAILS
clean:
	rm -f $(DATASET)/aln/*
	rm -f $(DATASET)/bed/*
	rm -f $(DATASET)/blocks/*
	rm -f $(DATASET)/*

cleanOnCluster:
	rm -f aln/*
	rm -f bed/*
	rm -f blocks/*
	rm -f log/*

compile:
	@echo use javaSource to compile the java source code
	@echo As long you haven't modified anything, there's no need to compile

# copy all local files to compute cluster
copyFiles: copyMake
	#@echo don\'t forget to run cleanCluster before running this
	rsync -ap --progress --delete ../fastaFragments/$(DATASET)/fasta/ $(CLUSTERHOST):$(CLUSTERDIR)/fasta/
	rsync -avp --delete ./javaSource/bin/ $(CLUSTERHOST):$(CLUSTERDIR)/bin/
	rsync -avp --delete --cvs-exclude ./javaSource/lib/ $(CLUSTERHOST):$(CLUSTERDIR)/lib/
	scp javaSource/run/enhSearch.xml $(CLUSTERHOST):$(CLUSTERDIR)/
	scp ../bin/lstOp $(CLUSTERHOST):$(CLUSTERDIR)/bin/
	scp contemplateToPsl $(CLUSTERHOST):$(CLUSTERDIR)/bin/

copyMake:
# copy makefile to compute cluster
	ssh $(CLUSTERHOST) mkdir -p $(CLUSTERDIR)/bed $(CLUSTERDIR)/blocks $(CLUSTERDIR)/aln
	scp makefile $(CLUSTERHOST):$(CLUSTERDIR)/
	scp ../config.mk $(CLUSTERHOST):$(CLUSTERDIR)/

copyOthers:
	# make a backup copy of liftUp files to cluster
	rsync --progress ../fastaFragments/$(DATASET)/* $(CLUSTERHOST):$(CLUSTERDIR)
	
toCluster: copyFiles copyMake

# RUN THIS ON YOUR CLUSTER with "make submitJobs"
submitJobs:
	mkdir -p jobs/ log/
	rm -f jobs/*
	echo fasta/*.fa | tr ' ' '\n' | sed -e 's/fasta\///g' | sed -e 's/.fa$$//g'  > inputFiles
	echo bed/*.bed  | tr ' ' '\n' | sed -e 's/bed\///g'   | sed -e 's/.bed$$//g' > doneFiles
	bin/lstOp remove inputFiles doneFiles > todoFiles
	splitFile todoFiles $(JOBCHUNKSIZE) jobs/chunk
	for i in jobs/*; do make CHUNK=`basename $$i` submitOneChunk; done
	@echo When the jobs have completed:
	@echo check the log directory with grep Exception \* for crashed jobs
	@echo if everything is OK, run "make fromCluster" on your local machine to 
	@echo copy the results back from the cluster to your own machine

# this gets called by submitJobs on the cluster, do not run this yourself
submitOneChunk:
	qsub -V -N $(CHUNK) -o log/\$$JOB_ID-\$$JOB_NAME.out -j y -b y -shell y -cwd -V 'export CLASSPATH=lib/nmica.jar:lib/biojava.jar:lib/bytecode.jar:bin; for j in `cat jobs/$(CHUNK)`; do echo processing $$j; $(JAVA) ModelRunner enhSearch.xml fasta/$$j.fa -bedOut bed/$$j.bed -blockListOut blocks/$$j.blocks -fastaOut aln/$$j.fa; echo OK; done'

# run this after your jobs have finished to check for errors
checkLogs:
	grep -H "in thread .main." log/*

# run this LOCALLY to get your files back from the cluster
fromCluster:
	mkdir -p $(DATASET)/bed
	rsync -ap --progress --delete  $(CLUSTERHOST):$(CLUSTERDIR)/bed/ $(DATASET)/bed/
	rsync -ap --progress --delete  $(CLUSTERHOST):$(CLUSTERDIR)/blocks/ $(DATASET)/blocks/
	rsync -ap --progress --delete  $(CLUSTERHOST):$(CLUSTERDIR)/aln/ $(DATASET)/aln/

# CONVERT BLOCKS TO PSL
# (THIS IS ONLY A LITTLE BIT OF SCRIPTING - EVEN CASEY COULD DO THIS)

# size of fasta seq is needed for psl output, used to check for errors by pslToBed
faSizes:
	mkdir -p $(DATASET)/temp
	faSize -detailed ../fastaFragments/$(DATASET)/fasta/*.fa > $(DATASET)/temp/fastaSizes.tab

# need to translate filename to fragment name
fileToFrag:
	head -vn1 $(DATASET)/bed/*.bed | gawk '/^==/ {filename=$$2} /^frag/ {OFS="\t"; fragname=$$1; print filename, fragname}' > $(DATASET)/temp/fileToFrag.tab

toPsl:
	contemplateToPsl $(QUERY) $(TARGET) $(DATASET)/temp/fastaSizes.tab $(DATASET)/blocks > $(DATASET)/raw.psl

liftUpPslToBed:
	# convert psl to bed
	pslToBed $(DATASET)/raw.psl $(DATASET)/raw.bed
	bedSort $(DATASET)/raw.bed $(DATASET)/raw.bed
	# liftup psl to genome
	liftUp -type=.bed8 $(DATASET)/constrainedElements.bed ../fastaFragments/$(DATASET)/liftUp.lft error $(DATASET)/raw.bed

blocksToPslToBed: faSizes fileToFrag toPsl liftUpPslToBed

# keep only longest features and remove exons
filterBed-longest:
	# filter to get only longest chains
	bedOverlapOnlyLongest --joinBlocks --minBlock=0 $(DATASET)/constrainedElements.bed $(DATASET)/constrainedElements.longest.bed

filterBed-exonOverlaps:
	bedRemoveExonOverlapBlocks $(DATASET)/constrainedElements.longest.bed geneAnnotations/refseq.$(TARGET).bed $(DATASET)/constrainedElements.longest.noExons.bed 
	#bedFilter/bedOverlapOnlyLongest --minBlock=0 $(DATASET)/constrainedElements.noExons.bed $(DATASET)/constrainedElements.noExons.longest.bed

# use filtered bed file to filter psl file
filterPsl:
	# get names of elements that passed filter and keep only these in the psl
	cut -f4 $(DATASET)/constrainedElements.longest.bed > /tmp/filteredNames.lst
	# filter psl file with bed names to keep only longest psls
	lstOp filter $(DATASET)/raw.psl /tmp/filteredNames.lst > $(DATASET)/longest.psl
	# lift psl to genome
	liftUp -nohead $(DATASET)/longest.liftUp.psl ../fastaFragments/$(DATASET)/liftUp.lft error $(DATASET)/longest.psl
	# sort it
	#sort -k 14,14 -k 16,16n  $(DATASET)/longest.liftUp.psl > $(DATASET)/longest.liftUp.sorted.psl 
	
	wc -l $(DATASET)/*.psl
	wc -l $(DATASET)/*.bed

# keep only features with >=2, >=3, >=4 blocks
filterBedByBlockCount:
	bedOverlapOnlyLongest --minBlock 2 $(DATASET)/constrainedElements.longest.noExons.bed $(DATASET)/constrainedElements.filtered.minBlocks2.bed
	bedOverlapOnlyLongest --minBlock 3 $(DATASET)/constrainedElements.longest.noExons.bed $(DATASET)/constrainedElements.filtered.minBlocks3.bed
	bedOverlapOnlyLongest --minBlock 4 $(DATASET)/constrainedElements.longest.noExons.bed $(DATASET)/constrainedElements.filtered.minBlocks4.bed
	bedSort $(DATASET)/constrainedElements.filtered.minBlocks2.bed $(DATASET)/constrainedElements.filtered.minBlocks2.bed
	bedSort $(DATASET)/constrainedElements.filtered.minBlocks3.bed $(DATASET)/constrainedElements.filtered.minBlocks3.bed
	bedSort $(DATASET)/constrainedElements.filtered.minBlocks4.bed $(DATASET)/constrainedElements.filtered.minBlocks4.bed

convertBlocks: blocksToPslToBed filterBed-longest filterBed-exonOverlaps filterPsl filterBedByBlockCount
	echo Number of chains through the various stages of the pipeline
	wc -l $(DATASET)/*.bed | sort -rn

loadTracks:
# basic stats of contemplate predictions
stats:
	>$(DATASET)/stats.txt
	for i in $(DATASET)/*.minBlocks2.bed; do \
		echo share of genome covered by syntenic regions | tee -a $(DATASET)/stats.txt; \
		featureBits dm3 ../engstromPipeline/$(DATASET)/dm3.*.bed  2>&1 | tee -a $(DATASET)/stats.txt; \
		echo Number of chains in $$i                  | tee -a $(DATASET)/stats.txt; \
		wc -l $$i                                     | tee -a $(DATASET)/stats.txt; \
		echo Number of bases covered by chains in $$i | tee -a $(DATASET)/stats.txt; \
		featureBits dm3 $$i 2>&1                      | tee -a $(DATASET)/stats.txt; \
		echo Number of bases covered by blocks in $$i | tee -a $(DATASET)/stats.txt; \
		bedToExons $$i /tmp/exons.bed; featureBits dm3 /tmp/exons.bed 2>&1                      | tee -a $(DATASET)/stats.txt; \
		echo Average length of chains in $$i          | tee -a $(DATASET)/stats.txt; \
		cat $$i   | bedlen  | cut -f5 | avg           | tee -a $(DATASET)/stats.txt; \
		echo Number of features                       | tee -a $(DATASET)/stats.txt; \
		wc -l $(DATASET)/*.bed                        | tee -a $(DATASET)/stats.txt; \
		echo Number of chains longer than 1000 bp     | tee -a $(DATASET)/stats.txt; \
		cat $$i | bedlen | cut -f5 | gawk '($$1>1000) {print}' | wc -l | tee -a $(DATASET)/stats.txt; \
	done
	@echo Number of blocks: | tee -a $(DATASET)/stats.txt
	bedToExons $(DATASET)/constrainedElements.filtered.minBlocks2.bed stdout | wc -l | tee -a $(DATASET)/stats.txt
	@echo Average blocklen: | tee $(DATASET)/stats.txt
	bedToExons $(DATASET)/constrainedElements.filtered.minBlocks2.bed stdout | bedlen | cut -f5 | gawk '{sum+=$$1} END {print (sum)/NR}' | tee -a $(DATASET)/stats.txt
	@echo Average number of blocks per chain: | tee $(DATASET)/stats.txt
	cat $(DATASET)/constrainedElements.filtered.minBlocks2.bed | gawk '{ blocks = split($$12, a, ","); print blocks}' |  gawk '{sum+=$$1} END {print (sum)/NR}' | tee -a $(DATASET)/stats.txt
	@echo Total coverage of blocks | tee $(DATASET)/stats.txt
	bedToExons $(DATASET)/constrainedElements.filtered.minBlocks2.bed /tmp/exons.bed
	featureBits dm3 /tmp/exons.bed 2>&1 | tee -a $(DATASET)/stats.txt
	
# -- TESTING --
ucscFromBergman:
	rsync -avzp bergman:trackDb_genome/drosophila/dm3/ ./ucsc/

test:
	rm test/test.bed test/test.blk
	$(CONTEMPLATE) javaSource/run/tll.fa -bedOut test/test.bed -blockListOut test/test.blk
	cat test/test.blk
	cat test/test.bed

test-big:
	for i in ../fastaFragments/chainToFasta/temp/*.fa; \
		do echo processing $$i; $(CONTEMPLATE) $$i -bedOut test/`basename $$i .fa`.bed;  \
	done

# -- RERUNNING CRASHED JOBS

findCrashedJobs:
	for i in log/*.out; do cat $$i | gawk '/processing/ {id = $$2} /Exception/ {print id}'; done > jobs/crashed
rerunCrashedJobs:
	make CHUNK=crashed submitOneChunk

# --- Search bug, without any psl files

shortGaps-blocks:
	for i in $(DATASET)/blocks/*.blocks; do echo filename `basename $$i .blocks` | cat - $$i | gawk '/^filename/ {filename=$$2; next} $$1!="#" {if (tstart1!=0 && qstart1!=0) {qstart2=$$3; tstart2=$$1; qdiff=qstart2-qstart1; tdiff=tstart2-tstart1; if (tdiff<10) {print filename, line1; print filename, $$0;} }; line1=$$0; tstart1=$$1; qstart1=$$3;}'; done > $(DATASET)/max10.blocks

shortGaps-bed:
	bedChain $(DATASET)/constrainedElements.filtered.bed --filter -d 10 | cut -f1-3 > $(DATASET)/max10.bed
	wc -l $(DATASET)/max10.bed

gapsAlignments:
	 mkdir -p $(DATASET)/gapAln
	 cat $(DATASET)/max20.blocks | gawk '{file=$$1; start=$$3; getline; end=$$2+1; print "faFragWithGaps","$(DATASET)/aln/"file".fa", start, end, "stdout > $(DATASET)/gapAln/"file".fa"}' > $(DATASET)/gapSequences.sh
	 sh $(DATASET)/gapSequences.sh

contemplateHelp:
	$(CONTEMPLATE)

# lift up fragments to chromosome coordinates
#liftUpBeds.old:
	#liftUp -type=.bed8 $(DATASET)/constrainedElements.fromBedFiles.bed ../fastaFragments/$(DATASET)/liftUp.lft error $(DATASET)/bed/*.bed

