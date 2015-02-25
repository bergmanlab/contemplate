include ../config.mk 

all: download run

compile:
	cd src/chainNetSynteny/
	make
	cd src/netChainSubsetRelaxed/
	make

download: getChains getGenomes prepareGenomes

run: runPipeline 

# ------

clean:
	rm -rf $(DATASET)

getQuery:
	rsync --progress rsync://hgdownload.cse.ucsc.edu/gbdb/$(QUERY)/$(QUERY).2bit fromUcsc/

getTarget:
	rsync --progress rsync://hgdownload.cse.ucsc.edu/gbdb/$(TARGET)/$(TARGET).2bit fromUcsc/

getGenomes: getQuery getTarget

getChains:
	mkdir -p $(QUERY)
	wget http://hgdownload.cse.ucsc.edu/goldenPath/$(TARGET)/vs$(QUERY_CAMELCASE)/$(TARGET).$(QUERY).all.chain.gz -O - | gunzip > $(QUERY)/dm3.$(QUERY).all.chain;

prepareGenomes:
	prepare_files.sh $(TARGET) $(QUERY)   # create reciprocal chain files

runPipeline: 
	mkdir -p $(DATASET)
	mksyn.sh $(DATASET) $(TARGET) $(QUERY) $(CHROMS) # generate .net and .synteny.bed files
	# add output in chain format
	chainFilter -t=$(CHROMS) $(QUERY)/$(TARGET).$(QUERY).recip-over.chain | chainNetSynteny -frag -longName -minGapShown=1 -minAli=10000 -breakingAli=10000 stdin $(QUERY)/$(TARGET).$(QUERY).net $(QUERY)/$(QUERY).$(TARGET).net $(DATASET)/$(TARGET).$(QUERY).synteny.chain /dev/null -outChain; \

#make TARGET=dm3 QUERY=$$1 QUERY_CAMELCASE=$$2 getQuery getChains prepareGenomes runPipeline; 
runAllFlies:
		for QUERYSTR in `cat allFlies.txt`; do \
			IFS="/"; \
			set $$QUERYSTR; \
			make TARGET=dm3 QUERY=$$1 QUERY_CAMELCASE=$$2 runPipeline; \
		done
