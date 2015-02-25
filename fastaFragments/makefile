include ../config.mk

run:
	rm -rf $(DATASET)
	mkdir -p $(DATASET)/fasta
	chainToFasta/chainToFasta -fragmentFile=$(DATASET)/fragments.bed -maxFragSize=$(MAXFRAGSIZE) -minOverlap=$(MINOVERLAP) -maxHole=$(MAXHOLE) ../engstromPipeline/$(TARGET).$(QUERY).synteny.chain ../engstromPipeline/fromUcsc/$(TARGET).2bit ../engstromPipeline/fromUcsc/$(QUERY).2bit $(DATASET)/fasta -liftT=$(DATASET)/liftUp.lft
	du -sch $(DATASET)/fasta

compile:
	cd chainToFasta; make

