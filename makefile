include config.mk

VERSION=0.2

SUPPL=paper/suppl
TABLES=paper/tables
FIG=paper/figures

all: suppl tables figures

suppl: suppl1 suppl2 suppl3 suppl4 suppl5 suppl6 suppl7
tables: table2 table3
figures: fig2-dotplot figx-histograms

init:
	mkdir -p $(FIG)
	mkdir -p $(TABLES)
	mkdir -p $(SUPPL)

clean:
	rm -f $(SUPPL)/*
	rm -f $(TABLES)/*

# create supplementary files and write to SUPPL
suppl1:
	zip -j $(SUPPL)/suppl1.engstromOutput.zip engstromPipeline/dm3.dp3.synteny.chain engstromPipeline/dm3.droVir3.synteny.chain

suppl2:
	cp contemplate/$(DATASET_DP)/constrainedElements.bed /tmp/contemplate.dm3dp3.bed
	cp contemplate/$(DATASET_DP)/constrainedElements.longest.bed /tmp/contemplate.dm3dp3.longest.bed
	cp contemplate/$(DATASET_DP)/constrainedElements.longest.noExons.bed /tmp/contemplate.dm3dp3.longest.nonCoding.bed
	cp contemplate/$(DATASET_DP)/constrainedElements.filtered.minBlocks2.bed /tmp/contemplate.dm3dp3.longest.nonCoding.atLeast2Blocks.bed
	zip -j $(SUPPL)/suppl2.contemplateOutput.dm3dp3.zip /tmp/contemplate.dm3dp3.bed /tmp/contemplate.dm3dp3.longest.bed /tmp/contemplate.dm3dp3.longest.nonCoding.bed /tmp/contemplate.dm3dp3.longest.nonCoding.atLeast2Blocks.bed 

suppl3:
	cp contemplate/$(DATASET_DV)/constrainedElements.bed /tmp/contemplate.dm3droVir3.bed
	cp contemplate/$(DATASET_DV)/constrainedElements.longest.bed /tmp/contemplate.dm3droVir3.longest.bed
	cp contemplate/$(DATASET_DV)/constrainedElements.longest.noExons.bed /tmp/contemplate.dm3droVir3.longest.nonCoding.bed
	cp contemplate/$(DATASET_DV)/constrainedElements.filtered.minBlocks2.bed /tmp/contemplate.dm3droVir3.longest.nonCoding.atLeast2Blocks.bed
	zip -j $(SUPPL)/suppl3.contemplateOutput.dm3droVir3.zip /tmp/contemplate.dm3droVir3.bed /tmp/contemplate.dm3droVir3.longest.bed /tmp/contemplate.dm3droVir3.longest.nonCoding.bed /tmp/contemplate.dm3droVir3.longest.nonCoding.atLeast2Blocks.bed 

suppl4:
	zip -j $(SUPPL)/suppl4.PierstorffChanData.zip comparePierstorff/$(DATASET_DP)/trim/*.bed comparePierstorff/bedBenchmark comparePierstorff/README.supplData.txt

suppl5:
	cp cbpAnnotation/cbpPeaks/modencodeHaeusslerComparison.xls $(SUPPL)/suppl5.cbp.modencodeVsThisStudy.xls

suppl6:
	mysql redfly -e 'select rc_id as redfly_reporterConstructId, state AS redfly_State, version as redfly_version, ReporterConstruct.name as redfly_reporterConstructName, ExpressionTerm.flybase_id as redfly_flybaseAnatomyId, ExpressionTerm.term as redfly_flybaseAnatomyTerm, stage AS mapped_stage from ReporterConstruct JOIN RC_has_ExprTerm USING (rc_id) JOIN ExpressionTerm USING (term_id) JOIN anatomyToStage USING (flybase_id) WHERE is_crm=1 ORDER BY rc_id;' > /tmp/redflyTable.tsv
	cp /tmp/redflyTable.tsv $(SUPPL)/suppl6.redflyStageAssignment.tsv
	ls -la $(SUPPL)

table2:
	mkdir -p $(TABLES)
	cp comparePierstorff/$(DATASET_DP)/results.tsv $(TABLES)/table2.noraBenchmark.tsv

table3:
	cp cbpAnnotation/stageBenchmark.xls $(TABLES)/table3.cbpStageBenchmark.tsv

fig2-dotplot:
	dottup -bsequence contemplate/training/loci/dm3/dpp_disk.fa -asequence contemplate/training/loci/droVir2/dpp_disk.fa -wordsize 15 -graph ps -gtitle "Dotplot of DPP locus" -gsubtitle "wordSize 15" -gytitle "D. melanogaster" -gxtitle "D. pseudoobscura" -gdirectory $(FIG) -send1 6000 -goutfile fig2.dv.dotplot
	dottup -bsequence contemplate/training/loci/dm3/dpp_disk.fa -asequence contemplate/training/loci/droVir2/dpp_disk.fa -wordsize 15 -graph ps -gtitle "Dotplot of DPP locus" -gsubtitle "wordSize 15" -gytitle "D. melanogaster" -gxtitle "D. pseudoobscura" -gdirectory $(FIG) -send1 6000 -goutfile fig2.dp.dotplot
	#evince $(FIG)/*.ps

figx-histograms:
	cp plotHistograms/lenHistograms.pdf $(FIG)/figx.histograms.pdf

# create source distribution package on max.smith
src:
	cd ..; tar cvfz contemplateGenome/contemplate-src.tar.gz --no-recursion --dereference `cat contemplateGenome/exportfiles`
	sudo mv contemplate-src.tar.gz /var/www/default/download/contemplate/
	scp /var/www/default/download/contemplate/contemplate-src.tar.gz  haeussle,contemplate-hmm@frs.sourceforge.net:/home/frs/project/c/co/contemplate-hmm/contemplate-src.$(VERSION).tgz

stats:
	echo PSEUDOOBSCURA STATS:
	cd contemplate; make -s stats DATASET=$(DATASET_DP)
	echo 
	echo
	echo VIRILIS STATS
	cd contemplate; make -s stats DATASET=$(DATASET_DV) 

# update current tree from max.smith
update:
	@if [ `hostname` != "max.smith.manchester.ac.uk" ]; then  \
		cd ..;  \
		rm contemplate-src.tar.gz -f;  \
		wget http://max.smith.man.ac.uk/download/contemplate/contemplate-src.tar.gz;  \
		tar xvfz contemplate-src.tar.gz; cd contemplateGenome;  \
	else  \
		echo stupid, max!; \
	fi

# LOAD TRACK 
# this loads tracks into a local UCSC browser 
# LOCAL MIRROR LOADING
loadTracks:
	# load stage assignment into local UCSC mirror
	grep -v track engstromPipeline/$(DATASET)/$(TARGET).$(QUERY).synteny.bed > /tmp/synteny.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)EngstromRegions /tmp/synteny.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)Frag fastaFragments/$(DATASET)/fragments.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)Raw contemplate/$(DATASET)/constrainedElements.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)Longest contemplate/$(DATASET)/constrainedElements.longest.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)NoExons contemplate/$(DATASET)/constrainedElements.longest.noExons.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)MinBlocks2 contemplate/$(DATASET)/constrainedElements.filtered.minBlocks2.bed
	hgLoadBed $(TARGET) cont$(QUERY_UP)CbpStages cbpAnnotation/$(DATASET)/contemplate.p300Stage.bed


# add public key from this machine to your clusterhost, no need to type passwords anymore
clusterPublicKey:
	ssh $(CLUSTERHOST) "mkdir -p ~/.ssh && cat - >> ~/.ssh/authorized_keys" < ~/.ssh/id_dsa.pub

# tar for casey
tar:
	find . -print | egrep 'notUsed|genomeBrowser' > /tmp/exclude.txt 
	tar cfzh /tmp/contemplateGenome.$(VERSION).tar . --exclude-from /tmp/exclude.txt 

