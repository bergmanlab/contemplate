include $(JKSRC)/inc/common.mk

L = -lm $(MYSQLLIBS) -lssl -lcrypto
MYLIBDIR = $(JKSRC)/lib/$(MACHTYPE)
MYLIBS =  $(MYLIBDIR)/jkhgap.a $(MYLIBDIR)/jkweb.a
HG_INC = -I$(JKSRC)/inc -I$(JKSRC)/hg/inc

O = chainToFasta.o

chainToFasta: $O $(MYLIBS)
	#${CC} ${COPT} ${CFLAGS} -o $(HOME)/bin/chainToFasta $O $(MYLIBS) $L
	#${STRIP} $(HOME)/bin/chainToFasta${EXE}
	${CC} ${COPT} ${CFLAGS} -o chainToFasta $O $(MYLIBS) $L
	${STRIP} chainToFasta

install:
	mv chainToFasta $(HOME)/bin/chainToFasta 

clean:
	rm -f $(O) 
	rm -rf temp/*
	rm -f test/*

test: chainToFasta
	rm -f temp/*
	chainToFasta -fragmentFile=temp/fragments.bed -maxFragSize=10000 -minOverlap=1000 test/test.chain ../../engstromPipeline/fromUcsc/dm3.2bit ../../engstromPipeline/fromUcsc/droVir3.2bit temp -lift=temp/test.lft -fragmentFile=test/test.bed -verbose=4

# negative strand
test2: chainToFasta
	rm -f temp/*
	chainToFasta -fragmentFile=temp/fragments.bed -maxFragSize=10000 -minOverlap=1000 test/test2.chain ../../engstromPipeline/fromUcsc/dm3.2bit ../../engstromPipeline/fromUcsc/droVir3.2bit temp -lift=temp/test.lft -fragmentFile=temp/test.bed -verbose=4

# pos strand
test3: chainToFasta
	rm -f temp/*
	chainToFasta -fragmentFile=temp/fragments.bed -maxFragSize=10000 -minOverlap=1000 test/test3.chain ../../engstromPipeline/fromUcsc/dm3.2bit ../../engstromPipeline/fromUcsc/droVir3.2bit temp -lift=temp/test.lft -fragmentFile=temp/test.bed -verbose=4

test-quiet: chainToFasta
	rm -f temp/*
	chainToFasta -fragmentFile=temp/fragments.bed -maxFragSize=10000 -minOverlap=1000 test.chain ../../engstromPipeline/fromUcsc/dm3.2bit ../../engstromPipeline/fromUcsc/droVir3.2bit temp -lift=temp/test.lft -fragmentFile=test/test.bed 

e:
	vim chainToFasta.c

