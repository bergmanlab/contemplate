
/*
 * Written by M Haeussler, 2010 for the contemplate package.
 * You can do whatever you want with this source code
 * Uses the UCSC Genome Browser library (http://hgdownload.cse.ucsc.edu/admin/jksrc.zip).
 */

#include "common.h"
#include "linefile.h"
#include "localmem.h"
#include "hash.h"
#include "options.h"
#include "dnaseq.h"
#include "dnautil.h"
#include "rbTree.h"
#include "chain.h"
#include "chainNet.h"
#include "portable.h"
#include "bed.h"
#include "twoBit.h"
#include "fa.h"


struct optionSpec optionSpecs[] =
{
    {"maxFragSize", OPTION_INT },
    {"fragmentFile", OPTION_STRING },
    {"minOverlap", OPTION_INT },
    {"maxHole", OPTION_INT },
    {"liftT", OPTION_STRING },
    {"liftQ", OPTION_STRING },
    {NULL, 0}
};

int maxFragSize = 15000;	/* maximum size of fasta fragments, 0 to deactivate splitting */
int minOverlap = 2000;	/* minimum size of overlap (if last fragment was long enough) */
int maxHole = 2500;     /* maximum size of gaps when going backwards */
char *fragmentFilename = 0;	/* filename of bedfile with fragments, 0 to deactivate */
char *liftSpecTFilename = 0;	/* filename of target liftUp file, 0 to deactivate*/
char *liftSpecQFilename = 0;	/* filename of query liftUp file, 0 to deactivate*/

void usage()
/* Explain usage and exit. */
{
    errAbort(
	"chainToFasta - extract fasta sequences for complete chains from 2bit files\n"
    "               create one fasta file per chain,\n"
    "               sequences can then be re-aligned with any other software\n"
	"usage:\n"
	"   chainToFasta in.chain target.2bit query.2bit outputDir\n"
	"where:\n"
	"   in.chain is the chain file \n"
	"   target.2bit is the target genome 2bit file\n"
	"   query.2bit is the query genome 2bit file\n"
    "   outputDir is the output directory for the fasta files\n"
    "   (one file will be created per chain, filename=chainId)\n"

	"options:\n"
	"   -maxFragSize=N         - Break fasta files if they are longer than x bp.\n"
    "                            Will only break at gaps. 0 to deactivate.\n"
    "                            Default %d.\n"
	"   -fragmentFile=bedFile  - Write fragments to a bed file for inspection\n"
	"   -liftT=liftUpFile       - Generate a liftUp-file to map\n"
    "                             target fragments back to the genome\n"
	"   -liftQ=liftUpFile       - Generate a liftUp-file to map query fragments\n"
    "                             back to the genome\n"
    "   -minOverlap=INT        - minimum Overlap of two tiles, by\n"
    "                            default %d\n"
	"   -verbose=INT           - generate debug info\n",
    maxFragSize, minOverlap);
}

void writeToBed(char* bedFilename, char* chrom, int start, int end, char* name, char* qName, int qSize, int qStart, int qEnd, char qStrand)
/* write to bedFilename */
{
    static FILE* bedFile = NULL;

    if (bedFilename==NULL) 
        return;

    if (bedFile==NULL) 
        bedFile = mustOpen(bedFilename, "w");

    char* buf = needMem(1000); 
    if (qStrand=='-') {
        int qStart2 = qSize - qEnd;
        qEnd   = qSize - qStart;
        qStart = qStart2;
    }
    sprintf(buf, "%s\t%d\t%d\t%s|%s:%d-%d(%c)\n", chrom, start, end, name, qName, qStart, qEnd, qStrand);
    mustWrite(bedFile, buf, strlen(buf));
    freeMem(buf);
}


char *makeLocString(char* db, char* chrom, int start, int end, char strand, int count)
/* Make a string of format "dm3-chr1-2000-2100" */
/* Make a string of format "dm3.fragmentXXX */
{
    char *s2;

	assert(start < 999999999); /* make sure pos is <= 9 digits long (i.e. < 1 gigabase) */
	assert(end < 999999999);   /* make sure pos is <= 9 digits long (i.e. < 1 gigabase) */
    s2 = needMem(strlen(db) + 1 + strlen(chrom) + 1 + 9 + 1 + 9 + 1+ 1+ 1); 
    //sprintf(s2, "%s %s:%d-%d %c", db, chrom, start, end, strand);
    //sprintf(s2, "%s-%s-%d-%d", db, chrom, start, end);
    strtok(db, ".");
    sprintf(s2, "fragment%d.%s", count, db);

    return s2;
}

struct chainFragList
{
    struct chainFragList *next;
    int tStart, tEnd, qStart, qEnd;
};


struct lastBlockStartList
{
    struct lastBlockStartList *next;
    struct lastBlockStartList *previous;
    int tStart, tEnd, qStart, qEnd;
};

struct chainFragList *splitChain(struct chain* chain, int maxFragSize) 
/* split chain into pieces of maxFragSize and return as a list of positions */
{
    struct cBlock *b;
    
    // a double linked list of all block starts, we need this to be able to walk backwards
    struct lastBlockStartList *lastBlockStarts = needMem(sizeof(struct lastBlockStartList));
    struct lastBlockStartList *blockStart1 = lastBlockStarts;
    lastBlockStarts->previous = NULL;

    // create a new fragment
    struct chainFragList *frag1 = needMem(sizeof(struct chainFragList));

    // frag is a pointer to the current fragment
    struct chainFragList *frag= frag1;
    frag->tStart=-1;

    // b is a pointer to the current block
    for (b = chain->blockList; b != NULL; b = b->next) 
    {
        verbose(4, "current block:    %d %d %d %d\n", b->tStart, b->tEnd, b->qStart, b->qEnd);
        verbose(4, "current fragment: %d %d %d %d\n", frag->tStart, frag->tEnd, frag->qStart, frag->qEnd);
        
        // if there is no data in the current fragment yet: store the first
        // start/end coordinates and stop immediately
        if (frag->tStart==-1) 
        {
            verbose(4, "first block, init'ing fragment and lastStarts\n");
            frag->tStart = b->tStart;
            frag->qStart = b->qStart;
            frag->tEnd = b->tEnd;
            frag->qEnd = b->qEnd;
            lastBlockStarts->tStart = b->tStart;
            lastBlockStarts->qStart = b->qStart;
            lastBlockStarts->tEnd = b->tEnd;
            lastBlockStarts->qEnd = b->qEnd;
            continue;
        }

        // append current block start to lastBlockStarts
        lastBlockStarts->next  = needMem(sizeof(struct lastBlockStartList));
        lastBlockStarts->next->previous = lastBlockStarts;
        lastBlockStarts = lastBlockStarts->next;
        lastBlockStarts->tStart = b->tStart;
        lastBlockStarts->qStart = b->qStart;
        lastBlockStarts->tEnd = b->tEnd;
        lastBlockStarts->qEnd = b->qEnd;
        verbose(4, "start appended to lastBlockStarts, tStart=%d\n", b->tStart);

        verbose(4, "target difference %d\n", b->tEnd - frag->tStart);
        verbose(4, "query difference %d\n", (b->qEnd - frag->qStart));

        if ( (maxFragSize==0) || ( ((b->tEnd - frag->tStart) < maxFragSize) && ( (b->qEnd - frag->qStart) < maxFragSize) ) )
        {
            // next fragment is not going to be too long for maxFragSize
            verbose(4, "next frag short enough -> extend fragment end, target diff %d\n", (b->tEnd - frag->tStart));
            verbose(4, "next frag short enough -> extend fragment end, query diff %d\n", (b->qEnd - frag->qStart));
            frag->tEnd = b->tEnd;
            frag->qEnd = b->qEnd;
        } 
        else 
        {
            frag->tEnd = b->tEnd;
            frag->qEnd = b->qEnd;

            verbose(4, "\nnew fragment\n");
            // search backwards in lastBlockStarts for a start which is far away enough
            struct lastBlockStartList *lb = lastBlockStarts;
            verbose(4, "overlap between frag and last Block: target %d, query %d\n", (b->tStart - lb->tStart), (b->qStart - lb->qStart));
            while ( ( ((b->tStart - lb->tStart) <= maxHole ) || (b->qStart - lb->qStart)  <= maxHole) && lb->previous!=NULL) {
                if ((lb->tEnd - lb->previous->tStart >= maxHole) || (lb->qEnd - lb->previous->qStart >= 2*maxHole))
                    break; // do not step over big holes that are > 2*maxHole, skip them
                lb = lb->previous;
                verbose(4, "Went back to %d, overlaps are %d/%d\n", lb->tStart, (b->tStart - lb->tStart), (b->qStart - lb->qStart));
            }
                
            // create a new fragment with these start positions
            struct chainFragList *newfrag = needMem(sizeof(struct chainFragList));
            newfrag->tStart = lb->tStart;
            newfrag->qStart = lb->qStart;
            verbose(4, "New fragment coordinates: tStart %d, tEnd %d\n", newfrag->tStart, newfrag->tEnd);
            
            // if the last block was very long, shorten it to minOverlap
            if (newfrag->tEnd - newfrag->tStart >= minOverlap) 
                newfrag->tStart = newfrag->tEnd - minOverlap;
            if (newfrag->qEnd - newfrag->qStart >= minOverlap) 
                newfrag->qStart = newfrag->qEnd - minOverlap;
           
            // to create only slightly overlapping fragments
            //newfrag->tStart = lastBlock->tStart;
            //newfrag->qStart = lastBlock->qStart;
            
            // to create non-overlapping fragments:
            //newfrag->tStart = b->tStart;
            //newfrag->qStart = b->qStart;

            newfrag->tEnd = b->tEnd;
            newfrag->qEnd = b->qEnd;

            // append newfrag to list of fragments
            frag->next = newfrag;
            frag = newfrag;
            frag->next=NULL;
        }

    }
    // free lastBlockStarts list at end of chain
    slFreeList(blockStart1);

    return frag1;
}

void writeLiftSpec(FILE* liftSpecFile, int start, char* idLine, int fragSize, char* name, int chromSize)
/*  LiftSpec
 *  is tab-delimited with each line of the form:
 *     offset oldName oldSize newName newSize
 *  example lines:
 *     0       test0   1000    tll_mel 3646
 *     1000    test1   1000    tll_mel 3646
 *     2000    test2   1000    tll_mel 3646
 *     3000    test3   646     tll_mel 3646
 *
 */
{
    if (liftSpecFile==NULL) 
        return;

    char* buf = needMem(strlen(idLine) + strlen(name) + (3*9) + 6); 
    sprintf(buf, "%d\t%s\t%d\t%s\t%d\n", start, idLine, fragSize, name, chromSize);
    mustWrite(liftSpecFile, buf, strlen(buf));
    freeMem(buf);
}

struct dnaSeq* twoBitReadStrandedSeqFrag(struct twoBitFile* twoBit, char* seqname, int seqLen, int start, int end, char strand) 
/* get qSeq, taking care to substract start/end from seqLen and revcomp sequence if strand is - */
{
    if (strand=='-') 
    {
        int end2  = seqLen - start;
        start     = seqLen - end;
        end       = end2;
    }
    struct dnaSeq *seq = twoBitReadSeqFrag(twoBit, seqname, start, end);
    if (strand=='-') {
        reverseComplement(seq->dna, seq->size);
    }
    return seq;
}

void chainToFasta(char *chainFilename, 
		     char *t2bitFilename, char *q2bitFilename,
		     char *outDirname)
/* chainToFasta - convert chains to one fasta file per chain */
{
    verbose(1, "maxFragSize %d, minOverlap %d, maxHole %d\n", maxFragSize, minOverlap, maxHole);
    struct twoBitFile* t2bit = twoBitOpen(t2bitFilename);
    char* tDb = basename(t2bitFilename);
    struct twoBitFile* q2bit = twoBitOpen(q2bitFilename);
    char* qDb = basename(q2bitFilename);

    FILE* liftSpecQFile = NULL;
    FILE* liftSpecTFile = NULL;

    if (liftSpecQFilename!=NULL) 
        liftSpecQFile = mustOpen(liftSpecQFilename, "w");

    if (liftSpecTFilename!=NULL) 
        liftSpecTFile = mustOpen(liftSpecTFilename, "w");

    struct lineFile *chainFile;
    struct chain *chain = NULL;

    /* Iterate over chains and extract sequences from 2bitFiles */
    verbose(1, "Reading, splitting chains from %s, generating fasta files...\n", chainFilename);
    chainFile = lineFileOpen(chainFilename, TRUE);
    int fileCount = 0;
    while ((chain = chainRead(chainFile)) != NULL)
    {
        struct chainFragList *chainFragments  = splitChain(chain, maxFragSize);
        struct chainFragList *chainFrag = chainFragments;

        while (chainFrag) {
            verbose(4, "writing chain frag, start %d, end %d, qstart %d, qend %d\n", chainFrag->tStart, chainFrag->tEnd, chainFrag->qStart, chainFrag->qEnd);

            char outFilename[1024];
            //sprintf(outFilename, "%s/%s-%d.fa", outDirname, chain->tName, chainFrag->tStart);
            sprintf(outFilename, "%s/%s-%d.fa", outDirname, chain->tName, chainFrag->tStart);

            FILE *outFile = mustOpen(outFilename, "w");

            char *tIdLine = makeLocString(tDb, chain->tName, chainFrag->tStart, chainFrag->tEnd, '+', fileCount);
            char *qIdLine = makeLocString(qDb, chain->qName, chainFrag->qStart, chainFrag->qEnd, chain->qStrand, fileCount);
            struct dnaSeq *tSeq = twoBitReadSeqFrag(t2bit, chain->tName, chainFrag->tStart, chainFrag->tEnd);
            verbose(5, "getting qseq\n");
            struct dnaSeq *qSeq = twoBitReadStrandedSeqFrag(q2bit, chain->qName, chain->qSize, chainFrag->qStart, chainFrag->qEnd, chain->qStrand);
            faWriteNext(outFile, tIdLine, tSeq->dna, tSeq->size);
            faWriteNext(outFile, qIdLine, qSeq->dna, qSeq->size);

            int tChromSize = twoBitSeqSize(t2bit, chain->tName);
            int tFragSize = chainFrag->tEnd - chainFrag->tStart;
            writeLiftSpec(liftSpecTFile, chainFrag->tStart, tIdLine, tFragSize, chain->tName, tChromSize);
            writeLiftSpec(liftSpecQFile, chainFrag->tStart, tIdLine, tFragSize, chain->tName, tChromSize);

            strtok(tIdLine, "."); // strip off species name from id line
            writeToBed(fragmentFilename, chain->tName, chainFrag->tStart,
                chainFrag->tEnd, tIdLine, chain->qName, chain->qSize, chainFrag->qStart, chainFrag->qEnd, chain->qStrand);

            fileCount++;

            freeMem(tIdLine);
            freeMem(qIdLine);
            carefulClose(&outFile);

            chainFrag = chainFrag->next;
        }

        slFreeList(&chainFragments);
    }

    verbose(1, "Done, generated %d fasta files\n", fileCount);
}


int main(int argc, char *argv[])
/* Process command line. */
{
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5)
	usage();
    maxFragSize = optionInt("maxFragSize", maxFragSize);
    minOverlap  = optionInt("minOverlap", minOverlap);
    maxHole     = optionInt("maxHole", maxHole);
    fragmentFilename = optionVal("fragmentFile", NULL);
    liftSpecTFilename = optionVal("liftT", NULL);
    liftSpecQFilename = optionVal("liftQ", NULL);
    chainToFasta(argv[1], argv[2], argv[3], argv[4]);

    return 0;
}

