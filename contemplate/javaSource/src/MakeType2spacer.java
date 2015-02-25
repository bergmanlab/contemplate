import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

import org.biojava.bio.symbol.IllegalSymbolException;

import org.biojava.bio.seq.DNATools;

import org.biojava.utils.ChangeVetoException;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dist.PairDistribution;

import java.io.PrintWriter;

/**
 * This class creates an instance of the Type 2 Spacer submodel 
 * described by the parameters supplied on the command line 
 * and then writes out an XML representation of the object to STDOUT.
 * A Type 2 Spacer submodel created and initialized by this class
 * can be used for subsequent training with TrainType2spacer.
 * <p>
 * Use as follows:
 * <pre>
 * java MakeType2spacer &lt;no of match/mismatch only stages&gt; &lt;no of mixed stages&gt; &lt;%GC seq 0&gt; &lt;%GC seq 1&gt; &lt;initial self transition weighting&gt;
 * </pre>
 * <p>
 * Match/mismatch stages emit only pairs of symbols. Mixed stages can also emit insert states in either seq 0 or seq 1.
 * The GC content of sequence 0 and 1 are specified as a percentage, while the probability of the self-transition
 * in the chain states is listed as a weighting (relative to 1 for a forward transition).
 * @author David Huen
 */
public class MakeType2spacer
{
    public static void main(String [] args)
        throws Exception
    {
        if (args.length != 5) {
            System.out.println("Usage: java MakeType2spacer <no of mismatch only stages> <no of mixed stages> <%GC seq 0> <%GC seq 1> <initial self transition weighting>");
            System.exit(1);
        }

        // create a simple Type2spacer
        int nMatchOnlyStages = Integer.parseInt(args[0]);
        int nMixedStages = Integer.parseInt(args[1]);
        Type2spacer thisType2spacer = new Type2spacer(nMatchOnlyStages, nMixedStages, DNATools.getDNA(), DNATools.getDNA());

        // create a writer
        PrintWriter pw = new PrintWriter(System.out);
        XMLWriter xw = new PrettyXMLWriter(pw);

        // create and set a seq 0 dist
        Distribution seq0Dist = makeDNADistribution(Double.parseDouble(args[2]));

        thisType2spacer.setSeq0Dist(seq0Dist);

        // create and set a seq 1 dist
        Distribution seq1Dist = makeDNADistribution(Double.parseDouble(args[3]));

        thisType2spacer.setSeq1Dist(seq1Dist);


        Distribution dna2Dist = new PairDistribution(seq0Dist, seq1Dist);

        thisType2spacer.setEmissionDist(dna2Dist);

        thisType2spacer.setModelParameters(Double.parseDouble(args[4]), 1.0, 1.0, 1.0);

        // now write out XML file describing this Type2spacer
        xw.printRaw("<?xml version=\"1.0\"?>\n");
        thisType2spacer.dumpParams(xw);

        pw.flush(); pw.close();
    }


    public static Distribution makeDNADistribution(double percentGC)
        throws IllegalSymbolException, ChangeVetoException
    {
        Distribution dist = new SimpleDistribution(DNATools.getDNA());
        double gc = 0.005 * percentGC;
        double at = 0.5 * (1.0 - gc - gc);
        dist.setWeight(DNATools.a(), at);
        dist.setWeight(DNATools.t(), at);
        dist.setWeight(DNATools.c(), gc);
        dist.setWeight(DNATools.g(), gc);

        return dist;
    }
}
