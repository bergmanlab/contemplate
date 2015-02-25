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
 * This class creates an instance of the Block submodel
 * described by the parameters supplied on the command line
 * and then writes out an XML representation of the object to STDOUT.
 * A Block submodel created and initialized by this class
 * can be used for subsequent training with TrainBlock.
 * <br><br>
 * Use as follows:<br>
 * <pre>
 * java MakeBlock &lt;no of match/mismatch states&gt; &lt;%GC seq 0&gt; &lt;%GC seq 1&gt; &lt;initial self transition weighting&gt;
 * </pre>
 * <p>
 * Match/mismatch stages are states that emit pairs of ungapped symbols. 
 * The GC content of sequence 0 and 1 are specified as a percentage, while the probability of the self-transition
 * in the chain states is listed as a weighting (relative to 1 for a forward transition).
 * @author David Huen
 */

public class MakeBlock
{
    public static void main(String [] args)
        throws Exception
    {
        if (args.length != 4) {
            System.out.println("Usage: java MakeBlock <no of states> <%GC seq 0> <%GC seq 1> <initial self transition weighting>");
            System.exit(1);
        }

        // create a simple chain
        int nStates = Integer.parseInt(args[0]);
        Block thisMatch = new Block(nStates, DNATools.getDNA(), DNATools.getDNA());

        // create a writer
        PrintWriter pw = new PrintWriter(System.out);
        XMLWriter xw = new PrettyXMLWriter(pw);

        // create and set a seq 0 dist
        Distribution seq0Dist = makeDNADistribution(Double.parseDouble(args[1]));
        Distribution seq1Dist = makeDNADistribution(Double.parseDouble(args[2]));

        Distribution dna2Dist = new PairDistribution(seq0Dist, seq1Dist);

        thisMatch.setEmissionDist(dna2Dist);

        thisMatch.setModelParameters(Double.parseDouble(args[3]), 1.0);

        // now write out XML file describing this chain
        xw.printRaw("<?xml version=\"1.0\"?>\n");
        thisMatch.dumpParams(xw);

        pw.flush(); pw.close();
    }


    private static Distribution makeDNADistribution(double percentGC)
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
