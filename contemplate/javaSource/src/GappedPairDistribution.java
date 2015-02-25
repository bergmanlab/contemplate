
import org.biojava.bio.BioError;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.ExtendedDistributionTrainer;
import org.biojava.bio.dist.GapDistribution;
import org.biojava.bio.dist.IgnoreCountsTrainer;

import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;

import org.biojava.bio.dist.PairDistribution;

import java.io.Serializable;
import java.util.List;
import java.util.Arrays;

/**
 * This class implements PairDistributions
 * in which one Distribution emits gaps only.
 * It is intended for use in MarkovModels and
 * differs from the PairDistribution in that
 * the addCount method explicitly expects a
 * gap in the Symbol it handles.
 *
 * @author David Huen
 */
public class GappedPairDistribution
    extends PairDistribution
{
    private final boolean seq0Gapped;
    private final Distribution dist;
    private final GapDistribution gapDist;

    GappedPairDistribution(GapDistribution gapDist, Distribution dist)
    {
        super(gapDist, dist);
        this.gapDist = gapDist;
        this.dist = dist;
        seq0Gapped = true;
    }

    GappedPairDistribution(Distribution dist, GapDistribution gapDist)
    {
        super(dist, gapDist);
        this.dist = dist;
        this.gapDist = gapDist;
        seq0Gapped = false;
    }

    public void registerWithTrainer(DistributionTrainerContext dtc) {
//        System.out.println("registering GappedPairTrainer");
        dtc.registerTrainer(this, new GappedPairTrainer());

        // also register the underlying distribution
        dist.registerWithTrainer(dtc);
    }

    private class GappedPairTrainer
        extends IgnoreCountsTrainer
        implements Serializable, ExtendedDistributionTrainer
    {
        public double getCount(DistributionTrainerContext dtc, AtomicSymbol as)
            throws IllegalSymbolException 
        {
            List symL = as.getSymbols();

            Symbol gap;
            Symbol lookup;
            if (seq0Gapped) {
                gap = (Symbol) symL.get(0);
                lookup = (Symbol) symL.get(1);
            }
            else {
                gap = (Symbol) symL.get(1);
                lookup = (Symbol) symL.get(0);
            }

            // we only emit gaps in one sequence and so there are no counts for anything else.
            if (gap != gapDist.getAlphabet().getGapSymbol()) return 0.0;

            dist.getAlphabet().validate(lookup);
            return (dtc.getCount(dist, lookup));     
        }

        public void addCount(
            DistributionTrainerContext dtc, AtomicSymbol sym, double times
            )
            throws IllegalSymbolException
        {
            addCount(dtc, (Symbol) sym, times);
        }

        public void addCount(
            DistributionTrainerContext dtc, Symbol sym, double times
            ) 
            throws IllegalSymbolException
        {
//            System.out.println("gappedpair (start): observed " + sym.getName());
            if (sym instanceof BasisSymbol) {
                List symL = ((BasisSymbol) sym).getSymbols();

                Symbol gap;
                Symbol lookup;

                if (seq0Gapped) {
                    gap = (Symbol) symL.get(0);
                    lookup = (Symbol) symL.get(1);
                }
                else {
                    gap = (Symbol) symL.get(1);
                    lookup = (Symbol) symL.get(0);
                }

                // I should not be getting this symbol but it is not illegal for this alphabet
//                System.out.println("gappedpair: observed " + lookup.getName());
//                if (gap != gapDist.getAlphabet().getGapSymbol()) return; /* commented out, for some reason, gap alfa is inconsistent */

                // increment count of the base Distribution
//                System.out.println("gappedpair: adding " + lookup.getName());
                dtc.addCount(dist, lookup, times);
            }
            else throw new IllegalSymbolException("A symbol supplied to PairDistribution must be a BasisSymbol.");
        }

        public Symbol sampleSymbol() 
        {
            try {
                Symbol [] sArray;

                if (seq0Gapped) {
                    sArray = new Symbol [] { gapDist.sampleSymbol(), dist.sampleSymbol() };
                }
                else {
                    sArray = new Symbol [] { dist.sampleSymbol(), gapDist.sampleSymbol() };
                }

                return getAlphabet().getSymbol(Arrays.asList(sArray));
            }
            catch (IllegalSymbolException ise) {
                throw new BioError("Couldn't sample symbol", ise);
            }
        }

        /**
         * it should not be necessary to implement a train() method as the underlying
         * Distribution is also independently registered with the context.
         */
/*
        public void train(DistributionTrainerContext dtc, double weight)
        {
            dtc.getTrainer(dist).train(dtc, weight);
        }
*/
    }
}

