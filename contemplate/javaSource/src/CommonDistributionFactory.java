/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

import java.util.ArrayList;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;

import org.biojava.bio.dist.AbstractDistribution;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTrainer;
import org.biojava.bio.dist.DistributionTrainerContext;

/**
 * This distribution is intended for use as a transition
 * Distribution for a State in a MarkovModel.  In models, it
 * is frequently desirable to have entire groups of States
 * have identical patterns of probabilities for self-transitions and
 * transitions of neighbouring States.  However, one
 * Distribution cannot be used for all these States
 * because the specific destination of the transitions of a State
 * are unique for that State.
 * <p>
 * This factory class allows you to define a Distribution that has
 * non-zero weights for a predefined number of Symbols
 * without actually specifying those Symbols.
 * Specific Symbols are be bound to those values
 * when obtaining a new dependent Distribution instance
 * from this class.  This new instance can then be supplied
 * to the setWeights method in MarkovModel to define the
 * transition probabilities from that State.
 * <p>
 * All the dependent Distributions act as if they were one
 * one Distribution.  Calling setWeight in one of these
 * for a Symbol will change the underlying weight for
 * the value associated with it as well as all the values
 * return by getWeight for all Symbols associated with that
 * value in the other dependent Distributions.
 * <p>
 * During training, the trainer accumulates counts from
 * all dependent Distributions and eventually updates the
 * single common Distribution which in turn results in all
 * the dependent Distributions being updated.
 *
 * @author David Huen
 * @since 1.4
 */

public class CommonDistributionFactory
{
    private int nNonzeroValues; // number of non-zero values in this Distribution
    private double [] probability; // probability values stored here

    private double [] counts; // used during training.
    private boolean countsTrained = true;

    /**
     * creates a Distribution with specified number of non-zero values.
     * @param nNonzeroValues required number of non-zero values in Distribution.
     */
    public CommonDistributionFactory(
        int nNonzeroValues
        )
    {
        this.nNonzeroValues = nNonzeroValues;
        probability = new double [nNonzeroValues];

        // initialise the probabilities
        for (int i=0; i<nNonzeroValues; i++) {
            probability[i] = 0.0;
        }

        counts = new double [nNonzeroValues];
    }

    /**
     * initialise the internal array to the specified values.
     * <p>
     * <b>This is a fairly dangerous method to use.  Be sure you know what you are doing!!!</b>
     */
    public void setProbabilities(double [] values)
        throws IllegalArgumentException
    {
        if (values.length == nNonzeroValues) {
            for (int i=0; i < nNonzeroValues; i++) {
                probability[i] = values[i];
            }
        }
        else
            throw new IllegalArgumentException("array size does not match expected number of non-zero values.");
    }

    /**
     * get the array that backs this distribution.
     */
    public double [] getProbabilities()
    {
        return probability;
    }

    /**
     * make a dependent Distribution with the specified mapping of Symbols to values.
     * @param alfa the FiniteAlphabet that the requested Distribution is over.
     * @param assocSymbols the mapping of Symbols to the specific values of the Distribution.
     */
    public Distribution getDependentDistribution(FiniteAlphabet alfa, AtomicSymbol [] assocSymbols)
        throws IllegalSymbolException
    {
        return new DependentDistribution(alfa, assocSymbols);
    }

    private class DependentDistribution
        extends AbstractDistribution
    {
        private ArrayList assocSymbols;
        private Distribution nullModel = null;
        private FiniteAlphabet alfa;

        private DependentDistribution(FiniteAlphabet alfa, AtomicSymbol [] assocSymbols)
            throws IllegalSymbolException
        {
            // validate number of Symbols
            if (assocSymbols.length != nNonzeroValues) throw new IllegalArgumentException("number of Symbols is incorrect.");

            // validate the symbols against the alphabet
            this.assocSymbols = new ArrayList(nNonzeroValues);

            for (int i=0; i < nNonzeroValues; i++) {
                if (!alfa.contains(assocSymbols[i]))
                    throw new IllegalSymbolException("the symbol specified is not in the alphabet for this Distribution.");

                this.assocSymbols.add(assocSymbols[i]);
            }

            this.alfa = alfa;
        }

        final public Alphabet getAlphabet() { return alfa; }
        final protected Distribution getNullModelImpl () { return nullModel; }

        /**
         * note this implementation permits the setting of a null model
         * for each individual dependent Distribution.
         */
        final protected void setNullModelImpl(Distribution nullModel)
        {
            this.nullModel = nullModel;
        }

        /**
         * this must map the specified Symbol back to the underlying
         * associated probability.
         */
        final protected double getWeightImpl(AtomicSymbol sym)
        {
            int idxSym;
            if ((idxSym = assocSymbols.indexOf(sym)) != -1) {
                // known Symbol, return value
                return probability[idxSym];
            }
            else
                return 0.0;
        }

        final protected void setWeightImpl(AtomicSymbol sym, double weight)
        {
            int idxSym = assocSymbols.indexOf(sym);
            if (idxSym != -1) {
                // known Symbol, return value
                probability[idxSym] = weight;
            }
        }

        public Distribution getNullModel() {
            return this.nullModel;
        }

        final public void registerWithTrainer(DistributionTrainerContext dtc)
        {
            dtc.registerTrainer(this, new Trainer());
        }

        /**
         * a private Trainer for this Distribution.
         * counts are only added when they are associated with a Symbol in this Distribution.
         */
        private class Trainer
            implements DistributionTrainer
        {
            public void addCount(DistributionTrainerContext dtc, AtomicSymbol sym, double times)
            {
                int idxSym = assocSymbols.indexOf(sym);
                if (idxSym != -1) {
//                    System.out.println("adding counts to " + sym.getName() + " at index " + idxSym);
                    counts[idxSym] += times;
                }
            }

            public double getCount(DistributionTrainerContext dtc, AtomicSymbol sym)
                throws IllegalSymbolException 
            {
                int idxSym = assocSymbols.indexOf(sym);
                if (idxSym != -1) {
                    return counts[idxSym];
                }
                else return 0.0;
            }

            public void clearCounts(DistributionTrainerContext dtc) 
            {
                for (int i=0; i< nNonzeroValues; i++) {
                    counts[i] = 0.0;
                }

                countsTrained = false;
            }


            public void train(DistributionTrainerContext dtc, double weight)
                throws ChangeVetoException 
            {
                if(!hasListeners())  {
                    trainImpl(dtc, weight);
                } else {
                    ChangeSupport changeSupport = getChangeSupport(Distribution.WEIGHTS);
                    synchronized(changeSupport) {
                        ChangeEvent ce = new ChangeEvent(
                        DependentDistribution.this,
                        Distribution.WEIGHTS
                        );
                        changeSupport.firePreChangeEvent(ce);
                        trainImpl(dtc, weight);
                        changeSupport.firePostChangeEvent(ce);
                    }
                }
            }

            protected void trainImpl(DistributionTrainerContext dtc, double weight) 
            {
                // we only compute the Distribution if it has not been done
                // already in this cycle.
                if (!countsTrained) {
                    double sum = 0.0;

                    for (int i=0; i < nNonzeroValues; i++) {
//                        System.out.println("training counts[" + i + "] = " + counts[i]);
                        counts[i] += weight;
                        sum += counts[i];
                    }

                    double invSum = 1.0/sum;

                    for (int i=0; i < nNonzeroValues; i++) {
//                        System.out.println("training probability[" + i + "] = " + probability[i]);
                        probability[i] = counts[i] * invSum;
                    }

                    countsTrained = true;
                }
            }
        }
    }
}

