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

//import java.util.List;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.GapDistribution;
//import org.biojava.bio.dist.PairDistribution;

import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.AlphabetManager;

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

import org.xml.sax.SAXException;

import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.NamedNodeMap;

/**
  A class that implements a multi-staged HMM that emits two
 * sequences of different length as a block of (mis)matched
 * nucleotides followed by successive gaps in the shorter sequence.
 * This class models long, unconstrained spacers between enhancers.
 * This class differs from Type1spacer in that spacers have non-zero 
 * length in both sequences.
 *
 * @author David Huen
 */
public class Type2spacer
    implements SubModelConstructor
{
    private String name = "Type2spacer";

    // alphabets for sequence 0 and 1
    private FiniteAlphabet seq0Alfa = null;
    private FiniteAlphabet seq1Alfa = null;

    // cross product alphabet for this model
    private FiniteAlphabet modelAlfa = null;

    // there is a common emission Distribution for all the States.
    private Distribution emissionDist = null;

    // Distribution of emission in sequence 0 when gaps are inserted in sequence 1
    private Distribution seq0Dist = null;
    // Distribution of emission in sequence 1 when gaps are inserted in sequence 0
    private Distribution seq1Dist = null;

    // Gapped pair distributions for states that insert gaps
    private Distribution insert0Dist = null;
    private Distribution insert1Dist = null;

    // parameters
    private int nMixedStages = 0;
    private int nMatchOnlyStages = 0;
    private int nExtraMixedStages = 0;
    private int nExtraMatchOnlyStages = 0;
    private int nTotalStates = 0;

    // distribution factories
    // type 1: mismatch, multidestination
    final private CommonDistributionFactory multiDestMismatchFact = new CommonDistributionFactory(4);
    // type 2: insert in sequence 0
    final private CommonDistributionFactory insert0Fact = new CommonDistributionFactory(2);
    // type 3: insert in sequence 1
    final private CommonDistributionFactory insert1Fact = new CommonDistributionFactory(2);
    // type 4: mismatch, single destination
    final private CommonDistributionFactory singleDestMismatchFact = new CommonDistributionFactory(2);

    final static int [] mismatchAdvance = { 1, 1};
    final static int [] insert0Advance = { 0, 1};
    final static int [] insert1Advance = { 1, 0};

    // XML string constants
    final static String SUBMODEL = "submodel";
    final static String NAME = "name";
    final static String TYPE = "type";
    final static String MISMATCH_TYPE2SPACER = "mismatch_type2spacer";
    final static String N_MIXED_STAGES = "no_mixed_stages";
    final static String N_MATCH_ONLY_STAGES = "no_match_only_stages";
    final static String SEQ_0_ALFA = "seq0Alfa";
    final static String SEQ_1_ALFA = "seq1Alfa";
    final static String TYPE_1_DIST = "type1dist";
    final static String TYPE_2_DIST = "type2dist";
    final static String TYPE_3_DIST = "type3dist";
    final static String TYPE_4_DIST = "type4dist";
    final static String EMISSION_DIST = "emission_dist";
    final static String SEQ_0_DIST = "seq0_dist";
    final static String SEQ_1_DIST = "seq1_dist";
    final static String SELF = "self";
    final static String MISMATCH = "mismatch";
    final static String INSERT0 = "insert0";
    final static String INSERT1 = "insert1";

    final static Maker maker = new SubModelConstructor.Maker()
    {
        public String getKey() { return MISMATCH_TYPE2SPACER; }
        public SubModelConstructor makeConstructor(Node node)
            throws IllegalArgumentException
        {
            return new Type2spacer(node);
        }
    };

    /**
     * Instantiate the Type2spacer from a DOM subtree.
     * @param paramNode root node of the subtree parameterising this object.
     */
    public Type2spacer(Node paramNode)
        throws IllegalArgumentException
    {
        // we expect a "mismatch Type2spacer" Element
        String tag;
        NamedNodeMap attrs;

        if ((paramNode.getNodeType() == Node.ELEMENT_NODE) &&
            (tag = paramNode.getNodeName()).equals(SUBMODEL) &&
            ((attrs = paramNode.getAttributes()) != null)
            ) {
            // pick up the attributes.
            Node attr;

            // verify that the node really describes this class
            if (((attr = attrs.getNamedItem(TYPE)) == null) ||
                !(attr.getNodeValue().equals(MISMATCH_TYPE2SPACER))
                ) {
                throw new IllegalArgumentException("This node is not appropriate for Block");
            }

            // pick up the name attribute
            if ((attr = attrs.getNamedItem(NAME)) != null) {
                name = attr.getNodeValue();
            }

            // pick up the no_mixed_stages attribute
            if ((attr = attrs.getNamedItem(N_MIXED_STAGES)) != null) {
                nMixedStages = Integer.parseInt(attr.getNodeValue());
            }

            // pick up the no_match_only_stages attribute
            if ((attr = attrs.getNamedItem(N_MATCH_ONLY_STAGES)) != null) {
                nMatchOnlyStages = Integer.parseInt(attr.getNodeValue());
            }

            // pick up the alphabets for the sequences
            if ((attr = attrs.getNamedItem(SEQ_0_ALFA)) != null) {
                seq0Alfa = (FiniteAlphabet) AlphabetManager.alphabetForName(attr.getNodeValue());
            }

            if ((attr = attrs.getNamedItem(SEQ_1_ALFA)) != null) {
                seq1Alfa = (FiniteAlphabet) AlphabetManager.alphabetForName(attr.getNodeValue());
            }

            initialise();

            // now iterate thru the child nodes to pick up the distribution vectors.
            NodeList children = paramNode.getChildNodes();
            for (int i=0; i<children.getLength(); i++) {
                Node thisNode = children.item(i);

                if (thisNode.getNodeType() == Node.ELEMENT_NODE) {
                    if ((tag = thisNode.getNodeName()).equals(TYPE_1_DIST)) {
                        attrs = thisNode.getAttributes();

                        double [] probs = new double[4];
                        for (int j=0; j < probs.length; j++) {
                            probs[j] = 0.0;
                        }

                        if ((attr = attrs.getNamedItem(SELF)) != null) {
                            probs[0] = Double.parseDouble(attr.getNodeValue());
                        }

                        if ((attr = attrs.getNamedItem(MISMATCH)) != null) {
                            probs[1] = Double.parseDouble(attr.getNodeValue());
                        }

                        if ((attr = attrs.getNamedItem(INSERT0)) != null) {
                            probs[2] = Double.parseDouble(attr.getNodeValue());
                        }

                        if ((attr = attrs.getNamedItem(INSERT1)) != null) {
                            probs[3] = Double.parseDouble(attr.getNodeValue());
                        }

                        multiDestMismatchFact.setProbabilities(probs);
                    }
                    else if (tag.equals(TYPE_2_DIST)) {
                        attrs = thisNode.getAttributes();

                        double [] probs = new double[2];
                        for (int j=0; j < probs.length; j++) {
                            probs[j] = 0.0;
                        }

                        if ((attr = attrs.getNamedItem(SELF)) != null) {
                            probs[0] = Double.parseDouble(attr.getNodeValue());
                        }

                        if ((attr = attrs.getNamedItem(INSERT0)) != null) {
                            probs[1] = Double.parseDouble(attr.getNodeValue());
                        } 
                        insert0Fact.setProbabilities(probs);
                    }
                    else if (tag.equals(TYPE_3_DIST)) {
                        attrs = thisNode.getAttributes();

                        double [] probs = new double[2];
                        for (int j=0; j < probs.length; j++) {
                            probs[j] = 0.0;
                        }

                        if ((attr = attrs.getNamedItem(SELF)) != null) {
                            probs[0] = Double.parseDouble(attr.getNodeValue());
                        }

                        if ((attr = attrs.getNamedItem(INSERT1)) != null) {
                            probs[1] = Double.parseDouble(attr.getNodeValue());
                        }
                        insert1Fact.setProbabilities(probs);
                    }
                    else if (tag.equals(TYPE_4_DIST)) {
                        attrs = thisNode.getAttributes();

                        double [] probs = new double[2];
                        for (int j=0; j < probs.length; j++) {
                            probs[j] = 0.0;
                        }

                        if ((attr = attrs.getNamedItem(SELF)) != null) {
                            probs[0] = Double.parseDouble(attr.getNodeValue());
                        }

                        if ((attr = attrs.getNamedItem(MISMATCH)) != null) {
                            probs[1] = Double.parseDouble(attr.getNodeValue());
                        }

                        singleDestMismatchFact.setProbabilities(probs);
                    }
                    else if (tag.equals(EMISSION_DIST)) {
                        setEmissionDist(loadDOMDistribution(thisNode));
                    }
                    else if (tag.equals(SEQ_0_DIST)) {
                        setSeq0Dist(loadDOMDistribution(thisNode));
                    }
                    else if (tag.equals(SEQ_1_DIST)) {
                        setSeq1Dist(loadDOMDistribution(thisNode));
                    }
//                    else System.out.println("tag is " + tag);
                }
            }

            if ((nMixedStages == 0) || (nMatchOnlyStages == 0)) throw new IllegalArgumentException("initialisation was unsuccessful. DOM ill-formed?");
        }
        else throw new IllegalArgumentException("An Element was expected.");
    }

    /**
     * Constructor requiring explicit specification of all parameters.
     *
     * @param nMixedStages number of states that have match/insert functionality.
     * @param nMatchOnlyStages number of stages that only have match functionality.
     * @param seq0Alfa FiniteAlphabet of sequence 0.
     * @param seq1Alfa FiniteAlphabet of sequence 1.
     */
    public Type2spacer(
        int nMixedStages,
        int nMatchOnlyStages,
        FiniteAlphabet seq0Alfa,
        FiniteAlphabet seq1Alfa
        )
        throws IllegalArgumentException
    {
        if (nMixedStages < 0) throw new IllegalArgumentException("there must be at least one mixed stage");
        if (nMatchOnlyStages < 0) throw new IllegalArgumentException("there must be at least one match only stage");

        this.nMixedStages = nMixedStages;
        this.nMatchOnlyStages = nMatchOnlyStages;

        this.seq0Alfa = seq0Alfa;
        this.seq1Alfa = seq1Alfa;

        initialise();
    }

    private void initialise()
    {
        // there is always a final mixed stage: it differs from the earlier
        // ones by having a mismatch state with only one destination (type 4).
        // insert states are either type 2 or type 3.
        // there may be other mixed stages.  These have type 1 mismatch transition Distributions
        // and type 2/3 insert states.
        nExtraMixedStages = nMixedStages - 1;

        // the final match only state has multiple destinations. This has a type 1 distribution.
        // there may be others but all those only have single destinations.  They have type 4 distributions.
        nExtraMatchOnlyStages = nMatchOnlyStages - 1;

        nTotalStates = nMatchOnlyStages + 3 * nMixedStages;

        if ((seq0Alfa != null) && (seq1Alfa != null)) {
            List aList = new ArrayList(); aList.add(seq0Alfa); aList.add(seq1Alfa);
            modelAlfa = (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(aList);
        }
    }

    private Distribution loadDOMDistribution(Node thisNode)
        throws IllegalArgumentException
    {
        // there should be a child node of tag "Distribution"
        NodeList kids = thisNode.getChildNodes();
        boolean notFound = true; int j=0;
        while (notFound && (j < kids.getLength())) {
            Node kid = kids.item(j); j++;
            if (kid.getNodeName().equals("Distribution")) {
                try {
                    return XMLDistTools.readDistFromDOM(kid);
                }
                catch (SAXException se) {
                    throw new IllegalArgumentException("SAXException occurred during parsing of a distribution. File malformed?");
                }
            }
        }
        throw new IllegalArgumentException("XML for Distribution expected but not found.");
    }

    public String getName() { return name; }
    public void setName(String name) { this.name = name; }
    public void setForTraining(boolean forTraining) {}

    /**
     * create a Model with the specified number of
     * match only and free stages.  The initial State
     * of the model is returned.
     * <p>
     * <b>NOTE: You have to set the transition probability of the source state yourself.</b>
     */
    public State insertSubModel(
        MarkovModel model,
        State source,
        State dest,
        Annotation userAnnotation
        )
        throws IllegalAlphabetException, IllegalSymbolException, ChangeVetoException
    {
        // check that the alphabets are compatible
        if (model.emissionAlphabet() != modelAlfa) throw new IllegalAlphabetException("incompatible alphabets!");

        int i = 0;
        State [] subModelStates = new State [nTotalStates];

        // the problem is that the state alphabet is changing
        // until the the model is completed and there is no
        // guarantee that Distribution objects can cope with
        // that, especially if they use some internal cache
        // that is indexed by an alphabet index. So we will create
        // and add all states with only the emission distribution.
        // In the next pass we add the transitions and the transition
        // distributions.

        /**************************************
         * make the required number of states *
         **************************************/
        State thisState;

        // mixed stages
        for (int j = 0; j < nMixedStages; j++) {
            subModelStates[i++] = thisState = new SimpleEmissionState(
                name + ":mixed_mismatch:" + Integer.toString(j),
                userAnnotation,
                mismatchAdvance,
                emissionDist
            );
            model.addState(thisState);

            subModelStates[i++] = thisState = new SimpleEmissionState(
                name + ":mixed_insert0:" + Integer.toString(j),
                userAnnotation,
                insert0Advance,
                insert0Dist
            );
            model.addState(thisState);

            subModelStates[i++] = thisState = new SimpleEmissionState(
                name + ":mixed_insert1:" + Integer.toString(j),
                userAnnotation,
                insert1Advance,
                insert1Dist
            );
            model.addState(thisState);
        }
        
        for (int j = 0; j < nMatchOnlyStages; j++) {
            subModelStates[i++] = thisState = new SimpleEmissionState(
                name + ":mismatch:" + Integer.toString(j),
                userAnnotation,
                mismatchAdvance,
                emissionDist
            );
            model.addState(thisState);
        }

        State mismatch;
        State insert0;
        State insert1;

        // now put together the topology
        // FINAL MULTIDESTINATION STAGE
        // type 4 mismatch state
        i = 0;
        mismatch = subModelStates[i++];
        model.createTransition(mismatch, mismatch);
        model.createTransition(mismatch, dest);
        AtomicSymbol [] dests = new AtomicSymbol[2];
        dests[0] = mismatch;
        dests[1] = dest;
        Distribution thisDist = singleDestMismatchFact.getDependentDistribution(model.transitionsFrom(mismatch), dests);
        model.setWeights(mismatch, thisDist);

        // insert 0 state
        insert0 = subModelStates[i++];
        model.createTransition(insert0, insert0);
        model.createTransition(insert0, dest);
        dests = new AtomicSymbol[2];
        dests[0] = insert0;
        dests[1] = dest;
        thisDist = insert0Fact.getDependentDistribution(model.transitionsFrom(insert0), dests);
        model.setWeights(insert0, thisDist);

        // insert 1 state
        insert1 = subModelStates[i++];
        model.createTransition(insert1, insert1);
        model.createTransition(insert1, dest);
        dests = new AtomicSymbol[2];
        dests[0] = insert1;
        dests[1] = dest;
        thisDist = insert1Fact.getDependentDistribution(model.transitionsFrom(insert1), dests);
        model.setWeights(insert1, thisDist);

        // ALL BUT LAST MULTIDESTINATION STAGES
        for (int j = 0; j < nExtraMixedStages; j++) {
            // mismatch, type 1
            thisState = subModelStates[i++];
            model.createTransition(thisState, thisState);
            model.createTransition(thisState, mismatch);
            model.createTransition(thisState, insert0);
            model.createTransition(thisState, insert1);
            dests = new AtomicSymbol[4];
            dests[0] = thisState;
            dests[1] = mismatch; mismatch = thisState;
            dests[2] = insert0;
            dests[3] = insert1;
            thisDist = multiDestMismatchFact.getDependentDistribution(model.transitionsFrom(thisState), dests);
            model.setWeights(thisState, thisDist);

            // insert0
            thisState = subModelStates[i++];
            model.createTransition(thisState, thisState);
            model.createTransition(thisState, insert0);
            dests = new AtomicSymbol[2];
            dests[0] = thisState;
            dests[1] = insert0; insert0 = thisState;
            thisDist = insert0Fact.getDependentDistribution(model.transitionsFrom(thisState), dests);
            model.setWeights(thisState, thisDist);

            // insert1
            thisState = subModelStates[i++];
            model.createTransition(thisState, thisState);
            model.createTransition(thisState, insert1);
            dests = new AtomicSymbol[2];
            dests[0] = thisState;
            dests[1] = insert1; insert1 = thisState;
            thisDist = insert1Fact.getDependentDistribution(model.transitionsFrom(thisState), dests);
            model.setWeights(thisState, thisDist);
        }

        // FINAL MISMATCH ONLY STATE
        // mismatch, type 1
        thisState = subModelStates[i++];
        model.createTransition(thisState, thisState);
        model.createTransition(thisState, mismatch);
        model.createTransition(thisState, insert0);
        model.createTransition(thisState, insert1);
        dests = new AtomicSymbol[4];
        dests[0] = thisState;
        dests[1] = mismatch; mismatch = thisState;
        dests[2] = insert0;
        dests[3] = insert1;
        thisDist = multiDestMismatchFact.getDependentDistribution(model.transitionsFrom(thisState), dests);
        model.setWeights(thisState, thisDist);

        // ALL BUT LAST MISMATCH STATES
        // mismatch, type 4
        for (int j = 0; j < nExtraMatchOnlyStages; j++) {
            thisState = subModelStates[i++];
            model.createTransition(thisState, thisState);
            model.createTransition(thisState, mismatch);
            dests = new AtomicSymbol[2];
            dests[0] = thisState;
            dests[1] = mismatch; mismatch = thisState;
            thisDist = singleDestMismatchFact.getDependentDistribution(model.transitionsFrom(thisState), dests);
            model.setWeights(thisState, thisDist); 
        }

        // now connect source to first stage
        model.createTransition(source, thisState);
        return thisState;
    }

    /**
     * set the cross product emission distribution for (mis)match states.
     */
    public void setEmissionDist(Distribution emissionDist)
    {
        // validate the alphabet
        if (emissionDist.getAlphabet() != modelAlfa) 
            throw new IllegalArgumentException("The argument is of alphabet " 
                + emissionDist.getAlphabet().getName() + " but the model requires " + modelAlfa.getName());

        this.emissionDist = emissionDist;
    }

    /**
     * set the emission distribution for sequence 0 when gaps are inserted in seq1.
     */
    public void setSeq0Dist(Distribution seq0Dist)
        throws IllegalArgumentException 
   {
        if (seq0Dist.getAlphabet() != seq0Alfa)
            throw new IllegalArgumentException("Alphabet " 
                + seq0Dist.getAlphabet().getName() + " was supplied but " + seq0Alfa.getName() + " was required.");

        this.seq0Dist = seq0Dist;
        insert1Dist = new GappedPairDistribution(seq0Dist, new GapDistribution(seq1Alfa));
    }

    /**
     * set the emission distribution for sequence 1 when gaps are inserted in seq1.
     */
    public void setSeq1Dist(Distribution seq1Dist)
        throws IllegalArgumentException
    {
        if (seq1Dist.getAlphabet() != seq1Alfa)
            throw new IllegalArgumentException("Alphabet "
                + seq1Dist.getAlphabet().getName() + " was supplied but " + seq1Alfa.getName() + " was required.");

        this.seq1Dist = seq1Dist;
        insert0Dist = new GappedPairDistribution(new GapDistribution(seq0Alfa), seq1Dist);
    }

    public Distribution getEmissionDist() { return emissionDist; }
    public Distribution getSeq0Alfa() { return seq0Dist; }
    public Distribution getSeq1Alfa() { return seq1Dist; }

    /**
     * use to initialise model parameters.
     */
    public void setModelParameters(
        double wSelf,
        double wMatch,
        double wInsert0,
        double wInsert1
        )
    {
        // compute the actual probabilities.
        double wTotal = wSelf + wMatch + wInsert0 + wInsert1;
        double pSelf = wSelf/wTotal;
        double pMatch = wMatch/wTotal;
        double pInsert0 = wInsert0/wTotal;
        double pInsert1 = wInsert1/wTotal;

        // when initialising, the pSelf of every state is kept
        // constant while the other parameters are varied as necessary.
        // type 1: multidestination mismatch
        double [] prob = new double [4];
        prob[0] = pSelf;
        prob[1] = pMatch;
        prob[2] = pInsert0;
        prob[3] = pInsert1;
        multiDestMismatchFact.setProbabilities(prob);

        // type 2: insert in sequence 0
        prob = new double [2];
        prob[0] = pSelf;
        prob[1] = 1.0 - pSelf;
        insert0Fact.setProbabilities(prob);

        // type 3: insert in sequence 0
        prob = new double [2];
        prob[0] = pSelf;
        prob[1] = 1.0 - pSelf;
        insert1Fact.setProbabilities(prob);

        // type 4: single destination mismatch
        prob = new double [2];
        prob[0] = pSelf;
        prob[1] = 1.0 - pSelf;
        singleDestMismatchFact.setProbabilities(prob);
    }

    public void dumpParams(XMLWriter xw)
        throws IOException
    {
        xw.openTag(SUBMODEL);
        xw.attribute(NAME, name);
        xw.attribute(TYPE, MISMATCH_TYPE2SPACER);
        xw.attribute(N_MIXED_STAGES, Integer.toString(nMixedStages));
        xw.attribute(N_MATCH_ONLY_STAGES, Integer.toString(nMatchOnlyStages));
        xw.attribute(SEQ_0_ALFA, seq0Alfa.getName());
        xw.attribute(SEQ_1_ALFA, seq1Alfa.getName());

        // write out common distribution vectors
        xw.openTag(TYPE_1_DIST);
        double [] probs = multiDestMismatchFact.getProbabilities();
        xw.attribute(SELF, Double.toString(probs[0]));
        xw.attribute(MISMATCH, Double.toString(probs[1]));
        xw.attribute(INSERT0, Double.toString(probs[2]));
        xw.attribute(INSERT1, Double.toString(probs[3]));
        xw.closeTag(TYPE_1_DIST);

        xw.openTag(TYPE_2_DIST);
        probs = insert0Fact.getProbabilities();
        xw.attribute(SELF, Double.toString(probs[0]));
        xw.attribute(INSERT0, Double.toString(probs[1]));
        xw.closeTag(TYPE_2_DIST);

        xw.openTag(TYPE_3_DIST);
        probs = insert1Fact.getProbabilities();
        xw.attribute(SELF, Double.toString(probs[0]));
        xw.attribute(INSERT1, Double.toString(probs[1]));
        xw.closeTag(TYPE_3_DIST);

        xw.openTag(TYPE_4_DIST);
        probs = singleDestMismatchFact.getProbabilities();
        xw.attribute(SELF, Double.toString(probs[0]));
        xw.attribute(MISMATCH, Double.toString(probs[1]));
        xw.closeTag(TYPE_4_DIST);

        // dump the emission Distribution
        if (emissionDist != null) {
            xw.openTag(EMISSION_DIST);
            XMLDistTools.writeDist2XML(emissionDist, xw);
            xw.closeTag(EMISSION_DIST);
        }

        if (seq0Dist != null) {
            xw.openTag(SEQ_0_DIST);
            XMLDistTools.writeDist2XML(seq0Dist, xw);
            xw.closeTag(SEQ_0_DIST);
        }

        if (seq0Dist != null) {
            xw.openTag(SEQ_1_DIST);
            XMLDistTools.writeDist2XML(seq1Dist, xw);
            xw.closeTag(SEQ_1_DIST);
        }

        xw.closeTag(SUBMODEL);
    }

    public static void registerMaker(Map map)
    {
        map.put(maker.getKey(), maker);
    }
}

