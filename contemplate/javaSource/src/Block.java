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

import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
//import org.biojava.bio.dist.GapDistribution;
//import org.biojava.bio.dist.PairDistribution;

import org.biojava.bio.dp.MarkovModel;
//import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
//import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.AlphabetManager;

import org.biojava.utils.ChangeVetoException;
//import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

import org.xml.sax.SAXException;

import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.NamedNodeMap;

/**
 * A class that implements a Block submodel that generates
 * ungapped alignments of conserved noncoding 
 * sequences with a negative binomial length distribution.
 */
public class Block
    implements SubModelConstructor
{
    private String name = "Block";

    private int nStates = 0; // no. of stages.
    private FiniteAlphabet seq0Alfa = null;
    private FiniteAlphabet seq1Alfa = null;

    // cross product alphabet for this model
    private FiniteAlphabet modelAlfa = null;

    // there is a common emission Distribution for all the States.
    private Distribution emissionDist = null;

    // transition Distribution factory
    final private CommonDistributionFactory transitionFact = new CommonDistributionFactory(2);

    final static int [] matchAdvance = { 1, 1};

    // XML string constants
    final static String SUBMODEL = "submodel";
    final static String NAME = "name";
    final static String TYPE = "type";
    final static String BLOCK = "block";
    final static String N_STATES = "no_states";
    final static String SEQ_0_ALFA = "seq0Alfa";
    final static String SEQ_1_ALFA = "seq1Alfa";
    final static String TRANSITION_DIST = "transDist";
    final static String EMISSION_DIST = "emission_dist";
    final static String SELF = "self";
    final static String NEXT = "next";

    final static Maker maker = new SubModelConstructor.Maker()
    {
        public String getKey() { return BLOCK; }
        public SubModelConstructor makeConstructor(Node node)
            throws IllegalArgumentException
        {
            return new Block(node);
        }
    };

    /**
     * Instantiate an instance from a XML DOM Node.
     * @param paramNode root node of the subtree parameterising this object.
     */
    public Block(Node paramNode)
        throws IllegalArgumentException
    {
        // we expect a "mismatch chain" Element
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
                !(attr.getNodeValue().equals(BLOCK))
                ) {
                throw new IllegalArgumentException("This node is not appropriate for Block");
            }

            // pick up the no_mixed_stages attribute
            if ((attr = attrs.getNamedItem(NAME)) != null) {
                name = attr.getNodeValue();
            }

            // pick up the no_mixed_stages attribute
            if ((attr = attrs.getNamedItem(N_STATES)) != null) {
                nStates = Integer.parseInt(attr.getNodeValue());
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
                    if ((tag = thisNode.getNodeName()).equals(TRANSITION_DIST)) {
                        attrs = thisNode.getAttributes();

                        double [] probs = new double[2];
                        for (int j=0; j < probs.length; j++) {
                            probs[j] = 0.0;
                        }

                        if ((attr = attrs.getNamedItem(SELF)) != null) {
                            probs[0] = Double.parseDouble(attr.getNodeValue());
//                            System.out.println("setting self to " + probs[0]);
                        }

                        if ((attr = attrs.getNamedItem(NEXT)) != null) {
                            probs[1] = Double.parseDouble(attr.getNodeValue());
//                            System.out.println("setting next to " + probs[1]);
                        }

                        transitionFact.setProbabilities(probs);
                    }
                    else if (tag.equals(EMISSION_DIST)) {
                        setEmissionDist(loadDOMDistribution(thisNode));
                    }
                }
            }

            if (nStates == 0) throw new IllegalArgumentException("initialisation was unsuccessful. DOM ill-formed?");
        }
        else throw new IllegalArgumentException("An Element was expected.");
    }

    public Block(int nStates, FiniteAlphabet seq0Alfa, FiniteAlphabet seq1Alfa)
    {
        this.nStates = nStates;
        this.seq0Alfa = seq0Alfa;
        this.seq1Alfa = seq1Alfa;

        initialise();
    }

    private void initialise()
    {
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

    public void setName(String name) { this.name = name; }
    public String getName() { return name; }
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
        State [] subModelStates = new State [nStates];

        /**************************************
         * make the required number of states *
         **************************************/
        State thisState = null;
        for (int j = 0; j < nStates; j++) {
            subModelStates[i++] = thisState = new SimpleEmissionState(
                name + ":match:" + Integer.toString(j),
                userAnnotation,
                matchAdvance,
                emissionDist
            );
            model.addState(thisState);
        }

        // set up topology
        State match = dest;
        i = 0;
        for (int j = 0; j < nStates; j++) {
            thisState = subModelStates[i++];
            model.createTransition(thisState, thisState);
            model.createTransition(thisState, match);
            AtomicSymbol [] dests = new AtomicSymbol[2];
            dests[0] = thisState;
            dests[1] = match;
            Distribution thisDist = transitionFact.getDependentDistribution(model.transitionsFrom(thisState), dests);
            model.setWeights(thisState, thisDist);
            match = thisState;
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
     * use to initialise model parameters.
     */
    public void setModelParameters(double wSelf, double wNext)
    {
        double wTotal = wSelf + wNext;
        double pSelf = wSelf/wTotal;
        double pMatch = 1.0 - pSelf;

        double [] prob = new double [2];
        prob[0] = pSelf;
        prob[1] = pMatch;
        transitionFact.setProbabilities(prob);
    }

    public void dumpParams(XMLWriter xw)
        throws IOException
    {
        xw.openTag(SUBMODEL);
        xw.attribute(NAME, name);
        xw.attribute(TYPE, BLOCK);
        xw.attribute(N_STATES, Integer.toString(nStates));
        xw.attribute(SEQ_0_ALFA, seq0Alfa.getName());
        xw.attribute(SEQ_1_ALFA, seq1Alfa.getName());

        xw.openTag(TRANSITION_DIST);
        double [] probs = transitionFact.getProbabilities();
        xw.attribute(SELF, Double.toString(probs[0]));
        xw.attribute(NEXT, Double.toString(probs[1]));
        xw.closeTag(TRANSITION_DIST);

        // dump the emission Distribution
        if (emissionDist != null) {
            xw.openTag(EMISSION_DIST);
            XMLDistTools.writeDist2XML(emissionDist, xw);
            xw.closeTag(EMISSION_DIST);
        }

        xw.closeTag(SUBMODEL);
    }


    public static void registerMaker(Map map)
    {
        map.put(maker.getKey(), maker);
    }
}

