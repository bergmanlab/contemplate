
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.SimpleDotState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
//import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.xml.XMLWriter;

import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * Implements an HMM with one Block submodel and
 * two Spacer submodels, one for constrained (Type1) spacers
 * and another for unconstrained (Type2) spacers.
 * <p>
 * Code is included to reinstantiate an instance
 * of this class from an XML description of its
 * topology and the statistical distributions that
 * specify its behaviour.
 *
 * @author David Huen
 * @since 1.4
 */
public class EnhSearch
    implements ModelConstructor
{
    private String myName = "EnhSearch";

    private SubModelConstructor matchConstructor;
    private SubModelConstructor constrainedConstructor;
    private SubModelConstructor unconstrainedConstructor;

    // left dot state transition dist
    final private CommonDistributionFactory leftFact = new CommonDistributionFactory(2);

    // right dot state transition dist
    final private CommonDistributionFactory rightFact = new CommonDistributionFactory(3);

    // magical state probabilities
    private double pMagicToLeft = 0.8;
    private double pMagicToRight = 0.2;

    private FiniteAlphabet alfa0 = null;
    private FiniteAlphabet alfa1 = null;
    private FiniteAlphabet modelAlfa = null;

    private static final String SUBMODEL = "submodel";
    private static final String MARKOV_MODEL = "markov_model";
    private static final String ENHSEARCH = "enhsearch";
    private static final String TYPE = "type";
    private static final String ALFA0 = "alpha0";
    private static final String ALFA1 = "alpha1";
    private static final String MAGIC = "magicalState";
    private static final String LEFT = "left";
    private static final String RIGHT = "right";
    private static final String LEFT_DOT_DIST = "left_dot";
    private static final String MATCH = "match";
    private static final String END = "end";
    private static final String RIGHT_DOT_DIST = "right_dot";
    private static final String CONSTRAINED = "constrained";
    private static final String UNCONSTRAINED = "unconstrained";
    private static final String MATCH_SUBMODEL = "match_submodel";
    private static final String CONSTRAINED_SPACER_SUBMODEL = "constrained_submodel";
    private static final String UNCONSTRAINED_SPACER_SUBMODEL = "unconstrained_submodel";

    /**
     * Static Instance that registers this class with
     * the ModelConstructorFactory.
     */
    final static Maker maker = new ModelConstructor.Maker()
    {
        public String getKey() { return ENHSEARCH; }
        public ModelConstructor makeConstructor(Node node)
            throws IllegalArgumentException
        {
            return new EnhSearch(node);
        }
    };

    /**
     * Constructor that instantiates an instance of this
     * class from an XML Node.
     * @param paramNode XML Node that describes the Model.
     */
    public EnhSearch(Node paramNode)
        throws IllegalArgumentException
    {
        // we expect a "mismatch chain" Element
        String tag;
        NamedNodeMap attrs;
        Node attr;

        if ((paramNode.getNodeType() == Node.ELEMENT_NODE) &&
            (tag = paramNode.getNodeName()).equals(MARKOV_MODEL) &&
            ((attrs = paramNode.getAttributes()) != null)
            ) {
            // check that this Node is really mine
            if (((attr = attrs.getNamedItem(TYPE)) != null) &&
                !(attr.getNodeValue().equals(ENHSEARCH)) ) {
                throw new IllegalArgumentException("A model of type " + attr.getNodeValue() + " when " + ENHSEARCH + "was expected");
            }

            // pick up the alphabet data
            if ((attr = attrs.getNamedItem(ALFA0)) != null) {
                alfa0 = (FiniteAlphabet) AlphabetManager.alphabetForName(attr.getNodeValue());
            }

            if ((attr = attrs.getNamedItem(ALFA1)) != null) {
                alfa1 = (FiniteAlphabet) AlphabetManager.alphabetForName(attr.getNodeValue());
            }

            // compute the model alphabet
            if ((alfa0 != null) && (alfa1 != null)) {
                modelAlfa = (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(
                    Arrays.asList(new Alphabet [] { alfa0, alfa1 })
                    );
            }
            else new IllegalArgumentException("Alphabets were not specified.");

            // now get child nodes and parse them
            NodeList children;
            if ((children = paramNode.getChildNodes()) != null) {
                for (int i = 0; i < children.getLength(); i++) {
                    Node child = children.item(i);

                    // we only want the Elements
                    String childTag;
                    if ((child.getNodeType() == Node.ELEMENT_NODE) && 
                        ((childTag = child.getNodeName()) != null) ) {
                        if  (childTag.equals(MAGIC)) {
                            // handle magical state attributes
                            Node childAttr;
                            NamedNodeMap childAttrs = child.getAttributes();
                            if (childAttrs != null) {
                                if ((childAttr = childAttrs.getNamedItem(LEFT)) != null) {
                                   pMagicToLeft = Double.parseDouble(childAttr.getNodeValue()); 
                                }

                                if ((childAttr = childAttrs.getNamedItem(RIGHT)) != null) {
                                   pMagicToRight = Double.parseDouble(childAttr.getNodeValue());
                                }
                            }
                        }
                        else if (childTag.equals(LEFT_DOT_DIST)) {
                            // handle magical state attributes
                            Node childAttr;
                            NamedNodeMap childAttrs = child.getAttributes();

                            if (childAttrs != null) {
                                double [] probs = new double [2];

                                if ((childAttr = childAttrs.getNamedItem(MATCH)) != null) {
                                   probs[0] = Double.parseDouble(childAttr.getNodeValue());
                                }

                                if ((childAttr = childAttrs.getNamedItem(END)) != null) {
                                   probs[1] = Double.parseDouble(childAttr.getNodeValue());
                                }

                                leftFact.setProbabilities(probs);
                            }
                        }
                        else if (childTag.equals(RIGHT_DOT_DIST)) {
                            // handle magical state attributes
                            Node childAttr;
                            NamedNodeMap childAttrs = child.getAttributes();

                            if (childAttrs != null) {
                                double [] probs = new double [3];

                                if ((childAttr = childAttrs.getNamedItem(CONSTRAINED)) != null) {
                                   probs[0] = Double.parseDouble(childAttr.getNodeValue());
                                }
                                if ((childAttr = childAttrs.getNamedItem(UNCONSTRAINED)) != null) {
                                   probs[1] = Double.parseDouble(childAttr.getNodeValue());
                                }

                                if ((childAttr = childAttrs.getNamedItem(END)) != null) {
                                   probs[2] = Double.parseDouble(childAttr.getNodeValue());
                                }

                                rightFact.setProbabilities(probs);
                            }
                        }
                        else if (childTag.equals(MATCH_SUBMODEL)) {
                            // get the contained Node and use to create constructor
                            Node submodelNode = child.getFirstChild();
                            boolean notFound = true;
                            while (notFound && (submodelNode != null)) {
                                if ((submodelNode.getNodeType() == Node.ELEMENT_NODE) &&
                                    (submodelNode.getNodeName().equals(SUBMODEL))
                                    ) {
                                    matchConstructor = SubModelConstructorFactory.getConstructor(submodelNode);
                                    notFound = false;
                                }
                                else {
                                    submodelNode = submodelNode.getNextSibling();
                                }
                            }
                        }
                        else if (childTag.equals(CONSTRAINED_SPACER_SUBMODEL)) {
                            // get the contained Node and use to create constructor
                            Node submodelNode = child.getFirstChild();
                            boolean notFound = true;
                            while (notFound && (submodelNode != null)) {
                                if ((submodelNode.getNodeType() == Node.ELEMENT_NODE) &&
                                    (submodelNode.getNodeName().equals(SUBMODEL))
                                    ) {
                                    constrainedConstructor = SubModelConstructorFactory.getConstructor(submodelNode);
                                    notFound = false;
                                }
                                else {
                                    submodelNode = submodelNode.getNextSibling();
                                }
                            }
                        }
                        else if (childTag.equals(UNCONSTRAINED_SPACER_SUBMODEL)) {
                            // get the contained Node and use to create constructor
                            Node submodelNode = child.getFirstChild();
                            boolean notFound = true;
                            while (notFound && (submodelNode != null)) {
                                if ((submodelNode.getNodeType() == Node.ELEMENT_NODE) &&
                                    (submodelNode.getNodeName().equals(SUBMODEL))
                                    ) {
                                    unconstrainedConstructor = SubModelConstructorFactory.getConstructor(submodelNode);
                                    notFound = false;
                                }
                                else {
                                    submodelNode = submodelNode.getNextSibling();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Constructor that builds a minimal instance of this class.
     * @param alfa0 FiniteAlphabet for sequence 0.
     * @param alfa1 FiniteAlphabet for sequence 1.
     */
    public EnhSearch(
        FiniteAlphabet alfa0,
        FiniteAlphabet alfa1
        )
    {
        this.alfa0 = alfa0;
        this.alfa1 = alfa1;
        modelAlfa = (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(
            Arrays.asList(new Alphabet [] { alfa0, alfa1 })
            );
    }

    /**
     * set the name of this object.
     */
    public void setName(String name) { myName = name; }

    /**
     * Initialise the match submodel from the specified node.
     * @param node XML node that defines the match submodel.
     */
    public void setMatchModel(Node node)
    {
        matchConstructor = SubModelConstructorFactory.getConstructor(node);
    }

    /**
     * Initialise the constrained spacer submodel from the specified node.
     * @param node XML node that defines the constrained submodel.
     */
    public void setConstrainedSpacer(Node node)
    {
        constrainedConstructor = SubModelConstructorFactory.getConstructor(node);
    }

    /**
     * Initialise the unconstrained spacer submodel from the specified node.
     * @param node XML node that defines the unconstrained submodel.
     */
    public void setUnconstrainedSpacer(Node node)
    {
        unconstrainedConstructor = SubModelConstructorFactory.getConstructor(node);
    }

    /**
     * Initialise the relative weights for transitions to the different
     * spacer chains and to termination.  The weights are normalised to
     * unity before the absolute weights are computed.
     * <p>
     * @param wConstrained weight of transition to constrained submodel.
     * @param wUnconstrained weight of transition to unconstrained submodel.
     * @param wEnd weight of transition to end state.
     */
    public void setModelParameters(double wConstrained, double wUnconstrained, double wEnd)
    {
        double wTotal = wConstrained + wUnconstrained + wEnd;

        double pConstrained = wConstrained/wTotal;
        double pUnconstrained = wUnconstrained/wTotal;
        double pEnd = wEnd/wTotal;

        // left dot state
        double [] probs = new double [2];
        probs[0] = 1.0 - pEnd;
        probs[1] = pEnd;
        leftFact.setProbabilities(probs);

        // right dot state
        probs = new double [3];
        probs[0] = pConstrained;
        probs[1] = pUnconstrained;
        probs[2] = pEnd;
        rightFact.setProbabilities(probs);
    }

    /**
     * Construct the complete model from the parameters specified earlier.
     */
    public MarkovModel constructModel(String name)
        throws ChangeVetoException, BioException
    {
        MarkovModel model = new SimpleMarkovModel(2, modelAlfa, name);

        // check that all submodel constructors are in place
        if ((matchConstructor == null) ||
            (constrainedConstructor == null) ||
            (unconstrainedConstructor == null)
            ) {
            throw new BioException("Factory is insufficiently parameterised.");
        }

        // set up model
        State left = new SimpleDotState(name + ":left");
        model.addState(left);
        State right = new SimpleDotState(name + ":right");
        model.addState(right);

        // set up the annotation objects
        Annotation matchAnn = new SmallAnnotation();
        Annotation constrainedAnn = new SmallAnnotation();
        Annotation unconstrainedAnn = new SmallAnnotation();
        StatePathOutput.setInteresting(matchAnn);
        StatePathOutput.setInteresting(constrainedAnn);
        StatePathOutput.setMatchClass(matchAnn);
        StatePathOutput.setConstrainedClass(constrainedAnn);
        StatePathOutput.setUnconstrainedClass(unconstrainedAnn);
        Annotation ignoreAnn = new SmallAnnotation();

        // set up submodels
        State matchStart = matchConstructor.insertSubModel(model, left, right, matchAnn);
        State constrainedStart = constrainedConstructor.insertSubModel(model, right, left, constrainedAnn);
        State unconstrainedStart = unconstrainedConstructor.insertSubModel(model, right, left, unconstrainedAnn);

        // define transitions
        model.createTransition(left, model.magicalState());
        State [] dests = new State [2];
        dests[0] = matchStart;
        dests[1] = model.magicalState();
        Distribution leftDist = leftFact.getDependentDistribution(model.transitionsFrom(left), dests);
        model.setWeights(left, leftDist);

        
        model.createTransition(right, model.magicalState());
        dests = new State [3];
        dests[0] = constrainedStart;
        dests[1] = unconstrainedStart;
        dests[2] = model.magicalState();
        Distribution rightDist = rightFact.getDependentDistribution(model.transitionsFrom(right), dests);
        model.setWeights(right, rightDist);

        model.createTransition(model.magicalState(), left);
        model.createTransition(model.magicalState(), right);
        Distribution dist = model.getWeights(model.magicalState());
        dist.setWeight(left, pMagicToLeft);
        dist.setWeight(right, pMagicToRight);
        
        return model;
    }

    /**
     * Write out the parameters required for reconstructing this instance in XML form.
     */
    public void dumpParams(XMLWriter xw)
        throws IOException
    {
        xw.openTag(MARKOV_MODEL);
        xw.attribute(TYPE, ENHSEARCH);
        xw.attribute(ALFA0, alfa0.getName());
        xw.attribute(ALFA1, alfa1.getName());

        xw.openTag(MAGIC);
        xw.attribute(LEFT, Double.toString(pMagicToLeft));
        xw.attribute(RIGHT, Double.toString(pMagicToRight));
        xw.closeTag(MAGIC);

        xw.openTag(LEFT_DOT_DIST);
        double [] probs = leftFact.getProbabilities();
        xw.attribute(MATCH, Double.toString(probs[0]));
        xw.attribute(END, Double.toString(probs[1]));
        xw.closeTag(LEFT_DOT_DIST);

        xw.openTag(RIGHT_DOT_DIST);
        probs = rightFact.getProbabilities();
        xw.attribute(CONSTRAINED, Double.toString(probs[0]));
        xw.attribute(UNCONSTRAINED, Double.toString(probs[1]));
        xw.attribute(END, Double.toString(probs[2]));
        xw.closeTag(RIGHT_DOT_DIST);

        xw.openTag(MATCH_SUBMODEL);
        matchConstructor.dumpParams(xw);
        xw.closeTag(MATCH_SUBMODEL);

        xw.openTag(CONSTRAINED_SPACER_SUBMODEL);
        constrainedConstructor.dumpParams(xw);
        xw.closeTag(CONSTRAINED_SPACER_SUBMODEL);

        xw.openTag(UNCONSTRAINED_SPACER_SUBMODEL);
        unconstrainedConstructor.dumpParams(xw);
        xw.closeTag(UNCONSTRAINED_SPACER_SUBMODEL);

        xw.closeTag(MARKOV_MODEL);
    }

    /**
     * Register this class with the ModelConstructorFactory.
     */
    public static void registerMaker(Map map)
    {
        map.put(maker.getKey(), maker);
    }
}
