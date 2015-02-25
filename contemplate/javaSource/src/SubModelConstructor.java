
/**
 * Interface implemented by factories that generate
 * SubModels.
 * @author David Huen
 */

import org.biojava.utils.xml.XMLWriter;
import org.biojava.utils.ChangeVetoException;
import org.biojava.bio.Annotation;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.w3c.dom.Node;

import java.io.IOException;
//import java.util.Map;

public interface SubModelConstructor
{

    /**
     * Interface implemented by clasess that make
     * instances that implement the SubModelFactory
     * interface.
     */
    public interface Maker
    {
        /**
         * Get the key associated with this ConstructorMaker.
         */
        public String getKey();

        /**
         * Make an instance of the SubModelFactory
         * that makes SubModels with the properties
         * defined by the node.
         */
        public SubModelConstructor makeConstructor(Node node)
            throws IllegalArgumentException;
    }

    /**
     * Get the name of this instance.
     */
    public String getName();

    /**
     * Set the name of this instance.
     */
    public void setName(String name);

    /**
     * Sets the class to create a training version of the
     * the submodel.  This is often a no-op method.
     */
    public void setForTraining(boolean forTraining);

    /**
     * Insert the submodel into the specified model
     * between the specified states.
     */
    public State insertSubModel(
        MarkovModel model,
        State source,
        State dest,
        Annotation userAnnotation
        )
        throws ChangeVetoException, IllegalSymbolException, IllegalAlphabetException;

    /**
     * Dump the parameters for this SubModelConstructor.
     */
    public void dumpParams(XMLWriter xw)
        throws IOException;
}
