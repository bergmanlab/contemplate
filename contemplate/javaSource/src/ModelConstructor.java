/**
 * Interface implemented by factories that generate
 * SubModels
 */

import org.biojava.utils.xml.XMLWriter;
import org.biojava.utils.ChangeVetoException;

import org.biojava.bio.BioException;
//import org.biojava.bio.Annotation;
import org.biojava.bio.dp.MarkovModel;
//import org.biojava.bio.dp.State;
//import org.biojava.bio.symbol.IllegalAlphabetException;
//import org.biojava.bio.symbol.IllegalSymbolException;
import org.w3c.dom.Node;

import java.io.IOException;
//import java.util.Map;

/**
 * Interface implemented by classes that construct models.
 *
 * @author David Huen
 */
public interface ModelConstructor
{

    /**
     * Interface implemented by objects that make
     * instances of this kind of SubModelFactory.
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
         * @param node DOM Node parameterising the instance.
         */
        public ModelConstructor makeConstructor(Node node)
            throws IllegalArgumentException;
    }

    /**
     * Insert the submodel into the specified model
     * between the specified states.
     * @param name name to be assigned to the instance
     * of the model that will be created by a call to
     * this method.
     */
    public MarkovModel constructModel(String name)
        throws ChangeVetoException, BioException;

    /**
     * Dump the parameters for this SubModelConstructor.
     * @param xw XMLWriter instance to use to dump the
     * XML output.
     */
    public void dumpParams(XMLWriter xw)
        throws IOException;
}
