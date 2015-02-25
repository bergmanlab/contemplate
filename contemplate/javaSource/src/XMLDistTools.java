
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.XMLDistributionReader;

import org.biojava.utils.xml.XMLWriter;

import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.FiniteAlphabet;

import org.biojava.bio.BioError;

import org.w3c.dom.Node;

import org.xml.sax.SAXException;

import java.io.IOException;
import java.util.Iterator;

public class XMLDistTools
{
    public static void writeDist2XML(Distribution d, XMLWriter xw) 
        throws IOException
    {
        xw.openTag("Distribution");
        xw.attribute("type", "Distribution");

        xw.openTag("alphabet");
        xw.attribute("name", d.getAlphabet().getName());
        xw.closeTag("alphabet");

        for (Iterator i = ((FiniteAlphabet) d.getAlphabet()).iterator();
            i.hasNext();) {
            Symbol sym = (Symbol) i.next();
            double weight = 0.0;

            try {
                weight = d.getWeight(sym);
            } catch (IllegalSymbolException ex) {
                throw new BioError("Distribution has been built with Illegal Symbols !?", ex);
            }

            xw.openTag("weight");
            xw.attribute("sym", sym.getName());
            xw.attribute("prob", Double.toString(weight));
            xw.closeTag("weight");

        }

        xw.closeTag("Distribution");

    } //end writeXML

    public static Distribution readDistFromDOM(Node root)
        throws SAXException
    {
        XMLDistributionReader handler = new XMLDistributionReader();

        DOM2SAXConverter.convertToSAX(root, handler);

        return handler.getDist();
    }
}

