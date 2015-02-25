
//import org.biojava.bio.dp.MarkovModel;
//import org.biojava.bio.dp.DP;
//import org.biojava.bio.dp.DPFactory;
//import org.biojava.bio.dp.StatePath;
//import org.biojava.bio.dp.ScoreType;
//import org.biojava.bio.dp.twohead.DPInterpreter;

import org.biojava.bio.seq.DNATools;
//import org.biojava.bio.seq.Sequence;
//import org.biojava.bio.seq.SequenceIterator;
//import org.biojava.bio.seq.io.SeqIOTools;

import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

import java.io.PrintWriter;
import java.io.File;
//import java.io.FileReader;
//import java.io.BufferedReader;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

/**
 * This class makes an instance of of the full EnhSearch model from XML definitions of each
 * of its constituent submodels.  It then writes out an XML definition of the
 * full model to STDOUT.  A full EnhSearch model can be used for either immediate
 * searching using the ModelRunner class or further training using the XX class.
 * <br><br>
 * Use as follows:<br>
 * <pre>
 * java MakeEnhSearch &lt;Block submodel file&gt; &lt;Type 1 spacer submodel file&gt; &lt;Type 2 spacer submodel file&gt;
 * </pre>
 * <p>
 * Each parameter specifies an XML file that parameterises a specific submodel in the full EnhSearch model.
 * @author David Huen
 */
public class MakeEnhSearch
{
    public static void main(String [] argv)
        throws Exception
    {
        if (argv.length != 3) {
            System.out.println("Usage: java MakeEnhSearch <match submodel file> <constrained spacer submodel file> <unconstrained spacer submodel file>");
            System.exit(1);
        }

        // load match submodel DOM
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document document = builder.parse( new File(argv[0]) );
        Node matchRoot = document.getDocumentElement();

        // load constrained spacer submodel DOM
        builder = factory.newDocumentBuilder();
        document = builder.parse( new File(argv[1]) );
        Node constrainedRoot = document.getDocumentElement();

        // load unconstrained spacer submodel DOM
        builder = factory.newDocumentBuilder();
        document = builder.parse( new File(argv[2]) );
        Node unconstrainedRoot = document.getDocumentElement();

        // create EnhSearch object
        EnhSearch enh = new EnhSearch(DNATools.getDNA(), DNATools.getDNA());

        enh.setMatchModel(matchRoot);
        enh.setConstrainedSpacer(constrainedRoot);
        enh.setUnconstrainedSpacer(unconstrainedRoot);
        enh.setModelParameters(17.0, 1.0, 0.1);

        // create a writer
        PrintWriter pw = new PrintWriter(System.out);
        XMLWriter xw = new PrettyXMLWriter(pw);

        xw.printRaw("<?xml version=\"1.0\"?>\n");
        enh.dumpParams(xw);
        pw.flush(); pw.close();
    }
}

