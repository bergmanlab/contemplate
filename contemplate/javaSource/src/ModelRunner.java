import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.twohead.DPInterpreter;

//import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

//import org.biojava.utils.xml.PrettyXMLWriter;
//import org.biojava.utils.xml.XMLWriter;

import java.io.*;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import net.derkholm.nmica.utils.CliTools;

/**
 * Class that finds optimal state path through a pair of sequences
 * under an EnhSearch model using the Viterbi algorithm
 * <br><br>
 * Use as follows:-<br>
 * <pre>
 * java ModelRunner &lt;enhSearch.xml&gt; &lt;seq 0&gt; &lt;seq 1&gt; -fastaOut &lt;gapped fasta format filename&gt; -blockListOut &lt;block format filename&gt; -laganOut &lt;lagan format filename&gt;
 * </pre>
 * <p> 
 *
 * @author David Huen, Thomas Down
 * @since 1.4
 */

public class ModelRunner
{
    private PrintStream laganOut;
    private PrintStream blockListOut;
    private PrintStream fastaOut;
    private PrintStream bedOut;

    
    public void setLaganOut(OutputStream os) {
        this.laganOut = new PrintStream(os);
    }
    
    public void setBlockListOut(OutputStream os) {
        this.blockListOut = new PrintStream(os);
    }
    
    public void setFastaOut(OutputStream os) {
        this.fastaOut = new PrintStream(os);
    }
    
    public void setBedOut(OutputStream os) {
        this.bedOut = new PrintStream(os);
    }

    public static void main(String [] argv)
        throws Exception
    {
        ModelRunner app = new ModelRunner();
        argv = CliTools.configureBean(app, argv);
        app.run(argv);
    }
    
    public void run(String[] argv)
        throws Exception
    {
        if (argv.length < 2 || argv.length > 4) {
            System.err.println("Usage: java ModelRunner <EnhSearch XML file> <seq0> <seq1> [<outputFormat>]");
            System.err.println("Supported output formats: BlockList, Lagan, Fasta");
            System.err.println("Example Usage: java ModelRunner <EnhSearch XML file> <seq0> <seq1> -blockListOut <blockFile> -laganOut <laganFile> -fastaOut <fastaFile> -bedOut <bedFile>");
            return;
        }
        
        File modelFile = new File(argv[0]);
        File seq0File = new File(argv[1]);

        // load EnhSearch model DOM
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document document = builder.parse(modelFile);
        Node modelRoot = document.getDocumentElement();

        // create ModelConstructor object
        ModelConstructor constructor = ModelConstructorFactory.getConstructor(modelRoot);

        MarkovModel model = constructor.constructModel("TestModel");

        // read in sequences
        BufferedReader br = new BufferedReader(new FileReader(seq0File));
        SequenceIterator iter = (SequenceIterator)SeqIOTools.fileToBiojava("fasta", "DNA", br);
        Sequence [] seqs = new Sequence [2];
        seqs[0] = iter.nextSequence();
        seqs[1] = iter.nextSequence();

        DPFactory fact = new DPFactory.DefaultFactory(new DPInterpreter.Maker());
        DP aligner = fact.createDP(model);

        StatePath result = aligner.viterbi(seqs, ScoreType.PROBABILITY);

        if (laganOut != null) {
            laganOut.println("# ContemplateLAGAN");
            StatePathOutput.printOutputLaganesque(new DualHeadStatePathTool(result), false, laganOut);
            laganOut.close();
        }
        if (blockListOut != null) {
            blockListOut.println("# ContemplateBlocks - format: source sequence, target start, target end, query start, query end, target blocks, query blocks");
            StatePathOutput.printBlocks(new DualHeadStatePathTool(result), seqs, blockListOut, false);
            blockListOut.close();
        }
        if (bedOut != null) {
            StatePathOutput.printBlocks(new DualHeadStatePathTool(result), seqs, bedOut, true);
            bedOut.close();
        }

        if (fastaOut != null) {
            StatePathOutput.printFasta(new DualHeadStatePathTool(result), fastaOut);
            fastaOut.close();
        }
        
    }
}