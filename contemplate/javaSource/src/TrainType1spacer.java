import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

import java.io.PrintWriter;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
//import java.io.PrintStream;

//import java.util.Iterator;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import javax.xml.parsers.DocumentBuilder; 
import javax.xml.parsers.DocumentBuilderFactory; 

import org.biojava.bio.BioException;
import org.biojava.bio.Annotation;

import org.biojava.bio.dist.Distribution;

import org.biojava.bio.dp.State;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.TrainingSet;
import org.biojava.bio.dp.SimpleTrainingSet;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.ViterbiTrainer;
import org.biojava.bio.dp.twohead.DPInterpreter;
//import org.biojava.bio.dp.XmlMarkovModel;

import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;

import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;

import org.biojava.bio.seq.io.SeqIOTools;

/**
 * This class is used to train Type1spacer objects using the Viterbi algorithm that were previously
 * created and parameterised by MakeType1spacer. On completion of training, the instance will write the parameters that resulted
 * from training to STDOUT.  These can be saved and used to reinstantiate a trained
 * model for further use.
 * <br><br>
 * Use as follows:
 * <pre>
 * java TrainType1spacer &lt;Type1spacer submodel XML filename&gt; &lt;training set filename&gt; &lt;min length for training&gt;
 * </pre>
 * <p>
 * The training set is a FASTA-formatted file in which successive pairs of sequences
 * form pairs for model training.
 * <p>
 * A lower limit can be set for sequence lengths that can be used for training.  This is
 * important in cases where the model cannot emit less than a certain number of
 * symbols because of its topology.  Attempting to train with such a set will
 * lead to an Exception being thrown and training will then fail.
 *
 * @author David Huen
 */
public class TrainType1spacer
{
    public static void main(String [] argv)
        throws Exception
    {
        if (argv.length != 3) {
            System.out.println("Usage: java TrainType1spacer <chain parameters filename> <training set filename> <min length for training>");
            System.exit(1);
        }

        // create a DOM factory
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();

        // create a document
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document document = builder.parse( new File(argv[0]) );

        // navigate to the root node of the chain parameters
        Node root = document.getDocumentElement();

        // create a simple chain
        Type1spacer thisType1spacer = new Type1spacer(root);

        // create a HMM
        MarkovModel model = new SimpleMarkovModel(2, AlphabetManager.generateCrossProductAlphaFromName("(DNA x DNA)"), "TestModel");

        // insert chain into model
        thisType1spacer.setForTraining(true);
        State start = thisType1spacer.insertSubModel(model, model.magicalState(), model.magicalState(), Annotation.EMPTY_ANNOTATION);
        Distribution sourceDist = model.getWeights(model.magicalState());
        sourceDist.setWeight(start, 1.0);

        // create the DP for this task
        DPFactory fact = new DPFactory.DefaultFactory(new DPInterpreter.Maker());
        DP dp = fact.createDP(model);

        // load a training set
        TrainingSet ts = loadTrainingSet(argv[1], Integer.parseInt(argv[2]));

        // create an instance of the trainer
        ViterbiTrainer vt = new ViterbiTrainer(dp);
        vt.train(ts, 1.0, new Stopper(0.00001, 30));

        // create a writer
        PrintWriter pw = new PrintWriter(System.out);
        XMLWriter xw = new PrettyXMLWriter(pw);
        xw.printRaw("<?xml version=\"1.0\"?>\n");
        thisType1spacer.dumpParams(xw);

        pw.flush(); pw.close();
    }

    public static TrainingSet loadTrainingSet(String filename, int minLength)
        throws FileNotFoundException, BioException
    {
        // create a 2-head training set for DNA
        FiniteAlphabet [] compAlfas = new FiniteAlphabet[2];
        compAlfas[0] = compAlfas[1] = DNATools.getDNA();
        SimpleTrainingSet chainTS = new SimpleTrainingSet(2, compAlfas);

        // fill in training set: file has sequential pairs of FASTA sequences
        BufferedReader br = new BufferedReader(new FileReader(filename));

        SequenceIterator iter = (SequenceIterator)SeqIOTools.fileToBiojava("fasta", "DNA", br);

        int count = 0;
        Sequence [] seqs = new Sequence[2];

        boolean lengthOK = true;

        while (iter.hasNext()) {
            Sequence thisSeq;
            seqs[count++] = thisSeq = iter.nextSequence();

            if (thisSeq.length() < minLength) lengthOK = false;

            // when we have a pair of sequences, add to training set
            if (count == 2) {
                if (lengthOK) {
//                    System.out.println("adding to training set: " + seqs[0].getName() + " " + seqs[1].getName());
                    chainTS.addTrainingCase(seqs);
                }
                count = 0;
                lengthOK = true;
            }
        }

        return chainTS;
    }
}
