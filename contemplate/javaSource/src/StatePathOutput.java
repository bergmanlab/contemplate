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

import java.util.AbstractCollection;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.xml.XMLWriter;

import org.biojava.bio.Annotation;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.GappedSymbolList;
//import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.DotState;
import org.biojava.bio.dp.EmissionState;

import com.sun.xml.internal.ws.util.StringUtils;
//import org.biojava.bio.dp.twohead.*;

import java.io.IOException;
import java.io.PrintStream;

/**
 *
 * @author David Huen
 */
public class StatePathOutput
{
	private static final int CONTEMPLATEFORMAT = 1;
	private static final int FORMAT_BEDSTARTS   = 2;
	private static final int FORMAT_BEDSIZES   = 3;
	private static final int FORMAT_BEDLENGTH   = 4;
	
	private static final Object OF_INTEREST = new Object();
    private static final SubModelClass MATCH_BLOCK = new SubModelClass();
    private static final SubModelClass CONSTRAINED_SPACER = new SubModelClass();
    private static final SubModelClass UNCONSTRAINED_SPACER = new SubModelClass();
    private static final String INTERESTING_STATE = "of_interest";
    private static final String SUBMODEL_CLASS = "submodel_class";


    /***************************************************************
     * Pretty-print the alignment that is defined by the StatePath *
     ***************************************************************/
    public static void printOutput(DualHeadStatePathTool path)
        throws IllegalSymbolException
    {
        int start, end;
        start = end = 0;
        boolean extendingBlock = false;

        // we want to find all high-homology blocks
        // and print them
        for (int i=1; i <= path.length(); i++) {
            // retrieve the state label
//            String stateLabel = path.getStateLabelAt(i).split(":",2)[0];
            State state = path.getStateAt(i);
//            System.out.println(i + " " + stateLabel + " " + seqSymbols.symbolAt(i)); // uncomment to dump path
//            if (i >9484) System.out.println(i + " " + stateLabel + " " + extendingBlock);

            if (extendingBlock) {
                // when extending, we continue to do so as long as the State is one
                // of the desired ones or that the current block is one that should
                // be ignored.
//                if (stateLabel.equals("ignore")) continue;
                if (state instanceof DotState) continue;
                if (wantedState(state)) {
                    end = i;
                    continue;
                }

                // extension terminated, print output
                if ((end-start) > 4) {
                    System.out.println();
                    System.out.println("new block");
                    displayBlock(path, start, end);
                }

                // reset
                start = end = i;
                extendingBlock = false;
            }
            else {
                if (wantedState(state)) {
                    extendingBlock = true;
                    start = end = i;
                }
            }
        }

        // deal with final loose end
        if (extendingBlock) {
            System.out.println();
            System.out.println("new block");
            displayBlock(path, start, end);
        }
    }
    
    public static void printBlocks(DualHeadStatePathTool path, Sequence[] seqs) 
        throws IllegalSymbolException
    {
        printBlocks(path, seqs, System.out, false);
    }
    
    public static void printBlocks(DualHeadStatePathTool path, Sequence[] seqs, PrintStream out, boolean bedFormat) 
        throws IllegalSymbolException
    {
        int blockStart = -1;
        for (int i = 1; i < path.length(); ++i) {
            State state = path.getStateAt(i);
            if (state instanceof EmissionState) {
                boolean isMatch = isMatchState(state);
                boolean isUnconstrained = isUnconstrainedState(state);
                
                if (isMatch && blockStart < 0) {
                    blockStart = i;
                } else if (isUnconstrained && blockStart > 0) {
                    printBlock(path, blockStart, i, seqs, out, bedFormat);
                    blockStart = -1;
                }
            }
        }
        if (blockStart > 0) {
            printBlock(path, blockStart, path.length(), seqs, out, bedFormat);
        }
    }
    
    private static void printBlock(DualHeadStatePathTool path, int firstMatch, int upperBound, Sequence[] seqs, PrintStream out, boolean bedFormat) 
        throws IllegalSymbolException
    {
        int lastMatch = upperBound;
        while (!isMatchState(path.getStateAt(lastMatch))) {
            --lastMatch;
        }
        // modified by Max to output bed format, also had to change getMatchBlockString & appendMatchBlock
        if (bedFormat) {
           int featStart = path.getComponentSourceCoord(firstMatch, 0) -1;
           int featEnd = path.getComponentSourceCoord(lastMatch, 0);
     	   out.println(
 			   seqs[0].getName()+"\t"+ 
 			   Integer.toString(featStart)+"\t"+ 
 			   Integer.toString(featEnd)+"\t"+
 			   seqs[0].getName()+"\t"+
 			   "0"+"\t"+
 			   "+"+"\t"+
 			   Integer.toString(featStart)+"\t"+
 			   Integer.toString(featEnd)+"\t"+
 			   "50,50,50\t"+
 			   getMatchBlockString(path, firstMatch, lastMatch, 0, featStart, FORMAT_BEDLENGTH)+"\t"+
               getMatchBlockString(path, firstMatch, lastMatch, 0, featStart, FORMAT_BEDSIZES)+"\t"+
               getMatchBlockString(path, firstMatch, lastMatch, 0, featStart, FORMAT_BEDSTARTS));
        } else {
            out.println(
            		seqs[0].getName()+"\t"+ 
                    path.getComponentSourceCoord(firstMatch, 0) + "\t" +
                    path.getComponentSourceCoord(lastMatch, 0) + "\t" +
                    path.getComponentSourceCoord(firstMatch, 1) + "\t" +
                    path.getComponentSourceCoord(lastMatch, 1) + "\t" +
                    getMatchBlockString(path, firstMatch, lastMatch, 0, 0, CONTEMPLATEFORMAT) + "\t" +
                    getMatchBlockString(path, firstMatch, lastMatch, 1, 0, CONTEMPLATEFORMAT)
                    );
	        }
    }
    
    private static String getMatchBlockString(DualHeadStatePathTool path, int firstMatch, int lastMatch, int seq, int featStart, int format) 
        throws IllegalSymbolException
    {
        StringBuffer sb = new StringBuffer();
        int mbStart = -1;
        for (int i = firstMatch; i <= lastMatch; ++i) {
            boolean isMatch = isMatchState(path.getStateAt(i));
            if (isMatch && mbStart < 0) {
                mbStart = i;
            } else if (!isMatch && mbStart >= 0) {
                appendMatchBlock(sb, path, mbStart, i - 1, seq, format, featStart);
                mbStart = -1;
            }
        }
        if (mbStart >=0 ) {
            appendMatchBlock(sb, path, mbStart, lastMatch, seq, format, featStart);
        }
        if (format==FORMAT_BEDLENGTH)
        	return Integer.toString(sb.toString().split(",").length);
        else
        	return sb.toString();
    }
    
    private static void appendMatchBlock(StringBuffer sb, DualHeadStatePathTool path, int firstMatch, int upperBound, int seq, int format, int featStart)
        throws IllegalSymbolException
    {
        int lastMatch = upperBound;
        while (!isMatchState(path.getStateAt(lastMatch))) {
            --lastMatch;
        }
        
        if (sb.length() > 0) {
        		sb.append(',');
        }
        
    	if (format==FORMAT_BEDSTARTS || format==FORMAT_BEDLENGTH) 
    		sb.append(path.getComponentSourceCoord(firstMatch, seq)-1-featStart);
    	else if (format==FORMAT_BEDSIZES) 
    		sb.append(path.getComponentSourceCoord(lastMatch, seq)-path.getComponentSourceCoord(firstMatch, seq)+1);
    	else {
	        sb.append(path.getComponentSourceCoord(firstMatch, seq));
	        sb.append('-');
	        sb.append(path.getComponentSourceCoord(lastMatch, seq));
    	}
    }
    
    /**
     * Print an alignment in a LAGAN-like format
     */
    
    public static void printOutputLaganesque(DualHeadStatePathTool path, boolean wantStates)
        throws IllegalSymbolException
    { 
        printOutputLaganesque(path, wantStates, System.out);
    }
        
    public static void printOutputLaganesque(DualHeadStatePathTool path, boolean wantStates, PrintStream out)
        throws IllegalSymbolException
    {
        int start = 1;
        int end = path.length();
        int width = 60;
        
        StringBuffer line1, line2, comp, ruler1, ruler2, states;

        line1 = new StringBuffer();
        line2 = new StringBuffer();
        ruler1 = new StringBuffer();
        ruler2 = new StringBuffer();
        comp = new StringBuffer();
        states = new StringBuffer();
        int rowStart = start;
        int rowEnd = start;

        int count=0;
        for (int i=start; i <= end; i++) {
            // get the symbols
//            BasisSymbol sym = (BasisSymbol) path.getAlignedSymbolsAt(i);

            // retrieve the state at this point
            State state = path.getStateAt(i);

            // DotStates do not contribute anything useful to this output
            if (state instanceof DotState) {
                if (count == 0) {
                    rowEnd = rowStart++;
                }
            }
            else {
                // get the paired symbols that make this symbol
                List syms = path.getAlignedSymbolsAt(i);
//                if (i > 2150) System.out.println("doing symbol " + i + " " + );

                char sym1  = DNATools.dnaToken(((Symbol) syms.get(0)));
                char sym2  = DNATools.dnaToken(((Symbol) syms.get(1)));

                // assemble output lines
                if (wantedState(state)) {
                    line1.append(("" + sym1).toUpperCase());
                    line2.append(("" + sym2).toUpperCase());
                } else {
                    line1.append(sym1);
                    line2.append(sym2);
                }
                
                if (sym1 != '-') {
                    int seq1pos = path.getComponentSourceCoord(i, 0);
                    if ((seq1pos % 10) == 0) {
                        String rString = "" + seq1pos;
                        int rulerSpaces = line1.length() - ruler1.length() - rString.length();
                        for (int s = 0; s < rulerSpaces; ++s) {
                            ruler1.append(' ');
                        }
                        ruler1.append(rString);
                    }
                }
                if (sym2 != '-') {
                    int seq2pos = path.getComponentSourceCoord(i, 1);
                    if ((seq2pos % 10) == 0) {
                        String rString = "" + seq2pos;
                        int rulerSpaces = line2.length() - ruler2.length() - rString.length();
                        for (int s = 0; s < rulerSpaces; ++s) {
                            ruler2.append(' ');
                        }
                        ruler2.append(rString);
                    }
                }

                if (sym1 == sym2) {
                    comp.append('|');
                }
                else {
                    comp.append(' ');
                }
                
                Object stateClass = state.getAnnotation().getProperty(SUBMODEL_CLASS);
                if (stateClass == CONSTRAINED_SPACER) {
                    states.append('C');
                } else if (stateClass == UNCONSTRAINED_SPACER) {
                    states.append('U');
                } else if (stateClass == MATCH_BLOCK) {
                    states.append('M');
                } else {
                    states.append(' ');
                }

                count++;
                rowEnd = i;
            }

            if (count >= width) {
                // get coordinates
                int seq1Start = path.getComponentSourceCoord(rowStart, 0);
                int seq2Start = path.getComponentSourceCoord(rowStart, 1);
                int seq1End = path.getComponentSourceCoord(rowEnd, 0);
                int seq2End = path.getComponentSourceCoord(rowEnd, 1);

                
                out.println("         " + ruler1);
                out.println("seq1     " + line1);
                if (wantStates) {
                    out.println("         " + states);
                } else {
                    out.println("         " + comp);
                }
                out.println("seq2     " + line2);
                out.println("         " + ruler2);
                out.println();

                // reset things
                rowStart = i+1;
                rowEnd = rowStart;
                line1 = new StringBuffer();
                line2 = new StringBuffer();
                comp = new StringBuffer();
                ruler1 = new StringBuffer();
                ruler2 = new StringBuffer();
                states = new StringBuffer();
                count = 0;
            }
        }

        if (count != 0) {
            out.println("         " + ruler1);
            out.println("seq1     " + line1);
            if (wantStates) {
                out.println("         " + states);
            } else {
                out.println("         " + comp);
            }
            out.println("seq2     " + line2);
            out.println("         " + ruler2);
            out.println();
        }
    }
    
    public static void printFasta(DualHeadStatePathTool path, PrintStream out)
        throws IOException, IllegalSymbolException
    {
        new FastaFormat().writeSequence(
            new SimpleSequence(
                path.getAlignmentComponentSymbolList(0),
                null,
                "seq0",
                Annotation.EMPTY_ANNOTATION
            ),
            out
        );
        new FastaFormat().writeSequence(
            new SimpleSequence(
                path.getAlignmentComponentSymbolList(1),
                null,
                "seq1",
                Annotation.EMPTY_ANNOTATION
            ),
            out
        );
    }

    public static void dumpStatePath(DualHeadStatePathTool path, int start, int end)
        throws IllegalSymbolException
    {
        if (end == -1) end = path.length();

        for (int i=start; i <= end; i++) {

            State state = path.getStateAt(i);
            String label = state.getName();

            if (!(state instanceof DotState)) {

                List syms = path.getAlignedSymbolsAt(i);

                if (syms != null) {
                    char sym1 = DNATools.dnaToken(((Symbol) syms.get(0)));
                    char sym2 = DNATools.dnaToken(((Symbol) syms.get(1)));

                    System.out.println(i + " " + label + " " + sym1 + " " + sym2);
                }
            }
            else {
                System.out.println(i + " " + label);
            }

        }
    }

    /**
     * Indicates that the State should form part of the output
     */
    private static boolean wantedState(State state)
    {
        Annotation ann;
        if ((ann = state.getAnnotation()) != null &&
            state.getAnnotation().containsProperty(INTERESTING_STATE) &&
            (state.getAnnotation().getProperty(INTERESTING_STATE) == OF_INTEREST) 
            ) {
            return true;
        }
        else return false;
    }

    private static boolean isConstrainedState(State state) {
        Annotation ann = state.getAnnotation();
        return ann.containsProperty(SUBMODEL_CLASS) && ann.getProperty(SUBMODEL_CLASS) == CONSTRAINED_SPACER;
    }
    
    private static boolean isUnconstrainedState(State state) {
        Annotation ann = state.getAnnotation();
        return ann.containsProperty(SUBMODEL_CLASS) && ann.getProperty(SUBMODEL_CLASS) == UNCONSTRAINED_SPACER;
    }
    
    private static boolean isMatchState(State state) {
        Annotation ann = state.getAnnotation();
        return ann.containsProperty(SUBMODEL_CLASS) && ann.getProperty(SUBMODEL_CLASS) == MATCH_BLOCK;
    }    
    
    private static void displayBlock(DualHeadStatePathTool path, int start, int end)
        throws IllegalSymbolException
    {
        int width = 50;
        StringBuffer line1, line2, comp;

        line1 = new StringBuffer();
        line2 = new StringBuffer();
        comp = new StringBuffer();
        int rowStart = start;
        int rowEnd = start;

        int count=0;
        for (int i=start; i <= end; i++) {
            // get the symbols
//            BasisSymbol sym = (BasisSymbol) path.getAlignedSymbolsAt(i);

            // retrieve the state at this point
            State state = path.getStateAt(i);

            // DotStates do not contribute anything useful to this output
            if (state instanceof DotState) {
                if (count == 0) {
                    rowEnd = rowStart++;
                }
            }
            else {
                // get the paired symbols that make this symbol
                List syms = path.getAlignedSymbolsAt(i);
//                if (i > 2150) System.out.println("doing symbol " + i + " " + );

                char sym1  = DNATools.dnaToken(((Symbol) syms.get(0)));
                char sym2  = DNATools.dnaToken(((Symbol) syms.get(1)));

                // assemble output lines
                line1.append(sym1);
                line2.append(sym2);

                if (sym1 == sym2) {
                    comp.append('|');
                }
                else {
                    comp.append(' ');
                }

                count++;
                rowEnd = i;
            }

            if (count >= width) {
                // get coordinates
                int seq1Start = path.getComponentSourceCoord(rowStart, 0);
                int seq2Start = path.getComponentSourceCoord(rowStart, 1);
                int seq1End = path.getComponentSourceCoord(rowEnd, 0);
                int seq2End = path.getComponentSourceCoord(rowEnd, 1);

                
                /* David's output format
                System.out.println(rightJustify(Integer.toString(seq1Start), 7) + " " + line1 + " " + seq1End);
                System.out.println(rightJustify(Integer.toString(rowStart), 7) + " " + comp  + " " + i);
                System.out.println(rightJustify(Integer.toString(seq2Start), 7) + " " + line2 + " " + seq2End);
                System.out.println();
                
                */
                
                System.out.println(rightJustify(Integer.toString(Math.abs(seq1Start)), 7) + " " + line1 + " " + absDown(seq1End));
                System.out.println("        " + comp);
                System.out.println(rightJustify(Integer.toString(Math.abs(seq2Start)), 7) + " " + line2 + " " + absDown(seq2End));
                System.out.println();

                // reset things
                rowStart = i+1;
                rowEnd = rowStart;
                line1 = new StringBuffer();
                line2 = new StringBuffer();
                comp = new StringBuffer();
                count = 0;
            }
        }

        if (count != 0) {
            // print this line
//            System.out.println("rowStart, rowEnd are " + rowStart + " " + rowEnd);
            System.out.println(
                rightJustify(Integer.toString(Math.abs(path.getComponentSourceCoord(rowStart,0))), 7)
                + " " + line1 + " " + absDown(path.getComponentSourceCoord(rowEnd, 0))
                );
            System.out.println(
                "        " + comp);
            System.out.println(
                rightJustify(Integer.toString(Math.abs(path.getComponentSourceCoord(rowStart,1))), 7)
                + " " + line2 + " " + absDown(path.getComponentSourceCoord(rowEnd, 1))
                );
        }
    }

    private static int absDown(int i) {
        if (i > 0) {
            return i;
        } else {
            return -i - 1;
        }
    }
    
    static String leftJustify(String orig, int length)
    {
        StringBuffer buffer;

        // check the number of chars to be padded on
        int padCount = length - orig.length();

        if (padCount > 0) {
            buffer = new StringBuffer(orig);

            for (int i=0; i < padCount; i++) {
                buffer.append(' ');
            }
            return buffer.toString();
        }
        else if (padCount < 0) {
            return orig.substring(0, length - 1);
        }
        else return orig;
    }

    static String rightJustify(String orig, int length)
    {
        StringBuffer buffer;

        // check the number of chars to be padded on
        int padCount = length - orig.length();

        if (padCount > 0) {
            buffer = new StringBuffer(length);

            for (int i=0; i < padCount; i++) {
                buffer.append(' ');
            }
            buffer.append(orig);

            return buffer.toString();
        }
        else if (padCount < 0) {
            return orig.substring(0, length - 1);
        }
        else return orig;
    }

    /***********************
     * dump matches as XML *
     ***********************/
    public static void generateMatchesAsXML(
        DualHeadStatePathTool path,
        XMLWriter xw,
        String seq0Name,
        String seq1Name
        )
        throws IllegalSymbolException, IOException
    {
        xw.printRaw("<?xml version=\"1.0\"?>\n");
        xw.openTag("alignment");
        xw.attribute("source", seq0Name);
        xw.attribute("target", seq1Name);

        // add code to put the sequence names in!!!!

        DualHeadSequenceBuilder builder = new DualHeadSequenceBuilder(path, seq0Name, seq1Name);

        for (int i=1; i <= path.length(); i++) {
            // add to current accumulator until done
            if (!builder.addStatePathAtIndex(i)) {

                printXML(builder, xw);

                // start a new accumulator
                builder = new DualHeadSequenceBuilder(path, seq0Name, seq1Name);
                builder.addStatePathAtIndex(i);
            }
        }

        printXML(builder, xw);

        xw.closeTag("alignment");
    }

    private static void printXML(
        DualHeadSequenceBuilder builder,
        XMLWriter xw
        )
        throws IllegalSymbolException, IOException
    {
        if (builder.hasSequence()) {
            if (builder.subModelClass() == MATCH_BLOCK) {
                builder.writeXML(xw);
            }
        }
    }

    /**
     * Generate FASTA files that contain sequences that can be used to train 
     * the alignment model.                                                  
     */
    public static void generateTrainingSequences(
        DualHeadStatePathTool path,
        PrintStream matchStream,
        PrintStream constrainedSpacerStream,
        PrintStream unconstrainedSpacerStream
        )
        throws IllegalSymbolException, IOException
    {
//        Alphabet modelAlfa = path.getAlignedSymbols().getAlphabet();
        // we collect emissions from all States that are marked
        // with a specific submodel class.  As long as the class
        // of the next State remains unchanged, we continue to
        // accumulate the States into the same SymbolLists.
        DualHeadSequenceBuilder builder = new DualHeadSequenceBuilder(path, "", "");

        for (int i=1; i <= path.length(); i++) {
            // add to current accumulator until done
            if (!builder.addStatePathAtIndex(i)) {

                printSequences(builder, matchStream, constrainedSpacerStream, unconstrainedSpacerStream);

                // start a new accumulator
                builder = new DualHeadSequenceBuilder(path, "", "");
                builder.addStatePathAtIndex(i);
            }
        }

        // handle last accumulator
        printSequences(builder, matchStream, constrainedSpacerStream, unconstrainedSpacerStream);
    }

    private static void printSequences(
        DualHeadSequenceBuilder builder,
        PrintStream matchStream,
        PrintStream constrainedSpacerStream,
        PrintStream unconstrainedSpacerStream
        )
        throws IllegalSymbolException, IOException
    {
        if (builder.hasSequence()) {
            if (builder.subModelClass() == MATCH_BLOCK) {
                builder.writeSequence(0, "", matchStream);
                builder.writeSequence(1, "", matchStream);
            }
            else if (builder.subModelClass() == CONSTRAINED_SPACER) {
                builder.writeSequence(0, "", constrainedSpacerStream);
                builder.writeSequence(1, "", constrainedSpacerStream);
            }
            else if (builder.subModelClass() == UNCONSTRAINED_SPACER) {
                builder.writeSequence(0, "", unconstrainedSpacerStream);
                builder.writeSequence(1, "", unconstrainedSpacerStream);
            }
            else System.out.println("unexpected SubModelClass!");
        }
    }

    private static class DualHeadSequenceBuilder
    {
        private DualHeadStatePathTool path;
        private SubModelClass myClass = null;
        private boolean extending = false;
        private final Alphabet alfa;
        private final List seqAlfa;
        private boolean hasSequence = false;
        private List [] seqSymList;
        private int [] min = new int [2];
        private int [] max = new int [2];
        private String [] seqName = new String [2];

        private DualHeadSequenceBuilder(
            DualHeadStatePathTool path,
            String seq0Name,
            String seq1Name
            )
        {
            this.path = path;
            seqName[0] = seq0Name;
            seqName[1] = seq1Name;

            // get the alphabet
            this.alfa = path.getAlignedSymbols().getAlphabet();
            seqAlfa = alfa.getAlphabets();
            seqSymList = new List [2];
            seqSymList[0] = new ArrayList();
            seqSymList[1] = new ArrayList();

            min[0] = min[1] = max[0] = max[1] = 0;
        }

        /**
         * @param index coordinate in StatePath frame.
         */
        private boolean addStatePathAtIndex(int index)
            throws IllegalSymbolException
        {
            // get the State
//            System.out.println("addStatePathAtIndex: " + index);
//            System.out.println("extending is " + extending);
            State myState = path.getStateAt(index);
            SubModelClass thisClass = getSubModelClass(myState);

            if (!extending) {
//                System.out.println("starting new block " + index);
                if (thisClass != null) {
                    // got a state of interest
                    myClass = thisClass;
                    extending = true;

//                    addSymbol(path.getAlignedSymbols().symbolAt(index));
                    addPathElement(index);
//                    int ungappedIndex = Math.abs(path.getStatePathUngappedCoord(index));
//                    System.out.println("ungappedIndex is " + ungappedIndex);
//                    min[0] = path.getComponentSourceCoord(ungappedIndex, 0);
//                    min[1] = path.getComponentSourceCoord(ungappedIndex, 1);

                    return true;
                }
                else return true;                
            }
            else {
                // in extension phase
                if (thisClass == myClass) {
                    // continue extending
                    addPathElement(index);

                    return true;
                }
                else return false;
            }
        }

        /**
         * @param index position in StatePath coordinate frame.
         */
        private void addPathElement(int index)
            throws IllegalSymbolException
        {
            // recover symbol
            Symbol sym = path.getAlignedSymbols().symbolAt(index);

            // validate alphabet
            if (!alfa.contains(sym)) throw new IllegalSymbolException(sym.getName() + " is not in alphabet " + alfa.getName());

            // the symbol can be a gap (from DotState), a BasisSymbol with a gap in one sequence, 
            // or a BasisSymbol with AtomicSymbols in both sequences.
            if ((sym == AlphabetManager.getGapSymbol()) ||
                (sym == alfa.getGapSymbol())
                ) {

            }
            else if (sym instanceof BasisSymbol) {
                List symList = ((BasisSymbol) sym).getSymbols();
//                System.out.println("processing " + index + " " + sym  );
                Symbol currSym;
                if (((currSym = (Symbol) symList.get(0)) instanceof AtomicSymbol) &&
                    ((Alphabet) seqAlfa.get(0)).contains(currSym)
                    ) {
//                    System.out.println("adding to 0 " + index);
                    max[0] = Math.abs(path.getComponentSourceCoord(index, 0));
                    if (min[0] == 0) min[0] = max[0];
                    hasSequence = true;
                }

                if (((currSym = (Symbol) symList.get(1)) instanceof AtomicSymbol) &&
                    ((Alphabet) seqAlfa.get(1)).contains(currSym)
                    ) {
//                    System.out.println("adding to 1 " + index);
                    max[1] = Math.abs(path.getComponentSourceCoord(index, 1));
                    if (min[1] == 0) min[1] = max[1];
                    hasSequence = true;
                }
//                System.out.println("(" + min[0] +"," + max[0] + ") (" + min[1] + "," + max[1]);
            }
        }
/*
        // gaps have -ve coordinates that are of the magnitude of the first base
        // after the gap.  We need to adjust these down when we want the first base before the gap.
        private int adjustGapCoord(int coord)
        {
            if (coord < 0) {
                return Math.abs(coord) - 1;
            } else
                return coord;
        }
*/
        private Object subModelClass() { return myClass; }

        private boolean hasSequence()
        {
            return hasSequence;
        }

        private void writeSequence(int index, String name, PrintStream ps)
            throws IllegalSymbolException, IOException
        {
            // do we have a sequence at all?
//            System.out.println("writeSequence " + index + " " + min[index] + " " + max[index]);
            if (max[index] > 0) {

                SymbolList symList = path.getSourceSymbolList(index).subList(min[index], max[index]);
                FastaFormat fmt = new FastaFormat();

                Sequence seq = new SimpleSequence(
                    symList,
                    "",
                    name,
                    Annotation.EMPTY_ANNOTATION
                    );

                fmt.writeSequence(seq, ps);
            }
            else {
                ps.println("> " + name);
                ps.println("");
            }
        }

        private void writeXML(XMLWriter xw)
            throws IllegalSymbolException, IOException
        {
            xw.openTag("match_block");
                xw.openTag("location");
                    xw.attribute("seq", seqName[0]);
                    xw.attribute("min", Integer.toString(min[0]));
                    xw.attribute("max", Integer.toString(max[0]));
                xw.closeTag("location");
                xw.openTag("location");
                    xw.attribute("seq", seqName[1]);
                    xw.attribute("min", Integer.toString(min[1]));
                    xw.attribute("max", Integer.toString(max[1]));
                xw.closeTag("location");

//            xw.attribute("min0", Integer.toString(min[0]));
//            xw.attribute("max0", Integer.toString(max[0]));
//            xw.attribute("min1", Integer.toString(min[1]));
//            xw.attribute("max1", Integer.toString(max[1]));
            xw.closeTag("match_block");
        }
    }

    private static class SubModelClass
    {
    }

    private static SubModelClass getSubModelClass(State state)
    {
        Annotation ann;
//        System.out.println("submodelclass ann " + state.getAnnotation());
        if (((ann = state.getAnnotation()) != null) &&
            (ann.containsProperty(SUBMODEL_CLASS))
            ) {
            return (SubModelClass) state.getAnnotation().getProperty(SUBMODEL_CLASS);
        }
        else return null;
    }

    /**
     * set a property within this Annotation to indicate
     * emissions from its associated State extend the region
     * of interest.
     */
    public static void setInterestingState(State state)
        throws IllegalArgumentException, ChangeVetoException
    {
        setInteresting(state.getAnnotation());
    }

    public static void setInteresting(Annotation ann)
        throws IllegalArgumentException, ChangeVetoException
    {
        ann.setProperty(INTERESTING_STATE, OF_INTEREST);
    }

    public static void setConstrainedClass(Annotation ann)
        throws IllegalArgumentException, ChangeVetoException
    {
        ann.setProperty(SUBMODEL_CLASS, CONSTRAINED_SPACER);
    }

    public static void setUnconstrainedClass(Annotation ann)
        throws IllegalArgumentException, ChangeVetoException
    {
        ann.setProperty(SUBMODEL_CLASS, UNCONSTRAINED_SPACER);
    }

    public static void setMatchClass(Annotation ann)
        throws IllegalArgumentException, ChangeVetoException
    {
        ann.setProperty(SUBMODEL_CLASS, MATCH_BLOCK);
    }

	
}


