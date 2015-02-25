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

import java.util.List;

import org.biojava.bio.dp.State;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.GappedSymbolList;
import org.biojava.bio.symbol.Alignment;
//import org.biojava.bio.Annotation;
//import org.biojava.utils.ChangeVetoException;
/**
 * Support class for analysing and printing StatePaths
 * @author David Huen
 */
public class DualHeadStatePathTool
{
    StatePath path;
    GappedSymbolList pairedSymList;
    Alignment alignedComponents;
    GappedSymbolList [] component = new GappedSymbolList[2];
    SymbolList stateLabels;

    public DualHeadStatePathTool(StatePath path)
    {
        this.path = path;

        // recover the statepath aligned sequences
        // This is a single GappedSymbolList with each symbol being
        // a List containing the aligned symbols at that position
        pairedSymList = (GappedSymbolList) path.symbolListForLabel(StatePath.SEQUENCE);

        // the pairedSymList is itself an Alignment object and the
        // component Sequences can be recovered by calling
        // symbolListForLabel to recover the individual GappedSymbolLists
        // that represent the source sequences.
        alignedComponents = (Alignment) pairedSymList.getSourceSymbolList();

        // I want to recover the component sequences
        List seqList = alignedComponents.getLabels();
        // System.out.println(seqList);

        // now recover the base SymbolLists
        component[0] = (GappedSymbolList) alignedComponents.symbolListForLabel(seqList.get(0));
        component[1] = (GappedSymbolList) alignedComponents.symbolListForLabel(seqList.get(1));

        // get state labels
        stateLabels = path.symbolListForLabel(StatePath.STATES);
    }

    public int length()
    {
        return pairedSymList.length();
    }

    public StatePath getStatePath()
    {
        return path;
    }

    public GappedSymbolList getAlignedSymbols()
    {
        return pairedSymList;
    }

    public List getAlignedSymbolsAt(int index)
    {
        return ((BasisSymbol) pairedSymList.symbolAt(index)).getSymbols();
    }

    /**
     * get State at given point in StatePath
     */
    public State getStateAt(int index)
    {
        return (State) stateLabels.symbolAt(index);
    }

    public String getStateLabelAt(int index)
    {
        return stateLabels.symbolAt(index).getName();
    }

    public int getStatePathUngappedCoord(int gapped)
    {
        return pairedSymList.viewToSource(gapped);
    }

    /**
     * @param gapped position in StatePath coordinates.
     */
    public int getComponentSourceCoord(int gapped, int compIndx)
        throws IllegalArgumentException
    {
        if ((compIndx != 0) && (compIndx != 1)) throw new IllegalArgumentException();
        return component[compIndx].viewToSource(getStatePathUngappedCoord(gapped));
    }

    public SymbolList getSourceSymbolList(int compIndx)
    {
        return component[compIndx].getSourceSymbolList();
    }
    
    public SymbolList getAlignmentComponentSymbolList(int compIndx) {
        return component[compIndx];
    }
}
