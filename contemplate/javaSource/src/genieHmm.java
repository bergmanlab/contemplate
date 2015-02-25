import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;


public class genieHmm {

	public static void main(String[] args) throws IllegalArgumentException, BioException, ChangeVetoException {
		 MarkovModel casino = createGenieModel();
    	 DP dp=DPFactory.DEFAULT.createDP(casino);
	     StatePath obs_rolls = dp.generate(300);
	     SymbolList roll_sequence = obs_rolls.symbolListForLabel(StatePath.SEQUENCE);
	     System.out.println(roll_sequence.seqString());
         
         SymbolList[] res_array = {roll_sequence};
         StatePath v = dp.viterbi(res_array, ScoreType.PROBABILITY);
         
         for(int i = 1; i <= v.length(); i++){
             System.out.println(v.symbolAt(StatePath.STATES, i).getName()+" == "+obs_rolls.symbolAt(i));
           }
	}


	public static MarkovModel createGenieModel() throws IllegalSymbolException, ChangeVetoException, IllegalAlphabetException {
		String[] ballColorNames = {"red", "green", "white"};
	    Symbol[] ballSymbols=new Symbol[ballColorNames.length];
	
	    //alphabet setup
	    SimpleAlphabet ballColors=new SimpleAlphabet();
	    ballColors.setName("DiceAlphabet");
	
	    for(int i=0;i<ballColorNames.length;i++) {
	        ballSymbols[i]= AlphabetManager.createSymbol(ballColorNames[i],Annotation.EMPTY_ANNOTATION);
	        ballColors.addSymbol(ballSymbols[i]);
	    }
	    
		// define emission symbols
	    int [] advance = { 1 };
		Distribution urn1Dist;
		Distribution urn2Dist;
	    urn1Dist = DistributionFactory.DEFAULT.createDistribution(ballColors);
	    urn2Dist = DistributionFactory.DEFAULT.createDistribution(ballColors);
		EmissionState urn1 = new SimpleEmissionState("Urn1", Annotation.EMPTY_ANNOTATION, advance, urn1Dist);
		EmissionState urn2 = new SimpleEmissionState("Urn2", Annotation.EMPTY_ANNOTATION, advance, urn2Dist);
		
	    // define model
		SimpleMarkovModel genieModel = new SimpleMarkovModel(1, ballColors, "genieModel");
        genieModel.addState(urn1);
        genieModel.addState(urn2);
        genieModel.createTransition(genieModel.magicalState(),urn1);
        genieModel.createTransition(genieModel.magicalState(),urn2);
        genieModel.createTransition(urn1,genieModel.magicalState());
        genieModel.createTransition(urn2,genieModel.magicalState());
        genieModel.createTransition(urn1,urn2);
        genieModel.createTransition(urn1,urn1);
        genieModel.createTransition(urn2,urn2);
        genieModel.createTransition(urn2,urn1);
        
        double[] urn1EmissionProbs = {0.1, 0.4, 0.5};
        double[] urn2EmissionProbs = {0.6, 0.3, 0.1};
        
        for(int i=0;i<ballColorNames.length;i++)   {
          urn1Dist.setWeight(ballSymbols[i],urn1EmissionProbs[i]);
          urn2Dist.setWeight(ballSymbols[i],urn2EmissionProbs[i]);
        }
     
        //set up transition scores.
        Distribution dist;

        dist = genieModel.getWeights(genieModel.magicalState());
        dist.setWeight(urn1, 0.6);
        dist.setWeight(urn2, 0.4);

        dist = genieModel.getWeights(urn1);
        dist.setWeight(urn1,               0.7);
        dist.setWeight(urn2,               0.6);
        dist.setWeight(genieModel.magicalState(), 0.01);

        dist = genieModel.getWeights(urn2);
        dist.setWeight(urn1,                 0.4);
        dist.setWeight(urn2,               0.6);
        dist.setWeight(genieModel.magicalState(), 0.01);  
        return genieModel;
	
	}
	
	
}

