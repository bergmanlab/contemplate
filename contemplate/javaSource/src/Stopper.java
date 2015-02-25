import org.biojava.bio.dp.TrainingAlgorithm;
import org.biojava.bio.dp.StoppingCriteria;

public class Stopper
    implements StoppingCriteria
{
    private double tol;
    private int maxCycles;

    public Stopper(double tol, int maxCycles) { this.tol = tol; this.maxCycles = maxCycles; }

    public boolean isTrainingComplete(TrainingAlgorithm ta)
    {
        if (ta.getCycle() <= maxCycles) {
//            System.out.println("score is " + ta.getCurrentScore());
            return (Math.abs(ta.getCurrentScore() - ta.getLastScore()) < tol);
        }
        else return true;
    }
}

