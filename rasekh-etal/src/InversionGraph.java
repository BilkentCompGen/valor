import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;


public class InversionGraph {
	public ArrayList<InversionNode> nodes;
	public final double LAMBDA;
	public final double GAMMA;
	public final int TABU;
	
	public InversionGraph (double lambda , double gamma, int tabu)
	{
		TABU = tabu;
		GAMMA = gamma;
		LAMBDA = lambda;
		nodes = new ArrayList<InversionNode>();
	}
	public int size()
	{
		return nodes.size();
	}
	public void initializeGraph(ArrayList<Inversion> splitClones) throws Exception
	{
		int edgeSize = 0;
		int nodeSize = splitClones.size();
		
		for (int i = 0; i < nodeSize; i++)
		{
			nodes.add(new InversionNode(splitClones.get(i)));
		}
		for (int i = 0; i < nodeSize; i++)
		{
			for (int j = i+1; j < nodeSize; j++)
			{
				if (splitClones.get(i).overlaps(splitClones.get(j)))
				{
					nodes.get(i).adj.add(nodes.get(j));
					nodes.get(j).adj.add(nodes.get(i));
					edgeSize++;
				}
			}
		}
		System.out.println("SplitCloneGraph made...");
		System.out.println("\t\tNodes:\t" + nodeSize + "\tEdges: " + edgeSize);
	}
	public static void sortByDegree(ArrayList<InversionNode> nodes)
	{
		// sort the splirClones according to support in decreasing order
		Collections.sort(nodes, new Comparator<InversionNode>() {
			@Override public int compare(InversionNode n1, InversionNode n2) {
			    return (int) (n1.degree() > n2.degree() ? 1 : 0);
			}
		});
	}
	public InversionClique findMaxBrokenClique (int seed, boolean sorted) {
		InversionClique foundClique;
		if (seed < 0 || seed >= nodes.size())
		{
			System.out.println("Seed out of boundary...");
		}
		// sort nodes if desired, may want to keep it in order of chrom/position
		if (sorted)
		{
			sortByDegree(nodes);
		}
		// reset everything to initial state
		for (InversionNode n : nodes)
		{
			n.dv = 0;
			n.tabu = 0;
		}
		// make a clique with initial seed
		// NOTICE: at the end clique might not contain this seed
		foundClique = new InversionClique(LAMBDA, nodes.get(seed));
		// mark all the adj of the seed to have dv = 1
		for (InversionNode n : foundClique.nodes.get(0).adj)
		{
			n.dv = 1;
		}
		// remove that node from the nodes
		nodes.remove(0);
		while(PlateauMove(foundClique));
		return foundClique;
	}
	public boolean removeNode(int i)
	{
		if(i >= nodes.size() || i < 0)
			return false;
		nodes.remove(i);
		return true;
	}
	public boolean PlateauMove (InversionClique Q)
	{		
		double max;
		boolean added, removed;
		double score;
		InversionNode best = null;
		//add
		max = -1;
		for (InversionNode n : nodes)
		{
			if (n.tabu > 0)
			{
				n.tabu --;
			}
			if (n.tabu == 0 && n.dv > GAMMA * (Q.VPrime - 1) && Q.EPrime + n.dv >= (LAMBDA * (Q.VPrime * (Q.VPrime-1)))/2.0)
			{
				score = (2.0 * (Q.EPrime + n.dv)) / (Q.VPrime * (Q.VPrime + 1));
				if (max < score)
				{
					best = n;
					max = score;
				}
			}
		}
		if (max == -1)
		{
			added = false;
		}
		else 
		{
			Q.addNode(best);
			nodes.remove(best);
			added = true;
			// update dv
			for (InversionNode n : best.adj)
			{
				n.dv++;
			}
			best.tabu = TABU;
		}
		//remove
		max = -1;
		for (InversionNode n : Q.nodes)
		{
			if (n.tabu > 0)
			{
				n.tabu --;
			}
			if (n.tabu == 0 && n.dv < GAMMA * (Q.VPrime - 1) && Q.EPrime - n.dv > (LAMBDA * (Q.VPrime - 2) * (Q.VPrime-1))/2.0)
			{
				score = (2.0 * (Q.EPrime - n.dv)) / ((Q.VPrime - 1) * (Q.VPrime -2));
				if (max < score)
				{
					best = n;
					max = score;
				}
			}
		}
		if (max == -1)
		{
			removed = false;
		}
		else 
		{
			Q.remNode(best);
			nodes.add(best);
			removed = true;
			// update dv
			for (InversionNode n : best.adj)
			{
				n.dv--;
			}
			best.tabu = TABU;
		}
		return (added | removed);
	}
	
}
