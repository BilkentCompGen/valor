import java.util.ArrayList;


public class InversionNode {

	Inversion self;
	boolean dead;
	ArrayList<InversionNode> adj;
	int dv;
	int tabu;
	public InversionNode(Inversion self) {
		this.self = self;
		adj = new ArrayList<InversionNode>();
		dv = 0;
		tabu = 0;
		dead = true;
	}
	public int degree ()
	{
		return adj.size();
	}
	public void calcDV(ArrayList<InversionNode> clique)
	{
		for (InversionNode e : clique)
		{
			if (adj.contains(e))
			{
				dv++;
			}
		}
	}
	public void kill()
	{
		dead = true;
	}

}
