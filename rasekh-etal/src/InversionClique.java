import java.util.ArrayList;


public class InversionClique {

	ArrayList <InversionNode> nodes;
	int EPrime;
	int VPrime;
	final double LAMBDA;
	public InversionClique (double lambda, InversionNode n) {
		LAMBDA = lambda;
		nodes = new ArrayList<InversionNode>();
		nodes.add(n);
		VPrime = 1;
		EPrime = 0;
	}
	public int addNode(InversionNode n)
	{
		nodes.add(n);
		VPrime++;
		for (InversionNode v : n.adj)
		{
			if (nodes.contains(v))
			{
				EPrime++;
			}
		}
		return EPrime;
	}
	public int remNode(InversionNode n)
	{
		nodes.remove(n);
		VPrime--;
		for (InversionNode v : n.adj)
		{
			if (nodes.contains(v))
			{
				EPrime--;
			}
		}
		return EPrime;
	}
	/*
	public int DV()
	{
		return (int) Math.ceil(LAMBDA*((VPrime+1)*VPrime) / 2.0) - EPrime;
	}
	public int EV()
	{
		return EPrime - (int) Math.ceil(LAMBDA*((VPrime-1)*(VPrime-2)) / 2.0) - EPrime;
	}
	public ArrayList<Node> Crit()
	{
		ArrayList<Node> crit = new ArrayList<Node>();
		for (Node v : nodes)
		{
			if (v.dv < Math.ceil(LAMBDA * VPrime))
			{
				crit.add(v);
			}
		}
		return crit;
	}
	public ArrayList<Node> RCrit()
	{
		ArrayList<Node> rcrit = new ArrayList<Node>();
		for (Node v : nodes)
		{
			if (v.dv - 1 < Math.ceil(LAMBDA * (VPrime-2)))
			{
				rcrit.add(v);
			}
		}
		return rcrit;
	}
	public ArrayList<Node> AddV(ArrayList<Node> allNodes)
	{
		ArrayList<Node> addv = new ArrayList<Node>();
		ArrayList<Node> crit = Crit();
		int dv = DV();
		for (Node v : allNodes)
		{
			boolean check = true;
			if (v.dv >= Math.max(LAMBDA*VPrime, dv))
			{
				for (Node e : crit)
				{
					if (!v.adj.contains(e))
					{
						check = false; 
						break;
					}
				}
				if (check)
				{
					addv.add(v);
				}
			}
		}
		return addv;
	}
	public ArrayList<Node> RemV()
	{
		ArrayList<Node> remv = new ArrayList<Node>();
		ArrayList<Node> rcrit = RCrit();
		int ev = EV();
		for (Node n : nodes)
		{
			if (n.dv <= ev)
			{
				boolean check = true;
				for (Node e : rcrit)
				{
					if (n.adj.contains(e))
					{
						check = false; 
						break;
					}
				}
				if (check)
				{
					remv.add(n);
				}
			}
		}
		return remv;
	}
	public boolean checkDV()
	{
		int sum = 0;
		for (Node n : nodes)
		{
			sum += n.dv;
		}
		return sum == EPrime*2;
	}
	public ArrayList<Node> PAdd(ArrayList<Node> V)
	{
		ArrayList<Node> padd = new ArrayList<Node>();
		for (Node v : V)
		{
			if (!nodes.contains(v) && v.dv > LAMBDA * (this.VPrime-1))
			{
				padd.add(v);
			}
		}
		return padd;
	}
	public ArrayList<Node> PCrit(Node w)
	{
		ArrayList<Node> pcrit = new ArrayList<Node>();
		for (Node v : nodes)
		{
			if (v.dv - 1 < LAMBDA * (this.VPrime-1) && !w.adj.contains(v))
			{
				pcrit.add(v);
			}
		}
		return pcrit;
	}
	public double rVPrimeW (Node w)
	{
		return EPrime + w.dv - LAMBDA * VPrime * (VPrime-1) / 2.0;
	}
	public ArrayList<Node> PRem(Node w)
	{
		ArrayList<Node> prem = new ArrayList<Node>();
		ArrayList<Node> pcrit = PCrit(w);
		boolean check;
		
		double rvprimew = rVPrimeW(w);
		for (Node v : nodes)
		{
			int dv = v.dv + (w.adj.contains(v) ? 1 : 0);
			if (dv <= rvprimew)
			{
				check = true;
				for (Node n : pcrit)
				{
					if (v.adj.contains(n))
					{
						check = false;
						break;
					}
				}
				check &= !v.adj.contains(w);
				if (check) 
				{
					prem.add(v);
				}
			}
		}
		return prem;
	}*/
}
