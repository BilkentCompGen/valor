import java.util.ArrayList;

public class InversionCluster {
	public final int CLIQUE_SIZE;
	public String chromosome;
	public ArrayList<Inversion> splitClones;
	public PairedDNAInterval breakPoint;
	public int minusMinusSupport;
	public int plusPlusSupport;
	
	/*
	 * CONSTRUCTOR
	 */
	public InversionCluster (Inversion splitClone, int cliqueSize) 
	{
		this.chromosome = splitClone.chromosome;
		this.CLIQUE_SIZE = cliqueSize;
		this.splitClones = new ArrayList<Inversion>();
		this.splitClones.add(splitClone);
		breakPoint = splitClone;
		
		this.minusMinusSupport = -1;
		this.plusPlusSupport = -1;
	}

	/*
	 * 
	 */
	public boolean addSplitClone(Inversion splitClone)
	{
		if (!this.isCompatible(splitClone) ||
			 this.splitClones.contains(splitClone))
		{
			return false;
		}
		this.splitClones.add(splitClone);
		if (breakPoint == null)
		{
			breakPoint = splitClone;
		}
		else
		{
			breakPoint.intersect(splitClone);
		}		
		
		return true;
	}
	/*
	 * 
	 */
	public boolean isCompatible (Inversion splitClone)
	{
		return (breakPoint == null) && (this.breakPoint.overlaps(splitClone));
	}
	public boolean isCompatible (InversionNode node)
	{
		return (breakPoint != null) && (this.breakPoint.overlaps(node.self));
	}
	/*
	 * 
	 */
	public boolean isCompatible (InversionCluster cluster)
	{
		return this.breakPoint.overlaps (cluster.breakPoint);
	}
	
	/*
	 * NOTE: this is before cluster
	 */
	public boolean merge (InversionCluster cluster)
	{
		InversionCluster that;
		that = cluster;
		if (!this.isCompatible(that)) {
			return false;
		}
		// can we use heuristic 1 ?
		for (Inversion sc : that.splitClones)
		{
			this.addSplitClone(sc);
		}
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString()
	{
        return(this.breakPoint + "\t"+
               this.CLIQUE_SIZE + "\t" +
			   this.plusPlusSupport + "\t" +
			   this.minusMinusSupport + "\t");


	}

	/*
	 * 
	 */
	public int optimizedSize()
	{
		int sum = 0;
		for (Inversion sc : splitClones)
		{
			if (!sc.covered)
			{
				sum++;
			}
		}
		return sum;
	}
	/*
	 * 
	 */
	public void cover()
	{
		for (Inversion sc : splitClones)
		{
			sc.covered = true;
		}
	}
	/*
	 * 
	 */
	/*public boolean setDelCount (ArrayList<DNAInterval> deletions)
	{
		DNAInterval inversion = new DNAInterval(this.chromosome, this.left.end, this.right.start);
		this.delCount = 0;
		for (DNAInterval del : deletions)
		{
			if (inversion.overlaps(del))
			{
				this.delCount++;
			}
		}
		return (this.delCount > 0 ? true : false);
	}*/
	/*
	 * 
	 */
	public boolean canPartition(Inversion splitClone) throws Exception
	{
		for (Inversion potential : this.splitClones)
		{
			if (potential.overlaps(splitClone))
			{
				return true;
			}
		}
		return false;
	}
	public void setReadSupport( ArrayList<PairedDNAInterval> plusPlus,
			ArrayList<PairedDNAInterval> minusMinus,
			int limit) throws Exception
	{
		PairedDNAInterval bp = new PairedDNAInterval(
				new DNAInterval(chromosome, Math.min(breakPoint.left.start, breakPoint.left.end) - limit, Math.max(breakPoint.left.start, breakPoint.left.end) + limit),
				new DNAInterval(chromosome, Math.min(breakPoint.right.start, breakPoint.right.end) - limit, Math.max(breakPoint.right.start, breakPoint.right.end) + limit));
		bp.setReadSupport(plusPlus);
		plusPlusSupport = bp.setReadSupport(plusPlus);
		minusMinusSupport = bp.setReadSupport(minusMinus);
	}
}
