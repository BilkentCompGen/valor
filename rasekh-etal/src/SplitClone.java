
public class SplitClone extends PairedDNAInterval {
	public String pool;
	
	public SplitClone (Region r1, Region r2) throws Exception
	{
		super (r1, r2);
		if (!r1.pool.equals(r2.pool))
		{
			throw new Exception("r1 and r2 are not from same pool @ PairedRegion");
		}
		pool = r1.pool;
	}
	/*
	 *  this.left-----------that.left                    this.right-----------that.right
	 *           <---GAP--->                                       <---GAP--->
	 */
	public boolean canBeSplitClone (SplitClone pairedInterval)
	{
		PairedDNAInterval that = pairedInterval;
		return (that.chromosome.equals(this.chromosome) &&		// on same chromosome
			this.left.start < that.left.start &&				// in order
			this.left.end < that.left.end &&				
			this.right.start < that.right.start &&
			this.right.end < that.right.end &&
			that.left.distance(this.left) > Config.INV_OVERLAP &&	// > overlap to ensure its not a deletion
			that.right.distance(this.right) > Config.INV_OVERLAP &&
			that.left.distance(this.left) < Config.INV_GAP &&		// not too far away
			that.right.distance(this.right) < Config.INV_GAP);
	}
	@Override
	public String toString()
	{
		return  left + "\t" +
				right + "\t" +
				pool;
	}

}
