import java.util.ArrayList;


public class PairedDNAInterval {
	DNAInterval left;
	DNAInterval right;
	String chromosome;
	
	public PairedDNAInterval(DNAInterval left, DNAInterval right) throws Exception {
		if (!left.chromosome.equals(right.chromosome))
			throw new Exception("Intervals are not on same chromosome @ PairedInterval");
		this.left = left;
		this.right = right;
		chromosome = left.chromosome;
	}
	/*
	 * 
	 */
	public boolean overlaps (PairedDNAInterval pairedInterval)
	{
		PairedDNAInterval that = pairedInterval;
		return (this.left.overlaps(that.left) && this.right.overlaps(that.right));
	}
	public boolean overlaps (PairedDNAInterval pairedInterval, final int limit)
	{
		PairedDNAInterval that = pairedInterval;
		return (this.left.overlaps(that.left, limit) && this.right.overlaps(that.right, limit));
	}
	public boolean isAlternate(PairedDNAInterval pairedInterval)
	{
		PairedDNAInterval that = pairedInterval;
		if (this.left.overlaps(that.left) || this.right.overlaps(that.right))
		{
			return true;
		}
		return false;
	}
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString()
	{
		return this.chromosome + "\t" + this.left.start + "\t" + this.left.end + "\t" + this.right.start + "\t" + this.right.end;
	}
	/*
	 * 
	 */
	public int overlaps (ArrayList<PairedDNAInterval> reads)
	{
		int sum = 0;
		int beg, last, middle = 0;
		// binary search the beginning assuming intervals is sorted
		beg = 0;
		last = reads.size() - 1;
		while (beg < last)
		{
			middle = (last + beg) / 2;
			if (reads.get(middle).left.end < this.left.start) // middle is before this
			{
				beg = middle + 1;
			}
			else //if (intervals.get(middle).start > this.end) // middle is after this
			{
				last = middle - 1;
			}
		}
		// from this beginning start checking, but end when there cannot be any more overlaps
		for (int i = middle; i < reads.size(); i++)
		{
			if (this.overlaps(reads.get(i)))
			{
				sum++;
			}
			if (this.left.end < reads.get(i).left.start)
			{
				// intervals passed each other
				break;
			}
		}
		return sum;
	}
	/*
	 * 
	 */
	public void intersect (PairedDNAInterval that)
	{
		this.left = that.left.intersect(this.left);
		this.right = that.right.intersect(this.right);
	}
	/*
	 * given the support reads (++ and --) calculate the support of the cluster
	 * First calculate the union of all split clones for this cluster
	 * Then search the read array lists for overlaps
	 */
	public int setReadSupport( ArrayList<PairedDNAInterval> reads)
	{
		return this.overlaps(reads);
	}
	
}
