public class DNAInterval {
	public String chromosome;
	public int start;
	public int end;
	/*
	 * CONSTRUCTOR
	 */
	public DNAInterval(String chrom, int start, int end)
	{
		this.chromosome = chrom;
		this.start = Math.min(start, end);
		this.end = Math.max(start, end);
	}
	/*
	 * 
	 */
	public double size()
	{
		return (double) (end - start);
	}
	/*
	 * 
	 */
	public boolean overlaps (DNAInterval interval)
	{
		DNAInterval that = interval;
		return (this.chromosome.equals(that.chromosome) &&
				((that.start >= this.start - Config.READ_SIZE && that.start <= this.end + Config.READ_SIZE) || 
				(that.end >= this.start - Config.READ_SIZE && that.end <= this.end + Config.READ_SIZE)));
	}
	public boolean overlaps (DNAInterval interval, final int limit)
	{
		DNAInterval that = interval;
		return (this.chromosome.equals(that.chromosome) &&
				((that.start >= this.start - limit - Config.READ_SIZE && that.start <= this.end + limit + Config.READ_SIZE) || 
				(that.end >= this.start - limit - Config.READ_SIZE && that.end <= this.end + limit + Config.READ_SIZE)));
	}

	/*
	 * 
	 */
	public int distance (DNAInterval interval)
	{
		DNAInterval that = interval;
		if (this.start < that.start)
		{
			return that.distance (this);
		}
		return this.start - that.end;
	}
	/*
	 * 
	 */
	public DNAInterval intersect (DNAInterval interval)
	{
		DNAInterval that = interval;
		if (!this.overlaps(that))
		{
			return null;
		}
		return new DNAInterval(chromosome, 
				Math.max(this.start, that.start), 
				Math.min(this.end, that.end));
	}
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString()
	{
		return chromosome + "\t" + start + "\t" + end;
	}
}
