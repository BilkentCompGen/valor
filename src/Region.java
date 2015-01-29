/*
 * a region is a mapped interval from the input bed file
 */
public class Region extends DNAInterval {
	public int pool;
	public double size;
	/*
	 * constructor
	 */
	public Region (String chromosome, int start, int end, double size, int pool)
	{
		super(chromosome, start, end);
		this.pool = pool;
		this.size = size;
	}
	public Region (String chromosome, int start, int end, int pool)
	{
		super(chromosome, start, end);
		this.pool = pool;
		this.size = size();
	}
	/*
	 * check if this region can be paired with r2
	 * look at the conditions
	 */
	public boolean canBePaired(Region r2)
	{
		Region r1 = this;
			
		return 
			// if they come from the same chromosome
			(r1.chromosome.equals(r2.chromosome)) &&
			// if they come from the same pool
			(r1.pool == r2.pool) &&
			// they should be smaller than normal
			(r1.size <= Config.NORMAL_SIZE) &&
			(r2.size <= Config.NORMAL_SIZE) &&
			// if the size is OK (|A|+|B| is normal)
			(r1.size + r2.size >= Config.LOW_CRITERIA) &&
			(r1.size + r2.size <= Config.UP_CRITERIA) &&
			// and the distance is between the desired inversion size
			(r2.start - r1.end >= Config.MIN_INVERSION_SIZE) &&
			(r2.start - r1.end <= Config.MAX_INVERSION_SIZE);
	}
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString ()
	{
		return chromosome + "\t" + start + " "  + end + " " + size + " "+ pool;
	}
}
