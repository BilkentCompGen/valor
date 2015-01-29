import java.util.ArrayList;

/*
 * 	A-----C				         B-----D
 *    gap	<--inversion size-->   gap
 */
public class Inversion extends PairedDNAInterval{

	public SplitClone AB;
	public SplitClone CD;
	public SplitClone breakpoint;
	public int plusPlusSupport;
	public int minusMinusSupport;
	String chromosome;
	boolean covered;
	/*
	 * constructor, support is initially zero, calculated later
	 */
	public Inversion(String chromosome, SplitClone AB, SplitClone CD) throws Exception
	{
		super ( new DNAInterval(chromosome, AB.left.end, CD.left.start),
				new DNAInterval(chromosome, AB.right.end, CD.right.start));
		this.AB = AB;
		this.CD = CD;
		this.plusPlusSupport = -1;
		this.minusMinusSupport = -1;
		this.chromosome = chromosome;
		covered = false;
		
	}
	public boolean overlaps(Inversion that) throws Exception
	{
		PairedDNAInterval bp1 = new PairedDNAInterval(
				new DNAInterval(this.chromosome, Math.min(this.AB.left.end, this.CD.left.start), Math.max(this.AB.left.end, this.CD.left.start)), 
				new DNAInterval(this.chromosome, Math.min(this.AB.right.end, this.CD.right.start), Math.max(this.AB.right.end, this.CD.right.start)));
		
		PairedDNAInterval bp2 = new PairedDNAInterval(
				new DNAInterval(that.chromosome, Math.min(that.AB.left.end, that.CD.left.start), Math.max(that.AB.left.end, that.CD.left.start)), 
				new DNAInterval(that.chromosome, Math.min(that.AB.right.end, that.CD.right.start), Math.max(that.AB.right.end, that.CD.right.start)));
		return bp1.overlaps(bp2);
	}
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString()
	{
		return  chromosome + "\t" +
				AB.left.end + "\t" +
				CD.left.start + "\t" +
				AB.right.end + "\t" +
				CD.right.start + "\t" +
				plusPlusSupport + "\t" +
				minusMinusSupport;
	}
	/*
	 * 
	 */
	public void setReadSupport( ArrayList<PairedDNAInterval> plusPlusReads,
									ArrayList<PairedDNAInterval> minusMinusReads)/*,
									ArrayList<DNAInterval> plusMinusReads)*/
	{
		plusPlusSupport = this.AB.setReadSupport(plusPlusReads);
		minusMinusSupport = this.CD.setReadSupport(minusMinusReads);
	}
}
