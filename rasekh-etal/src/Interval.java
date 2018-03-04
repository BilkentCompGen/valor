
public class Interval {
	public Interval (int start, int end)
	{
		this.start = start;
		this.end = end;
	}
	public int size()
	{
		return end - start + 1;
	}
	public int start;
	public int end;
	public String toString()
	{
		return start + "\t" + end;
	}

}
