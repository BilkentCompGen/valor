import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;
/*
 * class responsible for all IO operations, read file ect.
 */
public class IOTools {
	
	
	/*
	 * fixes the file to contain the set because the data is provided in such a format
	 * each one of the three sets has 8 letters (A to H) and each letter contains 12 digits
	 * each one is a pool, so for example set 1, char A, digit 1 is pool 1 and 
	 * set 3, char H, digit 12 is pool 288. Because of the probe sets, 12x8
	 * we have to fix this, the data is given in this format : 
	 * (given by the C code that we run for reading the BAM files)
	 *  cnt	set	char	digit	size	chromosome	startLoci	endLoci
	 */
	public static boolean fix (String regionFilename)
	{
		Scanner fileScan;
		PrintWriter fileWrite;
		String filename = regionFilename;
		String [] data;
		// new file name will have a -fixed at the end of it
		String newFilename = filename.substring(0, filename.lastIndexOf(".")) 
				+ "-fixed" 
				+ filename.substring(filename.lastIndexOf("."), filename.length());
		// try to open original file
		try
		{
			fileScan = new Scanner(new File(filename));
		}
		catch (Exception e)
		{
			System.out.println("ERROR OPENNING ORIGINAL FILE");
			System.out.println(e);
			return false;
		}
		// try to create the fixed version of the file
		try
		{
			fileWrite = new PrintWriter(new File(newFilename));
		}
		catch (Exception e)
		{
			System.out.println("ERROR CREATING FIXED FILE");
			System.out.println(e);
			fileScan.close();
			return false;
		}
		// read data line by line and calculate the pool, accordingly
		try {
			fileWrite.println(fileScan.nextLine().replace("\"", "") + "\tpool" );
			while (fileScan.hasNextLine())
			{
				String line = fileScan.nextLine().replace("\"", "");
				data = line.split("\t");
				// cnt	set	char	digit	size	chromosome	startLoci	endLoci
				// explain: read the method header
				int pool = (Integer.parseInt(data[1])-1) * 96 +  ((int) (data[2].charAt(0) - 'A') * 12) + Integer.parseInt(data[3]);
				fileWrite.println(line + "\t" + pool );
			}
			// close up so no strange errors occur
			fileScan.close();
			fileWrite.close();
			filename = newFilename;
			return true;
		}
		catch (Exception e)
		{
			fileScan.close();
			fileWrite.close();
			System.out.println("ERROR READING/WRITING FILE");
			System.out.println(e);
			return false;
		}
	}
	/*
	 * Read the filename file, assuming that it is fixed and contains the pool data
	 */
	public static ArrayList<Region> readRegions (String regionFilename, String chromosome)
	{
System.out.println("FILENAME:  " + regionFilename);
		ArrayList<Region> regionList = new ArrayList<Region>();
		Scanner fileScan;
		String filename = regionFilename;
		String [] line;
		// try to to open the file
		try {
			fileScan = new Scanner(new File(filename));
		}
		catch (Exception e)
		{
			System.out.println("ERROR OPENNING FIXED FILE");
			System.out.println(e);
			return null;
		}
		// read line by line and add the region to regionList
		try {
			fileScan.nextLine().replace("\"", "");
			while (fileScan.hasNextLine())
			{
				line = fileScan.nextLine().replace("\"", "").split("\t");
				// cnt	set	char	digit	size	chromosome	startLoci	endLoci pool
				if (line[0].equals(chromosome))
				{
					// String chromosome, long start, long end, String pool
					regionList.add(new Region(line[0],
							Integer.parseInt(line[1]), 
							Integer.parseInt(line[2]), 
							line[3]));
				}
			}
			fileScan.close();
			/*Collections.sort(regionList, new Comparator<Region>() {
				@Override public int compare(Region r1, Region r2) {
				    return (int) (r1.start - r2.start);
				}
			
			});*/
			return regionList;
		}
		catch (Exception e)
		{
			System.out.println("ERROR READING FIXED FILE");
			e.printStackTrace();
			fileScan.close();
			return null;
		}
	}
	/*
	 *  read the reads
	 */
	public static ArrayList<PairedDNAInterval> readBedFile (String bedFilename, String chromosome)
	{
		ArrayList<PairedDNAInterval> pairedIntervalList;
		Scanner fileScan;
		String [] line;
		int start1, start2, end1, end2;
		
		// try to to open the plus file
		try {
			fileScan = new Scanner(new File(bedFilename));
		}
		catch (Exception e)
		{
			System.out.println("ERROR OPENNING BED FILE ");
			e.printStackTrace();
			return null;
		}
		// set the readList
		pairedIntervalList = new ArrayList<PairedDNAInterval>();
		// read line by line and add the read to readList
		try {
			while (fileScan.hasNextLine())
			{
				line = fileScan.nextLine().replace("\"", "").split("\t");
				if (chromosome.equals("all") || line[0].equals(chromosome))
				{
					// String chromosome, long start, long end
					// add it to the readList
					start1 = Math.min(Integer.parseInt(line[1]), Integer.parseInt(line[2]));
					end1 = Math.max(Integer.parseInt(line[1]), Integer.parseInt(line[2]));
					start2 = Math.min(Integer.parseInt(line[3]), Integer.parseInt(line[4]));
					end2 = Math.min(Integer.parseInt(line[3]), Integer.parseInt(line[4]));
					if(start1 < start2) {
						pairedIntervalList.add(new PairedDNAInterval(new DNAInterval(chromosome, start1, end1), 
																	 new DNAInterval(chromosome, start2, end2)));
					}
				}
			}
			fileScan.close();
			/*Collections.sort(pairedIntervalList, new Comparator<PairedDNAInterval>() {
				@Override public int compare(PairedDNAInterval r1, PairedDNAInterval r2) {
					if (r1.left.start - r2.left.start == 0)
					{
						return (r1.right.start - r2.right.start);
					}
				    return (r1.left.start - r2.left.start);
				}
			
			});*/
			return pairedIntervalList;
		}
		catch (Exception e)
		{
			System.out.println("ERROR READING BED FILE ");
			e.printStackTrace();
			fileScan.close();
			return null;
		}
		
	}


	/*
	 * write inversion clusters
	 */
	public static boolean writeClusters(ArrayList<ArrayList<InversionCluster>> clusters)
	{
		PrintWriter chrom;
		// temp variables
		String chr, filename;
		boolean status;
		System.out.println("*************************************************");
		System.out.println("Creating Cluster files...");

		filename = "valor_inversion_calls.xls";
		try {
			chrom = new PrintWriter(new File(filename));
		} 
		catch (Exception e)
		{
			System.out.println("ERROR: COULD NOT OPEN OUTPUT FILE [" + filename + "]. @InputOutput.writeClustersByChrom");
			System.out.println(e.getMessage());
			return false;
		}
		chrom.println("chrom\tleft_start\tleft_end\tright_start\tright_end\tClique_size\tsupport++\tsupport--");
		try {
			for (int i = 0; i < clusters.size(); i++)
			{
				for (InversionCluster clu : clusters.get(i))
				{
					chrom.println(clu);
				}
			}
			status = true;
		}
		catch(Exception e)
		{
			System.out.println("ERROR: COULD NOT WRITE CLUSTER FILES @InputOutput.writeClustersByChrom\n");
			System.out.println(e.getMessage());
			status = false;
		}
		chrom.close();
		return status;
	}
	
	
}
