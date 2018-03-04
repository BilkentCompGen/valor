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
					// String chromosome, long start, long end, int pool
					regionList.add(new Region(line[0],
							Integer.parseInt(line[1]), 
							Integer.parseInt(line[2]), 
							Integer.parseInt(line[3])));
				}
			}
			fileScan.close();
			Collections.sort(regionList, new Comparator<Region>() {
				@Override public int compare(Region r1, Region r2) {
				    return (int) (r1.start - r2.start);
				}
			
			});
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
			Collections.sort(pairedIntervalList, new Comparator<PairedDNAInterval>() {
				@Override public int compare(PairedDNAInterval r1, PairedDNAInterval r2) {
					if (r1.left.start - r2.left.start == 0)
					{
						return (r1.right.start - r2.right.start);
					}
				    return (r1.left.start - r2.left.start);
				}
			
			});
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
	 * write bed files for pools
	 */
	public static boolean writePoolBeds(ArrayList<Inversion> inversions, String poolOutputDir)
	{
		if (!poolOutputDir.equals(""))
		{
			// printers
			PrintWriter [] bed;
			// temp variables
			String filename;
			int pool;
			final int AB, CD, POOL_COUNT;
			boolean status;
			// write the bed files for each pool
			System.out.println("*************************************************");
			System.out.println("Creating BED files by pool...");
			// a bed for each pool (total 288 pools) and 2 extra for all AB ones and all CD ones
			bed = new PrintWriter[290];
			// ab is index for PoolAB and cd is index for PoolCD
			POOL_COUNT = 288;
			AB = POOL_COUNT;
			CD = POOL_COUNT + 1;
			for (int i = 0; i < POOL_COUNT; i++)
			{
				filename = poolOutputDir + "/pool" + Integer.toString(i+1) + ".bed";
				try{
					bed[i] = new PrintWriter(new File(filename));
				} 
				catch (Exception e)
				{
					System.out.println("ERROR: COULD NOT MAKE POOL BED FILE [" + filename + "] @InputOutput.writePoolBeds");
					return false;
				}
			}
			filename = poolOutputDir + "/allPoolsAB.bed";
			try{
				bed[AB] = new PrintWriter(new File(filename));
			} 
			catch (Exception e)
			{
				System.out.println("ERROR: COULD NOT MAKE POOL BED FILE [" + filename + "] @InputOutput.writePoolBeds");
				return false;
			}
			filename = poolOutputDir + "/allPoolsCD.bed";
			try{
				bed[CD] = new PrintWriter(new File(filename));
			} 
			catch (Exception e)
			{
				System.out.println("ERROR: COULD NOT MAKE POOL BED FILE [" + filename + "] @InputOutput.writePoolBeds");
				return false;
			}
			// try to write them
			try {
				for (Inversion inv : inversions)
				{
					// each region (A B C or D) should be written to their own pool
					// and also A and B to "AB" and B and D to "CD"
					pool = inv.AB.pool;
					bed[pool-1].println(inv.AB);
					bed[AB].println(inv.AB);
					pool = inv.CD.pool;
					bed[pool-1].println(inv.CD);
					bed[CD].println(inv.CD);
				}
				status = true;
			}
			catch(Exception e)
			{
				System.out.println("ERROR: COULD NOT WRITE POOL BED FILES @InputOutput.writePoolBeds");
				status = false;
			}
			for (int i = 0; i < POOL_COUNT; i++)
			{
				bed[i].close();
			}
			bed[AB].close();
			bed[CD].close();
			return status;
		} // end of printing pool bed files
		System.out.println("POOL BED files not written - Not required by call method @ InputOutput.writePoolBeds");
		return false;
	}
	/*
	 * write split clones for chromosomes
	 */
	/*public boolean writeSplitClones(ArrayList<SplitClone> splitClones, String splitCloneOutputDir)
	{
		if (!splitCloneOutputDir.equals(""))
		{
			PrintWriter [] chrom;
			// temp variables
			String chr;
			int index, x, y, all;
			boolean status;
			System.out.println("*************************************************");
			System.out.println("Creating Split Clone files...");
			// a file for each chromosome (1-22 + X) and one for all chromosomes
			chrom = new PrintWriter[24];
			// x is for chromosome X and all is for all chromosomes
			x = 22;
			y = 23;
			all = 24;
			try {
				for (int i = 0; i < 22; i++)
				{
					chrom[i] = new PrintWriter(new File(splitCloneOutputDir + "/chr" + (i+1)) + ".sc");
				}
				chrom[x] = new PrintWriter(new File(splitCloneOutputDir + "/chrX.sc"));
				chrom[y] = new PrintWriter(new File(splitCloneOutputDir + "/chrY.sc"));
				chrom[all] = new PrintWriter(new File(splitCloneOutputDir + "/allChrom.sc"));
			} catch (Exception e)
			{
				System.out.println("ERROR: COULD NOT MAKE BED FILES.");
				return false;
			}
			try {
				for (SplitClone sc : splitClones)
				{
					// get the index according to the chromosome
					chr = sc.AB.left.chromosome;
					if (chr.indexOf('X') > 0)
					{
						index = x;
					}
					else
					{
						// don't forget to -1
						index = Integer.parseInt(chr.substring(3)) - 1;
					}
					chrom[index].println(sc);
					chrom[all].println(sc);
				}
				status = true;
			}
			catch(Exception e)
			{
				System.out.println("ERROR: COULD NOT WRITE BED FILES.");
				status = false;
			}
			// close them so thing don't get shitty
			for (int i = 0; i < 24; i++)
			{
				chrom[i].close();
			}
			return status;
		} // end of printing split clones
		System.out.println("Split Clone files not written - Not required by call method @ InputOutput.writeSplitClones");
		return false;
	}
	*/
	/*
	 * write split clones for chromosomes
	 */
	public static boolean writeClustersByChrom(ArrayList<ArrayList<InversionCluster>> clusters, String clusterOutputDir)
	{
		if (!clusterOutputDir.equals(""))
		{
			PrintWriter [] chrom;
			// temp variables
			String chr, filename;
			final int X, Y, ALL;
			boolean status;
			System.out.println("*************************************************");
			System.out.println("Creating Cluster files...");
			// a file for each chromosome (1-22 + X) and one for all chromosomes
			chrom = new PrintWriter[clusters.size() + 1];
			X = 22;
			Y = 23;
			ALL = 24;
			for (int i = 0; i < clusters.size(); i++)
			{
				if (i == X)
				{
					chr = "chrX";
				}
				else if (i == Y)
				{
					chr = "chrY";
				}
				else
				{
					chr = "chr" + Integer.toString(i+1);
				}
				
				filename = clusterOutputDir + chr + ".xls";
				try {
					chrom[i] = new PrintWriter(new File(filename));
				} 
				catch (Exception e)
				{
					System.out.println("ERROR: COULD NOT OPEN CLUSTER FILES [" + filename + "]. @InputOutput.writeClustersByChrom");
					System.out.println(e.getMessage());
					return false;
				}
			}
			filename = clusterOutputDir + "allChrom.xls";
			try {
				chrom[ALL] = new PrintWriter(new File(filename));
			} 
			catch (Exception e)
			{
				System.out.println("ERROR: COULD NOT OPEN CLUSTER FILES [" + filename + "]. @InputOutput.writeClustersByChrom");
				System.out.println(e.getMessage());
				return false;
			}
			
			try {
				for (int i = 0; i < clusters.size(); i++)
				{
					for (InversionCluster clu : clusters.get(i))
					{
						chrom[i].println(clu);
						chrom[ALL].println(clu);
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
			// close them so thing don't get shitty
			for (int i = 0; i < clusters.size(); i++)
			{
				chrom[i].close();
			}
			chrom[ALL].close();
			return status;
		} // end of printing split clones
		System.out.println("Cluster files not written - Not required by call method @InputOutput.writeClustersByChrom");
		return false;
	}
	
	
	/*
	 * write a bed file from 1000gdeletions stuff
	 */
	/*public static boolean writeDeletionBedFile(ArrayList<DNAInterval> deletions, String deletionOutputDir)
	{
		if (!deletionOutputDir.equals(""))
		{
			PrintWriter bedFile;
			// temp variables
			String filename = deletionOutputDir + "1000gDeletion.bed";
			boolean status;
			System.out.println("*************************************************");
			System.out.println("Creating deletion BED file from 1000genome...");
			// a file for each chromosome (1-22 + X) and one for all chromosomes
			try {
				bedFile = new PrintWriter(new File(filename));
			} 
			catch (Exception e)
			{
				System.out.println("ERROR: COULD NOT OPEN DELETION BED FILE @InputOutput.writeDeletionBedFile");
				return false;
			}
			try {
				for (DNAInterval del : deletions)
				{
					bedFile.println(del);
				}
				status = true;
			}
			catch(Exception e)
			{
				System.out.println("ERROR: COULD NOT DELETION BED FILE @InputOutput.writeDeletionBedFile");
				e.printStackTrace();
				status = false;
			}
			bedFile.close();
			return status;
		} // end of printing deletion bed file
		System.out.println("Deletion bed file not written - Not required by call method @InputOutput.writeDeletionBedFile");
		return false;
	}*/
	
}
