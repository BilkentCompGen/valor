
import java.util.Collections;
import java.util.Comparator;
import java.util.ArrayList;
/**
 * @author zee
 *
 */
public class InversionDetection {

	/**
	 * @param args
	 * 0. regionFilename
	 * 1. chromosome
	 * 2. readFilenamePlusPlus
	 * 3. readFilenameMinusMinus
	 * 4. clusterOutputDir			
	 * 5. poolOutputDir
	 * 
	 */
	/******INPUT DIRECTORIES****/
	public static String regionFilename;
	public static String chromosome;
	public static String readFilenamePlusPlus;	
	public static String readFilenameMinusMinus;

	public static void main(String[] args) throws Exception {

		/******************************/
		Config.print();
		/******************************/
		/******INPUT DIRECTORIES****/
		regionFilename = args[0];
		chromosome = args[1];
		readFilenamePlusPlus = args[2];	
		readFilenameMinusMinus = args[3];
		/******************************/
		ArrayList<Inversion> allInversions = new ArrayList<Inversion>();
		ArrayList<ArrayList<InversionCluster>> allClusterList = new ArrayList<ArrayList<InversionCluster>>();
		if (chromosome.equals("all") || chromosome.equals("any"))
		{
			for (int i = 1; i <= 22; i++)
			{
				chromosome = "chr" + i;
				runForAChromosome(chromosome, allClusterList, allInversions);
			}
			chromosome = "chrX";
			runForAChromosome(chromosome, allClusterList, allInversions);
			chromosome = "chrY";
			runForAChromosome(chromosome, allClusterList, allInversions);
		}
		else
		{
			runForAChromosome(chromosome, allClusterList, allInversions);
		}
		System.out.println("ALL INVERSIONS: " + allInversions.size());
		System.out.println("***************************************");
		System.out.println("Writing data to files...");
		IOTools.writeClusters(allClusterList);
		System.out.println("***************************************");
		System.out.println("Done!");
		System.gc();
		
	}
	public static void runForAChromosome(String chromosome, 
			ArrayList<ArrayList<InversionCluster>> allClusterList, 
			ArrayList<Inversion> allInversions) throws Exception
	{
		ArrayList<Region> regionList;
		ArrayList<SplitClone> pairedRegions;
		ArrayList<Inversion> Inversions;
		ArrayList<InversionCluster> clusterList = new ArrayList<InversionCluster>();
		ArrayList<PairedDNAInterval> minusMinusReads, plusPlusReads; //, plusMinusReads, deletionList;
		System.out.println(chromosome);

		System.out.println("Reading reads on " + chromosome + "...");
		minusMinusReads = IOTools.readBedFile(readFilenameMinusMinus, 
				 								  chromosome);
		plusPlusReads = IOTools.readBedFile (readFilenamePlusPlus, 
				 								 chromosome);
		System.out.println("#(++)\t" + plusPlusReads.size());
        System.out.println("#(--)\t" + minusMinusReads.size());
		
		//plusMinusReads = InputOutput.readBedFile(readFilenamePlusMinus.substring(0, readFilenamePlusMinus.lastIndexOf(".")) + chromosome.substring(3) + ".bed", 
		//										 chromosome);
		System.out.println("***************************************");
		regionList = IOTools.readRegions (regionFilename, chromosome);				
		System.out.println("Number of Regions: " + regionList.size());
		pairedRegions = createPairedRegions(regionList);
		System.out.println("Number of Paired Regions: " + pairedRegions.size());
		Inversions = createInversions(chromosome, pairedRegions);
		System.out.println("Number of Inversions: " + Inversions.size());
		// adding support information for inversions
		System.out.println("Updating support information for Inversions...");
		updateInversionSupport(Inversions, plusPlusReads, minusMinusReads);
		// removing unsupported inversions
		System.out.println("Removing unsupported inversions...");
		for (int i = 0; i < Inversions.size(); i++)
		{
			if (Inversions.get(i).minusMinusSupport == 0 || Inversions.get(i).plusPlusSupport == 0)
			{
				Inversions.remove(i);
				i--;
			}
		}
		allInversions.addAll(Inversions);
		System.out.println("Number of Inversions after removing non supported ones: " + Inversions.size());
		//clusterList = createClusters(inversions);
		
		clusterList = findClusters(Inversions);
		// updating cluster support
		System.out.println("Calculating inversion support for each cluster...");
		updateClusterSupport(clusterList, plusPlusReads, minusMinusReads);
		System.out.println("Number of Clusters: " + clusterList.size());
		Collections.sort(clusterList, new Comparator<InversionCluster>() {
                        @Override public int compare(InversionCluster clu1, InversionCluster clu2) {
                                return clu1.breakPoint.left.start - clu2.breakPoint.left.start;
                        }
                });

		allClusterList.add(clusterList);
		System.out.println();
		Inversions = null;
		clusterList = null;
		plusPlusReads = null;
		minusMinusReads = null;
		//plusMinusReads = null;
		System.gc();
	}
	/*
	 * 
	 */
	public static ArrayList<SplitClone> createPairedRegions (ArrayList<Region> regionList)
	{
		// VARIABLES
		ArrayList<SplitClone> pairedList;
		
		// PROGRAM CODE
		
		/**********************************************/
		// CREATE THE PAIRED REGIONS
		pairedList = new ArrayList<SplitClone>();
		try {
			for (Region r1 : regionList)
			{
				for (Region r2 : regionList)
				{
					if (r1.canBePaired(r2))
					{
						pairedList.add(new SplitClone(r1, r2));
					}				
				}
			}
		} catch (Exception e) {
			System.out.println(e);
		}
		return pairedList;
	}
	/*
	 * 
	 */
	public static ArrayList<Inversion> createInversions (String chromosome, ArrayList<SplitClone> pairedList)
	{
		// VARIABLES
		ArrayList<Inversion> splitCloneList;
		
		// PROGRAM CODE
		
		/**********************************************/
		// CREATE THE INVERSIONS
		splitCloneList = new ArrayList<Inversion>();
		try {
			for (SplitClone AB : pairedList)
			{
				for (SplitClone CD : pairedList)
				{
					if (AB.canBeSplitClone(CD))
					{
						splitCloneList.add(new Inversion(chromosome, AB, CD));
					}				
				}
			}
		} catch (Exception e) {
			System.out.println(e);
		}
		return splitCloneList;
	}

	public static void explain(int status)
	{
		System.out.println("1) To fix file enter 2 args: regionFilename*  \"fix\"*");
		/*System.out.println("2) regionFilename " + 				// 0 
							   "chromosome " + 				// 1
							   "readFilenamePlusPlus " +		// 2
							   "readFilenameMinusMinus " +		// 3
							   "readFilenamePlusMinus " +		// 4
							   "deletionFilename " +					// 5
							   "poolOutputDir " +				// 6
							   "clusterOutputDir " +			// 7
							   "");*/
		// TODO
		System.out.println("2) regionFilename " + 				// 0 
				   "chromosome " + 				// 1
				   "readFilenamePlusPlus " +		// 2
				   "readFilenameMinusMinus " +		// 3
				   "poolOutputDir " +				// 6
				   "clusterOutputDir " +			// 7
				   "");
		System.out.println("All fields are mandatory.");
		System.out.println("chromosome=all for all chromosomes");
		System.exit(status);
	}
	/*
	 * 
	 */
	public static void updateInversionSupport(ArrayList<Inversion> inversions, 
											   ArrayList<PairedDNAInterval> plusPlus, 
											   ArrayList<PairedDNAInterval> minusMinus)/*,
											   ArrayList<DNAInterval> plusMinus,
											   ArrayList<DNAInterval> deletions)*/
	{
		for (Inversion inv : inversions)
		{
			inv.setReadSupport(plusPlus, minusMinus);//, plusMinus);

		}
	}
	public static void updateClusterSupport(ArrayList<InversionCluster> clusters, 
			   ArrayList<PairedDNAInterval> plusPlus, 
			   ArrayList<PairedDNAInterval> minusMinus) throws Exception
	{
		for (InversionCluster clu : clusters)
		{
			clu.setReadSupport(plusPlus, minusMinus, Config.FRAG_SIZE*4) ;
		}
		// sort clusters by support
		Collections.sort(clusters, new Comparator<InversionCluster>() {
			@Override public int compare(InversionCluster clu1, InversionCluster clu2) {
				return (clu2.minusMinusSupport * clu2.plusPlusSupport - clu1.minusMinusSupport * clu1.plusPlusSupport);
			}
		});
		for (int i = 0; i < clusters.size(); i++)
		{
			for (int j = i+1; j < clusters.size(); j++)
			{
				if ((clusters.get(i).breakPoint.left.overlaps(clusters.get(j).breakPoint.left, Config.INV_GAP) || 
					 clusters.get(i).breakPoint.right.overlaps(clusters.get(j).breakPoint.right, Config.INV_GAP) ))
				// the sizes are almost the same
				//(Math.abs((clusters.get(i).breakPoint.right.end - clusters.get(i).breakPoint.left.start) - (clusters.get(j).breakPoint.right.end - clusters.get(j).breakPoint.left.start)) < Config.GAP)
				{
					clusters.remove(j);
					j--;
				}
			}
		}
               
	}
	public static ArrayList<InversionCluster> findClusters(ArrayList<Inversion> inversions) throws Exception
	{
		InversionClique Q;
		ArrayList<InversionCluster> clusterList = new ArrayList<InversionCluster>();
		InversionCluster clu;
		InversionGraph graph = new InversionGraph(Config.QCLIQUE_LAMBDA, 
													Config.QCLIQUE_GAMMA,
													3 + inversions.size() / 1000);
		graph.initializeGraph(inversions);
		System.out.println("Find cliques and clusters...");
		while (graph.size() > 0)
		{
			Q = graph.findMaxBrokenClique(0, true);
			if (Q != null && Q.VPrime > 0)
			{
				// make a cluster for it
				InversionGraph.sortByDegree(Q.nodes);
				clu = new InversionCluster(Q.nodes.get(0).self, Q.VPrime);
				for (InversionNode n : Q.nodes)
				{
					if (clu.isCompatible(n.self))
					{
						clu.addSplitClone(n.self);
					}
				}
				clusterList.add(clu);
				for (int i =0; i < graph.size(); i++)
				{
					if (clu.isCompatible(graph.nodes.get(i)))
					{
						clu.addSplitClone(graph.nodes.get(i).self);
						graph.nodes.remove(i);
						i--;
					}
				}
			}
			else
			{
				break;
			}
		} 
		return clusterList;
	}
}

