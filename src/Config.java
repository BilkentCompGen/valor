
public class Config {
    /***************READ INFOMATION****************************/
    public static final int READ_LENGTH = 100; // length of each pair end
    public static final int FRAG_SIZE = 1000; // max segment size (distance between paired end reads)

    /**************CLONE INFORMATION**************************/
    public static final int CLONE_MEAN = 150000;
    public static final int CLONE_STD_DEV = 40000;
    public static final int CLONE_MAX = CLONE_MEAN + 3 * CLONE_STD_DEV;
    public static final int CLONE_MIN = CLONE_MEAN - 3 * CLONE_STD_DEV; 
    /*************INVERSION INFORMATION****************************/
    public static final int INV_MIN_SIZE = 500000; // 500K
    public static final int INV_MAX_SIZE = 10000000; // 10M
    public static final int INV_GAP = CLONE_MEAN;
    public static final int INV_OVERLAP = -1 * CLONE_MEAN; // 1 clone size
    /*************GRAPH PROPERTIES****************************/
    public static final double QCLIQUE_LAMBDA = 0.5;
    public static final double QCLIQUE_GAMMA = 0.6;
    /***************INFER CLONES******************/
    public static final int WINDOW_SIZE = 2 * FRAG_SIZE; // min window size
    public static final double MIN_COVERAGE = 0.50; // min coverage of window size
    public static final int EXTENSION = FRAG_SIZE; // extension wing


    public static void print()
    {

        System.out.println("WINDOW_SIZE : " + WINDOW_SIZE);
		System.out.println("MIN_COVERAGE : " + MIN_COVERAGE);
		System.out.println("EXTENSION : " + EXTENSION);


        System.out.println("**************** RUN  DETAILS ********************");
		System.out.println("READ_LENGTH : " + READ_LENGTH);
        System.out.println("FRAG_SIZE : " + FRAG_SIZE);

        //System.out.println("CLONE_SIZE : " + CLONE_SIZE);
        System.out.println("CLONE_MEAN : " + CLONE_MEAN);
        System.out.println("CLONE_STD_DEV : " + CLONE_STD_DEV);
        System.out.println("CLONE_MAX : " + CLONE_MAX);
        System.out.println("CLONE_MIN : " + CLONE_MIN);

        System.out.println("INV_MIN_SIZE : " + INV_MIN_SIZE);
        System.out.println("INV_MAX_SIZE : " + INV_MAX_SIZE);
        System.out.println("INV_GAP : " + INV_GAP);
        System.out.println("INV_OVERLAP : " + INV_OVERLAP);

        System.out.println("QCLIQUE_LAMBDA : " + QCLIQUE_LAMBDA);
		System.out.println("QCLIQUE_GAMMA : " + QCLIQUE_GAMMA);


        System.out.println("*************************************************");  
    }
}
