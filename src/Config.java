
public class Config {
    /***************READ INFOMATION****************************/
    public static final int READ_DIST = 100; // length of each pair end
    public static final int READ_SIZE = 1000; // max segment size (distance between paired end reads)

    /**************CLONE INFORMATION**************************/
    public static final int NORMAL_SIZE = 150000;
    public static final double MEAN = 137000;
    public static final double STD_DEV = 40000;
    public static final double UP_CRITERIA = MEAN + 3 * STD_DEV;
    public static final double LOW_CRITERIA = MEAN - 3 * STD_DEV; 
    /*************INVERSION INFORMATION****************************/
    public static final int MIN_INVERSION_SIZE = 500000; // 500K
    public static final int MAX_INVERSION_SIZE = 10000000; // 10M
    public static final int GAP = 150000;
    public static final int OVERLAP = -2000; // 50K
    public static final int LIMIT = 1500; // max 
    /*************GRAPH PROPERTIES****************************/
    public static final double LAMBDA = 0.5;
    public static final double GAMMA = 0.6;
    /***************INFER CLONES******************/
    public static final int WINDOW = 6500; // min window size
    public static final double COVERAGE = 0.50; // min coverage of window size
    public static final int EXTENSION = 2000; // extension wing


    public static void print()
    {
        System.out.println("**************** RUN  DETAILS ********************");
		System.out.println("NORMAL SIZE : " + NORMAL_SIZE);
        System.out.println("MEAN : " + MEAN);
        System.out.println("STD_DEV : " + STD_DEV);
        System.out.println("UP_CRITERIA : " + UP_CRITERIA);
        System.out.println("LOW_CRITERIA : " + LOW_CRITERIA);
        System.out.println("MIN_INVERSION_SIZE : " + MIN_INVERSION_SIZE);
        System.out.println("MAX_INVERSION_SIZE : " + MAX_INVERSION_SIZE);
        System.out.println("GAP : " + GAP);
		System.out.println("OVERLAP : " + OVERLAP);
        System.out.println("*************************************************");  
    }
}
