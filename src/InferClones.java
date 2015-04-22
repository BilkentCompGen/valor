import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;


public class InferClones {

        /**
         * @param args
         */
        public static void main(String[] args) throws IOException {
                final int N = Config.WINDOW_SIZE;//Integer.parseInt(args[4], 10);
                final int OVERLAP = Config.EXTENSION;//Integer.parseInt(args[5], 10);
                final double COVERAGE = Config.MIN_COVERAGE;//Double.parseDouble(args[6]) / 100.0;

                String chr;
                String input = args[0]; // on each line: chr start end
                String output = args[1];
                String [] line;
                Scanner scan;
                PrintWriter write;
                Interval v, w, r;
                int covered, start, end, begin;
                //double percent;

                ArrayList<Interval> windows;
                ArrayList<Interval> regions;
                ArrayList<Interval> reads;
                //System.out.println(chr + "\t" + length);
                //System.out.println("READING READS...");
                write = new PrintWriter(new File(output));
                for (int cnt = 1; cnt <= 24; cnt++) {
                        if (cnt == 23) {
                                chr = "chrX";
                        }
			else if (cnt == 24) {
                                chr = "chrY";
                        }
                        else {
                                chr     = "chr" + Integer.toString(cnt);
                        }
                        System.out.println("chrom:\t" + chr);
                        windows = new ArrayList<Interval>();
                        regions = new ArrayList<Interval>();
                        reads = new ArrayList<Interval>();
                        scan = new Scanner(new File(input));

                        while (scan.hasNextLine())
                        {
                                line = scan.nextLine().split("\t");
                                if (chr.equals(line[0]))
                                {
                                        start = Integer.parseInt(line[1]);
                                        end   = Integer.parseInt(line[2]);
                                        if (start < end)
                                        {
                                                reads.add(new Interval(start, end));
                                        }
                                }
                        }
                        scan.close();
                        if (reads.size() == 0) continue;
                        //Sorting reads
                        //System.out.println("No. of reads:\t" + reads.size());
                        //System.out.println("SORTING READS...");
                        Collections.sort(reads, new Comparator<Interval>() {
                                @Override
                                public int compare(Interval  a, Interval  b)
                                {
                                    if (a.start < b.start) return -1;
                                    if (a.start > b.start) return 1;
                                    if (a.end < b.end) return -1;
                                    if (a.end > b.end) return 1;
                                    return 0;
                                }
                        });
                        begin = 0;
                        v = new Interval(0,0);
                        // making N-windows with more than 50% coverage
                        //System.out.println("MAKING >50% COVERAGE WINDOWS...");
                        //percent = 0.2;
                        //System.out.print("...");
                        int firstOverlap, lastOverlap;

                        for (int i = reads.get(0).start-N/2; i < reads.get(reads.size()-1).end+N/2; i+=Config.READ_LENGTH)
                        {
                                ////System.out.println((i+1)*100 / length);
                                /*if (1.0*i / (length-N) >= percent)
                                {
                                        //System.out.print(((int) (percent*100)) + "%...");
                                        percent += 0.2;
                                }*/
                                v = new Interval(i, i + N);
                                w = new Interval(v.start, v.start);
                                covered = 0;
                                firstOverlap = lastOverlap = -1;
                               for (int j = begin; j < reads.size(); j++)
                                {
                                        r = reads.get(j);
                                        if (r.end < v.start)
                                        {
                                                begin++;
                                        }
                                        else if (r.start <= v.end)
                                        {
                                                if (firstOverlap == -1)
                                                {
                                                        firstOverlap = r.start;
                                                }
                                                lastOverlap = r.end;
                                                // they overlap
                                                // r is in i
                                                /*if (r.start >= temp.start && r.end <= temp.end) {
                                                        // do nothing
                                                }
                                                // r extends i
                                                else */if (r.start <= w.end && r.end > w.end) {
                                                        // update end
                                                        w.end = r.end;
                                                }
                                                // i is in v but does not overlap with i
                                                else if (r.start > w.end)
                                                {
                                                        // update coverage
                                                        covered += w.size();
                                                        w.start = r.start;
                                                        w.end = Math.min(r.end, v.end);
                                                }
                                        }
                                        // break condition
                                        else if (r.start > v.end)
                                        {
                                                covered += w.size();
                                                break;
                                        }
                                }
                                double x = (1.0 * covered) / N;
                                if(x >= COVERAGE)
                                {
                                        windows.add(new Interval(Math.max(v.start, firstOverlap), Math.min(v.end, lastOverlap)));
                                }
                        }
                        System.gc();
                        //System.out.println("100%");
                        //System.out.println("No. of windows:\t" + windows.size());
                        //System.out.println("MERGING OVERLAPS...");
                        for (int i = 0; i < windows.size(); i++)
                        {
                                v = windows.get(i);
                                for (int j = i+1; j < windows.size(); j++)
                                {
                                        w = windows.get(j);
                                        // merge if overlaps
                                        if (w.start <= v.end + OVERLAP)
                                        {

                                                // extend temp
                                                v.end = w.end;
                                                // remove j
                                                windows.remove(j);
                                                j--;
                                        }
                                }
                                if (v.size() >= N)
                                        regions.add(v);
                        }
                        System.gc();
                        //System.out.println("No. of regions:\t" + windows.size());
                        //System.out.println("WRITING TO FILE...");

                        for (Interval i : regions)
                        {
                                write.println(chr + "\t" + i + "\t" + args[2]);
                        }
                }
                write.close();
        }
}

