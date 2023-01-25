import java.util.*;

public class RunPTAS2 {

    public static void main(String args[]){

        String data_set = args[0];
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
        
        int sample_size = Integer.parseInt(args[1]);
        int ROUNDS = 10;

        int[] k_vals = {10};

        for (int q = 0; q < k_vals.length; q++) {

        // Collect numbers here
        double[] pivotTimes = new double[ROUNDS];
        long[] pivotScores = new long[ROUNDS];

        double[] fillTimes = new double[ROUNDS];
        long[] fillScores = new long[ROUNDS];

        double[] randTimes = new double[ROUNDS];
        long[] randScores = new long[ROUNDS];

        double[] blendTimes = new double[ROUNDS];
        long[] blendScores = new long[ROUNDS];

        double[] hybridTimes = new double[ROUNDS];
        long[] hybridScores = new long[ROUNDS];


        double pivotTimeTotal = 0;
        double fillTimeTotal = 0;
        double randTimeTotal = 0;
        double blendTimeTotal = 0;
        double hybridTimeTotal = 0;
  
        int pivotNumClusters = 0;
        int fillNumClusters = 0;
        int randNumClusters = 0;
        int blendNumClusters = 0;
        int hybridNumClusters = 0;

        int largestPivotCluster = 0;
        int largestFillCluster = 0;
        int largestRandCluster = 0;
        int largestBlendCluster = 0;
        int largestHybridCluster = 0;

        long pivotScoreTotal = 0;
        long fillScoreTotal = 0;
        long randScoreTotal = 0;
	long blendScoreTotal = 0;
	long hybridScoreTotal = 0;

        int k = k_vals[q]; 

        System.out.println("Start: k = " + k);
        // System.out.println("Partition count for " + sample_size + ", " + k + ": " + PTAS.countP(sample_size, k));
        for (int j = 0; j < ROUNDS; j++) {

        ArrayList<ArrayList<Integer>> initial = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < k; i++) 
            initial.add(new ArrayList<Integer>());
        HashSet<Integer> tricky = new HashSet<Integer>();
        for (int i = 0; i < prob_matrix.size(); i++)
            tricky.add(i);

        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> pivot_result = NewPTAS.additiveKApproxAlg(prob_matrix, k, sample_size);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        pivotTimeTotal += (pivotTime / 1000.0);
        pivotTimes[j] = pivotTime;

        initial = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < k; i++) 
            initial.add(new ArrayList<Integer>());
        tricky = new HashSet<Integer>();
        for (int i = 0; i < prob_matrix.size(); i++)
            tricky.add(i);

        long fillStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fill_result = pivot_result; //  NewPTAS.RunPTAS(prob_matrix, k, 0.5, tricky, initial, 0);
        long fillTime = System.currentTimeMillis() - fillStart;
        fillTimeTotal += (fillTime / 1000.0);
        fillTimes[j] = fillTime;

        initial = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < k; i++) 
            initial.add(new ArrayList<Integer>());
        tricky = new HashSet<Integer>();
        for (int i = 0; i < prob_matrix.size(); i++)
            tricky.add(i);

        long randStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> rand_result = pivot_result; // NewPTAS.RunPTAS(prob_matrix, k, 0.25, tricky, initial, 0);
        long randTime = System.currentTimeMillis() - randStart;
        randTimeTotal += (randTime / 1000.0);
        randTimes[j] = randTime;

        long blendStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> blend_result = pivot_result;
        long blendTime = System.currentTimeMillis() - blendStart;
        blendTimeTotal += (blendTime / 1000.0); 
        blendTimes[j] = blendTime;

        long hybridStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fix = pivot_result; // minDisagree(prob_matrix, k, sample_size);
        long hybridTime = System.currentTimeMillis() - hybridStart;
        hybridTimeTotal += (hybridTime / 1000.0);
        hybridTimes[j] = hybridTime;

        pivotNumClusters += pivot_result.size();
        fillNumClusters += fill_result.size();
        randNumClusters += rand_result.size();
        blendNumClusters += blend_result.size();
        hybridNumClusters += fix.size();

        int max_pivot = 0;
        for (int i = 0; i < pivot_result.size(); i++) {
            if (pivot_result.get(i).size() > max_pivot)
                max_pivot = pivot_result.get(i).size();
        }
        largestPivotCluster += max_pivot;

        int max_fill = 0;
        for (int i = 0; i < fill_result.size(); i++) {
            if (fill_result.get(i).size() > max_fill)
                max_fill = fill_result.get(i).size();
        }
        largestFillCluster += max_fill;

        int max_rand = 0;
        for (int i = 0; i < rand_result.size(); i++) {
            if (rand_result.get(i).size() > max_rand)
                max_rand = rand_result.get(i).size();
        }
        largestRandCluster += max_rand;

        int max_blend = 0;
        for (int i = 0; i < blend_result.size(); i++) {
            if (blend_result.get(i).size() > max_blend)
                max_blend = blend_result.get(i).size();
        }
        largestBlendCluster += max_blend;
        
        int max_hybrid = 0;
        for (int i = 0; i < fix.size(); i++) {
            if (fix.get(i).size() > max_hybrid)
                max_hybrid = fix.get(i).size();
        }
        largestHybridCluster += max_hybrid;
        

        long pivot_score = Helper.quick_edit_dist(pivot_result, prob_matrix);
        pivotScoreTotal += pivot_score;
        pivotScores[j] = pivot_score;

        long fill_score = Helper.quick_edit_dist(fill_result, prob_matrix);
        fillScoreTotal += fill_score;
        fillScores[j] = fill_score;

        long rand_score = Helper.quick_edit_dist(rand_result, prob_matrix);
        randScoreTotal += rand_score;
        randScores[j] = rand_score;

        long blend_score = Helper.quick_edit_dist(blend_result, prob_matrix);
        blendScoreTotal += blend_score;
        blendScores[j] = blend_score;

        long hybrid_score = Helper.quick_edit_dist(fix, prob_matrix); 
        hybridScoreTotal += hybrid_score;
        hybridScores[j] = hybrid_score;

        }

        System.out.println("Finish");
        System.out.println();

        System.out.println("NewPTAS(10) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("NewPTAS(10) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotScores[i] + " ");
        System.out.println();
        System.out.println();
        /*
        System.out.println("NewPTAS(0.50) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(fillTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("NewPTAS(0.50) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(fillScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("NewPTAS(0.25) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(randTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("NewPTAS(0.25) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(randScores[i] + " ");
        System.out.println();
        System.out.println();
        
        System.out.println("MaxAgree(0.25) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MaxAgree(0.25) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendScores[i] + " ");
        System.out.println();
        System.out.println();
        
        System.out.println("MinDisagree times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MinDisagree scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridScores[i] + " ");
        System.out.println();
        System.out.println();
        */

        System.out.println("Average NewPTAS(10) time: " + pivotTimeTotal / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.50) time: " + fillTimeTotal / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.25) time: " + randTimeTotal / ((double) ROUNDS));
        //System.out.println("Average MaxAgree(0.25) time: " + blendTimeTotal / ((double) ROUNDS));
        //System.out.println("Average MinDisagree time: " + hybridTimeTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average NewPTAS(10) score: " + pivotScoreTotal / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.50) score: " + fillScoreTotal / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.25) score: " + randScoreTotal / ((double) ROUNDS));
	//System.out.println("Average MaxAgree(0.25) score: " + blendScoreTotal / ((double) ROUNDS));
	//System.out.println("Average MinDisagree score: " + hybridScoreTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average NewPTAS(10) num clusters: " + pivotNumClusters / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.50) num clusters: " + fillNumClusters / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.25) num clusters: " + randNumClusters / ((double) ROUNDS));
        //System.out.println("Average MaxAgree(0.25) num clusters: " + blendNumClusters / ((double) ROUNDS));
        //System.out.println("Average MinDisagree num clusters: " + hybridNumClusters / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average NewPTAS(10) max cluster size: " + largestPivotCluster / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.50) max cluster size: " + largestFillCluster / ((double) ROUNDS));
        //System.out.println("Average NewPTAS(0.25) max cluster size: " + largestRandCluster / ((double) ROUNDS));
        //System.out.println("Average MaxAgree(0.25) max cluster size: " + largestBlendCluster / ((double) ROUNDS));
        //System.out.println("Average MinDisagree max cluster size: " + largestHybridCluster / ((double) ROUNDS));     
        System.out.println();
        }




    }





}
