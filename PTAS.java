import java.util.*;
import java.io.*;

public class PTAS {

    // https://www.geeksforgeeks.org/count-number-of-ways-to-partition-a-set-into-k-subsets/
    public static int countP(int n, int k)
    {
        // Table to store results of subproblems
        int[][] dp = new int[n+1][k+1];
     
        // Base cases
        for (int i = 0; i <= n; i++)
            dp[i][0] = 0;
        for (int i = 0; i <= k; i++)
            dp[0][k] = 0;
     
        // Fill rest of the entries in dp[][]
        // in bottom up manner
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= k; j++)
                if (j == 1 || i == j)
                    dp[i][j] = 1;
                else
                    dp[i][j] = j * dp[i - 1][j] + dp[i - 1][j - 1];
         
        return dp[n][k];
     
    }

    private static ArrayList<ArrayList<Integer>> convertPartition(ArrayList<Integer> input, int[] partition, int k) {

        ArrayList<ArrayList<Integer>> clustering = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < k; i++) {
            ArrayList<Integer> cluster = new ArrayList<Integer>();
            clustering.add(cluster);
        }
        for (int i = 0; i < partition.length; i++) {
            //int val = partition[i] % k;
            //if (val < 0)
            //    val += k;           
            clustering.get(partition[i]).add(input.get(i));
        }       

        return clustering;
    }

    // https://www.informatik.uni-ulm.de/ni/Lehre/WS03/DMM/Software/partitions.pdf
    public static ArrayList<ArrayList<ArrayList<Integer>>> partition_k(ArrayList<Integer> input, int k) {

        ArrayList<ArrayList<ArrayList<Integer>>> clusterings = new ArrayList<ArrayList<ArrayList<Integer>>>();

        int num_nodes = input.size();
        // initial partition 
        int[] partition = new int[num_nodes];
        int[] max_sequence = new int[num_nodes];
        for (int i = num_nodes - k + 1; i < num_nodes; i++) {
            partition[i] = i - (num_nodes - k);
            max_sequence[i] = i - (num_nodes - k);
        }

        // iterate sequence

        // TODO: maybe don't store all possible clusterings, but generate one at a time
        boolean change = true;
        while (change) {
            // System.out.println(Arrays.toString(partition));
            // System.out.println(Arrays.toString(max_sequence));
            clusterings.add(convertPartition(input, partition, k));
            change = false;
            for (int i = num_nodes -1; i > 0; i--) {            
                if (partition[i] < k-1 && partition[i] <= max_sequence[i-1]) {
                    change = true;
                    partition[i] += 1;
                    max_sequence[i] = Math.max(max_sequence[i], partition[i]);
                    for (int j = i + 1; j <= num_nodes - (k - max_sequence[i]); j++) {
                        partition[j] = 0;
                        max_sequence[j] = max_sequence[i];
                    }
                    for (int j = num_nodes - (k - max_sequence[i]) + 1; j <= num_nodes - 1; j++) {
                        // order?
                        max_sequence[j] = k - (num_nodes - j);
                        partition[j] = max_sequence[j];                       

                    }
                    break;

                }

            }

        }

        return clusterings;

    }

    private static int intercluster_weight_rn(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<ArrayList<Integer>> prob_matrix) { 
        int total = 0;
        for (int i : cluster1) {
            for (int k = 0; k < cluster2.size(); k++) {
                int j = cluster2.get(k);
                if (prob_matrix.get(i).contains(j))
                    total += 1;
            }
        }
        return total; 
    }


    public static ArrayList<ArrayList<Integer>> maxAgree(ArrayList<ArrayList<Integer>> prob_matrix, int k, int sample_size, double epsilon) {

        // assume 0 < epsilon <= 1
        int num_nodes_graph = prob_matrix.size();
        int subset_num = (int) Math.ceil(4 / epsilon);
        int subset_size = (int) Math.ceil( ((double) num_nodes_graph) / subset_num);

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes_graph; i++)
            permutation.add(i);
        Collections.shuffle(permutation); // random ordering of input nodes

        // create random samples
        ArrayList<HashSet<Integer>> sample_sets = new ArrayList<HashSet<Integer>>();

        for (int i = 0; i < subset_num; i++) {
            HashSet<Integer> current_sample = new HashSet<Integer>();
            while (current_sample.size() < sample_size) {
                int nextVal = (int) (Math.random() * num_nodes_graph);
                if (!(nextVal >= i * subset_size && nextVal < (i+1) * subset_size)) // exclude one interval
                    current_sample.add(nextVal);
            }
            sample_sets.add(current_sample);
        }

        long best_score = num_nodes_graph * (num_nodes_graph - 1) / 2; // worst possible clustering score
        ArrayList<ArrayList<Integer>> best_clustering = null;

        for (int ii = 0; ii < sample_sets.size(); ii++) {
            // Generate all possible k clusterings of current sample
            HashSet<Integer> cur_sample = sample_sets.get(ii);
            ArrayList<Integer> inputs = new ArrayList<Integer>();
            inputs.addAll(cur_sample);

            // PARTITION STUFF 
            int num_nodes = inputs.size();
            int[] partition = new int[num_nodes];
            int[] max_sequence = new int[num_nodes];
            for (int i = num_nodes - k + 1; i < num_nodes; i++) {
                partition[i] = i - (num_nodes - k);
                max_sequence[i] = i - (num_nodes - k);
            }
    
            // iterate sequence
    
            // TODO: maybe don't store all possible clusterings, but generate one at a time
            boolean change = true;


            while(change) {
                // find better clustering...


                ArrayList<ArrayList<Integer>> sample_clustering = convertPartition(inputs, partition, k);
                ArrayList<ArrayList<Integer>> cluster_copy = new ArrayList<ArrayList<Integer>>(); // deep copy
                for (int a = 0; a < sample_clustering.size(); a++) {
                    ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                    cur_cluster.addAll(sample_clustering.get(a));
                    cluster_copy.add(cur_cluster);
                }

                // make clustering choice based on sample clustering ONLY
                for (int a = 0; a < num_nodes_graph; a++) {
                    int cur_node = permutation.get(a);
                    if (!cur_sample.contains(cur_node)) {
                        // ----------------------
                        int[] cluster1 = {cur_node};
                        int[] total_weights = new int[sample_clustering.size()];
                        int all_weights = 0;
                        for (int b = 0; b < sample_clustering.size(); b++){
                            int cur_icweight = intercluster_weight_rn(cluster1, sample_clustering.get(b), prob_matrix);
                            total_weights[b] = cur_icweight;
                            all_weights += cur_icweight;
                        }       
        
                        int cost_increase = all_weights + sample_clustering.get(0).size() -  2 * total_weights[0];
                        int cluster_assignment = 0;
                        for (int b = 1; b < sample_clustering.size(); b++) {
                            int cur_val = all_weights + sample_clustering.get(b).size() -  2 * total_weights[b];
                            if (cur_val < cost_increase) {
                                cost_increase = cur_val;
                                cluster_assignment = b;
                            }
                        }
                        cluster_copy.get(cluster_assignment).add(cur_node);
                        // ----------------------
                    }
                }

                // Check if newly formed clustering is better than current best
                long new_score = Helper.quick_edit_dist(cluster_copy, prob_matrix);
                if (new_score < best_score) {
                    best_score = new_score;
                    best_clustering = cluster_copy;
                }

                // PARTITION STUFF 
                change = false;
                for (int i = num_nodes -1; i > 0; i--) {            
                    if (partition[i] < k-1 && partition[i] <= max_sequence[i-1]) {
                        change = true;
                        partition[i] += 1;
                        max_sequence[i] = Math.max(max_sequence[i], partition[i]);
                        for (int j = i + 1; j <= num_nodes - (k - max_sequence[i]); j++) {
                            partition[j] = 0;
                            max_sequence[j] = max_sequence[i];
                        }
                        for (int j = num_nodes - (k - max_sequence[i]) + 1; j <= num_nodes - 1; j++) {
                            // order?
                            max_sequence[j] = k - (num_nodes - j);
                            partition[j] = max_sequence[j];                       
    
                        }
                        break;
    
                    }
    
                }


            }

        }

        // System.out.println("Best clustering cost: " + best_score);
        return best_clustering;
    }


    public static ArrayList<ArrayList<Integer>> minDisagree(ArrayList<ArrayList<Integer>> prob_matrix, int k, int sample_size) {

        int num_nodes_graph = prob_matrix.size();
        // int DIVIDING_FACTOR = 3;

        if (num_nodes_graph < k) // base case 
            k = num_nodes_graph;

        HashSet<Integer> current_sample = new HashSet<Integer>();
        if (num_nodes_graph > sample_size) {
            while (current_sample.size() < sample_size) {
                int nextVal = (int) (Math.random() * num_nodes_graph);
                current_sample.add(nextVal);
            }
        } else {
            for (int i = 0; i < num_nodes_graph; i++)
                current_sample.add(i);
        }

        ArrayList<ArrayList<Integer>> best_clustering = maxAgree(prob_matrix, k, sample_size, 1.0);
        long best_score = Helper.quick_edit_dist(best_clustering, prob_matrix); // num_nodes * (num_nodes - 1) / 2; // worst possible clustering score
        

        // get all possible clusterings
        ArrayList<Integer> inputs = new ArrayList<Integer>();
        inputs.addAll(current_sample);
        // ArrayList<ArrayList<ArrayList<Integer>>> sample_clusterings = partition_k(inputs, k);

        int num_nodes = inputs.size();
        int[] partition = new int[num_nodes];
        int[] max_sequence = new int[num_nodes];
        for (int i = num_nodes - k + 1; i < num_nodes; i++) {
            partition[i] = i - (num_nodes - k);
            max_sequence[i] = i - (num_nodes - k);
        }

        
        boolean change = true;
        while(change) {
            // find better clustering...
            ArrayList<ArrayList<Integer>> sample_clustering = convertPartition(inputs, partition, k);
            ArrayList<ArrayList<Integer>> cluster_copy = new ArrayList<ArrayList<Integer>>(); // deep copy
            for (int a = 0; a < sample_clustering.size(); a++) {
                ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                cur_cluster.addAll(sample_clustering.get(a));
                cluster_copy.add(cur_cluster);
            }

            // make clustering choice based on sample clustering ONLY
            for (int a = 0; a < num_nodes_graph; a++) {
                int cur_node = a;
                if (!current_sample.contains(cur_node)) {
                    // ----------------------
                    int[] cluster1 = {cur_node};
                    int[] total_weights = new int[sample_clustering.size()];
                    int all_weights = 0;
                    for (int b = 0; b < sample_clustering.size(); b++){
                        int cur_icweight = intercluster_weight_rn(cluster1, sample_clustering.get(b), prob_matrix);
                        total_weights[b] = cur_icweight;
                        all_weights += cur_icweight;
                    }       
        
                    int cost_increase = all_weights + sample_clustering.get(0).size() -  2 * total_weights[0];
                    int cluster_assignment = 0;
                    for (int b = 1; b < sample_clustering.size(); b++) {
                        int cur_val = all_weights + sample_clustering.get(b).size() -  2 * total_weights[b];
                        if (cur_val < cost_increase) {
                            cost_increase = cur_val;
                            cluster_assignment = b;
                        }
                    }
                    cluster_copy.get(cluster_assignment).add(cur_node);
                    // ----------------------
                }
            }

            /* SKIP RECURSIVE STEP FOR FASTER RUNTIME
            // Scan for small clusters
            int threshold = (int) Math.ceil( num_nodes / (2.0 * k));
            ArrayList<Integer> smallClusters = new ArrayList<Integer>();
            for (int a = 0; a < cluster_copy.size(); a++) {
                if (cluster_copy.get(a).size() < threshold) {
                    smallClusters.add(a);
                }
            }

            
            if (smallClusters.size() > 0 && num_nodes > sample_size) { // base case: don't recurse if no nodes were added to the sample
                ArrayList<Integer> smallElements = new ArrayList<Integer>();
                // remove small clusters from current clustering
                for (int a = smallClusters.size() - 1; a >= 0; a--) {
                    smallElements.addAll(cluster_copy.get(smallClusters.get(a)));
                    cluster_copy.remove(smallClusters.get(a));
                }

                // make recursive function call
                HashMap<Integer, Integer> relabelling = new HashMap<Integer, Integer>();
                for (int a = 0; a < smallElements.size(); a++)
                    relabelling.put(smallElements.get(a), a);

                ArrayList<ArrayList<Integer>> cur_prob_matrix = new ArrayList<ArrayList<Integer>>();
                for (int i = 0; i < smallElements.size(); i++) {
                    ArrayList<Integer> new_edges = new ArrayList<Integer>();
                    ArrayList<Integer> cur_edges = prob_matrix.get(smallElements.get(i));
                    for (int a = 0; a < cur_edges.size(); a++){
                        int cur_edge = cur_edges.get(a);
                        if (smallElements.contains(cur_edge))
                            new_edges.add(relabelling.get(cur_edge));
                    }
                    cur_prob_matrix.add(new_edges);
                }

                
                ArrayList<ArrayList<Integer>> adjusted_clustering = minDisagree(cur_prob_matrix, smallClusters.size(), sample_size);

                // replace node labels
                for (int a = 0; a < adjusted_clustering.size(); a++) {
                    ArrayList<Integer> adj_cluster = adjusted_clustering.get(a);
                    for (int i = 0; i < adj_cluster.size(); i++)
                        adj_cluster.set(i, smallElements.get(adj_cluster.get(i)));
                    cluster_copy.add(adj_cluster);
                } 
            }
            */

            // Check if newly formed clustering is better than current best
            long new_score = Helper.quick_edit_dist(cluster_copy, prob_matrix);
            if (new_score <= best_score) {
                best_score = new_score;
                best_clustering = cluster_copy;
            }


            // PARTITION STUFF
            change = false;
            for (int i = num_nodes -1; i > 0; i--) {            
                if (partition[i] < k-1 && partition[i] <= max_sequence[i-1]) {
                    change = true;
                    partition[i] += 1;
                    max_sequence[i] = Math.max(max_sequence[i], partition[i]);
                    for (int j = i + 1; j <= num_nodes - (k - max_sequence[i]); j++) {
                        partition[j] = 0;
                        max_sequence[j] = max_sequence[i];
                    }
                    for (int j = num_nodes - (k - max_sequence[i]) + 1; j <= num_nodes - 1; j++) {
                        // order?
                        max_sequence[j] = k - (num_nodes - j);
                        partition[j] = max_sequence[j];                       

                    }
                    break;

                }

            }

        }

        // System.out.println("Best clustering cost: " + best_score);
        return best_clustering;

    }



    public static void main(String args[]){

        String data_set = "cor_gym"; // args[0];
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
        
        int sample_size = 10; // Integer.parseInt(args[1]);
        //int k = 3;
        //double epsilon = 0.5;
        


        
        /* TEST
        ArrayList<Integer> input =  new ArrayList<Integer>();
        for (int i = 0; i < size; i++)
            input.add(i);
        // {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        // 
        ArrayList<ArrayList<ArrayList<Integer>>> clusterings = partition_k(input, k);
        System.out.println(clusterings.size());
        */

        /*
        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> maxResult = maxAgree(prob_matrix, k, sample_size, epsilon);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        System.out.println("MaxAgree score: " + Helper.quick_edit_dist(maxResult, prob_matrix));
        System.out.println("MaxAgree time: " + pivotTime / 1000.0);

        pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> minResult = minDisagree(prob_matrix, k, sample_size);
        pivotTime = System.currentTimeMillis() - pivotStart;
        System.out.println("MinDisagree score: " + Helper.quick_edit_dist(minResult, prob_matrix));
        System.out.println("MinDisagree time: " + pivotTime / 1000.0);
        */

        int ROUNDS = 2;

        int[] k_vals = {3};

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
        System.out.println("Partition count for " + sample_size + ", " + k + ": " + countP(sample_size, k));
        for (int j = 0; j < ROUNDS; j++) {

        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> pivot_result = minDisagree(prob_matrix, k, sample_size);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        pivotTimeTotal += (pivotTime / 1000.0);
        pivotTimes[j] = pivotTime;

        long fillStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fill_result = pivot_result;
        long fillTime = System.currentTimeMillis() - fillStart;
        fillTimeTotal += (fillTime / 1000.0);
        fillTimes[j] = fillTime;

        long randStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> rand_result = pivot_result;
        long randTime = System.currentTimeMillis() - randStart;
        randTimeTotal += (randTime / 1000.0);
        randTimes[j] = randTime;

        long blendStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> blend_result = pivot_result;
        long blendTime = System.currentTimeMillis() - blendStart;
        blendTimeTotal += (blendTime / 1000.0); 
        blendTimes[j] = blendTime;

        long hybridStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fix = pivot_result;
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

        System.out.println("MaxAgree(1.0) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MaxAgree(1.0) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MaxAgree(0.75) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(fillTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MaxAgree(0.75) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(fillScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MaxAgree(0.5) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(randTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("MaxAgree(0.5) scores: ");
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

        System.out.println("Average MaxAgree(1.0) time: " + pivotTimeTotal / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.75) time: " + fillTimeTotal / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.5) time: " + randTimeTotal / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.25) time: " + blendTimeTotal / ((double) ROUNDS));
        System.out.println("Average MinDisagree time: " + hybridTimeTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average MaxAgree(1.0) score: " + pivotScoreTotal / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.75) score: " + fillScoreTotal / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.5) score: " + randScoreTotal / ((double) ROUNDS));
	System.out.println("Average MaxAgree(0.25) score: " + blendScoreTotal / ((double) ROUNDS));
	System.out.println("Average MinDisagree score: " + hybridScoreTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average MaxAgree(1.0) num clusters: " + pivotNumClusters / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.75) num clusters: " + fillNumClusters / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.5) num clusters: " + randNumClusters / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.25) num clusters: " + blendNumClusters / ((double) ROUNDS));
        System.out.println("Average MinDisagree num clusters: " + hybridNumClusters / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average MaxAgree(1.0) max cluster size: " + largestPivotCluster / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.75) max cluster size: " + largestFillCluster / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.5) max cluster size: " + largestRandCluster / ((double) ROUNDS));
        System.out.println("Average MaxAgree(0.25) max cluster size: " + largestBlendCluster / ((double) ROUNDS));
        System.out.println("Average MinDisagree max cluster size: " + largestHybridCluster / ((double) ROUNDS));     
        System.out.println();
        }




    }





}
