import java.util.*;


public class NewPTAS {

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

    private static ArrayList<ArrayList<Integer>> removeEmpty(ArrayList<ArrayList<Integer>> cur_clustering) {
        ArrayList<Integer> empty = new ArrayList<Integer>();
        for (int i = 0; i < cur_clustering.size(); i++) {
            if (cur_clustering.get(i).size() == 0) 
                empty.add(i);
        }
        if (empty.size() != 0) {
            for (int i = empty.size() - 1; i >= 0; i--) {
                cur_clustering.remove(empty.get(i));
            }
        }
        return cur_clustering;

    }

    private static boolean containsDuplicates(ArrayList<ArrayList<Integer>> cur_clustering) {
        boolean flag = false;
        HashSet<Integer> discovered = new HashSet<Integer>();
        for (int i = 0; i < cur_clustering.size(); i++) {
            ArrayList<Integer> cur_cluster = cur_clustering.get(i);
            for (int j = 0; j < cur_cluster.size(); j++) {
                if (discovered.contains(cur_cluster.get(j)))
                    return true;
                discovered.add(cur_cluster.get(j));
            }

        }
        System.out.println("Num nodes discovered = " + discovered.size());  

        return flag;


    }

    public static ArrayList<ArrayList<Integer>> additiveKApproxAlg(ArrayList<ArrayList<Integer>> prob_matrix, int k, int sample_size) {
        /* 
         * Steps:
         * 1. take sample of t <= 1 / (epsilon^2) currently unassigned nodes
         * 2. Find best cut among t nodes
         * 3. Plug in best cut to current clustering assignment
         * 4. Greedily assign remaining nodes in random order (i.e ConstrainedVote) to fill out assignment         * 
         */

        // assume 0 < epsilon < 1
        // System.out.println(sample_size);
        // if (sample_size > 20)
        //    sample_size = 20.0; // just to make sure sample size doesn't blow up
        int num_nodes_graph = prob_matrix.size();

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes_graph; i++)
            permutation.add(i);
        Collections.shuffle(permutation); // random ordering of input nodes

        long best_score = num_nodes_graph * (num_nodes_graph - 1) / 2; // worst possible clustering score
        ArrayList<ArrayList<Integer>> best_clustering = null;

        // HashSet<Integer> settled_nodes = new HashSet<Integer>();    

        // Generate all possible cuts of sample
        ArrayList<Integer> input = new ArrayList<Integer>();
        int index = 0;
        while (input.size() < sample_size && index < num_nodes_graph) {
            input.add(permutation.get(index));
            index += 1;
        }

        // ArrayList<ArrayList<ArrayList<Integer>>> sample_clusterings = partition_k(input, k); // find k partitions
        int num_nodes = input.size();
        int[] partition = new int[num_nodes];
        int[] max_sequence = new int[num_nodes];
        for (int i = num_nodes - k + 1; i < num_nodes; i++) {
            partition[i] = i - (num_nodes - k);
            max_sequence[i] = i - (num_nodes - k);
        }

        // generate one k partition at a time
        boolean change = true;
        while (change) {
            // System.out.println(Arrays.toString(partition));
            // System.out.println(Arrays.toString(max_sequence));
            ArrayList<ArrayList<Integer>> sample_clustering = convertPartition(input, partition, k);
            change = false;

            // First, find best location of cut within current clustering        

            // Second, finish clustering via constrained Vote
            // ROUND 1

            HashMap<Integer, Integer> settled_nodes1 = new HashMap<Integer, Integer>();
            PriorityQueue<Pair> cluster_sizes1 = new PriorityQueue<Pair>();

            for (int a = 0; a < sample_clustering.size(); a++) {
                ArrayList<Integer> cur_cluster = sample_clustering.get(a);
                for (int b = 0; b < cur_cluster.size(); b++)
                    settled_nodes1.put(cur_cluster.get(b), a);
                Pair cur_cluster_info = new Pair(cur_cluster.size(), a);
                cluster_sizes1.add(cur_cluster_info);
            }

            int cur_index = 0;
            while (settled_nodes1.size() < num_nodes_graph) {
                while (settled_nodes1.containsKey(permutation.get(cur_index)))
                    cur_index++;
                Integer cur_node = (int) permutation.get(cur_index);

                double[] max_cluster = DNode.min_prob_edit_dist_choice_rn(cur_node, sample_clustering, k,
                        cluster_sizes1, settled_nodes1, prob_matrix);

                int new_index = (int) max_cluster[1];
                Pair cluster_info = new Pair(sample_clustering.get(new_index).size(), new_index);
                cluster_sizes1.remove(cluster_info);
                // System.out.println(worked);

                sample_clustering.get(new_index).add(cur_node);
                settled_nodes1.put(cur_node, new_index);

                cluster_info.incrementKey();
                cluster_sizes1.add(cluster_info);
            }

            // Finally, check if newly formed clustering is better than current best
            sample_clustering = removeEmpty(sample_clustering);
            long new_score = Helper.quick_edit_dist(sample_clustering, prob_matrix);
            if (new_score < best_score || best_clustering == null) {
                best_score = new_score;
                best_clustering = sample_clustering;
            }

            // K-Partition stuff
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


    public static ArrayList<ArrayList<Integer>> additiveApproxAlg(ArrayList<ArrayList<Integer>> prob_matrix, ArrayList<ArrayList<Integer>> cur_clustering, int k, double epsilon) {
        /* 
         * Steps:
         * 1. take sample of t <= 1 / (epsilon^2) currently unassigned nodes
         * 2. Find best cut among t nodes
         * 3. Plug in best cut to current clustering assignment
         * 4. Greedily assign remaining nodes in random order (i.e ConstrainedVote) to fill out assignment         * 
         */


        // assume 0 < epsilon < 1
        double sample_size = 1 / (Math.pow(epsilon, 2));
        // System.out.println(sample_size);
        // if (sample_size > 20)
        //    sample_size = 20.0; // just to make sure sample size doesn't blow up
        int num_nodes = prob_matrix.size();

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        Collections.shuffle(permutation); // random ordering of input nodes

        long best_score = num_nodes * (num_nodes - 1) / 2; // worst possible clustering score
        ArrayList<ArrayList<Integer>> best_clustering = null;

        HashSet<Integer> settled_nodes = new HashSet<Integer>();
        // sample nodes that are currently unassigned 
        for (int a = 0; a < cur_clustering.size(); a++) {
            ArrayList<Integer> cur_cluster = cur_clustering.get(a);
            for (int b = 0; b < cur_cluster.size(); b++)
                settled_nodes.add(cur_cluster.get(b));
        }        


        // Generate all possible cuts of sample
        ArrayList<Integer> inputs = new ArrayList<Integer>();
        int index = 0;
        while (inputs.size() < sample_size && index < num_nodes) {
            if (!settled_nodes.contains(permutation.get(index)))
                inputs.add(permutation.get(index));
            index += 1;
        }

        ArrayList<ArrayList<ArrayList<Integer>>> sample_clusterings = partition_k(inputs, 2); // 2 for cut
        for (int i = 0; i < sample_clusterings.size(); i++) {

            // First, find best location of cut within current clustering
            ArrayList<ArrayList<Integer>> sample_clustering = sample_clusterings.get(i);
            
            for (int j = 0; j < k; j++) {
                for (int jj = j + 1; jj < k; jj++) {
                    // copy current clustering 
                    ArrayList<ArrayList<Integer>> cluster_copy = new ArrayList<ArrayList<Integer>>(); // deep copy
                    ArrayList<ArrayList<Integer>> cluster_copy2 = new ArrayList<ArrayList<Integer>>(); // deep copy

                    for (int a = 0; a < cur_clustering.size(); a++) {
                        ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                        cur_cluster.addAll(cur_clustering.get(a));
                        cluster_copy.add(cur_cluster);

                        cur_cluster = new ArrayList<Integer>();
                        cur_cluster.addAll(cur_clustering.get(a));
                        cluster_copy2.add(cur_cluster);               
                    }

                    cluster_copy.get(j).addAll(sample_clustering.get(0));
                    cluster_copy.get(jj).addAll(sample_clustering.get(1));

                    cluster_copy2.get(jj).addAll(sample_clustering.get(0));
                    cluster_copy2.get(j).addAll(sample_clustering.get(1));

                    // Second, finish clustering via constrained Vote
                    // ROUND 1

                    HashMap<Integer, Integer> settled_nodes1 = new HashMap<Integer, Integer>();
                    PriorityQueue<Pair> cluster_sizes1 = new PriorityQueue<Pair>();

                    for (int a = 0; a < cluster_copy.size(); a++) {
                        ArrayList<Integer> cur_cluster = cluster_copy.get(a);
                        for (int b = 0; b < cur_cluster.size(); b++)
                            settled_nodes1.put(cur_cluster.get(b), a);
                        Pair cur_cluster_info = new Pair(cur_cluster.size(), a);
                        cluster_sizes1.add(cur_cluster_info);
                    }

                    int cur_index = 0;
                    while (settled_nodes1.size() < num_nodes) {
                        while (settled_nodes1.containsKey(permutation.get(cur_index)))
                            cur_index++;
                        Integer cur_node = (int) permutation.get(cur_index);
                        
                        double[] max_cluster = DNode.min_prob_edit_dist_choice_rn(cur_node, cluster_copy, k, cluster_sizes1, settled_nodes1, prob_matrix);
                
                        int new_index = (int) max_cluster[1];
                        Pair cluster_info = new Pair(cluster_copy.get(new_index).size(), new_index);
                        cluster_sizes1.remove(cluster_info);
                        // System.out.println(worked);
            
                        cluster_copy.get(new_index).add(cur_node);
                        settled_nodes1.put(cur_node, new_index);
            
                        cluster_info.incrementKey();
                        cluster_sizes1.add(cluster_info);
                    }

                    // Finally, check if newly formed clustering is better than current best
                    cluster_copy = removeEmpty(cluster_copy);
                    long new_score = Helper.quick_edit_dist(cluster_copy, prob_matrix);
                    if (new_score < best_score || best_clustering == null) {
                        best_score = new_score;
                        best_clustering = cluster_copy;
                    }

                    if (settled_nodes.isEmpty())
                        break;

                    // ROUND 2

                    HashMap<Integer, Integer> settled_nodes2 = new HashMap<Integer, Integer>();
                    PriorityQueue<Pair> cluster_sizes2 = new PriorityQueue<Pair>();

                    for (int a = 0; a < cluster_copy2.size(); a++) {
                        ArrayList<Integer> cur_cluster = cluster_copy2.get(a);
                        for (int b = 0; b < cur_cluster.size(); b++)
                            settled_nodes2.put(cur_cluster.get(b), a);
                        Pair cur_cluster_info = new Pair(cur_cluster.size(), a);
                        cluster_sizes2.add(cur_cluster_info);
                    }

                    cur_index = 0;
                    while (settled_nodes2.size() < num_nodes) {
                        while (settled_nodes2.containsKey(permutation.get(cur_index)))
                            cur_index++;
                        Integer cur_node = (int) permutation.get(cur_index);
                        
                        double[] max_cluster = DNode.min_prob_edit_dist_choice_rn(cur_node, cluster_copy2, k, cluster_sizes2, settled_nodes2, prob_matrix);
                
                        int new_index = (int) max_cluster[1];
                        Pair cluster_info = new Pair(cluster_copy.get(new_index).size(), new_index);
                        cluster_sizes2.remove(cluster_info);
                        // System.out.println(worked);
            
                        cluster_copy2.get(new_index).add(cur_node);
                        settled_nodes2.put(cur_node, new_index);
            
                        cluster_info.incrementKey();
                        cluster_sizes2.add(cluster_info);
                    }

                    cluster_copy2 = removeEmpty(cluster_copy2);
                    new_score = Helper.quick_edit_dist(cluster_copy2, prob_matrix);
                    if (new_score < best_score) {
                        best_score = new_score;
                        best_clustering = cluster_copy2;
                    }
                }
                if (settled_nodes.isEmpty())
                    break;
            }

        }


        // System.out.println("Best clustering cost: " + best_score);
        return best_clustering;
    }

    private static int function_p(int u, int v, int i, int j, ArrayList<ArrayList<Integer>> prob_matrix) {
        /* computes p_{u, v} (i, j) */

        int val1 = 0;
        if (!prob_matrix.get(u).contains(v))
            val1 = 1;

        int val2 = 0;
        if (i != j)
            val2 = 1;

        int combined_val = Math.abs(val1 - val2);

        return combined_val;
    }

    private static int function_b(HashMap<Integer, Integer> x, int v, int i, ArrayList<ArrayList<Integer>> prob_matrix) {
        int total = 0;
        for (int u = 0; u < prob_matrix.size(); u++) {
            if (u != v) 
                total += function_p(u, v, x.get(u), i, prob_matrix);
        } 
        return total;
    }


    public static ArrayList<ArrayList<Integer>> RunPTAS(ArrayList<ArrayList<Integer>> prob_matrix, int k, double epsilon, HashSet<Integer> tricky, ArrayList<ArrayList<Integer>> cur_clustering, int depth  ) {

        // assume incoming clustering already has k clusters defined 
        // "tricky" variables are currently unassigned 
        // assume delta = 1 for correlation clustering 

        int num_nodes = prob_matrix.size();

        // Line 1: call additive approx alg    
        ArrayList<ArrayList<Integer>> best_clustering = additiveApproxAlg(prob_matrix, cur_clustering, k, (epsilon) / (1 + epsilon));
        long best_score = Helper.quick_edit_dist(best_clustering, prob_matrix);

        // Lines 2 - 3
        if (best_score >= (Math.pow(tricky.size(), 2)) / (31104 * Math.pow(k, 3)) || depth >= k + 1) // constant 
            return best_clustering; 

        // System.out.println("Continue");
        
        // Line 5
        double s = 186624 * Math.pow(k, 4) * Math.log(1440 * Math.pow(k, 3)) / 2.0;
        if (s > 10)
            s = 10.0; // to make sure sample size doesn't blow up

        // Line 6
        HashSet<Integer> sample = new HashSet<Integer>();
        if (s >= tricky.size())
             sample = tricky;
        else {
            // sample s items from tricky
            ArrayList<Integer> permutation = new ArrayList<Integer>();
            permutation.addAll(tricky);
            Collections.shuffle(permutation);
            for (int i = 0; i < s; i++)
                sample.add(permutation.get(i));
        }

        // Line 7 -- just add clusters in order provided for now
        ArrayList<Integer> inputs = new ArrayList<Integer>();
        inputs.addAll(sample);
        ArrayList<ArrayList<ArrayList<Integer>>> sample_clusterings = partition_k(inputs, k);
        for (int a = 0; a < sample_clusterings.size(); a++) {
            ArrayList<ArrayList<Integer>> cur_partition = sample_clusterings.get(a);
            HashMap<Integer, Integer> partition_map = new HashMap<Integer, Integer>();
            for (int i = 0; i < k; i++) {
                ArrayList<Integer> cur_cluster = cur_partition.get(i);
                for (int j = 0; j < cur_cluster.size(); j++)
                    partition_map.put(cur_cluster.get(j), i);
            }

            // Line 8 -- assume "i" in the paper ranges over clusters
            HashMap<Integer, Integer> x_one = new HashMap<Integer, Integer>(); // map v to cluster assignment
            // non-tricky variables are assigned to place in cur_clustering

            for (int v : tricky) {
                // need to find arg min b_hat to assign v 
                double cur_min = tricky.size() * num_nodes;
                int cur_loc = 0; 

                for (int i = 0; i < k; i++) {

                    // first sum 
                    int sum_1 = 0;
                    for (int j = 0; j < s; j++) {
                        // p_{v_j, v}(x*_{v_j}, i)
                        int v_j = inputs.get(j);
                        int vj_loc = partition_map.get(v_j);
                        sum_1 += function_p(v, v_j, vj_loc, i, prob_matrix);
                    }
                    double part1 = (tricky.size() / s) * sum_1;

                    int sum_2 = 0;
                    for (int j = 0; j < cur_clustering.size(); j++) {
                        ArrayList<Integer> cur_cluster = cur_clustering.get(j);
                        for (int jj = 0; jj < cur_cluster.size(); jj++) 
                            sum_2 += function_p(v, cur_cluster.get(jj), j, i, prob_matrix);                   
                    }
                    double total_score = part1 + sum_2;

                    if (total_score < cur_min) {
                        cur_min = total_score;
                        cur_loc = i;
                    }              

                }
                // Line 9
                x_one.put(v, cur_loc);
            }

            // add assignments of non-tricky variables to x_one
            for (int i = 0; i < cur_clustering.size(); i++) {
                ArrayList<Integer> cur_cluster = cur_clustering.get(i);
                for (int j = 0; j < cur_cluster.size(); j++)
                    x_one.put(cur_cluster.get(j), i);
            }

            // Line 10
            HashMap<Integer, Integer> x_two = new HashMap<Integer, Integer>(); // map v to cluster assignment
            for (int v : tricky) {
                double cur_min = num_nodes;
                int cur_loc = 0; 

                for (int i = 0; i < k; i++) {
                    int total = function_b(x_one, v, i, prob_matrix);

                    if (total < cur_min) {
                        cur_min = total;
                        cur_loc = i;
                    }
                }
                x_two.put(v, cur_loc);

            }

            // Line 11
            HashSet<Integer> big_C = new HashSet<Integer>();
            for (int v : tricky) {
                boolean flag = true;
                for(int j = 0; j < num_nodes; j++) {
                    if (j != x_two.get(v)) {
                        if (function_b(x_one, v, x_two.get(v), prob_matrix) >= function_b(x_one, v, j, prob_matrix) - (tricky.size() / (12 * k))) {
                            flag = false;
                            break;
                        }                    
                    }
                }

                if (flag)
                    big_C.add(v);

            }

            // Line 12
            HashSet<Integer> T_prime = new HashSet<Integer>();
            for (int v : tricky) {
                if (!big_C.contains(v))
                    T_prime.add(v);
            }

            // Line 13
            // first make copy of cur_clustering
            ArrayList<ArrayList<Integer>> new_clustering = new ArrayList<ArrayList<Integer>>();
            for (int i = 0; i < cur_clustering.size(); i++) {
                ArrayList<Integer> cur_cluster = cur_clustering.get(i);
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                for (int j = 0; j < cur_cluster.size(); j++)
                    new_cluster.add(cur_cluster.get(j));
                new_clustering.add(new_cluster);
            }
            // now add variables from big_C
            for (int v : big_C) 
                new_clustering.get(x_two.get(v)).add(v);

            // Line 14 
            if (!T_prime.isEmpty()) {
                new_clustering = RunPTAS(prob_matrix, k, epsilon, T_prime, new_clustering, depth + 1);
            }

            long new_score = Helper.quick_edit_dist(new_clustering, prob_matrix);
            if (new_score < best_score) {
                best_score = new_score;
                best_clustering = new_clustering;
            }

        }

        // Line 16
        // System.out.println("Best clustering cost: " + best_score);
        return best_clustering;

    }



    public static void main(String args[]){

        String data_set = "cor_gym";
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
        
        // set up initial values
        int k = 5;
        int sample_size = 10;
        // double epsilon = 0.4;

        ArrayList<ArrayList<Integer>> initial = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < k; i++) 
            initial.add(new ArrayList<Integer>());
        HashSet<Integer> tricky = new HashSet<Integer>();
        for (int i = 0; i < prob_matrix.size(); i++)
            tricky.add(i);
        
        System.out.println("Starting algorithm");
        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> result = additiveKApproxAlg(prob_matrix, k, sample_size); // RunPTAS(prob_matrix, k, epsilon, tricky, initial, 0);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        System.out.println("Finished algorithm");
        System.out.println("Elapsed time: " + (pivotTime / 1000.0));
        System.out.println("Clustering cost: " + Helper.quick_edit_dist(result, prob_matrix));
        System.out.println("Number of clusters: " + result.size());
        // System.out.println("Duplicates = " + containsDuplicates(result));


    }





}
