import java.util.*;

public class KLocalSearch {


    // LOCAL SEARCH IMPROVEMENTS

    private static double[] min_cost_choice_ls(int cur_node, ArrayList<Integer> cluster_sizes, PriorityQueue<Pair> cluster_sizes_heap, HashMap<Integer, Integer> cluster_map, ArrayList<ArrayList<Integer>> prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cluster_sizes.size();         
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);

        // Calculate positive connections to established clusters 
        HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
        int total_connections = 0;
        for (int i = 0; i < cur_edges.size(); i++) {
            int cur_edge = cur_edges.get(i);
            int cluster_label = cluster_map.get(cur_edge);
            if (pos_connections.keySet().contains(cluster_label)) 
                pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
            else
                pos_connections.put(cluster_label, 1);
            total_connections += 1;

        }

        // Figure out best clustering choice 
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i : pos_connections.keySet()) { // i is the cluster label now
            if (i != cluster_map.get(cur_node)) {    
                Double[] cur_val = {total_connections + cluster_sizes.get(i) -  (2.0 * pos_connections.get(i)), (double) i};
                cost_increase.add(cur_val);
            } else {  // decrease cluster size by 1 to account for cur_node
                Double[] cur_val = {total_connections + (cluster_sizes.get(i) - 1) -  (2.0 * pos_connections.get(i)), (double) i};
                cost_increase.add(cur_val);
            }
        }

        Pair min_size_cluster = cluster_sizes_heap.peek();
        if (!pos_connections.keySet().contains(min_size_cluster.getValue())) {
            // add one more choice to consider

            if (min_size_cluster.getValue() != cluster_map.get(cur_node)) { 
                Double[] cur_val = {(double) (min_size_cluster.getKey() + total_connections), (double) min_size_cluster.getValue()};
                cost_increase.add(cur_val); 
            } else {
                Double[] cur_val = {(double) (min_size_cluster.getKey() - 1 + total_connections), (double) min_size_cluster.getValue()};
                cost_increase.add(cur_val);
            }
        }

        Double[] last_val = {(double) total_connections, (double) num_clusters};
        Double[] least_cost = DNode.find_min(cost_increase);

        double[] r_val = {1, least_cost[1]}; // will always be one of the current clusters 
        return r_val;

    }

    public static ArrayList<ArrayList<Integer>> initializeClustering(int num_nodes, int k) {
        ArrayList<ArrayList<Integer>> cur_clustering = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < k; i++) {
            ArrayList<Integer> new_cluster = new ArrayList<Integer>();
            cur_clustering.add(new_cluster);
        }
        for (int i = 0; i < num_nodes; i++) {
            int nextVal = (int) (Math.random() * k);
            cur_clustering.get(nextVal).add(i);
        }
        return cur_clustering;

    }

    public static ArrayList<ArrayList<Integer>> k_local_search_network(ArrayList<ArrayList<Integer>> prob_matrix, int k, double threshold) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
        ArrayList<ArrayList<Integer>> cur_clustering = initializeClustering(num_nodes, k);        

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        PriorityQueue<Pair> cluster_sizes_heap = new PriorityQueue<Pair>(); // to track min cluster size
        for (int i = 0; i < cur_clustering.size(); i++) {
            cluster_sizes.add(cur_clustering.get(i).size());
            Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
            cluster_sizes_heap.add(cluster_info);
        }

        int empty_moves = 0;
        int tolerance = (int) Math.ceil(threshold * num_nodes);    

        while(empty_moves < tolerance) { // while not converged...
            Integer cur_node = (int) (Math.random() * num_nodes);           
            double[] max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_sizes_heap, cluster_map, prob_matrix);
    
            if (max_cluster[0] >= 0.5) {    
                // add to existing cluster
                int old_index = cluster_map.get(cur_node);
                int new_index = (int) max_cluster[1];

                if (old_index != new_index) {

                    cluster_map.put(cur_node, new_index);
                    cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                    cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                    // add cur_node to the clustering
                    Pair cluster_info = new Pair(cur_clustering.get(new_index).size(), new_index);
                    cluster_sizes_heap.remove(cluster_info);
                    Pair cluster_info2 = new Pair(cur_clustering.get(old_index).size(), old_index);
                    cluster_sizes_heap.remove(cluster_info2);

                    cur_clustering.get(old_index).remove(cur_node);
                    cur_clustering.get(new_index).add(cur_node);
                    cluster_info.incrementKey();
                    cluster_sizes_heap.add(cluster_info);
                    cluster_info2.decrementKey();
                    cluster_sizes_heap.add(cluster_info2);

                    // empty_moves = 0; // reset?
                    empty_moves += 1; // NEW: just count total number of node calculations made

                } else {
                    empty_moves += 1;
                }
            }
            
        }

        // remove any empty clusters
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

    // VERSION 2: permute nodes, make complete passes until convergence or early stopping

    public static ArrayList<ArrayList<Integer>> k_local_search_network2(ArrayList<ArrayList<Integer>> prob_matrix, int k, double threshold) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
        ArrayList<ArrayList<Integer>> cur_clustering = initializeClustering(num_nodes, k);   
        long cur_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
        
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list); // get permutation of node set 

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        PriorityQueue<Pair> cluster_sizes_heap = new PriorityQueue<Pair>(); // to track min cluster size
        for (int i = 0; i < cur_clustering.size(); i++) {
            cluster_sizes.add(cur_clustering.get(i).size());
            Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
            cluster_sizes_heap.add(cluster_info);
        }

        int iterations = 0;

        while(true) { // while not converged...
            if (threshold >= 0 && iterations >= threshold) // use threshold = -1 to wait until convergence
                break; 

            for (int i = 0; i < num_nodes; i++) {

                Integer cur_node = ecs_list.get(i);        
                double[] max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_sizes_heap, cluster_map, prob_matrix);
    
                if (max_cluster[0] >= 0.5) {    
                    // add to existing cluster
                    int old_index = cluster_map.get(cur_node);
                    int new_index = (int) max_cluster[1];

                    if (old_index != new_index) {

                        cluster_map.put(cur_node, new_index);
                        cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                        cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                        // add cur_node to the clustering
                        Pair cluster_info = new Pair(cur_clustering.get(new_index).size(), new_index);
                        cluster_sizes_heap.remove(cluster_info);
                        Pair cluster_info2 = new Pair(cur_clustering.get(old_index).size(), old_index);
                        cluster_sizes_heap.remove(cluster_info2);

                        cur_clustering.get(old_index).remove(cur_node);
                        cur_clustering.get(new_index).add(cur_node);
                        cluster_info.incrementKey();
                        cluster_sizes_heap.add(cluster_info);
                        cluster_info2.decrementKey();
                        cluster_sizes_heap.add(cluster_info2);

                    } 
                }
            }
            iterations += 1; 
            long new_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
            if (new_score >= cur_score)
                break; // convergence, even if num iterations not met yet
            cur_score = new_score;
            
            
        }

        // remove any empty clusters
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

    public static ArrayList<ArrayList<Integer>> k_local_search_network3(ArrayList<ArrayList<Integer>> prob_matrix, ArrayList<ArrayList<Integer>> cur_clustering, int k, double threshold) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  

            Input outside clustering to be improved
        */
        int num_nodes = prob_matrix.size();
        long cur_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
        
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list); // get permutation of node set 

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        PriorityQueue<Pair> cluster_sizes_heap = new PriorityQueue<Pair>(); // to track min cluster size
        for (int i = 0; i < cur_clustering.size(); i++) {
            cluster_sizes.add(cur_clustering.get(i).size());
            Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
            cluster_sizes_heap.add(cluster_info);
        }

        int iterations = 0;

        while(true) { // while not converged...
            if (threshold >= 0 && iterations >= threshold) // use threshold = -1 to wait until convergence
                break; 

            for (int i = 0; i < num_nodes; i++) {

                Integer cur_node = ecs_list.get(i);        
                double[] max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_sizes_heap, cluster_map, prob_matrix);
    
                if (max_cluster[0] >= 0.5) {    
                    // add to existing cluster
                    int old_index = cluster_map.get(cur_node);
                    int new_index = (int) max_cluster[1];

                    if (old_index != new_index) {

                        cluster_map.put(cur_node, new_index);
                        cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                        cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                        // add cur_node to the clustering
                        Pair cluster_info = new Pair(cur_clustering.get(new_index).size(), new_index);
                        cluster_sizes_heap.remove(cluster_info);
                        Pair cluster_info2 = new Pair(cur_clustering.get(old_index).size(), old_index);
                        cluster_sizes_heap.remove(cluster_info2);

                        cur_clustering.get(old_index).remove(cur_node);
                        cur_clustering.get(new_index).add(cur_node);
                        cluster_info.incrementKey();
                        cluster_sizes_heap.add(cluster_info);
                        cluster_info2.decrementKey();
                        cluster_sizes_heap.add(cluster_info2);

                    } 
                }
            }
            iterations += 1; 
            long new_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
            if (new_score >= cur_score)
                break; // convergence, even if num iterations not met yet
            cur_score = new_score;
            
            
        }

        // remove any empty clusters
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

    public static ArrayList<ArrayList<Integer>> k_local_search_network_timed(ArrayList<ArrayList<Integer>> prob_matrix, ArrayList<ArrayList<Integer>> cur_clustering, int k, long startTime, long limit) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  

            Input outside clustering to be improved
        */
        int num_nodes = prob_matrix.size();
        long currentTime; 

        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list); // get permutation of node set 

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        PriorityQueue<Pair> cluster_sizes_heap = new PriorityQueue<Pair>(); // to track min cluster size
        for (int i = 0; i < cur_clustering.size(); i++) {
            cluster_sizes.add(cur_clustering.get(i).size());
            Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
            cluster_sizes_heap.add(cluster_info);
        }

        boolean flag = true;

        while(flag) { // while not converged...

            for (int i = 0; i < num_nodes; i++) {

                currentTime = System.currentTimeMillis() - startTime;
                if (currentTime > limit) {
                    flag = false;
                    break;
                }

                Integer cur_node = ecs_list.get(i);        
                double[] max_cluster = min_cost_choice_ls(cur_node, cluster_sizes, cluster_sizes_heap, cluster_map, prob_matrix);
    
                if (max_cluster[0] >= 0.5) {    
                    // add to existing cluster
                    int old_index = cluster_map.get(cur_node);
                    int new_index = (int) max_cluster[1];

                    if (old_index != new_index) {

                        cluster_map.put(cur_node, new_index);
                        cluster_sizes.set(old_index, cluster_sizes.get(old_index) - 1);
                        cluster_sizes.set(new_index, cluster_sizes.get(new_index) + 1);

                        // add cur_node to the clustering
                        Pair cluster_info = new Pair(cur_clustering.get(new_index).size(), new_index);
                        cluster_sizes_heap.remove(cluster_info);
                        Pair cluster_info2 = new Pair(cur_clustering.get(old_index).size(), old_index);
                        cluster_sizes_heap.remove(cluster_info2);

                        cur_clustering.get(old_index).remove(cur_node);
                        cur_clustering.get(new_index).add(cur_node);
                        cluster_info.incrementKey();
                        cluster_sizes_heap.add(cluster_info);
                        cluster_info2.decrementKey();
                        cluster_sizes_heap.add(cluster_info2);

                    } 
                }


            }
            /*
            long new_score = Helper.quick_edit_dist(cur_clustering, prob_matrix);
            if (new_score >= cur_score)
                break; // convergence, even if num iterations not met yet
            cur_score = new_score;
            */

            
        }

        // remove any empty clusters
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


    public static void main(String args[]) {

        String data_set = args[0];
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);


        int ROUNDS = 10;

        int[] k_vals = {5}; // {10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000};

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
        for (int j = 0; j < ROUNDS; j++) {

        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> pivot_result = k_local_search_network2(prob_matrix, k, 0.0);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        pivotTimeTotal += (pivotTime / 1000.0);
        pivotTimes[j] = pivotTime;

        long fillStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fill_result = k_local_search_network2(prob_matrix, k, 1.0);
        long fillTime = System.currentTimeMillis() - fillStart;
        fillTimeTotal += (fillTime / 1000.0);
        fillTimes[j] = fillTime;

        long randStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> rand_result = k_local_search_network2(prob_matrix, k, 5.0);
        long randTime = System.currentTimeMillis() - randStart;
        randTimeTotal += (randTime / 1000.0);
        randTimes[j] = randTime;

        long blendStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> blend_result = k_local_search_network2(prob_matrix, k, 10.0);
        long blendTime = System.currentTimeMillis() - blendStart;
        blendTimeTotal += (blendTime / 1000.0); 
        blendTimes[j] = blendTime;

        long hybridStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fix = k_local_search_network2(prob_matrix, k, -1);
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

        System.out.println("LocalSearch(0) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(0) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(1) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(fillTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(1) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(fillScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(5) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(randTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(5) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(randScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(10) times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(10) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(-1)  times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("LocalSearch(-1) scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridScores[i] + " ");
        System.out.println();
        System.out.println();

        System.out.println("Average LocalSearch(0) time: " + pivotTimeTotal / ((double) ROUNDS));
        System.out.println("Average LocalSearch(1) time: " + fillTimeTotal / ((double) ROUNDS));
        System.out.println("Average LocalSearch(5) time: " + randTimeTotal / ((double) ROUNDS));
        System.out.println("Average LocalSearch(10) time: " + blendTimeTotal / ((double) ROUNDS));
        System.out.println("Average LocalSearch(-1) time: " + hybridTimeTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average LocalSearch(0) score: " + pivotScoreTotal / ((double) ROUNDS));
        System.out.println("Average LocalSearch(1) score: " + fillScoreTotal / ((double) ROUNDS));
        System.out.println("Average LocalSearch(5) score: " + randScoreTotal / ((double) ROUNDS));
	System.out.println("Average LocalSearch(10) score: " + blendScoreTotal / ((double) ROUNDS));
	System.out.println("Average LocalSearch(-1) score: " + hybridScoreTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average LocalSearch(0) num clusters: " + pivotNumClusters / ((double) ROUNDS));
        System.out.println("Average LocalSearch(1) num clusters: " + fillNumClusters / ((double) ROUNDS));
        System.out.println("Average LocalSearch(5) num clusters: " + randNumClusters / ((double) ROUNDS));
        System.out.println("Average LocalSearch(10) num clusters: " + blendNumClusters / ((double) ROUNDS));
        System.out.println("Average LocalSearch(-1) num clusters: " + hybridNumClusters / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average LocalSearch(0) max cluster size: " + largestPivotCluster / ((double) ROUNDS));
        System.out.println("Average LocalSearch(1) max cluster size: " + largestFillCluster / ((double) ROUNDS));
        System.out.println("Average LocalSearch(5) max cluster size: " + largestRandCluster / ((double) ROUNDS));
        System.out.println("Average LocalSearch(10) max cluster size: " + largestBlendCluster / ((double) ROUNDS));
        System.out.println("Average LocalSearch(-1) max cluster size: " + largestHybridCluster / ((double) ROUNDS));     
        System.out.println();
        }



    }


}
