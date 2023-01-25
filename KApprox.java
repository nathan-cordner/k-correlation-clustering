import java.util.*;


public class KApprox {

   
    public static ArrayList<ArrayList<Integer>> TwoCluster(ArrayList<ArrayList<Integer>> prob_matrix, int i) {

        int num_nodes = prob_matrix.size();

        ArrayList<Integer> pivotCluster = new ArrayList<Integer>();
        pivotCluster.add(i);

        ArrayList<Integer> neighbors = prob_matrix.get(i);
        for (int j = 0; j < neighbors.size(); j++)
            pivotCluster.add(neighbors.get(j));

        ArrayList<Integer> otherCluster = new ArrayList<Integer>();
        for (int j = 0; j < num_nodes; j++)
            otherCluster.add(j);
        otherCluster.removeAll(pivotCluster);

        ArrayList<ArrayList<Integer>> cur_clustering = new ArrayList<ArrayList<Integer>>();
        cur_clustering.add(pivotCluster);
        if (!otherCluster.isEmpty()) // only non-empty clusters allowed :)
            cur_clustering.add(otherCluster);

        return cur_clustering;

    }

    public static ArrayList<ArrayList<Integer>> BestTwoCluster(ArrayList<ArrayList<Integer>> prob_matrix) {

        int num_nodes = prob_matrix.size();
        if (num_nodes == 1) { // base case
            ArrayList<ArrayList<Integer>> best_clustering = new ArrayList<ArrayList<Integer>>();
            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(0);
            best_clustering.add(cur_cluster);
            best_clustering.add(new ArrayList<Integer>());
            return best_clustering;

        }

        ArrayList<ArrayList<Integer>> best_clustering = null;
        long best_score = num_nodes * (num_nodes - 1) / 2;


        for (int i = 0; i < num_nodes; i++) {

            ArrayList<ArrayList<Integer>> cur_clustering = TwoCluster(prob_matrix, i);

            // pick best 
            long cur_cost =  Helper.quick_edit_dist(cur_clustering, prob_matrix);
            if (cur_cost < best_score) {
                best_score = cur_cost;
                best_clustering = cur_clustering;
            }            
        }

        return best_clustering; 

    }

    // --- LOCAL SEARCH METHODS ---

    private static double[] min_cost_choice_ls(int cur_node, ArrayList<Integer> cluster_sizes, PriorityQueue<Pair> cluster_sizes_heap, HashMap<Integer, Integer> cluster_map, ArrayList<ArrayList<Integer>> prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
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
        // double cur_cluster_cost = -1;
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

        Double[] least_cost = DNode.find_min(cost_increase);

        double[] r_val = {least_cost[0], least_cost[1]}; // will always be one of the current clusters 
        //if (least_cost[0] == cur_cluster_cost)
        //    r_val[1] = (double) cluster_map.get(cur_node); // prefer to stay in current cluster if possible

        return r_val;

    }

    public static ArrayList<ArrayList<Integer>> k_local_search_network(ArrayList<ArrayList<Integer>> prob_matrix, ArrayList<ArrayList<Integer>> cur_clustering) {
        /*
            Make one pass through given clustering, using "best of one-element moves" strategy
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();

        // get map of <node, cluster_label> pairs
        HashMap<Integer, Integer> cluster_map = Helper.get_clustering_map(cur_clustering);
        
        ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
        PriorityQueue<Pair> cluster_sizes_heap = new PriorityQueue<Pair>(); // to track min cluster size
        for (int i = 0; i < cur_clustering.size(); i++) {
            cluster_sizes.add(cur_clustering.get(i).size());
            Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
            cluster_sizes_heap.add(cluster_info);
        }

        long cur_clustering_cost = Helper.quick_edit_dist(cur_clustering, prob_matrix);

        while(true) { // while not converged...

            // Select best move
            Integer cur_node = -1;
            double[] max_cluster = {num_nodes * 2.0, -1.0};        

            for (int i = 0; i < num_nodes; i++) {           
                double[] cur_value = min_cost_choice_ls(i, cluster_sizes, cluster_sizes_heap, cluster_map, prob_matrix);
                if (cur_value[0] < max_cluster[0]) {
                    int old_index = cluster_map.get(i);
                    int new_index = (int) cur_value[1];
                    if (old_index != new_index) {
                        // update best move
                        cur_node = i;
                        max_cluster[0] = cur_value[0];
                        max_cluster[1] = cur_value[1];

                    }                 
                }
            }

            if (cur_node == -1)
                break; // LS hath converged

            // add to existing cluster
            int old_index = cluster_map.get(cur_node);
            int new_index = (int) max_cluster[1];

            // System.out.println("moving " + cur_node + " from " + old_index + " to " + new_index);

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

            long new_clustering_cost = Helper.quick_edit_dist(cur_clustering, prob_matrix);
            if (new_clustering_cost < cur_clustering_cost)
                cur_clustering_cost = new_clustering_cost;
            else
                break; // no new improvement            
            
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

    public static ArrayList<ArrayList<Integer>> KApproxAlg(ArrayList<ArrayList<Integer>> prob_matrix, int k) {

        int num_nodes = prob_matrix.size();
        ArrayList<ArrayList<ArrayList<Integer>>> all_clusterings = new ArrayList<ArrayList<ArrayList<Integer>>>();

        for (int a = 0; a < num_nodes; a++) {
            // System.out.println("Creating clustering " + (a + 1));

            ArrayList<ArrayList<Integer>> cur_clustering = TwoCluster(prob_matrix, a);

            int lower = 1; 

            while (cur_clustering.size() > lower && cur_clustering.size() < k) { // avoid empty case
                ArrayList<Integer> cluster = cur_clustering.get(lower);
                cur_clustering.remove(lower);
                int c_size = cluster.size();

                HashMap<Integer, Integer> relabelling = new HashMap<Integer, Integer>();
                for (int i = 0; i < c_size; i++)
                    relabelling.put(cluster.get(i), i);
    
                ArrayList<ArrayList<Integer>> cur_prob_matrix = new ArrayList<ArrayList<Integer>>();
                for (int i = 0; i < c_size; i++) {
                    ArrayList<Integer> new_edges = new ArrayList<Integer>();
                    ArrayList<Integer> cur_edges = prob_matrix.get(cluster.get(i));
                    for (int j = 0; j < cur_edges.size(); j++){
                        int cur_edge = cur_edges.get(j);
                        if (cluster.contains(cur_edge))
                            new_edges.add(relabelling.get(cur_edge));
                    }
                    cur_prob_matrix.add(new_edges);
                }
            
                ArrayList<ArrayList<Integer>> adjusted_clustering = BestTwoCluster(cur_prob_matrix);
    
                // replace node labels
                for (int j = 0; j < adjusted_clustering.size(); j++) {
                    ArrayList<Integer> adj_cluster = adjusted_clustering.get(j);
                    for (int i = 0; i < adj_cluster.size(); i++)
                        adj_cluster.set(i, cluster.get(adj_cluster.get(i)));
                    if (!adj_cluster.isEmpty())
                        cur_clustering.add(adj_cluster);
                }                 
                

                lower += 1;
            }

            all_clusterings.add(cur_clustering);
        
        }

        // TODO: Save all generated clusterings, then run LS

        ArrayList<ArrayList<Integer>> best_clustering = null;
        long best_score = num_nodes * (num_nodes - 1) / 2;
        
        for (int i = 0; i < all_clusterings.size(); i++) {
            ArrayList<ArrayList<Integer>> cur_clustering = all_clusterings.get(i);

            // System.out.println("LS on clustering " + (i+1));
            cur_clustering = k_local_search_network(prob_matrix, cur_clustering);

            // pick best 
            long cur_cost =  Helper.quick_edit_dist(cur_clustering, prob_matrix);
            if (cur_cost < best_score) {
                best_score = cur_cost;
                best_clustering = cur_clustering;
            }    
        }

        return best_clustering; 

    }

   


    // --- TEST DRIVER ---

    public static void main (String args[]){

        String data_set = "cor_cora";
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
        
        // set up initial values
        int k = 5;

        
        System.out.println("Starting algorithm");
        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> result = KApproxAlg(prob_matrix, k);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        System.out.println("Finished algorithm");
        System.out.println("Elapsed time: " + (pivotTime / 1000.0));
        System.out.println("Clustering cost: " + Helper.quick_edit_dist(result, prob_matrix));
        System.out.println("Number of clusters: " + result.size());

    }

}