import java.util.*;

public class DNode {

    // ECS Reverse Compare
    private static class ECSComparator implements Comparator<Double[]> {  
        public int compare (Double[] a1, Double[] a2) {  

            if (a1[0] == a2[0]) {   
                return 0;
            } else if (a1[0] > a2[0]) {  
                return -1;
            } else if (a1[0] < a2[0]) {
                return 1;
            }
            return 0;  
        }  
    }  

    // Deterministic Node
    private static double sum(ArrayList<Double> input_list) {
        double my_sum = 0;
        for (double i : input_list)
            my_sum += i;
        return my_sum;

    }

    private static double sum(double[] input_list) {
        double my_sum = 0;
        for (double i : input_list)
            my_sum += i;
        return my_sum;

    }

    public static Double[] find_min(ArrayList<Double[]> input_list) {

        double cur_min = input_list.get(0)[0];
        int cur_min_index = 0;

        for (int i = 0; i < input_list.size(); i++) {
            double cand = input_list.get(i)[0];
            if (cand < cur_min) {
                cur_min = cand;
                cur_min_index = i;
            }
        }
        return input_list.get(cur_min_index);

    }

    private static double intercluster_weight(int[] cluster1, ArrayList<Integer> cluster2, double[][] prob_matrix) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        for (int i : cluster1) {
            for (int j : cluster2) {
                cluster_weights.add(prob_matrix[i][j]);
            }
        }        

        return sum(cluster_weights) / ( (double) cluster_weights.size()); 
    }
    
    public static ArrayList<Double[]> node_expected_cluster_size(double[][] prob_matrix, int num_nodes) {
        // rank nodes by probability sums 
        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>();

        for (int i = 0; i < num_nodes; i++) {

            double[] cur_row = prob_matrix[i];

            Double[] cur_val = {sum(cur_row), (double) i};
            ecs_list.add(cur_val);
            // System.out.println(cur_val[0]);

        }

        ecs_list.sort(new ECSComparator());
        
        return ecs_list;
    }
    
    private static double[] min_prob_edit_dist_choice(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, double[][] prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size(); 
        
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight(cluster1, cur_clustering.get(i), prob_matrix);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster descreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
    }
    
    public static ArrayList<ArrayList<Integer>> adapted_node(double[][] prob_matrix) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        int num_nodes = prob_matrix.length; 
    
        // top-v():
        ArrayList<Double[]> ecs_list = node_expected_cluster_size(prob_matrix, num_nodes);
        ArrayList<Integer> settled_nodes = new ArrayList<Integer>();
        Integer cur_node = (int) Math.floor(ecs_list.get(0)[1]); // should be integer-valued anyway...
        settled_nodes.add(cur_node);

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
    
        int node_counter = 1;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) Math.floor(ecs_list.get(node_counter)[1]);
            node_counter += 1;
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice(cur_node, cur_clustering, prob_matrix);
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                cur_clustering.add(new_cluster);
            }
    
            settled_nodes.add(cur_node);
        }
            
        return cur_clustering;

    }    

    // --- 0/1 NETWORK ---

    public static ArrayList<Double[]> node_expected_cluster_size_network(ArrayList<ArrayList<Integer>> prob_matrix, int num_nodes) {
        // rank nodes by probability sums 

        ArrayList<Double[]> ecs_list = new ArrayList<Double[]>(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            Double[] cur_val = {(double) prob_matrix.get(i).size(), (double) i};
            ecs_list.add(cur_val);
        }

        ecs_list.sort(new ECSComparator());
        
        return ecs_list;
    }

    private static double intercluster_weight(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<ArrayList<Integer>> prob_matrix, double REPS_PROB) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        int num_reps = (int) Math.ceil(REPS_PROB * cluster2.size());
        for (int i : cluster1) {
            for (int k = 0; k < num_reps; k++) {
                int j = cluster2.get(k);
                double edge_weight = 0.0;
                if (prob_matrix.get(i).contains(j))
                    edge_weight = 1.0;
                cluster_weights.add(edge_weight);
            }
        }      

        return sum(cluster_weights) / ( (double) num_reps); 

        // return sum(cluster_weights) / ( (double) cluster_weights.size()); 
    }
    
    private static double[] min_prob_edit_dist_choice(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, HashMap<Integer, Integer> settled_nodes, ArrayList<ArrayList<Integer>> prob_matrix, double REPS_PROB) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size(); 
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();

        if (settled_nodes.size() < cur_edges.size()) { // DO IT THE OLD WAY

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight(cluster1, cur_clustering.get(i), prob_matrix, REPS_PROB);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster descreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};

        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
        } else { // TAKE ADVANTAGE OF NEIGHBORHOOD ORACLE TO REDUCE "QUERIES" 

            // Calculate positive connections to established clusters 
            HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
            int total_connections = 0;
            for (int i = 0; i < cur_edges.size(); i++) {
                int cur_edge = cur_edges.get(i);
                if (settled_nodes.keySet().contains(cur_edge)) {
                    int cluster_label = settled_nodes.get(cur_edge);
                    if (pos_connections.keySet().contains(cluster_label)) 
                        pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
                    else
                        pos_connections.put(cluster_label, 1);
                    total_connections += 1;
                }
            }

            // Figure out best clustering choice 
            ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
            for (int i : pos_connections.keySet()) { // i is the cluster label now              

                Double[] cur_val = {total_connections + cur_clustering.get(i).size() -  (2.0 * pos_connections.get(i)), (double) i};
                cost_increase.add(cur_val);

            }

            Double[] last_val = {(double) total_connections, (double) num_clusters};
            cost_increase.add(last_val); // cost of opening a new cluster
            Double[] least_cost = find_min(cost_increase);

            if (least_cost[1] < num_clusters) {
                double[] r_val = {1, least_cost[1]};
                return r_val;
            }
            else {
                double[] r_val = {0, least_cost[1]};
                return r_val; // force to open new cluster
            }
        }
    }

    public static ArrayList<ArrayList<Integer>> adapted_node_network(ArrayList<ArrayList<Integer>> prob_matrix, int MAX_REPS) {
        /*
            My approach: order nodes by expected cluster size
            -- settle nodes in that order
            -- iterate over all edges between current node and the settled set
            -- only query those edges that we can't yet infer
            -- break out of the loop if we find a match for the current node
        */
        int num_nodes = prob_matrix.size();
        MAX_REPS = Math.min(MAX_REPS, num_nodes); // don't want max representatives to exceed cur cluster size
    
        // top-v():
        ArrayList<Double[]> ecs_list = node_expected_cluster_size_network(prob_matrix, num_nodes);
        // ArrayList<Integer> settled_nodes = new ArrayList<Integer>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();
        Integer cur_node = (int) Math.floor(ecs_list.get(0)[1]); // should be integer-valued anyway...
        settled_nodes.put(cur_node, 0);

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
    
        int node_counter = 1;
        double REPS_PROB = 1.0;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) Math.floor(ecs_list.get(node_counter)[1]);
            node_counter += 1;

            if (MAX_REPS < settled_nodes.size())
                REPS_PROB = ((double) MAX_REPS) / settled_nodes.size();
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice(cur_node, cur_clustering, settled_nodes, prob_matrix, REPS_PROB);
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                settled_nodes.put(cur_node, cur_clustering.size());
                cur_clustering.add(new_cluster);
            }
    
            // settled_nodes.add(cur_node);
        }
            
        return cur_clustering;

    }   

    // --- RANDOM NODE WITH NEIGHBORHOOD ORACLE, USING AT MOST K CLUSTERS --- // 


    private static double intercluster_weight_rn(int[] cluster1, ArrayList<Integer> cluster2, ArrayList<ArrayList<Integer>> prob_matrix) { 
        // find average value of edges between clusters
        ArrayList<Double> cluster_weights = new ArrayList<Double>();
        for (int i : cluster1) {
            for (int k = 0; k < cluster2.size(); k++) {
                int j = cluster2.get(k);
                double edge_weight = 0.0;
                if (prob_matrix.get(i).contains(j))
                    edge_weight = 1.0;
                cluster_weights.add(edge_weight);
            }
        }      

        // return sum(cluster_weights) / ( (double) cluster2.size()); 

        return sum(cluster_weights) / ( (double) cluster_weights.size()); 
    }
    
    public static double[] min_prob_edit_dist_choice_rn(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, int k, PriorityQueue<Pair> cluster_sizes, HashMap<Integer, Integer> settled_nodes, ArrayList<ArrayList<Integer>> prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size();         
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);
        if (settled_nodes.size() < cur_edges.size()) { // DO IT THE OLD WAY

        for (int i = 0; i < num_clusters; i++){

            int[] cluster1 = {cur_node};
            double cur_icweight = intercluster_weight_rn(cluster1, cur_clustering.get(i), prob_matrix);
            Double[] cur_val = {cur_icweight, (double) i};
            avg_weights.add(cur_val);

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size(); 
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster decreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};
        if (cur_clustering.size() < k)
            cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
        } else { // TAKE ADVANTAGE OF NEIGHBORHOOD ORACLE TO REDUCE "QUERIES" 

            // Calculate positive connections to established clusters 
            HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
            int total_connections = 0;
            for (int i = 0; i < cur_edges.size(); i++) {
                int cur_edge = cur_edges.get(i);
                if (settled_nodes.keySet().contains(cur_edge)) {
                    int cluster_label = settled_nodes.get(cur_edge);
                    if (pos_connections.keySet().contains(cluster_label)) 
                        pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
                    else
                        pos_connections.put(cluster_label, 1);
                    total_connections += 1;
                }
            }

            // Figure out best clustering choice 
            ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
            for (int i : pos_connections.keySet()) { // i is the cluster label now              

                Double[] cur_val = {total_connections + cur_clustering.get(i).size() -  (2.0 * pos_connections.get(i)), (double) i};
                cost_increase.add(cur_val);

            }
            if (cur_clustering.size() == k) {
                Pair min_size_cluster = cluster_sizes.peek();
                if (!pos_connections.keySet().contains(min_size_cluster.getValue())) {
                    // add one more choice to consider
                    Double[] cur_val = {(double) (min_size_cluster.getKey() + total_connections), (double) min_size_cluster.getValue()};
                    cost_increase.add(cur_val); 
                }
            }

            Double[] last_val = {(double) total_connections, (double) num_clusters};
            if (cur_clustering.size() < k)
                cost_increase.add(last_val); // cost of opening a new cluster // REMOVING FOR K CLUSTER
            Double[] least_cost = find_min(cost_increase);

            if (least_cost[1] < num_clusters) {
                double[] r_val = {1, least_cost[1]};
                return r_val;
            }
            else {
                double[] r_val = {0, least_cost[1]};
                return r_val; // force to open new cluster
            }

        }
    }

    public static ArrayList<ArrayList<Integer>> k_random_node_network(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*
            RandomNode ordering // TODO: use same order as corresponding Pivot clustering?
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
    
        // top-v():
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list);

        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();
        Integer cur_node = (int) ecs_list.get(0); 
        settled_nodes.put(cur_node, 0); // first node is assigned to cluster 0 

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);

        // TODO: make finding min cluster more efficient
        PriorityQueue<Pair> cluster_sizes = new PriorityQueue<Pair>(); // to track min cluster size
        // Pair cur_cluster_info = new Pair(first_cluster.size(), cur_clustering.size() - 1);
        // cluster_sizes.add(cur_cluster_info);
    
        int node_counter = 1;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) ecs_list.get(node_counter);
            node_counter += 1;
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice_rn(cur_node, cur_clustering, k, cluster_sizes, settled_nodes, prob_matrix);
    
            if (max_cluster[0] >= 0.5 && cur_clustering.size() == k) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                Pair cluster_info = new Pair(cur_clustering.get(new_index).size(), new_index);
                cluster_sizes.remove(cluster_info);

                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);

                cluster_info.incrementKey();
                cluster_sizes.add(cluster_info);

            } else if (max_cluster[0] >= 0.5) {
                // no need to update heap just yet
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                settled_nodes.put(cur_node, cur_clustering.size());
                cur_clustering.add(new_cluster);

                if (cur_clustering.size() == k) {
                    // initialize heap
                    for (int i = 0; i < k; i++) {
                        Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
                        cluster_sizes.add(cluster_info);
                    }
                }
            }
    
            
        }
            
        return cur_clustering;

    }


    public static ArrayList<ArrayList<Integer>> exactly_k_random_node_network(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*
            RandomNode ordering // TODO: use same order as corresponding Pivot clustering?
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
    
        // top-v():
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list);

        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();
        Integer cur_node = (int) ecs_list.get(0); 
        settled_nodes.put(cur_node, 0); // first node is assigned to cluster 0 

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);

        // TODO: make finding min cluster more efficient
        PriorityQueue<Pair> cluster_sizes = new PriorityQueue<Pair>(); // to track min cluster size
        // Pair cur_cluster_info = new Pair(first_cluster.size(), cur_clustering.size() - 1);
        // cluster_sizes.add(cur_cluster_info);
    
        int node_counter = 1;
        
        while (settled_nodes.size() < num_nodes) {

            if ((num_nodes - settled_nodes.size()) + cur_clustering.size() == k) {
                // put all remaining nodes in singleton clusters
                while (settled_nodes.size() < num_nodes) {

                    cur_node = (int) ecs_list.get(node_counter);
                    node_counter += 1;

                    ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                    new_cluster.add(cur_node);
                    settled_nodes.put(cur_node, cur_clustering.size());
                    cur_clustering.add(new_cluster);
                }
                break;
            }

            cur_node = (int) ecs_list.get(node_counter);
            node_counter += 1;

            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = min_prob_edit_dist_choice_rn(cur_node, cur_clustering, k, cluster_sizes, settled_nodes, prob_matrix);
    
            if (max_cluster[0] >= 0.5 && cur_clustering.size() == k) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                Pair cluster_info = new Pair(cur_clustering.get(new_index).size(), new_index);
                cluster_sizes.remove(cluster_info);

                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);

                cluster_info.incrementKey();
                cluster_sizes.add(cluster_info);

            } else if (max_cluster[0] >= 0.5) {
                // no need to update heap just yet
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                settled_nodes.put(cur_node, cur_clustering.size());
                cur_clustering.add(new_cluster);

                if (cur_clustering.size() == k) {
                    // initialize heap
                    for (int i = 0; i < k; i++) {
                        Pair cluster_info = new Pair(cur_clustering.get(i).size(), i);
                        cluster_sizes.add(cluster_info);
                    }
                }
            }
    
            
        }
            
        return cur_clustering;

    }

    // --- RANDOM NODE WITH NEIGHBORHOOD ORACLE, USING AT MAX CLUSTER SIZE K --- // 


    public static double[] max_k_min_prob_edit_dist_choice_rn(int cur_node, ArrayList<ArrayList<Integer>> cur_clustering, int k, HashMap<Integer, Integer> settled_nodes, ArrayList<ArrayList<Integer>> prob_matrix) {
        // compute clustering choice that minimizes increase of prob edit dist
        int num_clusters = cur_clustering.size();         
        ArrayList<Double[]> avg_weights = new ArrayList<Double[]>();
        ArrayList<Integer> cur_edges = prob_matrix.get(cur_node);
        if (settled_nodes.size() < cur_edges.size()) { // DO IT THE OLD WAY

        for (int i = 0; i < num_clusters; i++){

            if (cur_clustering.get(i).size() < k) { // ONLY CONSIDER CLUSTERS THAT HAVE SPACE TO ADD NODES
                int[] cluster1 = {cur_node};
                double cur_icweight = intercluster_weight_rn(cluster1, cur_clustering.get(i), prob_matrix);
                Double[] cur_val = {cur_icweight, (double) i};
                avg_weights.add(cur_val);
            } else {
                // ADD PLACEHOLDER SO THAT CLUSTER WON'T BE CHOSEN
                Double[] cur_val = {2.0 * k, (double) i};
                avg_weights.add(cur_val);
            }

        }
        
        double[] total_weights = new double[num_clusters];
        for (int i = 0; i < num_clusters; i++) {
            total_weights[i] = avg_weights.get(i)[0] * cur_clustering.get(i).size();
        }
        
        double all_weights = sum(total_weights);
    
        // figure out which cluster decreases total score the least, or if we should create a new cluster
        ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
        for (int i = 0; i < num_clusters; i++) {

            Double[] cur_val = {all_weights + cur_clustering.get(i).size() -  2 * total_weights[i], (double) i};
            cost_increase.add(cur_val);

        }

        Double[] last_val = {all_weights, (double) num_clusters};
        cost_increase.add(last_val); // cost of opening a new cluster
        Double[] least_cost = find_min(cost_increase);

        if (least_cost[1] < num_clusters) {
            double[] r_val = {1, least_cost[1]};
            return r_val;
        }
        else {
            double[] r_val = {0, least_cost[1]};
            return r_val; // force to open new cluster
        }
        } else { // TAKE ADVANTAGE OF NEIGHBORHOOD ORACLE TO REDUCE "QUERIES" 

            // Calculate positive connections to established clusters 
            HashMap<Integer, Integer> pos_connections = new HashMap<Integer, Integer>();
            int total_connections = 0;
            for (int i = 0; i < cur_edges.size(); i++) {
                int cur_edge = cur_edges.get(i);
                if (settled_nodes.keySet().contains(cur_edge)) {
                    int cluster_label = settled_nodes.get(cur_edge);
                    if (pos_connections.keySet().contains(cluster_label)) 
                        pos_connections.put(cluster_label, pos_connections.get(cluster_label) + 1);
                    else
                        pos_connections.put(cluster_label, 1);
                    total_connections += 1;
                }
            }

            // Figure out best clustering choice 
            ArrayList<Double[]> cost_increase = new ArrayList<Double[]>();
            for (int i : pos_connections.keySet()) { // i is the cluster label now              
                if(cur_clustering.get(i).size() < k) {
                    Double[] cur_val = {total_connections + cur_clustering.get(i).size() -  (2.0 * pos_connections.get(i)), (double) i};
                    cost_increase.add(cur_val);
                }
            }

            Double[] last_val = {(double) total_connections, (double) num_clusters};
            cost_increase.add(last_val); // cost of opening a new cluster
            Double[] least_cost = find_min(cost_increase);

            if (least_cost[1] < num_clusters) {
                double[] r_val = {1, least_cost[1]};
                return r_val;
            }
            else {
                double[] r_val = {0, least_cost[1]};
                return r_val; // force to open new cluster
            }

        }
    }

    public static ArrayList<ArrayList<Integer>> max_k_random_node_network(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*
            RandomNode ordering // TODO: use same order as corresponding Pivot clustering?
            Assume 0/1 graph; Make use of Neighborhood oracle provided  
        */
        int num_nodes = prob_matrix.size();
    
        // top-v():
        ArrayList<Integer> ecs_list = new ArrayList<Integer>(num_nodes);
        for (int i = 0; i < num_nodes; i++)
            ecs_list.add(i);
        Collections.shuffle(ecs_list);

        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();
        Integer cur_node = (int) ecs_list.get(0); 
        settled_nodes.put(cur_node, 0); // first node is assigned to cluster 0 

        // clustering is incomplete until algorithm finishes
        ArrayList<ArrayList<Integer> > cur_clustering = new ArrayList<ArrayList<Integer> >();
        ArrayList<Integer> first_cluster = new ArrayList<Integer>();
        first_cluster.add(cur_node);
     
        cur_clustering.add(first_cluster);
        
   
        int node_counter = 1;
        
        while (settled_nodes.size() < num_nodes) {
            cur_node = (int) ecs_list.get(node_counter);
            node_counter += 1;
            
            // get all cluster weights between cur_node and cur_clusters
            // cluster_weights = [[intercluster_weight([cur_node], cur_clustering[i], prob_matrix), i] for i in range(len(cur_clustering))]
            // max_cluster = max(cluster_weights)
            double[] max_cluster = max_k_min_prob_edit_dist_choice_rn(cur_node, cur_clustering, k, settled_nodes, prob_matrix);
    
            if (max_cluster[0] >= 0.5) {
                // add cur_node to the clustering
                int new_index = (int) max_cluster[1];
                cur_clustering.get(new_index).add(cur_node);
                settled_nodes.put(cur_node, new_index);
            }
            else {
                ArrayList<Integer> new_cluster = new ArrayList<Integer>();
                new_cluster.add(cur_node);
                settled_nodes.put(cur_node, cur_clustering.size());
                cur_clustering.add(new_cluster);
            }
    
            
        }
            
        return cur_clustering;

    }

  

    
}
