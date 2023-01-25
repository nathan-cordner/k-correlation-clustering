import java.util.*;
import java.io.*;
import java.util.concurrent.ThreadLocalRandom;

public class PKwik {

    // --- UNCONSTRAINED PIVOT ---

    // LARGE NETWORK: READ AS ADJACENCY LIST
    public static ArrayList<ArrayList<Integer>> pKwikClustering(ArrayList<ArrayList<Integer>> prob_matrix) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        // int loop1 = 0;
        // int loop2 = 0;
        // long start = 0;

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
      

            // start = System.currentTimeMillis();
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                }
            }
            // loop1 += System.currentTimeMillis() - start;

            // start = System.currentTimeMillis();
            available_nodes.removeAll(cur_cluster);
            // for (int i = 1; i < cur_cluster.size(); i++) 
            //     available_nodes.remove(cur_cluster.get(i));
            // loop2 += System.currentTimeMillis() - start;
            clusters.add(cur_cluster);            
        }
        // System.out.println("Loop 1 finished in: " + loop1 / 1000.0 + " s");
        // System.out.println("Loop 2 finished in: " + loop2 / 1000.0 + " s");
            
        return clusters;
    }

    // INPUT: probability matrix
    public static ArrayList<ArrayList<Integer>> pKwikClusteringProb(ArrayList<ArrayList<Double[]>> prob_matrix) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            for (Double[] values : prob_matrix.get(cur_node)) {
                double temp = values[0];
                int neighbor = (int) temp;
                double prob = values[1];

                if (prob >= 0.5 && available_nodes.contains(neighbor))
                    cur_cluster.add(neighbor);
            }

            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);            
        }
            
        return clusters;
    }

    // --- AT MOST K CLUSTERS ---


    // LARGE NETWORK: READ AS ADJACENCY LIST
    public static ArrayList<ArrayList<Integer>> k_pivotBlend(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*
            "BLEND METHOD": form k pivot clusters, then proceed to RN

            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
            settled_nodes.put(cur_node, clusters.size());
      
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                    settled_nodes.put(neighbor, clusters.size());
                }
            }

            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);
            
            if (clusters.size() == k)
                break;
        }

        PriorityQueue<Pair> cluster_sizes = new PriorityQueue<Pair>();
        if (settled_nodes.size() < num_nodes) {
           // Build heap for node algorithm
           for (int i = 0; i < clusters.size(); i++) {
                Pair cur_cluster_info = new Pair(clusters.get(i).size(), i);
                cluster_sizes.add(cur_cluster_info);
           }
        }

        // finish adding nodes by RandomNode method 

        while (settled_nodes.size() < num_nodes) {
            while (settled_nodes.containsKey(permutation.get(cur_index)))
                cur_index++;
            Integer cur_node = (int) permutation.get(cur_index);
            
            double[] max_cluster = DNode.min_prob_edit_dist_choice_rn(cur_node, clusters, k, cluster_sizes, settled_nodes, prob_matrix);
    
            int new_index = (int) max_cluster[1];
            Pair cluster_info = new Pair(clusters.get(new_index).size(), new_index);
            cluster_sizes.remove(cluster_info);
            // System.out.println(worked);

            clusters.get(new_index).add(cur_node);
            settled_nodes.put(cur_node, new_index);

            cluster_info.incrementKey();
            cluster_sizes.add(cluster_info);
            
        }
            
        return clusters;
    }

    public static ArrayList<ArrayList<Integer>> k_pivot(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*

            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node

            Method: merge new pivot clusters once k are formed (alway add new cluster to current smallest cluster)
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();

        int cur_index = 0;

        // Priority Queue to track pivot cluster sizes once k clusters are formed
        PriorityQueue<Pair> cluster_sizes = new PriorityQueue<Pair>();

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
            settled_nodes.put(cur_node, clusters.size());
      

            // start = System.currentTimeMillis();
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                    settled_nodes.put(neighbor, clusters.size());
                }
            }
            // loop1 += System.currentTimeMillis() - start;

            // start = System.currentTimeMillis();
            available_nodes.removeAll(cur_cluster);
            // for (int i = 1; i < cur_cluster.size(); i++) 
            //     available_nodes.remove(cur_cluster.get(i));
            // loop2 += System.currentTimeMillis() - start;
            if (clusters.size() < k) {
                clusters.add(cur_cluster);
                if (clusters.size() == k) {
                    // build heap
                    for (int i = 0; i < clusters.size(); i++) {
                        Pair cur_cluster_info = new Pair(clusters.get(i).size(), i);
                        cluster_sizes.add(cur_cluster_info);
                    }
                } 
            }
            else {           

                // add to smallest
                Pair cur_min = cluster_sizes.poll();
                int smallest_index = cur_min.getValue();
                clusters.get(smallest_index).addAll(cur_cluster);
                cur_min.updateKey(clusters.get(smallest_index).size());
                cluster_sizes.add(cur_min);

            }
        }

            
        return clusters;
    }

    public static ArrayList<ArrayList<Integer>> exactly_k_pivot(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*

            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node

            Method: merge new pivot clusters once k are formed (alway add new cluster to current smallest cluster)
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();

        int cur_index = 0;

        // Priority Queue to track pivot cluster sizes once k clusters are formed
        PriorityQueue<Pair> cluster_sizes = new PriorityQueue<Pair>();

        while (!available_nodes.isEmpty()) {

            if (clusters.size() + available_nodes.size() == k) { // put remaining nodes into singletons
                for (int my_node : available_nodes) {
                    ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                    cur_cluster.add(my_node);
                    clusters.add(cur_cluster);
                }
                return clusters;
            }

            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
            settled_nodes.put(cur_node, clusters.size());
      

            // start = System.currentTimeMillis();
            int added_nodes = 1;
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                    settled_nodes.put(neighbor, clusters.size());
                    added_nodes += 1;

                    if (clusters.size() + 1 + (available_nodes.size() - added_nodes) == k) {
                        clusters.add(cur_cluster);
                        available_nodes.removeAll(cur_cluster);
                        // add remaining nodes into singletons

                        for (int my_node : available_nodes) {
                            cur_cluster = new ArrayList<Integer>();
                            cur_cluster.add(my_node);
                            clusters.add(cur_cluster);
                        }
                        return clusters;

                    }

                }


            }
            // loop1 += System.currentTimeMillis() - start;

            // start = System.currentTimeMillis();
            available_nodes.removeAll(cur_cluster);
            // for (int i = 1; i < cur_cluster.size(); i++) 
            //     available_nodes.remove(cur_cluster.get(i));
            // loop2 += System.currentTimeMillis() - start;
            if (clusters.size() < k) {
                clusters.add(cur_cluster);
                if (clusters.size() == k) {
                    // build heap
                    for (int i = 0; i < clusters.size(); i++) {
                        Pair cur_cluster_info = new Pair(clusters.get(i).size(), i);
                        cluster_sizes.add(cur_cluster_info);
                    }
                } 
            }
            else {           

                // add to smallest
                Pair cur_min = cluster_sizes.poll();
                int smallest_index = cur_min.getValue();
                clusters.get(smallest_index).addAll(cur_cluster);
                cur_min.updateKey(clusters.get(smallest_index).size());
                cluster_sizes.add(cur_min);

            }
        }

            
        return clusters;
    }


    public static ArrayList<ArrayList<Integer>> k_pivot_fill(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*

            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node

            Method: add pivot clusters to 0 to k-1, repeating in order
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();

        int cur_index = 0;
        int cluster_count = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
            settled_nodes.put(cur_node, clusters.size());
      
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                    settled_nodes.put(neighbor, clusters.size());
                }
            }

            available_nodes.removeAll(cur_cluster);

            if (clusters.size() < k) 
                clusters.add(cur_cluster);
            else {
                clusters.get(cluster_count).addAll(cur_cluster);
                cluster_count++;
                cluster_count = cluster_count % k;

            }           

        }

            
        return clusters;
    }

    public static ArrayList<ArrayList<Integer>> k_pivot_rand(ArrayList<ArrayList<Integer>> prob_matrix, int k) {
        /*

            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node

            Method: add new cluster to randomly chosen existing cluster
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();
        HashMap<Integer, Integer> settled_nodes = new HashMap<Integer, Integer>();

        int cur_index = 0;
        int cluster_count = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
            settled_nodes.put(cur_node, clusters.size());
      
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                    settled_nodes.put(neighbor, clusters.size());
                }
            }

            available_nodes.removeAll(cur_cluster);

            if (clusters.size() < k) 
                clusters.add(cur_cluster);
            else // choose random existing cluster to add to
                clusters.get(ThreadLocalRandom.current().nextInt(0, k)).addAll(cur_cluster);         

        }

            
        return clusters;
    }


    // --- MAX CLUSTER SIZE K ---

    // LARGE NETWORK: READ AS ADJACENCY LIST
    public static ArrayList<ArrayList<Integer>> maxKPivot(ArrayList<ArrayList<Integer>> prob_matrix, int max_size) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            ArrayList<Integer> neighbors = prob_matrix.get(cur_node);
            Collections.shuffle(neighbors); // add neighbors in random order up to max size

            for (int neighbor : neighbors) {
                if (available_nodes.contains(neighbor) && cur_cluster.size() < max_size) {
                    cur_cluster.add(neighbor);
                }
            }

            available_nodes.removeAll(cur_cluster);

            clusters.add(cur_cluster);            
        }
            
        return clusters;
    }


    
}
