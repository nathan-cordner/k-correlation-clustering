Data Format:
* Data files are stored in Data/[name of data set]/graph.txt
  * All data sets are either publicly available or provided here
  * We provide data set samples to show the format for the various experiments 
* Delimiter is hard-coded in driver files
* For K Correlation Clustering: 
  * The first line of the file must contain the total number of nodes (anything that follows it on the line will be ignored)
  * Rest of file lists positive edges as [node1] [node2]

To Compile: javac *.java
To Run: java [DriverName] [data set folder name]
* additional heap space may be needed for some experiments; increase the max heap size with the -Xmx flag

Drivers
-------

Hard coded parameters:
* delimiter for data set (use "," for Facebook data sets and "\\s" for all others)
* number of Pivot rounds, number of attributes used, etc. 

RunKApprox.java
* Input: positive edge list
* Method: Il'ev and Navrotskaya K-approx algorithm

RunKLocalSearch.java
* Input: positive edge list
* Method: Initialize random K clustering and perform local search

RunPTAS.java
* Input: positive edge list
* Method: Giotis and Guruswami MaxAgree and MinDisagree PTAS
* Additional command line argument: sample_size

RunPTAS2.java
* Input: positive edge list
* Method: Karpinski and Schudy PTAS
* Additional command line argument: sample_size

RunKCorrelation.java
* Input: positive edge list
* Method: "neighborhood oracle" for K-Pivot, Blend, and K-Vote

RunKBlend.java, RunKPivotLS.java, RunKVoteLS.java
* Input: positive edge list
* Method: Blend, Pivot, and Vote with additional LS rounds (e.g. 1, 5, 10, *)

RunTimedBlend.java, RunTimedLS.java, RunTimedPivot.java, RunTimedVote.java
* Input: positive edge list
* Method: Blend, LS, Pivot, and Vote with additional LS operations until time constraint reached

Code Files
----------

DNode.java
* Implementations of Vote

Helper.java
* Helper functions for reading data sets

Pair.java
* Custom data structure used for heap implementation

PKwik.java
* Pivot algorithm implementations

KApprox.java
* K-approx algorithm implementation

PTAS.java
* Giotis and Guruswami PTAS implementation

NewPTAS.java
* Karpinski and Schudy PTAS implementation

KLocalSearch.java
* LocalSearch implementations


Plots
-----

k_cc_plots.ipynb
-- Jupyter Notebook for data plots
