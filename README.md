# Simulated Annealing Heuristic for Virtual Network Reallocation 

## Dependencies

The implementation uses boost C++ libraries. Boost can be downloaded from [here](http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.bz2). After downloading to a local directory (say `/usr/local/lib`) extract the compressed archive.
```
tar xvf boost_1_60_0.tar.bz2
```
Then add the boost library path (in this example: `/usr/local/lib/boost_1_60_0`) in the `Makefile` as follows:
```
INCLUDE_PATHS = /usr/local/lib/boost_1_60_0
```

## File Organization
  * io.h: Utility functions for I/O.
  * datastructure.h: Contains necessary data structure definitions.
  * util.h(.cc): Contains several utility functions.
  * graph_util.h: Contains utility functions for Graph manipulation.
  * sa_util.h: Contains utility functions for Simulated Annealing.
  * simulated_annealing.h: Contains functions for Simulated Annealing algorithm.
  * vnr.cc: Contains the main function.
## How to run
```
$ make
$ ./vne_reallocation --case_directory=<case_directory> --max_iterations=<max_iterations> --iterations_per_temperature=<iterations_per_temperature>
```
`--max_iterations` and `iterations_per_temperature` are optional parameters for
setting the maximum number simulated annealing iterations and the number of
iterations to perform for each temperature value, respectively. The default
values are 300 and 150, respectively. 

Each test case directory should have the following format:
```
case5
├── sn.txt
├── optimization_para.txt
└── vnr
    ├── vn0.txt
    ├── vn0.txt.nmap
    ├── vn0.txt.semap
    ├── vn1.txt
    ├── vn1.txt.nmap
    ├── vn1.txt.semap
    ├── vnloc0.txt
    └── vnloc1.txt
```

Description of these files are as follows:
  * sn.txt = Specification of a physical network
  * optimization_para.txt = Parameters for optimization
  * vn$i.txt = Specification of the i-th virtual network request
  * vn$i.txt.nmap = The given node mapping of the i-th VN
  * vn$i.txt.semap = The given link mapping of the i-th VN
  * vnloc$i = Location constraint for the i-th VN
  
## Input file format

A topology file contains the list of edges. Each line contains a description of
an edge in a comma separated value (CSV) format. Format of a line is as follows:
```
<LinkID>,<SourceNodeId>,<DestinationNodeId>,<PeerID>,<Cost>,<Bandwidth>,<Latency>
```
Where,
  * LinkID = Id of a link. Ignored for both physical and virtual topology.
  * SourceNodeId = 0-based node index of the source of the link
  * DestinationNodeId = 0-based node index of the destination of the link
  * PeerID = Ignored
  * Cost = Cost of provisioning unit bandwidth on this link. Cost is ignored for
           virtual links.
  * Bandwidth = Residual bandwidth of a physical link. In case of virtual link,
                this is the bandwidth requirement
  * Delay = Latency of a physical link. In case of virtual link, this is the
            maximum delay requirement for the virtual link. (Not used)

A location constraint file contains as many lines as the number of virtual
nodes. Each line is a comma separated list of values. The first value indicates
the id of a virtual node followed by the ids of physical nodes where this
virtual node can be mapped.

The node mapping file for a VN contains as many lines as the number of virtual
nodes in the corresponding virtual network. Each line has the following format:
```
<virtual_node_id> <mapped_physical_node_id>
```

The link mapping file contains multiple lines, where each line corresponds to a
physical link and a virtual link mapped onto that physical link. Each line is
formatted as follows:
```
<plink_endpoint_0> <plink_endpoint_1> <mapped_vlink_endpoint0> <mapped_vlink_endpoint_1>
```

The optimization_para.txt file contains the following lines:
```
Goal Utilizatino = x%
alpha = <alpha>
beta = <beta>
```
Goal utilization is the utilization threshold for determining if a physical link is bottleneck or not. alpha and beta are the weights of bandwidth cost and bottleneck link cost in the objective function, respectively.

Please see the provided files as an example for better clarification.

Note: Nodes are numberded from `0 ... (n - 1)` in a network with `n` nodes.

## Output Files

Currently the solver prints output to the standard output and writes them to
the following output files inside `vnr` directory: 

* prev_cost = Cost of embedding before reoptimization.
* prev_bnecks = Number of bottleneck links before reoptimization.
* prev_bw_cost = Bandwidth cost before reoptimization.
* prev_max_plink_util = Maximum physical link utilization before reoptimization.
* sol_time = Execution time (in seconds)
* new_cost = Cost of the reoptimized embedding according to the cost function.
* new_bnecks = Number of bottleneck links after reoptimization
* new_bw_cost = Bandwidth cost after reoptimization
* new_max_plink_util = Maximum plink utilization after reoptimization
* vn*.node_remap = node mapping of that VN in the reoptimized embedding
* vn*.edge_remap = edge mapping of that VN in the reoptimized embedding
