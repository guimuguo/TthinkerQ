# A General-Purpose Framework for Answering Subgraph Queries in a Big Graph

Given a big graph G, a subgraph query Q finds a set of subgraphs of G that satisfy certain querying conditions. Applications of subgraph queries include finding a dense community containing designated members to organize an event, and querying a knowledge graph by subgraph matching. 

Inspired by the recently emerged think-like-a-task (TLAT) parallel paradigm for solving compute-intensive problems offline, we proposes a novel programming framework, called T-thinkerQ, for answering online subgraph queries in parallel following the TLAT paradigm. The design of T-thinkerQ inherits the benefits of prior TLAT systems such as high task concurrency, bounded memory consumption, and straggler elimination by a timeout strategy for load balancing. 



## Program Checklist
- **The `system` folder:** it contains the code for our T-thinkerQ engine, which is a task-based general-purpose framework for writing parallel subgraph query programs. In the folder, `worker.h` is the main thread that creates other computing threads (aka. compers) to work on tasks. 

In T-thinkerQ system, queries are submitted by client programs into an inter-process message queue *Q<sub>ipc<sub>*, and the worker of the server program periodically probes *Q<sub>ipc<sub>* and appends its queries to a query-waiting queue *Q<sub>wait<sub>* where the compers may then fetch queries for processing when they become available. T-thinkerQ utilizes a novel active task-queue list to ensure the fairness that queries are answered in the received order: a later query is processed only if earlier queries cannot saturate the available computation resources. 

<p align="center">
  <!-- <img src="imgs/System_Architecture.png"/> -->
  <img align="center" src="https://github.com/guimuguo/TthinkerQ/blob/main/img/System_Architecture.png"/>
</p>

To track query progress so that users are timely notified when a query completes, T-thinkerQ also adopts a novel lineage-based design that keeps track of how subtasks are generated by straggler tasks for divide-and-conquer processing. It associates each task object *t* with a task progress object *t.prog* of type *task_prog*. When a task *t* (and its progress object) is created, the child-task counter is initialized as 0. Whenever a child task is created by *t*, the counter is incremented. When *t* finishes its call of *compute(.)*, we set the complete-flag of its progress object. Also, when a child-task finishes, the counter is decremented. When the counter is decremented to 0 and complete-flag is also set, we delete *t*’s progress object and remove its entry from *T<sub>prog</sub>*, and update the progress object of its parent task in *T<sub>prog</sub>* by decrementing its child-task counter; if the counter becomes 0, the same process propagates to the parent of this parent task, and this propagation may continue all the way up to the root progress object in which case q is considered to have finished its current phase.

<p align="center">
  <!-- <img src="imgs/Lineage_Tracking.png"/> -->
  <img align="center" src="https://github.com/guimuguo/TthinkerQ/blob/main/img/Lineage_Tracking.png" />
</p>

### We use four kinds of subgraph queries to demonstrate the programming friendliness of T-thinkerQ as well as its excellent CPU-scalability.

- **The `app_kernel_ol` folder:** Given a minimum degree threshold *γ* ∈ [0, 1], a *γ*-quasi-clique is a subgraph *g* = (*V<sub>g<sub>* , *E<sub>g<sub>*) where each vertex *v* connect to at least *γ* fraction of the other vertices in *g*. Given a input vertex set S = {v<sub>1</sub>, v<sub>2</sub>, ...}, parameter γ and minsize, this is the application code for mining maximal γ-quasi-cliques, which include all vertices in the vertex set S and has size > minsize.

- **The `maximal_check` folder:** This is the postprocessing step, used to remove non-maximal quasi-cliques from the output of `app_kernel_ol`.

<p align="center">
  <!-- <img src="imgs/qc.png"/> -->
  <img align="center" src="https://github.com/guimuguo/TthinkerQ/blob/main/img/qc.png" width="350" height="250"/>
</p>
  
- **The `app_scs` folder:** This query finds a subgraph with the largest min-degree among all connected subgraph *g* that contain the query vertex *v<sub>q<sub>* and have *l*≤|*v<sub>q<sub>*|≤h.

<p align="center">
  <!-- <img src="imgs/scs.png"/> -->
  <img align="center" src="https://github.com/guimuguo/TthinkerQ/blob/main/img/scs.png" width="450" height="300"/>
</p>
  
- **The `app_hpcycle` folder:** Given two distinct vertices *s* and *t* in *G*, and a hop constraint *k*, this query outputs all paths from *s* to *t* with length at most *k*.

<p align="center">
  <!-- <img src="imgs/hpcycle.png"/> -->
  <img align="center" src="https://github.com/guimuguo/TthinkerQ/blob/main/img/hpcycle.png" width="350" height="200"/>
</p>
  
  
- **The `app_gmatch` folder:** This application sketches Ullmann’s recursive algorithm for subgraph matching, which, given a query graph *G<sub>q<sub>*, retrieves all subgraphs of a data graph *G* that are isomorphic to *G<sub>q<sub>*.

<p align="center">
  <!-- <img src="imgs/gm.png"/> -->
  <img align="center" src="https://github.com/guimuguo/TthinkerQ/blob/main/img/gm.png" width="450" height="200"/>
</p>
 
## Compilation
In each folder, `app_kernel_ol`, `maximal_check`, `app_scs`, `app_hpcycle`, `app_gmatch`, `client/batchFile_version` `client/console_version`, there is a Makefile. Just enter each folder and use the command `make` to compile, and a program named `run` will be generated.

## Execution
** Maximal Quasi-Clique Query**
  1. Start the Quasi-clique mining server:
 
      Run the program in the `app_kernel_ol` folder: `app_kernel_ol/run [input_data] [thread_num] [time_split_threshold] [basic_min_size] [basic_gamma] [bigtask_threshold]`, where: 
        - input_data: input graph file where the *i*-th row records the adjacency list of Vertex *i*
        - thread_num: number of threads. We also call each computing thread a comper
        - time_split_threshold: timeout duration threshold. A task running longer than the threshold will decompose into subtasks 
        - basic_min_size: basic minimum size threshold for all queries; each returned result should have at least so many vertices. All the following query's min_size should be no less than this basic_min_size
        - basic_gamma: basic user-specified minimum degree-ratio threshold for all queries. All the following query's gamma should be no less than this basic_gamma
        - bigtask_threshold: bigtask differentiated threshold. For Maximal Quasi-Clique Query, tasks whose candidate set's size larger than this threshold will be regarded as bigTask

        Example: `app_kernel_ol/run input_graph 32 1 30 0.85 200`

  2. Queries submission:
    
  (1)console_version. To submit the query, users run the program in the `client/console_version` folder: `client/console_version/run`. Then type your query in the prompt: `[min_size] [gamma] [Vertex1_ID] [Vertex2_ID] ...`
      
  Example: `client/console_version/run 31 0.9 0 1`
  
  (2)batchFile_version. To submit a batch of queries once, users need to run the program in the `client/batchFile_version` folder: `client/batchFile_version/run`. User need put all the queries in the batch_file. The run the porgram: `client/batchFile_version/run batch_file`
        
  3. Postprocessing:
      - The output for query_i will be in `ol_out\query_i_QueryString`
      - Each thread (Comper *i*) will write the results it finds to a file named *i*
      - Aggregate all quasi-cliques outputs into one file: `cat ol_out\query<sub>i<sub>_QueryString\* > results`
      - Remove non-maximal quasi-cliques: `maximal_check/quasiCliques results max_results`


## Demo
Click [here](https://www.youtube.com/watch?v=7doA--qe11U) for a video demo. It then runs the quasi-clique mining program to find all maximal results including the query's vertices.
  
## Requirements

* C++11
* [OpenMP](https://www.openmp.org/)

<!-- ## Contributors
* **Guimu Guo (guimuguo@uab.edu)**
* **Da Yan (yanda@uab.edu)**
* **Lyuheng Yuan (lyuan@uab.edu)**
* **Saugat Adhikari (saugat@uab.edu)**

The authors are affiliated with the Department of Computer Science, University of Alabama at Birmingham -->
