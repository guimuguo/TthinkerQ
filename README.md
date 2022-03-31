# A General-Purpose Framework for Answering Subgraph Queries in a Big Graph

Given a big graph G, a subgraph query Q finds a set of subgraphs of G that satisfy certain querying conditions. Applications of subgraph queries include finding a dense community containing designated members to organize an event, and querying a knowledge graph by subgraph matching. 

Inspired by the recently emerged think-like-a-task (TLAT) parallel paradigm for solving compute-intensive problems offline, we proposes a novel programming framework, called T-thinkerQ, for answering online subgraph queries in parallel following the TLAT paradigm. The design of T-thinkerQ inherits the benefits of prior TLAT systems such as high task concurrency, bounded memory consumption, and straggler elimination by a timeout strategy for load balancing. 

In addition, T-thinkerQ utilizes a novel active task-queue list to ensure the fairness that queries are answered in the received order: a later query is processed only if earlier queries cannot saturate the available computation resources. To track query progress so that users are timely notified when a query completes, T-thinkerQ also adopts a novel lineage-based design that keeps track of how subtasks are generated by straggler tasks for divide-and-conquer processing. We use four kinds of subgraph queries to demonstrate the programming friendliness of T-thinkerQ as well as its excellent CPU-scalability.



## Program Checklist
- **The `system` folder:** it contains the code for our T-thinkerQ engine, which is a task-based general-purpose framework for writing parallel subgraph query programs. In the folder, `worker.h` is the main thread that creates other computing threads (aka. compers) to work on tasks. 

- **The `app_kernel_ol` folder:** given a input vertex set S = {v<sub>1</sub>, v<sub>2</sub>, ...}, parameter γ and minsize, this is the application code for mining maximal γ-quasi-cliques, which include all vertices un S and has size > minsize.


- **The `maximal_check` folder:** This is the postprocessing step, used to remove non-maximal quasi-cliques from the output of `app_kernel_ol`.

## Compilation
In each folder, `app_kernel_ol` and `maximal_check`, there is a Makefile. Just enter each folder and use the command `make` to compile, and a program named `run` will be generated.

## Execution
**Workflow A: to Mine Maximal Quasi-Cliques **
  1. Quasi-clique mining:
 
      Run the program in the `app_kernel_ol` folder: `app_kernel_ol/run [input_data] [thread_num] [out-gamma] [in-gamma] [min_size] [time_split_threshold]`, where: 
        - input_data: input graph file where the *i*-th row records the adjacency list of Vertex *i*
        - thread_num: number of threads. We also call each computing thread a comper
        - out-gamma: user-specified minimum outdegree-ratio threshold
        - in-gamma: user-specified minimum indegree-ratio threshold
        - min_size: minimum size threshold; each returned result should have at least so many vertices
        - time_split_threshold: timeout duration threshold. A task running longer than the threshold will decompose into subtasks 

        Example: `app_qc/run input_graph 5 0.8 0.9 10 5`

  2. Postprocessing:
      - Each thread (Comper *i*) will write the results it finds to a file `output_i`
      - Aggregate all quasi-cliques outputs into one file: `cat output_* > results`
      - Remove non-maximal quasi-cliques: `maximal_check/quasiCliques results max_results`


## Demo
Click [here](https://colab.research.google.com/drive/1Cn0cB9uZ8uOtlPbAfTWw9g0NM9qBrkxC?usp=sharing) for a demo on Google Colab. The notebook first clones the repo and download the [Google-Web](https://snap.stanford.edu/data/web-Google.html) dataset. It then runs the quasi-clique mining program to find maximal results. Finally, it plots the first and second largest quasi-cliques.

## Requirements

* C++11
* [OpenMP](https://www.openmp.org/)

<!-- ## Contributors
* **Guimu Guo (guimuguo@uab.edu)**
* **Da Yan (yanda@uab.edu)**
* **Lyuheng Yuan (lyuan@uab.edu)**
* **Saugat Adhikari (saugat@uab.edu)**

The authors are affiliated with the Department of Computer Science, University of Alabama at Birmingham -->
