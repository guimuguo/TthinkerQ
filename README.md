# Maximal Directed Quasi-Clique Mining in a Large Graph

This project develops a parallel solution for mining maximal quasi-cliques on a large directed graph such as web graphs and social networks.
We generalize the concept of quasi-cliques to directed graphs by proposing (γ<sub>1</sub>, γ<sub>2</sub>)-quasi-cliques which have density requirements in both inbound and outbound directions of each vertex in a quasi-clique. An efficient recursive algorithm is proposed to and maximal (γ<sub>1</sub>, γ<sub>2</sub>)-quasi-cliques which integrates many effective pruning rules that are validated by ablation studies. Our algorithm is also fully compatible with the recent parallel compute-intensive graph mining paradigm G-thinker, and can thus scale up and scale out effectively.

We target a **single-machine multi-core** environment, since a distributed cluster is not always readily available to an average end user while most modern computers are multi-core. Our system achieved a close to ideal speedup ratio as illustrated in the figure below for the [Bitcoin](http://konect.cc/networks/soc-sign-bitcoinotc/) dataset.

<p align="center">
  <!-- <img src="imgs/img3.PNG" width="450" height="300" /> -->
  <img align="center" src="https://github.com/guimuguo/Tthinker_DQC/blob/main/imgs/img3.PNG" width="450" height="300" />
</p>




## Program Checklist
- **The `system` folder:** it contains the code for our T-thinker engine, which is a task-based general-purpose framework for writing parallel programs. In the folder, `worker.h` is the main thread that creates other computing threads (aka. compers) to work on tasks. When task queues are near empty, T-thinker will generate new tasks from data items to refill the queues; while if too many tasks are created (e.g. due to decomposing a big task), tasks will be spilled to local disk to keep memory bounded, and these tasks will be loaded back for processing when task queues have space. The figure below shows the tuned system parameters for our task queues:

<p align="center">
<img src="https://github.com/guimuguo/Tthinker_DQC/blob/main/imgs/img2.PNG" width="450" height="300" />
</p>

- **The `app_qc` folder:** this is the application code for mining maximal directed quasi-cliques, which runs on top of T-thinker. The figure below shows an example of the second largest quasi-clique for the [Google-Web](https://snap.stanford.edu/data/web-Google.html) dataset (directed network of hyperlinks between web pages) found by our application code.

<p align="center">
<img src="https://github.com/guimuguo/Tthinker_DQC/blob/main/imgs/img1.PNG" width="450" height="300" />
</p>


- **The `maximal_check` folder:** This is the postprocessing step, used to remove non-maximal quasi-cliques from the output of `app_qc`.

## Compilation
In each folder, `app_qc` and `maximal_check`, there is a Makefile. Just enter each folder and use the command `make` to compile, and a program named `run` will be generated.

## Execution
**Workflow A: to Mine Maximal Quasi-Cliques Directly**
  1. Quasi-clique mining:
 
      Run the program in the `app_qc` folder: `app_qc/run [input_data] [thread_num] [out-gamma] [in-gamma] [min_size] [time_split_threshold]`, where: 
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

Click [here](https://colab.research.google.com/drive/1V4LOIcPLcaiwQUkGZPGEgUt4tcXa3et3?usp=sharing) for another demo on Google Colab. The notebook first clones the repo and download the [Citation Graph](https://lfs.aminer.cn/lab-datasets/soinf/citation-raw.txt) dataset. It then runs the quasi-clique mining program to find maximal results for both (0.75, 0.75)-Quasi-Clique and (0.6, 0.4)-Quasi-Clique to show versatility of our method.

## Requirements

* C++11
* [OpenMP](https://www.openmp.org/)

<!-- ## Contributors
* **Guimu Guo (guimuguo@uab.edu)**
* **Da Yan (yanda@uab.edu)**
* **Lyuheng Yuan (lyuan@uab.edu)**
* **Jalal Khalil (jalalk@uab.edu)**

The authors are affiliated with the Department of Computer Science, University of Alabama at Birmingham -->
