# SH-assembly
SH-assembly is a bioinformatic pipeline for de novo genome assembly.
SH-assembly is a de Bruijn-graph based method.
It counts the k-mers using CQF-deNoise followed by contracting unbranched paths in DBG graph to unitigs.
The unitig graph then will be simplified using the graph simplification module in [Minia 3](https://github.com/GATB/minia).
In the original Minia 3 program, we can not directly call its simplication module.
Here we provide the modification upon Minia 3 (v3 git commit n3eb6f54), using which we can directly call the simplification module.

## External libraries
- z
- bz2
- boost
- tbb

## Installation
    git clone https://github.com/Christina-hshi/SH-assembly.git
    #in project root directory, 
    mkdir release
    cd release
    cmake ..
    #this will make all targets, including a customed version of CQF-deNoise, and assembly program.
    make 

## User mannual
The basic workflow of running SH-assembly is 
    1. run "CQF-deNoise" to counting kmers
    2. run "Contiger" to building unitig graph 
    3. run "minia" to simplify unitig graph, and produce a list of contigs.
For each program, run the program without arguments for usage instructions. 

## Contact
- hshi@cse.cuhk.edu.hk


