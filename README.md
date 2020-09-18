# SH-assembly
SH-assembly is a pipeline for de novo genome assembly.
It uses a de Bruijn graph (DBG)-based method.
It counts the k-mers using CQF-deNoise followed by contracting unbranched paths in DBG to unitigs.
After that, the unitig graph is simplified using the graph simplification module in [Minia](https://github.com/GATB/minia).
In the original Minia program, we can not directly call its simplication module.
Here we provide a customized version of Minia (v3 git commit n3eb6f54), using which we can directly call its simplification module.

## External libraries
- **z**    (at least 1.2.3.5)
- **bz2**
- **boost**
- **tbb**

## Installation
    git clone https://github.com/Christina-hshi/SH-assembly.git
    #go to project root directory
    cd SH-assembly
    mkdir release
    cd release
    cmake ..
    #this will make all targets, including a customized version of CQF-deNoise, and assembly program.
    make
    cd ..

    #install a customized version of minia, where the simplification module can be run as an idenpendent step.
    git clone https://github.com/Christina-hshi/Minia.git
    cd Minia
    git submodule init
    git submodule update --remote
    mkdir release    
    cd release
    cmake ..
    make -j8

## User manual
The basic workflow to run SH-assembly is 

1. run **CQF-deNoise** to count kmers.
2. run **Contiger** to build unitig graph. 
3. run **Minia** to simplify unitig graph, and produce contigs.

For each program, run the program without arguments for usage instructions. Following is the detailed manual for each step.
#### step 1: run "CQF-deNoise" to count k-mers
```
./CQF-deNoise  <options>
Options:
  -h [ --help ]            print help messages
  -k arg                   k-mer size
  -n [ --trueKmer ] arg    number of unique true k-mers
  -N arg                   total number of k-mers
  -e [ --alpha ] arg (=-1) average base error rate, when specified, the 
                           <errorProfile> is ignored
  --errorProfile arg       error profile in a file, each line with error rate 
                           for the corresponding base, e.g. the error rate of 
                           the second base is specified in the second line
  --fr arg (=0)            tolerable rate of true k-mers being wrongly removed,
                           default: 1/<trueKmer>
  --deNoise arg (=-1)      number of rounds of deNoise, when specified, the 
                           <fr> is ignored
  --endDeNoise             call deNoise after processing all the k-mers (not 
                           counted into the <deNoise>)
  -t arg (=16)             number of threads
  -f [ --format ] arg      format of the input: g(gzip); b(bzip2); f(plain 
                           fastq)
  -i [ --input ] arg       a file containing a list of read file name(s), should
                           be in the same directory as the fastq file(s)
  -o [ --output ] arg      output file name
```
Users are required to specify either ```<alpha>``` or ```<errorProfile>```. 
If you have no idea about the error rate or the error profile of the sequencing techniques used to generate your data, you can use [ntCard](https://github.com/bcgsc/ntCard.git) to first estimate the k-mer frequency histogram. 
[NtCard](https://github.com/bcgsc/ntCard.git) is quite efficient and it produces a lot of useful statistics, including the total number of k-mers (```<N>```), number of k-mers with different occurrence counts. 
And then by taking k-mers with low counts (e.g. singletons or k-mers with counts <=2, depending on the sequencing depth) as potential false k-mers and the rest of the k-mers as true k-mers, you can estimate an average base error rate by computing ```1 - (T/N)^(1/k)```, where ```T``` is the total number of true k-mers (please notice that ```T``` is not ```<n>```). 
If you don't know the number of unique true k-mers (```<n>```), you can use the number of k-mers with enough occurrence counts (e.g. >=2 or >=3) as an approximation, or you can use the size of the genome as the approximation of the number of unique true k-mers (```<n>```). 

Following is an example of how to set the parameters to run CQF-deNoise based on the k-mer frequency histogram estimated by ntCard to count the 47-mers in a C.elegans data set.

The output of the ntCard (run with '-k47') contains 
```
F1      16506371070 #this is the total number of k-mers <N>
F0      1810841770  #this is the total number of unique k-mers
f1      1665561610  #this is the number of unique k-mers occuring once
f2      26122317    #this is the number of unique k-mers occuring twice
f3      6172811
f<x>      ....
```

Based on the output of the ntCard, we can set 
```
N = F1 = 16506371070
n = F0 - f1 - f2 = 119157843 (which is similar to the size of the C.elegans genome, ~100 Mb)
e = 1 - ((F1 - f1 - 2*f2)/F1)^(1/k) = 0.00234
```
Please notice that in this example, we consider k-mers with counts <= 2 as potential false k-mers, since the sequencing depth of our data is very high (>100x). If you are using a data set of low sequencing depth, we recommmed you consider only singleton k-mers as potential false k-mers to estimate the values of parameters as stated above.

For this example, we can run
```
./CQF-deNoise -k 47 -N 16506371070 -n 119157843 -e 0.00234 -i ReadFiles.txt -o k47.cqf
```

CQF-deNoise can save a lot of space especially when the sequencing depth is high or there are a huge number of false k-mers. When the data set is of low sequencing depth or contains only a few number of false k-mers, we recommend you run CQF-deNoise without removing any k-mers by setting ```--deNoise=0```.

#### step 2: run "Contiger" to build unitig graph
```
./Contiger  <options>
Options::
  -h [ --help ]                         print help messages
  -k arg                                k-mer size
  -i [ --input ] arg                    a file containing a list of read file 
                                        name(s), should be absolute address if 
                                        the files are not in the running 
                                        directory
  -f [ --format ] arg (=f)              format of the input: g(gzip); b(bzip2);
                                        f(plain fastq)
  -c [ --cqf ] arg                      the counting quotient filter built with
                                        the same 'k'
  -s [ --abundance_min ] arg (=2)       minimum coverage of k-mers used to 
                                        extend the assembly
  -x [ --solid_abundance_min ] arg (=2) minimum coverage of a solid k-mer to 
                                        start the assembly
  -X [ --solid_abundance_max ] arg (=1000000)
                                        maximum coverage of a solid k-mer to 
                                        start the assembly
  -t arg (=16)                          number of threads
  -o [ --output ] arg (=unitigs.fa)     output contig file name (fasta)
```

For example, we can run
```
./Contiger -k 47 -i ReadFiles.txt -c k47.cqf -o unitigs.fa
```

#### step 3: run "Minia" to simplify unitig graph and produce contigs
Detailed manual to run Minia can be found [here](https://github.com/GATB/minia.git). In our customized version of Minia, we introduce an extra parameter "-unitig". When specified, it indicates that the input is unitig, and the graph simplification module will be called.

For example, we can run
```
./minia -kmer-size 47 -unitig -in unitigs.fa
```

## Contact
- hshi@cse.cuhk.edu.hk
- kevinyip@cse.cuhk.edu.hk


