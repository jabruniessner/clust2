### Clust 2: A clustering algorithm for [SDA7](https://mcm.h-its.org/sda/doc/doc_sda7/index.html)

This project was intially written to do complexes clustering for SDA. SDA still contains some dependencies.
When cloned the prjects complilation is as simple as typing 

```
make
```

into the terminal inside this projects folder. The compilation is currently only setup for the gcc compiler. 
Once the program is built, it consists of three parts: The 

```
clust
```
executable is the main part of the program. It can be invoked with the following signature

```
./program_name complexes_file pdb_file num_of_cluster method [ncomplexes]
```

The input paramters have the following meaning

complexes_file: This file is the complexes file as it is output by the SDA simulation. It needs to be in ASCII format, which 
is specified [here](https://mcm.h-its.org/sda/doc/doc_sda7/complexes_file.html). In previous versions, the complexes file needed to be brought into a 
"format suported by the clust program" using the read_record -format option in SDA (for reference: see [here](https://mcm.h-its.org/sda/doc/doc_sda7/tools.html)) 
This step is no longer necessary. The step was originally necessary, because the old clust program was written for an older version of SDA.

pdb_file: The pdb file of solute 2 

method: The method used for clustering. Currently implemented methods are 'a' average, 's' single 'm' maximum and 'c' centroid linkage. For reference see [here](https://en.wikipedia.org/wiki/Hierarchical_clustering)

ncomplexes(optional): The first ncomplexes to be used found in the complexes file (defaults to all)

The program creates two files. The first is the Treecuting.out file. It contains a list of all complexes and the cluster of the num_of_cluster representation clusters to which it belongs.

Second there is the cluster_(method)_linkage.out file, where method is any of the following keywords (average, single, maximum, centroid).
It contains 5 columns. The header is given to be 

```
#Node CycleNo    Item1    Item2   Distance
```

Where Node is the number of clusters after this cycle. CycleNo is the number of the cycle, Item1 and Item2 refer to the clusters that were merged in this cycle.
The numbering is for complexes (clusters with a single member complex) {0, ..., ncomplexes-1} and summarized cluster are numbered {-1, ..., -ncomplexes}. The last 
value is the distance between the two cluster that are merged.


This outputfile can then be used as an input for the script aux/complete_scoring.awk. This is an awk script that extends the output file by 4 columns:
increments, avg. ClSize, max. ClSize and min ClSize.
