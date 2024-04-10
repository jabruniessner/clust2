### Clust 2: A clustering algorithm for [SDA7](https://mcm.h-its.org/sda/doc/doc_sda7/index.html)

This project was initially written to do complexes clustering for SDA. The program is based on the 
[C clustering library](http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm).
When cloned the project's compilation is as simple as typing.

```
make
```

into the terminal inside this project's folder. The compilation is currently only setup for the gcc compiler. 
Once the program is built, it consists of three parts: The 

```
clust
```
executable is the main part of the program. It can be invoked with the following signature

```
./program_name complexes_file pdb_file num_of_cluster method [ncomplexes]
```

The input parameters have the following meaning

complexes_file: This file is the complexes file as it is output by the SDA simulation. It needs to be in ASCII format, which 
is specified [here](https://mcm.h-its.org/sda/doc/doc_sda7/complexes_file.html). In previous versions, the complexes file needed to be brought into a 
"format supported by the clust program" using the read_record -format option in SDA (for reference: see [here](https://mcm.h-its.org/sda/doc/doc_sda7/tools.html)) 
This step is no longer necessary. The step was originally necessary because the old clust program was written for an older version of SDA.

pdb_file: The pdb file of solute 2 

method: The method used for clustering. Currently implemented methods are 'a' average, 's' single 'm' maximum and 'c' centroid linkage. For reference see [here](https://en.wikipedia.org/wiki/Hierarchical_clustering)

ncomplexes(optional): The first ncomplexes to be used are found in the complexes file (defaults to all)

The program creates the cluster_(method)_linkage.out file, where the method is any of the following keywords (average, single, maximum, centroid).
It contains 5 columns. The header is given to be 

```
#Node CycleNo    Item1    Item2   Distance
```

Where Node is the number of clusters after this cycle. CycleNo is the number of the cycle, Item1 and Item2 refer to the clusters merged in this cycle.
The numbering is for complexes (clusters with a single member complex) {0, ..., ncomplexes-1} and summarized clusters are numbered {-1, ..., -ncomplexes}. The last 
value is the distance between the two clusters that are merged (which is also the minimal distance between any two clusters in the previous cycle).


This output file can then be used as an input for the script aux/complete_scoring.awk. This is an awk script that extends the output file by 4 columns:
increments, avg. ClSize, max. ClSize and min ClSize. increments is the percentage of the increase in distance between the first cycle and the last cycle in percent. 
I.e.: Let $d_i$ be the distance in the $i^{th}$ cycle, then the increment is given by

$$increment_i = \frac{d_i-d_{i-1}}{d_n-d_1}.$$

avg. ClSize is the average size of clusters in that cycle, which is given by the number of complexes divided by the number of clusters in that cycle.
max. ClSize is the maximum size of all clusters and min. ClSize is the minimum size of all clusters. The header after this step looks like this

```
#Node CycleNo    Item1    Item2   Distance  increments   avg. ClSize   max. ClSize   min. ClSize
```

The output can be saved to a file with the `>` in bash. I.e.

```
path/to/aux/complete_scoring.awk path/to/clustering_(method)_linkage.out > path/to/output_file.out
```

Beware that this script is dependent on numpy. Versions, with which it has been tested are numpy 1.24.3<=version<=1.26.4.
Its documentation can be found by either callig the script without arguments or by calling the script with the "--help" or "-h" flag in the argument.

In the last step, it is possible to plot the increments and the distances against the number of clusters. This can quickly be done with gnuplot and the aux/plot_cycles.gnu script.
It can be invoked with the command 

```
gnuplot -c path/to/aux/plot_cycles.gnu filein fileout "title name" lineno 
```


Where filein should be replaced by the output of the previous command `path/to/output_file.out`, fileout is the filename where you would like to save the plot(will be saved in pdf format),
"title name" is the title for the figure and lineno is the number of lines containing data. (Or the maximum number of clusters you had during the clustering.)

Remark: For a detailed introduction to gnuplot see [gnuplot website](http://www.gnuplot.info/) or this [book](https://alogus.com/g5script/gnuplot5/) or this [book](https://github.com/tianxiao/gnuplot-study/blob/master/gnuplot/Gnuplot%20in%20Action.pdf).

A typical plot from this script looks like this


![plot](https://github.com/jabruniessner/clust2/assets/95258260/30c8e879-83ca-47b7-a544-685eacf2cca2)


Based on this information it is possibel to choose how many clusters can be considered distinct. Based on this, it is typical to perform a tree cut.
The cut prgram can be invoked via

```
path/to/cut [-old] cluster_(method)_linkage.out complexes_file pdb_file num_of_cluster 
```
This program has three output files.

The first is the Treecuting.out file. It contains a list of all complexes and the cluster of the num_of_cluster representation clusters to which it belongs.

Second, there is the cluster_out.ascii it contains the same header as the complexes file and the entries of the representatives of the clusters. This file can then be used to create
corresponding PDB files, for instance with the script generateFortComplexesPdbs-SDA7.py, which is part of the SDA7 package (found in the aux directory of SDA). 

The [] around -old indicate that this flag is optional. If the -old option is specified the clusters in are ordered by ClSize, otherwise they are ordered by ClFSize. 
Ordering them by ClSize is an unadvisable thing to do, but was added for backwards compatibility.

With this data, it is possible to write the members of the clusters to files. This can be done with the 

``` 
path/to/aux/print_cluster_members.py
```

An example can be found in the example directory

Have fun using!
