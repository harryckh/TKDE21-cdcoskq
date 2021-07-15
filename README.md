# TKDE21-cdcoskq

Notes
=======================

	1. This code was used for the empirical study of the TKDE 2021 paper 
	"Cost-Aware and Distance-Constrained Collective Spatial Keyword Query".
	2. This code is developed by Harry Kai-Ho Chan (kai-ho@ruc.dk).
	3. This code is written in C/C++.
	4. This code runs under Unix/Linux.
	5. In case that you encounter any problems when using this code,
	please figure out the problem by yourself 
	(The code in fact is easy to read and you can modify it for your own purpose).


Usage
=======================

Step 1: Compile the source code

    make

Step 2: Run the code

    ./CDCoSKQ

Step 3: Input Files (Format of those files can be found in Appendix A)
By default, we have file "config.txt" for configuration

Step 4: Collect the querying results and running statistics 
[you can ignore this step if you don't want to collect the information of
querying results and running statistics]

the querying results are stored in "result.txt"
which format is explained in Appendix D.

the running statistics are stored in "stat.txt"
which format is explained in Appendix E.


Appendix A. The format of config.txt
============================

    <Cost function indicator>
    <Weight function indicator>
    <Algorithm indicator> 
    <Method indicator>
    <percentage of B>
    <# of dimensions>
    <# of objects>
    <Location file>
    <# of keywords>
    <Keyword file>
    <Weight option>
    <IR-tree file>
    <# of query keywords>
    <query set size>
    <Percentile lower bound>
    <Percentile upper bound>
    <Random seed>

Explanation of the content in config.txt
-----------------------

    <Cost function indicator>
      =   1: dist_{MaxSum}
      =   2: dist_{Dia}

    <Weight function indicator>
      =   1: cost_{Max}
      =   2: cost_{Sum}

    <Algorithm indicator> 
      = 1: CD-Exact
      = 2: CD-Appro
      = 3: Appro-Adapt
      = 4: Exact-Baseline

    <Method indicator> 
      if Algorithm indicator = 1: = 1 
      if Algorithm indicator = 2: = 2
      if Algorithm indicator = 3: 
        = 1: Caoo-Appro
        = 2: Long-Appro
      if Algorithm indicator = 4: = 1

    <# of dimensions>
      : the number of dimensions of spatial space.

    <# of objects>
      : the number of spatial objects.

    <Location file>
      : the file containing the locations of the spatial objects,
    which format is explained in Appendix II.

    <# of keywords>
      : the total number of all possible keywords contained by the objects.

    <Keyword file>
      : the file containing the keywords of the spatial objects,
    which format is explained in Appendix III.

    <Weight option>
      = 0 : generate weight by a normal distribution
      = 1 : weight of each object in 2nd column of the location file (refer to Appendix II)

    <IR-tree option>
      = 0: the IR-tree has been built (which will be built and stored in <IR-tree file>)
      = 1: the IR-tree has been built (which is stored in <IR-tree file file>)

    <IR-tree file>
      : the file for storing a new (or an existing) IR-tree.

    <# of query keywords>
      : the number of keywords in the query (i.e., |q.\psi|).

    <query set size>
      : the number of queries that will be performed for the same setting, 
      the average statistics based on which will be used.

    <Percentile lower bound>
      : the percentile lower bound that is used for generating queries.

    <Percentile upper bound>
      : the percentile upper bound that is used for generating queries.

    <Random seed>
      : the random seed for generating the queries. Use the same value to make sure same set of queries is used when comparing different algorithms.

    (See file config.txt in the folder for example)

Appendix B. The format of "Location file"
============================
Depending on the "Weight option",
We read two different format for location file,

    (1) "Weight option" = 1 :
    ------------------------
    <object ID1>, <weight>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    <object ID2>, <weight>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    <object ID3>, <weight>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    ...
    <object IDn>, <weight>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    ------------------------

    (2) "Weight option" = 0 :
    ------------------------
    <object ID1>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    <object ID2>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    <object ID3>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    ...
    <object IDn>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
    ------------------------

Note that
	n = # of objects
	m = # of dimensions

(See file running-loc in the folder for example)

Appendix C. The format of "Keyword file"
=============================

    <object ID1>, <1st keyword>, <2nd keyword>, ...
    <object ID2>, <1st keyword>, <2nd keyword>, ...
    <object ID3>, <1st keyword>, <2nd keyword>, ...
    ...
    <object IDn>, <1st keyword>, <2nd keyword>, ...

(See file running-doc in the folder for example)

Appendix D. The format of <result.txt>
=============================

    
    Query #1:
    Keywords: <Keywords of the query>
    Locations: <1st coordinate>	<2nd coordinate> ...

    <# of objects in the solution>

    <object ID1>: (<weight of this objec>) <1st keyword of this object> <2nd keyword of this object> ...
    <object ID2>: (<weight of this objec>) <1st keyword of this object> <2nd keyword of this object> ...
    ...
    <object IDk>: (<weight of this objec>) <1st keyword of this object> <2nd keyword of this object> ...

    Dist: <distance function value of the solution>
    B: <Budget for the query>
    Weight: <Cost function value of the solution>
    Time: <Execution time of the query>

    Query #2:
    same format as for Query #1.

    ...

    Query #t:
    same format as for Query #1

Note that 
	k = the number of objects in the solution
	t = query set size

(See file result.txt in the folder for example)

Appendix E. The format of "stat.txt"
=============================
    <Average weight function value>
    <Average distance function value>

    <the time of building the IR-tree>
    <the average time of performing a query>

    <Memory usage>
    <IR-tree memory usage>

    <Average cost ratio $\alpha$> (only for approximation algorithms)
    <Average distance ratio $\beta$> (only for approximation algorithms)

(See file stat.txt in the folder for example)



