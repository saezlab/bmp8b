Python module `csim`
====================
A module for calculation various similarity metrics between disjunct clusterings.

This module implements all methods collected by Wagner & Wagner here:
http://www.cs.ucsb.edu/~veronika/MAE/wagner07comparingclusterings.pdf

For description and comparison of the methods and their properties see the paper above.

Installation
------------

.. code:: bash
    
    pip install git+https://bitbucket.org/deeenes/csim.git


Requirements
------------

Requires numpy and scipy, and optionally igraph.
Should work under Linux, other Unix and Windows, both in Python 2 and 3.

Quick start
-----------

.. code:: python
        
        import csim

        # instantiate the main class:
        cs = csim.ClusteringSet()

        # read clusterings data from multiple files:
        cs.add_file_list(['clustering1.csv', 'clustering2.csv'], 
                         param = {'header': True, 'sep': ';', 'cols': [1]})

        # obtain a similarity matrix using the igraph implementation of
        # Danon adjusted mutual information:
        sim = cs.danon_mi_igraph()

Authors
-------
Dénes Türei -- denes@ebi.ac.uk (feedback, bug reports)

License
-------
GPLv3
