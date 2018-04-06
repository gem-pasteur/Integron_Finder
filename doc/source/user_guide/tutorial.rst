.. IntegronFinder - Detection of Integron in DNA sequences

.. _tutorial:

********
Tutorial
********

We assume here that the program is :ref:`installed <install>`.

Basic use
=========
.. note::
   The different options will be shown separately, but they can be used
   alltogether unless otherwise stated.

You can see all available options with::

    integron_finder -h

You can go to directory containing your sequence, or specify the path to that
sequence and call::

    integron_finder mysequence.fst

or::

    integron_finder path/to/mysequence.fst

It will perform a search, and outputs the results in a directory called
``Results_Integron_Finder_mysequence``. Within this directory, you can find:

- mysequence.integrons
   A tabular file with the annotations of the different elements
- mysequence.gbk
   A GenBank file with the sequence annotated with the same annotations from
   the previous file.
   generated on;y if ``--gbk`` option is set.
- mysequence_X.pdf
   For each complete integron, a simple graphic of the region is depicted
   generated only if ``--pdf`` option is set
- other
   A folder containing outputs of the different step in the program. It includes
   notably the protein file in fasta (mysequence.prt).
   This directory is available only if option ``--keep-tmp`` is set.

.. _local_max:

Thorough local detection
------------------------

This option allows a more sensitive search. It will be slower if integrons are
found, but will be as fast if nothing is detected.

.. code-block:: bash

    integron_finder mysequence.fst --local-max

.. _func_annot:

Functional annotation
---------------------

This option allows to annotate cassettes given HMM profiles. As Resfams database
is distributed, to annotate antibiotic resistance genes, just use::

    integron_finder mysequence.fst --func-annot

IntegronFinder will look in the directory
``Integron_Finder-x.x/data/Functional_annotation`` and use all ``.hmm`` files
available to annotate. By default, there is only ``Resfams.hmm``, but one can
add any other HMM file here. Alternatively, if one wants to use a database which
is present elsewhere on the user's computer without copying it into that
directory, one can specify the following option ::

    integron_finder mysequence.fst --path_func_annot bank_hmm

where ``bank_hmm`` is a file containing one absolute path to a hmm file per
line, and you can comment out a line ::

  ~/Downloads/Integron_Finder-x.x/data/Functional_annotation/Resfams.hmm
  ~/Documents/Data/Pfam-A.hmm
  # ~/Documents/Data/Pfam-B.hmm

Here, annotation will be made using Pfam-A et Resfams, but not Pfam-B. If a
protein is hit by 2 different profiles, the one with the best e-value will be kept.

.. _parallel:

Parallelization
---------------

The time limiting part are HMMER and INFERNAL. So IntegronFinder does not have
parallel implementation (yet?), but the user can set the number of CPU used by HMMER and
INFERNAL::

  integron_finder mysequence.fst --cpu 4

Default is 1.

Topology
--------

By default, IntegronFinder assumes that

    * your replicon is considered as **circular** if there is **only one replicon** in the input file.
    * your replicons are considered as **linear** if there are **several replicons** in the input file.

However, you can change this default behavior and specify the default topology with options
``--circ`` or ``--lin``::

    integron_finder --lin mylinearsequence.fst
    integron_finder --circ mycircularsequence.fst


If you have multiple replicon in the input file with different topologies you can specify a topology for each
replicon by providing a topology file.
The syntax for the topology file is simple:

    * one topology by line
    * one line start by the seqid followed by 'circ' or 'lin' for circular or linear topologies.

example::

    seq_id_1 circ
    seq_id_2 lin

You can also mix the options ``--circ`` or ``--lin`` with option ``--topology-file``::

    integron_finder --circ --topology-file path/to/topofile mysequences.fst

In the example above the default topology is set to *circular*.
The replicons specified in topofile supersede the default topology.


.. warning::
    However, if the replicon is smaller than ``4 x dt``
    (where ``dt`` is the distance threshold, so 4kb by default), the replicon is considered linear
    to avoid clustering problem.
    The topology used to searching integron is report in the *\*.integrons file*


.. _advance:

Advanced options
================

.. _distance_threshold:

Clustering of elements
----------------------

*attC* sites are clustered together if they are on the same strand and if they
are less than 4 kb apart. To cluster an array of *attC* sites and an integron
integrase, they also must be less than 4 kb apart. This value has been
empirically estimated and is consistent with previous observations showing that
biggest gene cassettes are about 2 kb long. This value of 4 kb can be modify
though::

    integron_finder mysequence.fst --distance-thresh 10000

or, equivalently::

    integron_finder mysequence.fst -dt 10000

This sets the threshold for clustering to 10 kb.

.. note::
    The option ``--outdir`` allows you to chose the location of the Results folder (``Results_Integron_Finder_mysequence``).
    If this folder already exists, IntegronFinder will not re-run analyses already done, except functional annotation.
    It allows you to re-run rapidly IntegronFinder with a different ``--distance-thresh`` value.
    Functional annotation needs to re-run each time because depending on the aggregation parameters,
    the proteins associated with an integron might change.


*attC* evalue
-------------

The default evalue is 1. Sometimes, degenerated *attC* sites can have a evalue
above 1 and one may want to increase this value to have a better sensitivity,
to the cost of a much higher false positive rate.

::

    integron_finder mysequence.fst --evalue-attc 5

Palindromes
-----------

*attC* sites are more or less palindromic sequences, and sometimes, a single
*attC* site can be detected on the 2 strands. By default, the one with the
highest evalue is discarded, but you can choose to keep them with the following
option::

    integron_finder mysequence.fst --keep-palindromes

Keep intermediate results
-------------------------

Integrons finder needs some intermediate results, It includes notably the protein file in fasta (mysequence.prt).
A folder containing these outputs is generated for each replicon and have name ``other_<replicon_id>``
This directory is remove at the end. You can keep this directory to see analyse each ``integron_finder`` steps
with the option ``--keep-tmp``.


Verbosity of outputs
--------------------

You can control the verbosity of the outputs with the options ``-v`` or ``-q`` to
respectively increase or decrease the verbosity.
These options are cumulative ``-vv`` or ``-qqq``.
