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
- mysequence_X.pdf
   For each complete integron, a simple graphic of the region is depicted
- other
   A folder containing outputs of the different step in the program. It includes
   notably the protein file in fasta (mysequence.prt).

.. _local_max:

Thorough local detection
------------------------

This option allows a more sensitive search. It will be slower if integrons are
found, but will be as fast if nothing is detected::

    integron_finder mysequence.fst --local_max

.. _func_annot:

Functional annotation
---------------------

This option allows to annotate cassettes given HMM profiles. As Resfams database
is distributed, to annotate antibiotic resistance genes, just use::

    integron_finder mysequence.fst --func_annot

IntegronFinder will look in the directory
``Integron_Finder-x.x/data/Functional_annotation`` and use all ``.hmm`` files
available to annotate. By default, there is only ``Resfams.hmm``, but one can
add any other HMM file here. Alternatively, if one wants to use a database which
is present elsewhere on the user's computer without copying it into that
directory, one can specify the following option::

    integron_finder mysequence.fst --path_func_annot bank_hmm

where ``bank_hmm`` is a file containing one absolute path to a hmm file per
line, and you can comment out a line::

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

.. _advance:

Advanced use
============

.. _distance_threshold:

Clustering of elements
----------------------

*attC* sites are clustered together if they are on the same strand and if they
are less than 4 kb apart. To cluster an array of *attC* sites and an integron
integrase, they also must be less than 4 kb apart. This value has been
empirically estimated and is consistent with previous observations showing that
biggest gene cassettes are about 2 kb long. This value of 4 kb can be modify
though::

    integron_finder mysequence.fst --distance_thresh 10000

or, equivalently::

    integron_finder mysequence.fst -dt 10000

This sets the threshold for clustering to 10 kb.

*attC* evalue
-------------

The default evalue is 1. Sometimes, degenerated *attC* sites can have a evalue
above 1 and one may want to increase this value to have a better sensitivity,
to the cost of a much higher false positive rate.

::

    integron_finder mysequence.fst --evalue_attc 5

Circularity
-----------

By default, IntegronFinder assumes replicon to be circular. However, if they
aren't, or if it's PCR fragments or contigs, you can specify that it's a linear
fragment::

    integron_finder mylinearsequence.fst --linear

However, if ``--linear`` is not used and the replicon is smaller than ``4 x dt``
(where ``dt`` is the distance threshold, so 4kb by default), the replicon is
considered linear to avoid clustering problem

Palindromes
-----------

*attC* sites are more or less palindromic sequences, and sometimes, a single
*attC* site can be detected on the 2 strands. By default, the one with the
highest evalue is discarded, but you can choose to keep them with the following
option::

    integron_finder mysequence.fst --keep_palindromes
