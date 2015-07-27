.. IntegronFinder - Detection of Integron in DNA sequences

.. _quickstart:

**********
Quickstart
**********

We assume here that the program is :ref:`installed <install>`.

Basic use
============

You can go to directory containing your sequence, or specify the path to that
sequence and call::

    integron_finder mysequence.fst

or::

    integron_finder path/to/mysequence.fst

It will perform a rapid search, and outputs the results in a directory called
``Results_Integron_Finder_mysequence``. Within this directory, you can find:

- mysequence.integrons
   A tabular file with the annotations of the different element
- mysequence.gbk
   A GenBank file with the sequence annotated with the same annotations from
   the previous file.
- mysequence_X.pdf
   For each complete integron, a simple graphic of the region is depicted

.. _eagle_eyes:

Thorough local detection
------------------------

This option allows a more sensitive search. It will be slower if integrons are
found, but will be as fast if nothing is detected::

    integron_finder mysequence.fst --eagle_eyes

.. _func_annot:

Function annotation
-------------------

This option allows to annotate genes given HMM profiles. As Resfams database is
distributed, to annotate antibiotic resistance genes, just use::

    integron_finder mysequence.fst --func_annot

IntegronFinder will look in the directory
``Integron_Finder-x.x/data/Functional_annotation`` and use all ``.hmm`` files
available to annotate. By default, there is only ``Resfams.hmm``, but one can
add any other HMM file here. Alternatively, if one wants to use a database which
is present elsewhere on the user's computer without copying it into that
directory, one can specify the following option::

    integron_finder mysequence.fst --dir_func_annot path/to/DB/HMM/

where ``path/to/DB/HMM/`` contains one or many ``.hmm`` (or ``.HMM``) files.

.. _parallel:

Parallelization
---------------

The time limiting part are HMMER and INFERNAL. So IntegronFinder has parallel
implementation, but the user can set the number of CPU used by HMMER and
INFERNAL::

  integron_finder mysequence.fst --cpu 4

Default is 1.

To start IntegronFinder on many nucleotide sequences, one can use "manual"
parallelization by calling multiple times IntegronFinder in ``bash``.

.. _advance:

Advanced use
============

.. _distance_threshold:

Clustering of elements
----------------------
