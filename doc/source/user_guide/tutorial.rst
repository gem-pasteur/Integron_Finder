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

The time limiting part are HMMER and INFERNAL.
So if you have to analyze one or few replicon the user can set the number of CPU used by HMMER and INFERNAL::

  integron_finder mysequence.fst --cpu 4

Default is 1.


If you want to deal with lot of replicons (from 10 to more than thoushand) we provide a workflow to parallelize
the execution by the data. This mean that we cut the data input in chunks (by default of one replicon) then execute
integron_finder in parallel on each replicon (the number of parallel tasks can be limitted) then agregate the results
in one golbal summary.
The workflow use the `nextflow <https://www.nextflow.io/>`_ framework and can be run on a single machine or a cluster.

So to run parallel_integron_finder.nf you have to install nextflow first, and  :ref:`integron_finder <install>`.
Then we provide 2 files parallel_integron_finder.nf which is the workflow itself in nextflow syntax and nextflow.profile.config
which is a configuration file to execute the workflow.
The wokflow file should not be modified. Whereas the profile must be adapted to the local architecture.

The file nextflow.profile.config provide two profile a standard and cluster profile

.. warning::

    On Ubuntu Bionic Beaver (18.04) The default java is not suitable to run nextflow
    So you have to install another jvm

        sudo add-apt-repository ppa:webupd8team/java
        sudo apt-get update
        sudo apt-get install oracle-java8-installer

    for more details see: https://medium.com/coderscorner/installing-oracle-java-8-in-ubuntu-16-10-845507b13343

    so now install nextflow.
    If you have  capsule error like ::

        CAPSULE EXCEPTION: Error resolving dependencies. while processing attribute Allow-Snapshots: false (for stack trace, run with -Dcapsule.log=verbose)
        Unable to initialize nextflow environment

    install nextflow as follow (change the nextflow version with the last release) ::

        wget -O nextflow http://www.nextflow.io/releases/v0.30.2/nextflow-0.30.2-all
        chmod 777 nextflow

    for more details see: https://github.com/nextflow-io/nextflow/issues/770#issuecomment-400384617


standard profile
""""""""""""""""

You can specify the number of tasks in parallel by setting the queueSize value ::

    standard {
            executor {
                name = 'local'
                queueSize = 7
            }
            process{
                container = 'Singularity/integron_finder.simg'
                executor = 'local'
                $integron_finder{
                    errorStrategy = 'ignore'
                    cpu=params.cpu
                }
            }
     }

This profile is to work with the integron_finder singularity image.
If you don't want to use it but you prefer to use an installed version you can remove the `container` line as following. ::

    standard {
            executor {
                name = 'local'
                queueSize = 7
            }
            process{
                executor = 'local'
                $integron_finder{
                    errorStrategy = 'ignore'
                    cpu=params.cpu
                }
            }
     }


All options available in non parallel version are also available for the parallel one.
A typical command line will be: ::

    ./parallel_integron_finder.nf -profile standard --replicons all_coli.fst --circ  --out E_Coli_all

.. note::
    the option starting with one dash are for nextflow, whereas the options starting by two dashes are for integron_finder

if you execute this line 2 directory will be created one named `work` containing lot of subdirectories this for all jobs
launch by nextflow and a directory named `Results_Integron_Finder_E_Coli_all` which contain the final results:

    * integron_report.html
    * integron_timeline.html
    * integron_trace.txt
    * Results_Integron_Finder_E_Coli_all

:integron_report.html: is an HTML execution report: a single document which includes many useful metrics about
    a workflow execution. For further details see https://www.nextflow.io/docs/latest/tracing.html#execution-report

:integron_timeline.html: is an HTML timeline for all processes executed in your pipeline.
    For further details see https://www.nextflow.io/docs/latest/tracing.html#timeline-report

:integron_trace.txt: creates an execution tracing file that contains some useful information about
    each process executed in your pipeline script, including: submission time, start time, completion time,
    cpu and memory used. For further details see https://www.nextflow.io/docs/latest/tracing.html#trace-report

:Results_Integron_Finder_E_Coli_all: contains the actual results as in non parallel version.


cluster profile
"""""""""""""""

The cluster profile is intented to work on a cluster managed by SLURM.
If You cluster is managed by an other drm change executor name by the right value
(see `nextflow supported cluster <https://www.nextflow.io/docs/latest/executor.html>`_ )

You can also managed

* The number of task in parallel with the `executor.queueSize` parameter (here 500).
  If you remove this line, the system will send in parallel as many jobs as there are replicons in your data set.
* The queue with `process.queue` parameter (here common,dedicated)
* and some options specific to your cluster management systems with `process.clusterOptions` parameter ::


    cluster {
        executor {
            name = 'slurm'
            queueSize = 500
        }

        process{
            container = 'Singularity/integron_finder.simg'
            executor = 'slurm'
            queue= 'common,dedicated'
            clusterOptions = '--qos=fast'
            $integron_finder{
                cpu=params.cpu
            }
        }
    }

    singularity{
        enabled = true
        runOptions = '-B /pasteur'
        autoMounts = false
    }



The profile above is intended to work with singularity.
If you want to work with an installed version of `integron_finder`
remove the `process.container` line and the `singularity block`. ::

   cluster {
        executor {
            name = 'slurm'
            queueSize = 500
        }

        process{
            executor = 'slurm'
            queue= 'common,dedicated'
            clusterOptions = '--qos=fast'
            $integron_finder{
                cpu=params.cpu
            }
        }
    }


To run the parallel version on cluster, for instance on a cluster managed by slurm,
I can launch the main nextflow process in one slot. The parallelization and the submission on the other slots
is made by nextflow itself so the command line look like: ::

    sbatch --qos fast -p common nextflow run  parallel_integron_finder.nf -profile cluster --replicons all_coli.fst --cpu 2 --local-max --gbk --circ --out E_Coli_all


The results will be the same as describe in local execution.




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
