.. IntegronFinder - Detection of Integron in DNA sequences

.. _changes:

************
What's new ?
************

.. _changesV2:

V2.0.6
======

- fix issue about pseudo genes https://github.com/gem-pasteur/Integron_Finder/issues/118
- use pyproject.toml for packaging

V2.0.5
======

- fix **major** bug introduce in 2.0.3 in the fix of pandas deprecation warning.
  In some conditions with --local-max with CALIN the number of CALIN is wrong.
- fix issue `Wrong coordinates in results <https://github.com/gem-pasteur/Integron_Finder/issues/114>`_
  If an attC site is find on the very first position and the model is truncated the computation of the position
  was wrong.
- add support of new gembase format

V2.0.3
======

**DO NOT USE** this version (see V2.0.5)

- Improve the support of gembase format (V1 and V2)
- fix deprecation warning, in `biopython`, `pandas` and `python` to be compliant with latest version of these libraries.
- fix `crashes on a multiple contigs file, using the gembase format, in the middle of the contigs being processed. <https://github.com/gem-pasteur/Integron_Finder/issues/103>`_
- fix `IntegronFinder bugs when there is a "space" character in the path to the genome file <https://github.com/gem-pasteur/Integron_Finder/issues/99>`_


In Version 2.0
==============

Here are the major changes between versions 1.x and 2.0.
Essentially, it has be designed such as it becomes easier to find integrons with high confidence in huge datasets
(but it works also for small datasets).

- IntegronFinder now accepts multifasta files as input.
- Only three files are created by default (see output section for details about the other possible output files):

    - A file with all integrons and their elements detected in all sequences in the input file.
    - A summary file with the number and type of integrons per sequence.
    - A file with standard output

- IntegronFinder can be run in parallel easily with a provided Nextflow script that is (almost) ready to use.
- We diversify the installation methods, so it can be easily deployed on a variety of machine. Notably,
  we built a docker container which will allow a smooth installation on clusters (via docker or singularity).
- CALINs are now reported when they have at least 2 *attC* sites (instead of 1 before).
  This value can be changed by the user with `--calin-threshold x`
- Promoters and attI sites are not detected by default to increase speed
- It is now easy to obtain multiple alignments of detected attC sites
- Improve the documentation, especially on the developer part so anyone can contribute.
- Add unit (or non regression) tests.
