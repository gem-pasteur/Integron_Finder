.. IntegronFinder - Detection of Integron in DNA sequences

.. _changes:

************
What's new ?
************

.. _changesV2:

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
