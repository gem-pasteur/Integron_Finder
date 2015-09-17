.. IntegronFinder - Detection of Integron in DNA sequences
   documentation master file, created by
   sphinx-quickstart on Mon Jul 27 15:07:43 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to IntegronFinder's documentation!
==========================================

IntegronFinder is a program that detects integrons in DNA sequences.
The program is available on a webserver (link), or by command line
(`IntegronFinder on github`_).

First, IntegronFinder annotates the DNA sequence's CDS with Prodigal.

Second, IntegronFinder detects independently integron integrase and *attC*
recombination sites. The Integron integrase is detected by using the intersection
of two HMM profiles:

- one specific of tyrosine-recombinase (PF00589)
- one specific of the patch III domain of the integron integrase.

The *attC* recombination site is detected with a covariance model (CM) which,
models the secondary structure in addition to the few conserved sequence
positions.

Third, the results are intergrated, and IntegronFinder distinguishes 3 types of
elements:

- complete integron
    Integron with integron integrase nearby *attC* site(s)
- In0 element
    Integron integrase only, without any *attC* site nearby
- attC0 element
    *attC* sites only, without integron integrase nearby.
    A rule of thumb to avoid false positive is to filter out singleton of
    *attC* site.

IntegronFinder can also annotate gene cassettes (CDS nearby *attC* sites) using
Resfams, a database of HMM profiles aiming at annotating antibiotic resistance
genes. This database is provided but the user can add any other HMM profiles
database of its own interest.

When available, IntegronFinder annotates the promoters and attI sites by pattern
matching.

.. image:: _static/pipeline.*
     :width: 400px
     :align: center
     :alt: IntegronFinder Pipeline

.. _`IntegronFinder on github`: https://github.com/gem-pasteur/Integron_Finder


.. toctree::
   :maxdepth: 2

   installation
   tutorial
   mobyle

..
  Indices and tables
  ==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`
