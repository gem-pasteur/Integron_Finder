.. IntegronFinder - Detection of Integron in DNA sequences
   documentation master file, created by
   sphinx-quickstart on Mon Jul 27 15:07:43 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to IntegronFinder's documentation!
==========================================

IntegronFinder is a program that detects integrons in DNA sequences.
The program is available on a webserver :ref:`(Mobyle) <mobyle>`, or by command line (`IntegronFinder on github`_).

You already read the :ref:`paper <references>` and want to install it ? Click :ref:`here <install>`

Integrons are major genetic element, notorious for their major implication in the spread of antibiotic resistance genes. More generally,  integrons are gene-capturing device, whose broader evolutionary role remains poorly understood. IntegronFinder is able to detect with high accuracy integron in DNA sequences. It is accurate because it combines the use of HMM profiles for the detection of the essential protein, the site-specific integron integrase, and the use of Covariance Models for the detection of the recombination site, the *attC* site.

|integron schema|

**How does it work ?**

- First, IntegronFinder annotates the DNA sequence's CDS with Prodigal.

- Second, IntegronFinder detects independently integron integrase and *attC*
  recombination sites. The Integron integrase is detected by using the intersection
  of two HMM profiles:

  - one specific of tyrosine-recombinase (PF00589)
  - one specific of the integron integrase, near the patch III domain of tyrosine recombinases.

The *attC* recombination site is detected with a covariance model (CM), which
models the secondary structure in addition to the few conserved sequence
positions.


- Third, the results are integrated, and IntegronFinder distinguishes 3 types of
  elements:

  - complete integron
      Integron with integron integrase nearby *attC* site(s)
  - In0 element
      Integron integrase only, without any *attC* site nearby
  - CALIN element
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
     :align: middle
     :alt: IntegronFinder Pipeline

.. |integron schema| image:: _static/schema.*
      :align: middle
      :width: 300px
      :alt: Integron Schema


.. _`IntegronFinder on github`: https://github.com/gem-pasteur/Integron_Finder


.. toctree::
   :maxdepth: 2

   installation
   tutorial
   mobyle
   references

..
  Indices and tables
  ==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`
