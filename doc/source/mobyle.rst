.. IntegronFinder - Detection of Integron in DNA sequences

.. _mobyle:

******
Mobyle
******

You can access IntegronFinder online, on the `Mobyle server of the Pasteur institute`_

.. note::
  This links to the development version of Mobyle, only accessible from Pasteur.
  If you're outside Pasteur, please provide me your IP so you can access it.

How to use it
=============

1. Copy your sequence or upload it in the appropriate field.
2. Select the option you want
3. Click on Run

If you want more options:

3. Click on advanced options
4. Select the options you want
5. Click on Run

You can see the role of the different functions in the :ref:`tutorial <tutorial>` page,
or by clicking on the |red question mark| in the corresponding field.

After submitting your job, you might need to enter your email.

Results
=======

Once the job is finished, you have a results page. It contains:

- integron_finder.out:
    Log of the run. It tells you how many integrons have been found of each types along with the number of *attC* sites.

- le pdf : replicon_X.pdf:
    Simple representation of one or more complete integrons found.
    The representation is very basic and a better representation can be
    obtained from the GenBank file and a software (eg Geneious) to represent it.
- replicon.gbk:
    The GenBank file of the input sequence with the annotation corresponding to
    the elements found (integrase, *attC*, promoter, attI, etc...)
- replicon.integrons:
    A tabular file listing all the elements and their caracteristics.

Finally, you have the sequence of the replicon you input and the command line used.

For each of the aforementionned file, you can save them by clicking on the save
button |savebutton|.



.. _`Mobyle server of the Pasteur institute`: http://mobyle-dev.pasteur.fr/cgi-bin/portal.py#forms::integron_finder
.. |red question mark| image:: _static/questionmark.png
   :height: 13
   :width: 13
.. |savebutton| image:: _static/mobyle_save.png
