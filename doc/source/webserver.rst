.. IntegronFinder - Detection of Integron in DNA sequences

.. _webserver:

##########
web server
##########

.. _galaxy:


Galaxy
******

You can access IntegronFinder online, on the `Galaxy server of the Pasteur institute`_


How to use it
=============

Registration on the Galaxy server of the Pasteur institute is not required to use the tool. Yet, if you wish to keep your history, we recommend you to register.

1. Upload your sequence with **Get Data - Upload File** in the menu on the left 
2. Select your file in the **Replicon file** list of Integron Finder
3. Select the options you want
4. Click on **Execute**

If you want more options:

3. Select **Show** on advanced parameters
4. Select the options you want
5. Click on **Execute**

You can see the role of the different functions in the :ref:`tutorial <tutorial>` page.

Results
=======

Once the job is finished, you get your results on right panel. All files contain the log of the run which tells you how many integrons have been found for each types along with the number of *attC* sites per type. There are 4 different outputs created:

- **Raw results archive**:
    An archive containing all raw results.
- **Integrons annotations**:
    A tabular file listing all the elements and their caracteristics.
- **GenBank**:
    The GenBank file of the input sequence with the annotation corresponding to
    the elements found (integrase, *attC*, promoter, attI, etc...).
- **Graphics**:
    Simple representation of one or more complete integrons found.
    The representation is very basic and a better representation can be
    obtained from the GenBank file and a software (eg Geneious) to represent it.


For each of the aforementionned files, you can save them by clicking on the download button.

.. _`Galaxy server of the Pasteur institute`: https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F1.5.1

******
Mobyle
******

You can access IntegronFinder online, on the `Mobyle server of the Pasteur institute`_


How to use it
=============

1. Copy your sequence or upload it in the appropriate field.
2. Select the options you want
3. Click on Run

If you want more options:

3. Click on advanced options (instead of Run)
4. Select the options you want
5. Click on Run

You can see the role of the different functions in the :ref:`tutorial <tutorial>` page,
or by clicking on the |red question mark| in the corresponding field.

After submitting your job, you may need to enter your email.

Results
=======

Once the job is finished, you have a result page, which contains:

- integron_finder.out:
    Log of the run. It tells you how many integrons have been found for each types along with the number of *attC* sites per type.

- **Schema of complete integron(s)** : replicon_X.pdf
    Simple representation of one or more complete integrons found.
    The representation is very basic and a better representation can be
    obtained from the GenBank file and a software (eg Geneious) to represent it.
- **annotated sequence** : replicon.gbk
    The GenBank file of the input sequence with the annotation corresponding to
    the elements found (integrase, *attC*, promoter, attI, etc...).
- **putative integrons** : replicon.integrons
    A tabular file listing all the elements and their caracteristics.

Finally, you have your initial sequence of the replicon and the command line used.

For each of the aforementionned files, you can save them by clicking on the save
button |savebutton|.



.. _`Mobyle server of the Pasteur institute`: http://mobyle.pasteur.fr/cgi-bin/portal.py#forms::integron_finder
.. |red question mark| image:: _static/questionmark.png
   :height: 13
   :width: 13
.. |savebutton| image:: _static/mobyle_save.png