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

Registration on the Galaxy server of the Pasteur institute is not required to use the tool.
Yet, if you wish to keep your history, we recommend you to register.

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

Once the job is finished, you get your results on right panel.
All files contain the log of the run which tells you how many integrons have been found for each types
along with the number of *attC* sites per type. There are 4 different outputs created:

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


For each of the aforementioned files, you can save them by clicking on the download button.

.. _`Galaxy server of the Pasteur institute`: https://galaxy.pasteur.fr/tool_runner?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F2.0%2Bgalaxy0


