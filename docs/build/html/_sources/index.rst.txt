.. pyTEnrich documentation master file, created by
   sphinx-quickstart on Wed Dec 30 11:29:02 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: images/CLIFF_logo_design.png

.. toctree::
   :maxdepth: 1
   :caption: Contents: 
   
   usage/installation.rst
   usage/execution.rst
   usage/detailmethods.rst
   usage/analysis_ex.rst
   source/pyTEnrich.rst
   
Overview of the methods
=======================

**CLIFF** integrates cell subtype abundance and expression with bulk drug sensitivity to infer the cell subtypes’ drug susceptibility. CLIFF relies on a new mathematical model using cell-subtype- associated latent variables, which we resolve with an Expectation-Maximization algorithm. 

**Mathematical Modelling of Cell-Type-Specific Drug Sensitivity** 

.. image:: images/cliff_schemes1.png
.. image:: images/cliff_schemes2.png
.. image:: images/cliff_schemes3.png
.. image:: images/cliff_schemes4.png
.. image:: images/cliff_schemes5.png
.. image:: images/cliff_schemes6.png
.. image:: images/cliff_schemes7.png
.. image:: images/cliff_schemes8.png
.. image:: images/cliff_schemes9.png
.. image:: images/cliff_schemes10.png
.. image:: images/cliff_schemes11.png

:ref:`detailmethods`
