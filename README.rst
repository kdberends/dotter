======
Dotter
======
|Build|
|Cover|
|Lic|
|Doc|
|Binder|


This repository contains a set of eco-hydraulic tools to analyse (small) streams. It is aimed at (regional) water managers to aid in getting insight, decision making and visualisation.

The dottertools are distributed as a Python package.

**UNDER DEVELOPMENT**
This package is currently under active development and has very little in the way of testing. You may well have problems running it. Use at your own risk. 

Getting started
===========

The easiest way to get started is to try out the interactive notebooks on binder. 

You can also visit our documentation. 


Installation
===========


The current development branch of dotter can be installed from GitHub using ``pip``:

::

    pip install git+https://github.com/kdberends/dotter


Dependencies
============

Dotter is tested on Python 3.6 and depends on NumPy,
SciPy, Pandas, and Matplotlib (see ``requirements.txt`` for version
information).

Additionally, some notebooks make use of interactive widgets. To be able to correctly run the notebooks, you need to have ipywidgets installed and enabled; see `here <https://ipywidgets.readthedocs.io/en/stable/user_install.html>`__

To compile the sphinx documentation you need `pandoc <http://pandoc.org/>`_ installed as well. 

License
=======

`GPL License, Version
3.0 <https://github.com/kdberends/dotter/blob/master/LICENSE.txt>`__


References
====

This project has been set up using PyScaffold 2.5.8. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.



.. |Cover| image:: https://img.shields.io/coveralls/github/kdberends/dotter.svg
   :target: https://coveralls.io/github/kdberends/dotter?branch=master

.. |Build| image:: https://img.shields.io/travis/kdberends/dotter.svg
    :target: https://travis-ci.org/kdberends/dotter

.. |Lic| image:: https://img.shields.io/github/license/kdberends/dotter.svg
   :target: https://github.com/kdberends/dotter/blob/master/LICENSE.txt

.. |Doc| image:: https://img.shields.io/readthedocs/dotter.svg
   :target: http://dotter.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |Binder| image:: https://mybinder.org/badge.svg 
   :target: https://mybinder.org/v2/gh/kdberends/dotter/master?filepath=examples/notebooks
