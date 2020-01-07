Introduction
============

The **SAGE**-analysis package was developed to handle the data products
produced by the **SAGE** model, available `here`_.

.. _here: https://github.com/sage-home/sage-model

Why Did We Create A New Repo?
-----------------------------

**SAGE** is an extremely modular and flexible semi-analytic model. The base
model (presented in `Croton et al., 2016`_) has been adjusted and altered to
answer a number of science questions on a variety of topics including galactic HI and angular
momentum properties in `DARK SAGE`_, the Epoch of Reionization in `RSAGE`_, and
the star formation history and galactic dust in `DUSTY SAGE` _.

Due to the wide array of science that **SAGE** can cover and the number of
models that spawned from its development, there has been a need to develop a
framework to ingest in data from (ideally) any **SAGE** variant. This repo
represents such an effort.  It represents a series of modules intended to
provide the easy ingestion and analysis of the base **SAGE** data, whilst
serving as a template for analysing any other **SAGE** flavours.

.. _Croton et al., 2016: https://arxiv.org/abs/1601.04709
.. _DARK SAGE: https://arxiv.org/abs/1605.00647
.. _RSAGE: https://arxiv.org/abs/1902.01611
.. _DUSTY SAGE: https://arxiv.org/abs/1902.01611

Advantages of the Package
-------------------------

* Easy analysis and plotting of multiple different **SAGE** models.  For
  example, comparing **SAGE** models with/without supernova feedback.
* Memory efficient analysis of **SAGE** output. All calculations are performed
  using only a single output file at a time, ensuring no extra memory overhead
  associated with opening many files.
* Support for the user to implement their own functions to analysis + plotting.
* Template for creating custom data classes to ingest any arbitrary **SAGE**
  data output. Useful if you're looking to develop using **SAGE**.
