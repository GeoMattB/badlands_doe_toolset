# Badlands-DoE-Toolset
I built this toolset to undertake the analysis fo the results from the Badlands models I built for my Masters at the University of Sydney.

A set of tools that combines the setup and analysis of Badlands using Design of Experiments (DoEgen) python software.
* Template to setup and run the DoEgen design of experiments (DoE) python software that will generate the design for multiple experiments.
* A utility to take the experiment design outputs from DoEgen and build the equivalent (multiple) Badlands experiment configuration files
* Scripts to run the multiple badlands models simulataneously using python multiprocessing module.
* A set of post-processing and evaluation tools to compare and evaluate the multiple models that result from the experiments.
* Addtional utilities to allow further investigation of Badlands models.

see 
https://badlands.readthedocs.io/
and
https://github.com/sebhaan/DoEgen


This repository is also available on pip, to install:

pip install badlands-doe-toolset


An example on the functionality and description of how to use this is available here:

https://github.com/GeoMattB/Badlands_DoE_synrift


The DoEgen component of this research was supported by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney.

It is an extension/implementaton of :

https://github.com/sebhaan/DoEgen

https://github.com/Sydney-Informatics-Hub/DoEgen-Geo
 
