# mTORevaluation
Scripts for mTOR evaluation of the overlap between the mTOR pathway and the output of NLP systems.

Detailed information can be found : http://sbnlp.github.io/mtorevaluation

Material describing the evaluation:

* [Extracting Biological Pathway Models From NLP Event Representations](http://aclweb.org/anthology/W/W16/W16-2916.pdf) - Proceedings of the 15th Workshop on Biomedical Natural Language Processing, Berlin, Germany, 2016, pages 119-127. Association for Computational Linguistics

* [Supplementary Materials](https://sbnlp.github.io/mTOR-evaluation/supp.pdf)

It might be difficult to run the software out of the box. If you have any trouble, don't hesitate to contact us.

### Prerequisites Software

* [libSBML with python bindings](See http://sbml.org/Software/libSBML/docs/python-api/)

* [simstring](http://www.chokkan.org/software/simstring/)

* [fuzzywuzzy](https://pypi.python.org/pypi/fuzzywuzzy)

* [networkx](https://networkx.github.io)

* [numpy](http://www.numpy.org) (only used for some statistics)

* [pandas](http://pandas.pydata.org) (only used for some statistics)


Most of these prerequisites can be installed on Ubuntu for python2.7 using the following commands
* sudo pip2 install libsbml fuzzywuzzy numpy networkx pandas

Installing simstring is a bit tricky. It needs to be compiled from scratch and on Ubuntu (14.04 at least) has to be patched. If you have trouble with simstring, please contact us.

### Prerequisites Data

The scripts relies on data such as SBO data and gene name compilations for simstring. This data is available from [here](https://www.dropbox.com/s/18b6nwg1wxpk2qk/SBNLP-mTorevaluation-sbo-gene_map-gene_list.zip). On unix systems you can try to run
* wget -qO- https://www.dropbox.com/s/18b6nwg1wxpk2qk/SBNLP-mTorevaluation-sbo-gene_map-gene_list.zip | bsdtar -xvf-

In the git directory directly next to these files. Expect to have the following files: sbo.pickle, gene_list.simstring and gene_map.pickle.

The data for the mTOR pathway is not public and cannot be made available through this repository. If you have questions how to retrieve the data, do not hesitate to contact us.

### Scripts

* cd2sbml.py - [CellDesigner](http://www.celldesigner.org) to SBML format conversion. CellDesigner uses a custom SBML format. This script convertes the CellDesigner customizations into standard SBML.
* annotate_sbml_entrez.py - Script for annotating species in an SBML file with Entrez Gene identifiers.
* networkx_analysis.py - Helpers for network overlap analysis, visualization etc
* networkx_analysis_mtor.py - Main script for running network overlap analysis and visualizing results


### Usage
0. (optional) If your NLP output is in standoff format, convert to SBML using https://github.com/sbnlp/mTORevaluation
1. Convert the CellDesigner SBML to pure SBML using cd2sbml.py
2. Annotate all SBML files with Entrez Gene identifiers using annotate_sbml_entrez.py
3. Run network analysis using networkx_analysis_mtor.py

### Tested Platforms

This software has only been tested and should work on *nix systems (tested on Ubuntu 14.04, Mac OSX)

### Contact

michael dot spranger at gmail.com AND sucheendrak at gmail dot com

### Copyright

Copyright [2015-2017] [Michael Spranger, Sucheendra K. Palaniappan, Samik Ghosh]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
