# DeepMetaPSICOV 1.0
### Deep residual neural networks for protein contact prediction

### NB This repo is still under development and is not yet ready for production use. We estimate that it will be ready in a few weeks once all testing is complete.

Shaun M. Kandathil, Joe G. Greener and David T. Jones

University College London

Requirements:
-------------
- Bash shell
- Perl5 (tested on 5.22.1 and 5.16.3)
- C and C++ compilers (tested with GCC 4.4.2, 4.8.5, and 5.4.0)
- Python 2 or 3 (preferably miniconda/anaconda, as this makes the PyTorch install much easier)
- The following Python modules:
  - PyTorch 0.3.1 
  
- Third-party programs:
  - HH-suite v3.0+ and a recent UniClust30 database (for making alignments; skip if you will only use pre-made alignments)
  - CCMpred v0.1.0 (github commit id XXXX)
  - FreeContact 1.0.21
  - Legacy PSI-BLAST 2.2.26 (executable `blastpgp`) and a suitable non-redundant database, e.g. Uniref90, formatted using `formatdb` (needed to generate PSIPRED and SOLVPRED inputs)

All other required programs written by our group are now bundled in this repo and do not need to be installed separately.

On some distributions, the C++ compiler is a separate add-on package and may not be installed by default. For example, on CentOS you will need to `yum install` packages `gcc` AND `gcc-c++`.

### Installing PyTorch using conda
Use `conda install -c pytorch pytorch=0.3.1`

### GPU support
If you are installing PyTorch using `conda`, `conda` should automatically detect a usable GPU and install an appropriate version with GPU support.
NB a GPU is not necessary for predicting contacts, but for faster runtimes we recommend using one.

Setup and testing:
------------------

TODO

Running:
--------

TODO

Citing:
-------
If you find DeepMetaPSICOV useful, please cite our paper (details to be added).
