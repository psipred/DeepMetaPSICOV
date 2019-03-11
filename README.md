# DeepMetaPSICOV 1.0
### Deep residual neural networks for protein contact prediction

Shaun M. Kandathil and David T. Jones

University College London

Requirements:
-------------
- Bash shell
- Perl5 (tested on 5.22.X and 5.16.3)
- C and C++ compilers (tested with GCC 4.8.5, 5.4.0)
- Python 2 or 3 (preferably miniconda/anaconda, as this makes the PyTorch install much easier)
- The following Python modules:
  - PyTorch 0.3.0+
  
- Third-party programs:
  - HH-suite v3.0+ and a recent UniClust30 database (for making alignments)
  - CCMpred v0.1.0 (github commit id XXXX)
  - FreeContact 1.0.21
  - Legacy PSI-BLAST 2.2.26 (executable blastpgp) and a suitable non-redundant database, e.g. Uniref90
  

On some distributions, the C++ compiler is a separate add-on package and may not be installed by default. For example, on CentOS you will need to `yum install` packages `gcc` AND `gcc-c++`.

### GPU support
If you are installing PyTorch from conda, conda should automatically detect a usable GPU and install an appropriate version with GPU support.
Although not necessary for predicting contacts, for faster runtimes we recommend using a GPU.

Setup and testing:
------------------



Running:
--------


Citing:
-------
If you find DeepMetaPSICOV useful, please cite our paper (details to be added).
