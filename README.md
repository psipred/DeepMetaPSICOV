# DeepMetaPSICOV 1.0
### Deep residual neural networks for protein contact prediction

Shaun M. Kandathil, Joe G. Greener and David T. Jones

University College London

Requirements:
-------------
- Bash shell
- C and C++ compilers (tested with GCC 4.4.2, 4.8.5, and 5.4.0)
- Python 2 or 3 (preferably miniconda/anaconda, as this makes the PyTorch install much easier)
- The following Python modules:
  - PyTorch 0.3.1 
  
- Third-party programs:
  - HH-suite v3.0+ and a recent UniClust30 database (for making alignments; skip if you will only use pre-made alignments)
  - CCMpred v0.1.0
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

### Build bundled programs
`cd src/; make; make install`

### Specify paths to external dependencies
Edit `run_DMP.sh` to indicate the paths to the third-party programs listed above, as well as other variables such as the number of threads to use for various programs. User-editable variables are demarcated by comment lines. We do not recommend changing anything outside this region unless you know what you are doing.

Testing TODO

Running:
--------
Run `/path/to/run_DMP.sh -h` to see the available options. DMP runs a number of programs to generate input features; their outputs are stored in a number of intermediate files. By default, DMP will attempt to reuse any files with the correct filenames (this is useful for debugging and allows you to 'continue' a failed run). You can force regeneration of intermediate files with the `--force` option.

At a minimum, you must provide the path to a FASTA-formatted target sequence in order to run DMP. There are a few different ways to run it:

### From sequence only (requires legacy BLAST and HHblits):
`/path/to/run_DMP.sh -i input.fasta`

### From sequence and pre-made alignment in PSICOV format (requires legacy BLAST only):
`/path/to/run_DMP.sh -i input.fasta -a input.aln`

### From sequence, and PSSM in legacy BLAST makemat format (requires HHblits only):
`/path/to/run_DMP.sh -i input.fasta -m input.mtx`

### From sequence, pre-made alignment, and PSSM in legacy BLAST makemat format (does not require BLAST or HHblits):
`/path/to/run_DMP.sh -i input.fasta -a input.aln -m input.mtx`

Citing:
-------
If you find DeepMetaPSICOV useful, please cite our paper at bioRxiv: https://www.biorxiv.org/content/10.1101/586800v2
