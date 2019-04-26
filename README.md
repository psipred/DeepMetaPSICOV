# DeepMetaPSICOV 1.0.0
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
  - CCMpred v0.1.0 (Need this exact version; available [here](http://bioinfadmin.cs.ucl.ac.uk/downloads/ccmpred-0.1.0/CCMpred-0.1.0.tar.gz))
  - FreeContact 1.0.21 (available [here](https://rostlab.org/owiki/index.php/FreeContact))
  - Legacy BLAST 2.2.26 (executables `blastpgp` and `makemat`) and a suitable non-redundant database, e.g. Uniref90, formatted using `formatdb` (needed to generate PSIPRED and SOLVPRED inputs)

All other required programs written by our group are now bundled in this repo and do not need to be installed separately.

On some distributions, the C++ compiler is a separate add-on package and may not be installed by default. For example, on CentOS you will need to `yum install` packages `gcc` AND `gcc-c++`.

### Installing PyTorch using conda
Use `conda install -c pytorch pytorch=0.3.1`

[comment]: # (### GPU support)
[comment]: # (If you are installing PyTorch using `conda`, `conda` should automatically detect a usable GPU and install an appropriate version with GPU support.)
[comment]: # (NB a GPU is not necessary for predicting contacts, but for faster runtimes we recommend using one.)

Setup and testing:
------------------

### Build bundled programs
`cd src/; make; make install`

### Specify paths to external dependencies
Edit `run_DMP.sh` to indicate the paths to the third-party programs listed above, as well as other variables such as the number of threads to use for various programs. User-editable variables are demarcated by comment lines. We do not recommend changing anything outside this region unless you know what you are doing.

### Testing
`cd test; ./testDMP.sh`

The script will use the configuration you have provided in `run_DMP.sh` and run a test contact prediction.
*NB:* the test script will not run PSI-BLAST or HHBlits; the script runs only the remaining parts of the DMP pipeline (using running Option 4 below).

Different versions of OSs, compilers etc. can lead to differing contact scores (as well as outputs from the feature generation programs), so we only test the ranking of the top-L predicted contacts against a reference output.

Running:
--------
Run `/path/to/run_DMP.sh -h` to see the available options. DMP runs a number of programs to generate input features; their outputs are stored in a number of intermediate files. By default, DMP will attempt to reuse any files with the correct filenames (this is useful for debugging and allows you to 'continue' a failed run). You can force regeneration of intermediate files with the `--force` option.

At a minimum, you must provide the path to a FASTA-formatted target sequence in order to run DMP. There are a few different ways to run it:

### Option 1: From sequence only (requires legacy BLAST and HHblits):
`/path/to/run_DMP.sh -i input.fasta`

### Option 2: From sequence and pre-made alignment in PSICOV format (requires legacy BLAST only):
`/path/to/run_DMP.sh -i input.fasta -a input.aln`

### Option 3: From sequence, and PSSM in legacy BLAST makemat format (requires HHblits only):
`/path/to/run_DMP.sh -i input.fasta -m input.mtx`

### Option 4: From sequence, pre-made alignment, and PSSM in legacy BLAST makemat format (does not require BLAST or HHblits):
`/path/to/run_DMP.sh -i input.fasta -a input.aln -m input.mtx`

Citing:
-------
If you find DeepMetaPSICOV useful, please cite our paper at bioRxiv: https://www.biorxiv.org/content/10.1101/586800v2

FAQs:
-----
### Do I have to use the exact versions of the programs that you mention?
Yes. DMP was trained on data output by specific versions of the feature generation programs, and you need to use the same versions during inference.

### Your paper mentions a bug in one of the feature generation programs; is this release affected?
The version of `alnstats` in this repository is not affected by the bug in question. Since we verified that the training of DMP did not suffer from the bug, we are releasing the bug-free version for inference. If for any reason you'd like the buggy version of `alnstats`, do get in touch.

### The version of CCMpred you use is ancient!
We know. Keen users will also have spotted that we use a number of input features in common with MetaPSICOV. We wanted to assess whether we improve over MetaPSICOV, DeepCov etc. using exactly the same training data, where possible. We are currently working on the next version of DMP, which will use the latest version of CCMpred, among other changes.

### The version of PyTorch you use is ancient!
The development of the DMP v1 began when PyTorch 0.3.0 was current.
By the end of CASP13, v1.0 was being prepared for release.
The DMP models were trained on v0.3.0 and although it is possible to read in and run the trained models using PyTorch v0.4-1.0+, we find that there are occasionally significant differences in the contact scores output.
In order to keep things as close as possible to the version we ran in CASP13, we are recommending that users use PyTorch 0.3.0 or 0.3.1 for DMP v1.0. We intend to upgrade the version of PyTorch used in the next release.
