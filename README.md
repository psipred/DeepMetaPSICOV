# DeepMetaPSICOV 1.1.0
### Deep residual neural networks for protein contact prediction

Shaun M. Kandathil, Joe G. Greener and David T. Jones

University College London

Changes in version 1.1.0
------------------------

- Changes for compatibility with PyTorch 0.4.0 and above but reading in trained model saved with 0.3.0. The 'reference' version using PyTorch 0.3 is still provided.
- Minor bugfixes in test script.
- Updated documentation.

Requirements:
-------------
- Bash shell
- C and C++ compilers (tested with GCC 4.4.2, 4.8.5, and 5.4.0)
- Python 2 or 3 (preferably miniconda/anaconda, as this makes the PyTorch install much easier; tested on miniconda Python3)
- The following Python modules:
  - PyTorch >= 0.4.0 (tested on 1.1.0)
  
- Third-party programs:
  - HH-suite v3.0+ and a recent UniClust30 database (for making alignments; skip if you will only use pre-made alignments)
  - CCMpred v0.1.0 (Need this exact version; available [here](http://bioinfadmin.cs.ucl.ac.uk/downloads/ccmpred-0.1.0/CCMpred-0.1.0.tar.gz))
  - FreeContact 1.0.21 (available [here](https://rostlab.org/owiki/index.php/FreeContact))
  - Legacy BLAST 2.2.26 (executables `blastpgp` and `makemat`) and a suitable non-redundant database, e.g. Uniref90, formatted using `formatdb` (needed to generate PSIPRED and SOLVPRED inputs)

All other required programs written by our group are now bundled in this repo and do not need to be installed separately.

On some distributions, the C++ compiler is a separate add-on package and may not be installed by default. For example, on CentOS you will need to `yum install` packages `gcc` AND `gcc-c++`.

### Installing PyTorch
For most `conda` users, `conda install -c pytorch pytorch` should work. Alternatively, visit https://pytorch.org/get-started/locally/

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
If you find DeepMetaPSICOV useful, please cite our paper in _Proteins_: https://onlinelibrary.wiley.com/doi/full/10.1002/prot.25779

FAQs:
-----
### Do I have to use the exact versions of the programs that you mention?
Yes. DMP was trained on data output by specific versions of the feature generation programs, and you need to use the same versions during inference.

### Your paper mentions a bug in one of the feature generation programs; is this release affected?
The version of `alnstats` in this repository is not affected by the bug in question. Since we verified that the training of DMP did not suffer from the bug, we are releasing the bug-free version for inference. If for any reason you'd like the buggy version of `alnstats`, do get in touch.

### The version of CCMpred you use is ancient!
We know. Keen users will also have spotted that we use a number of input features in common with MetaPSICOV. We wanted to assess whether we improve over MetaPSICOV, DeepCov etc. using exactly the same training data, where possible. We are currently working on the next version of DMP, which will use the latest version of CCMpred, among other changes.
