#!/bin/bash

# run DMP using example fasta, aln and mtx.

../run_DMP.sh -i 1guuA.fasta -a 1guuA.aln -m 1guuA.mtx

if [ $? == 0 ]; then
    echo 'DMP test ran OK. DMP is ready to use.'
else
    echo 'DMP test failed; please raise an issue on the GitHub page or email psipred@cs.ucl.ac.uk.'
fi

# cleanup
for ext in ss ss2 solv colstats pairstats ccmpred evfold psicov deepmetapsicov.21c deepmetapsicov.map deepmetapsicov.fix deepmetapsicov.con
do
    rm -f 1guuA.${ext}
done
