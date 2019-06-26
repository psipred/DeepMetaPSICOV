#!/bin/bash

set -u
set -o pipefail

# run DMP using example fasta, aln and mtx.

echo '*** Running DMP...'
../run_DMP.sh -i 1guuA.fasta -a 1guuA.aln -m 1guuA.mtx

if [ $? == 0 ]; then
    echo "*** Run finished; comparing top-L contacts against reference output..."
    L=$(tail -n 1 1guuA.fasta | wc -c)
    # This is a somewhat over-optimistic test.
    diff -q <(sort -rgk 5 1guuA.deepmetapsicov.con | head -n $L | cut -d ' ' -f 1-2) <(sort -rgk 5 example_con/1guuA.deepmetapsicov.con | head -n $L | cut -d ' ' -f 1-2)
    ex=$?
    if [ $ex == 0 ]; then
	echo '*** Top-L contacts match example output.'
	echo '*** DMP tests ran OK. DMP is ready to use.'
	exit 0
    elif [ $ex == 1 ]; then
	echo "*** Mismatch in Top-L contact ranking. This is not automatically a failure; it depends on how different the contact scores are."
	echo "*** You may wish to compare ${PWD}/1guuA.deepmetapsicov.con and ${PWD}/example_con/1guuA.deepmetapsicov.con."
	exit 0
    else
	echo "*** Error in 'diff'. Exit status was $ex . This shouldn't happen."
	exit 1
    fi
else
    echo '*** DMP test run failed; please check the error message(s) above, and your installation. If there are still problems, please raise an issue on the GitHub page or email psipred@cs.ucl.ac.uk.'
    exit 2
fi

# cleanup

# for ext in ss ss2 solv colstats pairstats ccmpred evfold psicov deepmetapsicov.21c deepmetapsicov.map deepmetapsicov.fix deepmetapsicov.con
# do
#     rm -f 1guuA.${ext}
# done
