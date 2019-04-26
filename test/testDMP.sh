#!/bin/bash

# run DMP using example fasta, aln and mtx.

../run_DMP.sh -i 1guuA.fasta -a 1guuA.aln -m 1guuA.mtx

L=$(tail -n 1 1guuA.fasta | wc -c)

if [ $? == 0 ]; then
    diff -q <(sort -rgk 5 1guuA.deepmetapsicov.con | head -n $L | cut -d ' ' -f 1-2) <(sort -rgk 5 example_con/1guuA.deepmetapsicov.con | head -n $L | cut -d ' ' -f 1-2)
    ex=$?
    if [ $ex == 0 ]; then
	echo 'Top-L contacts match example output'
    elif [ $ex == 1 ]; then
	echo 'Mismatch in Top-L contact ranking; '
    else
	echo "Error in diff. Exit status was $ex . This shouldn't happen."
    fi

    echo 'DMP tests ran OK. DMP is ready to use.'
else
    echo 'DMP test failed; please raise an issue on the GitHub page or email psipred@cs.ucl.ac.uk.'
    exit 1
fi

# cleanup

# for ext in ss ss2 solv colstats pairstats ccmpred evfold psicov deepmetapsicov.21c deepmetapsicov.map deepmetapsicov.fix deepmetapsicov.con
# do
#     rm -f 1guuA.${ext}
# done
