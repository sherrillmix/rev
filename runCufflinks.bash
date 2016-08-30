#!/bin/bash

baseDir=work/align
tmpDir=/tmp
echo "Running cuffdiff"
diffDir=work/diff
~/installs/cufflinks/cuffdiff --max-bundle-frags 100000000 -o $diffDir -p 16 --min-reps-for-js-test 1 -v ~/installs/star/index/knownGene.gtf work/align/control_1.bam,work/align/control_1.bam work/align/rev_1.bam,work/align/rev_2.bam

