============= version 2.01 =============
(1) Change the output format while option '-j' set to be 0 in single-end reads correction. Empty reads will be discarded in output file.

============= version 2.00 =============
(1) Add the parallel reads mode in KmerFreq, to increase speed.
(2) Add "FAST" method in error correction.
(3) Revise the branch and bound trie(BBT) to increase efficiency in correction.
(4) Add some parameters for trim strategy.
(5) Change the parallel reads mode in Corrector, to increase speed and reduce memory usage.
(6) Support larger kmer size, two sets program: *_AR (for small kmer <= 17bp) and *_HA (for big kmer > 17bp).
(7) Support space-kmer in KmerFreq_AR and Duo-kmer(consective and space kmer) in Corrector_AR.
(8) Support input compressed read file in *gz format.
(9) Add the quality control in correction, output *.QC.xls file.
(10) Support output fastq format, and change the output pair format in Pair-End file, replace the merge pair file.
(11) Support compressed ouput(*.gz).

============= version 1.00 =============
The first released version.