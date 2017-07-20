## Seed-and-Wobble-Python-Translation

### Original Notes sent to Professor
Please note that code was written on various machines and various different editors, and the indentation between different files,
and sometimes even between different functions in SAndW.py and seedAndWobbleModsLim.py is not uniform. This has not given me any run-time issues,
but will be cleaned up in future.

Version 0 (Not included on Git-hub) is my first working translation of the PERL seed and wobble into python [perhaps with slight edits].
It is more procedural than object oriented.
Contains SAndW.py - the main program - and seedAndWobbleModsLim.py - the modules required for SAndW.py
Takes about 280 seconds to run on testing data*

Version 1 is a tidied up version. The code is better documented, hopefully easier to read, and slightly more efficient.
Optional flags have been added to: 1) Allow the user to view and change the running parameters, 
2) Disable automatic quitting on [some] I/O errors ,and 3) Calculate resource usage - time taken and, in future, memory usage.
More work can probably be done on both cleaning up the code and the efficiency. The module file has hardly been changed.
Contains SAndW.py - the main program, and seedAndWobbleModsLim.py - the modules required for SAndW.py
Takes about 260 seconds to run on testing data

Version 2 is an object oriented version. The code is clearer, but I have not yet been able to optimise it for efficiency.
Again, more work can be done, and the modules file has only been changed insofar as it relates to the added objects.
Contains SAndW.py - the main program, seedAndWobbleModsLim.py - the modules required for SAndW.py, and SWclasses.py - the classes 
which are used in the programs.  
Takes about 330 seconds to run on testing data - hopefully this can be improved upon.

Command-line for testing is as follows: 

Version 0:

python SAndW.py proteinA_v1_combinatorial.txt 8 patterns_8of10.txt patterns_4x44k_all_8mer.txt Test

Version 1 + 2: 

python SAndW.py proteinA_v1_combinatorial.txt 8 patterns_8of10.txt patterns_4x44k_all_8mer.txt Test F T T

The files used for testing are available in the TestFiles directory.
Tested on an E-Machines laptop, intel core I5-430M 2.27 GHz * 4, with 4 GB ram. Tested in Ubuntu 14.04, Python v 2.7.6

### Additional notes and brief explanation of seed and wobble methodology
The seed and wobble methodology is used to find likely binding sites for transcription factors. These are proteins which bind to the DNA and control by either allowing, enhancing ordampening the transcription/translation of Genes further down the strem. These binding sites are fairly difficult to define accurately, so a lot of work has been done trying to  find fairly accurate, sensitive and specific representations. A common one is a Position Weight Matrix (PWM), which will create a matrix of a certain length, and give the likelihood of a certain base (A,C, G or T) being in that position for a binding site.

Improvements in DNA sequencing technologies has allowed scientists to create Protein Binding Micro-arrays (PBMs) which contain a combinatorial sequencing where all possible bases between given lengths (such as from 8 to 14 base-pairs long) will occur a certain amount of times. When the transcription factor of interest is applied to these PBM's, they can be tested, and will glow at different arrays. By measuring all of the Microarrays in which a certain Base-pair sequence (known as a kmer) occurs, one can get an enrichment score for that Kmer, and using that, one can calculate likely PWMs for the transcription factor. 

There is some pre-processing required to transform the scores taken from a PBM experiment into a format usable by seed and wobble, but, since we did not have access to our own PBM expiriments, and there would be no reason to make changes in this processing in any case, this was not translated from Perl into Python.

How the seed and wobble algorithm works: 

Seed and wobble takes as input a combinatorial list of the intensities of all the microarrays and their intensities, ranked from highest to lowest, and a list of canddate seed patterns. These seed patterns are simply a string of "1"s and "."s, where a "1" stands for a defined base, and a "." is a wildcard, that could be any base. For each seed pattern, each possible Kmer Pattern (Naively, it's4^8 possibilities, for a kmer of length 8, with or without wildcards, but it's actually somehwat less, since there are reverse compliments - the DNA double helix can be "read" in either direction, where we substitute "T"s for "A"s and "C"s for "G"s when reading in reverse.)

The Enrichment of ALL of these Kmers are calculated. Rather than averaging the intensity value of the microarrays  in which the selected Kmers appear, the enrichment is counted on a basis of how highly the intensities of their microarrays Rank, giving a ranking of between 1 and minus 1 - for example, if the microarrays containing  the kmer are the topped ranked intensity microarrays, the enrichment score would be 1, while if they were all the bottom ranks, the enrichment score would be minus 1, and if the rankings were clustered around tghe central intensity value, the score would be 0). This is done to avoid inaccuracies inherent in the PBM expirement The enrichment scores for all kmers which meat the seed patterns are calculated, and then the wobble step is used to calculate PWMs for the top 3 (by default) Kmers. Basically, from the candidate KMer - each position is wobbled - i.e. the enrichment is calculated for the defined positions, for each possible base in that position. The relative enrichment of these Kmers (4 possible for each position) is used to calculate a percentage likelihood for that base to be in that position. 
