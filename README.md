# RNADigestCID
RNADigestCID is a package that can calculate post-RNase digestion RNA sequences and mass. It can further calculate RNA fragment sequences and mass generated by collision-induced dissociation (CID) in mass spectrometry. To run RNADigestCID, 1) install R then download "seqinr", "stringr" and "rstudioapi" packages; 2) Create a blank txt file named "FileDirect.txt" in the work directory; 3) Export path to the work directory.

# Calculate post-RNase digested RNA sequences and mass
Download and save RNA sequences in fasta format.
Paste fasta file path to the first line of FileDirect.txt.
In terminal, run: Rscript RNase( ).R
The output csv file will be saved in the work directory.

# Calculate CID fragmentation RNA sequences and mass
Paste post-RNase digested sequence csv file path to the second line of FileDirect.txt.
In terminal, run: Rscript CID_( ).R
The csv file with CID fragmentation RNA sequences and mass will be saved in the work diectory.

# Calculate CID fragments mass after base loss
Paste CID fragmentation RNA sequence csv file path to the third line of FileDirect.txt.
In terminal, run: Rscript CID_a-B.R
The csv file with CID base loss RNA sequences and mass will be saved in the work directory.

