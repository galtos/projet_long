# projet_long
Code for the creation,training of a deep learning model for the prediction
of interface residue given a protein.
The environnement is in the .yml file

an exemple for the prediction of the interface residues given the .rsa file
given by Naccess, the pssm file genereted by psiblast, aaindex file:

python3 projet_long -test -pdb "../data/test_file/9rsa0A0B_2.pdb"
                          -pssm "../data/test_file/9rsa0A0B_2.pssm"
                          -rsa "../data/test_file/9rsa0A0B_2.rsa"
                          
