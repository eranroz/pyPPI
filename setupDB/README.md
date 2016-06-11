setupDB script creates a protein database in SQL:
A database of interesting properties of the proteins based on scripts of this library.

This should be easy to use script for invoking the most important scripts of the library and store them in DB
for easy retrieve.

How to use:
1. Create a folder and place there some file with list of PDBs to analyze. (each in new line)
2. The program will create the following directory structure in the same directory:
    ./pdbs/ - list of pdbs downloaded from PDB data bank
    ./results/ - some results of the analysis scripts
    ./debug/ - raw results and pymol scripts
    