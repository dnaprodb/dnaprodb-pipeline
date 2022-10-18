# DNAproDB Processing Pipeline
The DNAproDB processing pipeline is built as a collection of python scripts that take as input a PDB entry (in mmcif or PDB file format) and produce a single JSON file as output containing all information about the DNA and protein structural characteristics and the interactions between protein residues and nucleotides in the complex. Additional annotations from several other databases are used as well. A schematic of the pipeline is shown below.


It is recommended to use the mmcif file format, as the header info is parsed and additonal annotations can be extracted from it, which would be missing for the PDB file format.

## Install Instructions
To install, it is first recommended to create a new python enviroment via conda:

```
conda create -n dnaprodb python=3.8
conda activate dnaprodb
```

Next, clone this repository and install python dependencies via pip

```
git clone https://github.com/jaredsagendorf/dnaprodb-back.git
cd dnaprodb-back
pip install -r python-requirements.txt
```

Dependencies
reduce - 
hbplus - 
Curves -
3DNA -
pdb2pqr - 
dssr - 
dssp - 
SNAP - 
