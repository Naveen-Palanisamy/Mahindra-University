Peptide Docking using ADCP v1.1
===============================

**Activate adcpsuite environment**  
```sh
conda activate adcpsuite  
```
  
**Download pdb file:**  
https://files.rcsb.org/download/3Q47.pdb  
  
**In Chimera:**  
Clear Water
Save receptor and peptide pdb files  
  
**Add H:**
Use reduce command in terminal to add hydrogens
```sh
reduce 3q47_rec.pdb > 3q47_recH.pdb  
reduce 3q47_pep.pdb > 3q47_pepH.pdb  
```
  
**Preprare Molecules:**  
```sh
prepare_receptor -r 3q47_recH.pdb  
prepare_ligand -l 3q47_pepH.pdb
```
  
**Create target File:**  
```sh
agfr -r 3q47_recH.pdbqt -l 3q47_pepH.pdbqt -asv 1.1 -o 3q47 
```
  
**Display target file with agfrgui**    
```sh
agfrgui
```
  
**run docking with ADCP**  
```sh
adcp -T 3q47.trg -s "npisdvd" -N 2 -n 2000 -w testDir -o testNSAAandD -O
```
  
**minimize 2 solutions after docking**  
```sh
adcp -T 3q47.trg -s "npisdvd" -N 2 -n 2000 -w testDir -o testNSAAandD -O -nmin 2
```
    
**Display files**  
```sh
chimera 3q47.pdb testDir/testNSAAandD_omm_rescored_out.pdb
```
  
**Using Non-standard AA and D-type**  
***Files***  
```sh
cat  $CONDA_PREFIX/lib/python3.7/site-packages/ADCP/data/rotamers/SWISS_ROTAMER_LIST.dat
```
https://www.swisssidechain.ch/data/family_table.pdf  

```sh
adcp -T 3q47.trg -s "n<HGA>pisdvd" -N 5 -n 2000 -L swiss -w testDir -o testNSAAandD -O -nmin 2
```
   
**try something random**  
```
adcp -T 3q47.trg -s "n<ADAM>pisdvd" -N 5 -n 2000 -L swiss -w testDir -o testNSAAandD -O -nmin 2
```
  
**More details:**  
https://ccsb.scripps.edu/adcpv11/documentation/  
