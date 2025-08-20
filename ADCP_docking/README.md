ADCP v1.1 installation
======================
This repositiory contains installation scripts for ADCP v1.1  
**Using Conda Environment**  
Use adcpsuite_temp_conda.sh script  
**Commands**  
```sh
USERNAME@SYSTEM:~$ chmod +x adcpsuite_temp_conda.sh
USERNAME@SYSTEM:~$ ./adcpsuite_temp_conda.sh
```
In case conda/minoconda is not installed it will prompt a messege to install conda.  
In case conda/miniconda is installed but base environment is not activated, the code will ask you to activate it.  
**Command to activate conda base environment:** 
```sh
USERNAME@SYSTEM:~$ conda activate
(base) USERNAME@SYSTEM:~$ ## the standrad appearance of activated conda base environment
```

**Use ADCP**
After successful installation of ADCP v1.1. You should activate "adcpsuite" environment to use ADCP and related packages.
```sh
(base) USERNAME@SYSTEM:~$ conda activate adcpsuite
(adcpsuite) USERNAME@SYSTEM:~$  ## Appearance for activated adcpsuite environment
```
**Run adcp command**
```sh
(adcpsuite) USERNAME@SYSTEM:~$ adcp
```
The (trucated) output should be similar to
```
usage: usage: python runADCP.py -s GaRyMiChEL -T rec.trg -w WorkFolder -o output

AutoDock CrankPep V1.1.21

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -s SEQUENCE, --sequence SEQUENCE
                        initialize peptide from sequence, lower case for coil
                        and UPPER case for helix
  -p PARTITION, --partition PARTITION
                        partition for starting from a mixture of helix/coil
                        conformation, percentage(helix)=partition/100 note
                        this option will overwrite the CaSe in sequence
.
.
.
```
**MORE DETAILS ON RUNNING ADCP BASED DOCKING CAN BE FOUND FROM THE CCSB SITE**  
https://ccsb.scripps.edu/adcpv11/documentation/ 
