ADCP v1.1 installation
This repositiory contains installation scripts for ADCP v1.1
Using Conda Environment
Use adcpsuite_temp_conda.sh script
Commands

USERNAME@SYSTEM:~$ chmod +x adcpsuite_temp_conda.sh
USERNAME@SYSTEM:~$ ./adcpsuite_temp_conda.sh
In case conda/minoconda is not installed it will prompt a messege to install conda.
In case conda/miniconda is installed but base environment is not activated, the code will ask you to activate it.
Command to activate conda base environment:

USERNAME@SYSTEM:~$ conda activate
(base) USERNAME@SYSTEM:~$ ## the standrad appearance of activated conda base environment
Use ADCP After successful installation of ADCP v1.1. You should activate "adcpsuite" environment to use ADCP and related packages.

(base) USERNAME@SYSTEM:~$ conda activate adcpsuite
(adcpsuite) USERNAME@SYSTEM:~$  ## Appearance for activated adcpsuite environment
Run adcp command

(adcpsuite) USERNAME@SYSTEM:~$ adcp
The (trucated) output should be similar to

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
We have also provided a tester file that checks commands required for basic docking
Download the the tester.sh file from this directory and run is using this command in activated adcpsuite environment:

(adcpsuite) USERNAME@SYSTEM:~$ chmod +x tester.sh
(adcpsuite) USERNAME@SYSTEM:~$ ./tester.sh
This command may take a few minutes to complete. After completion, it will generate a log file that can be analyzed to find errors.

cat test_ADCP\log.out
MORE DETAILS ON RUNNING ADCP BASED DOCKING CAN BE FOUND FROM THE CCSB SITE
https://ccsb.scripps.edu/adcpv11/documentation/
