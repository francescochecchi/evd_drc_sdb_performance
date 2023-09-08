## Performance of a safe and dignified burial intervention during an Ebola epidemic in the eastern Democratic Republic of Congo, 2018-2019
Warsame et al. (2021) 

## Notes on data and R analysis code
17 December 2022  
Francesco Checchi  
Department of Infectious Disease Epidemiology  
Faculty of Epidemiology and Population Health 
London School of Hygiene and Tropical Medicine  
Francesco.checchi@lshtm.ac.uk  

### Donor
R2HC / ELRHA.

### Background on the study
This repository contains data and R scripts required to replicate an analysis of a Safe and Dignified Burial (SDB) service for suspect or confirmed Ebola Virus Disease (EVD) cases during an epidemic in eastern Democratic Republic of Congo (DRC). The epidemic lasted from 2018 to 2020, but data for this analysis only include 2018 and 2019. The analysis was conducted by the London School of Hygiene and Tropical Medicine (www.lshtm.ac.uk) at the request of the International Federation of Red Cross and Red Crescent Societies (IFRC), who supported the SDB service and oversaw collection of data on individual instances of SDB service delivery.

### Datasets
All datasets required to replicate the analysis are included in the file <evd_drc_sdb_performance_datasets_pub.xlsx>. The file includes a ‘dictionary’ tab.

### R analysis script
The script `evd_drc_sdb_performance_pub.R` loads necessary R packages, sets parameters (these can be modified by the user), loads datasets, manages data, runs analyses and outputs some tables and figures. Only this script needs to be run to replicate the analysis. Computation time should be <5 min on a standard laptop.

### Requirements to replicate the analysis
The dataset and R script files must be stored on the same folder, where output files will be saved. The directory for reading files is set automatically while the R script is run.
It is recommended to run the code from the open-source RStudio interface (https://rstudio.com/ ).
