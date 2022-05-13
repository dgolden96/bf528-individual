# Roles Selected: Project 2 Analyst, Project 3 Analyst

This repository contains Daniel Golden's R scripts used in the completion of the analyst roles for projects 2 and 3 for the individual project final of BF528.

In this project, I have completed the tasks associated with the analyst role in projects 2 and 3, analyzing either instructor-provided data or programmer-generated data associated with data analysis in the recreation of the results reported by O'Meara et al. and Wang et al.

# Repository Contents

## O'Meara et al. analyst script
All that is required to generate the plots and data for the analyst role in project 2 is to open the following .R file in Rstudio, substitute the filepath used in the script to read in the gene_exp.diff file with the filepath for that file locally stored on the user's computer (or the following filepath if Rstudio is being run on the SCC: "/project/bf528/project_2/data/P4_vs_P7_cuffdiff_out/gene_exp.diff"), and run the script. 
```
BF528_final_omeara_analyst.R
```

## Wang et al. analyst script
All that is required to generate the plots and data for the analyst role in project 3 is to open the following .R file in Rstudio, substitute the filepaths used in the script to read in the RMA expression matrix and toxgroup 6 metadata with the filepaths for those files locally stored on the user's computer (or the following filepaths if Rstudio is being run on the SCC: "/project/bf528/project_3/samples/liver-normalization-rma.txt" and "/project/bf528/project_3/groups/group_6_mic_info.csv"), and run the script.
```
BF528_final_wang_analyst.R
```
