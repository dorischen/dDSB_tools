# dDSB_tools

Essential scripts from our publication “**Spo11 generates gaps through concerted cuts at sites of topological stress**” (Prieler, S., Chen, D., Huang, L. et al., *Nature (2021)*, https://doi.org/10.1038/s41586-021-03632-x) are listed.

The R scripts are provided as stand-alone scripts to offer high flexibility and user-friendliness.
They consist of five major sections:
1.  **PACKAGES**
2.  **SUB-FUNCTIONS**
3.  **USER INPUT**
4.  **DEFAULTS**
5.  **FUNCTION**

Users need to install the packages (libraries, 1) and run the sub-functions (2) only once per session. Script options are changed in the “USER INPUT” part (3); each time changes are made, this part needs to be run before executing the default settings (4) and the function (5). Default settings contain e.g. output file suffixes or genome data details.

It is highly recommended to use *RStudio* (https://www.rstudio.com/) to execute these scripts. Sample data displaying the required input formats are provided (and referred to in section 3). Please note that these scripts were originally highly customised to the purposes of the published analyses and will be updated for more general use in future. Furthermore, scripts for creating certain data visualisations (like arc plots, profile plus annotation plots) are being converted into more generic versions and will potentially be published as R package.

For further questions, requests or comments please contact Doris Chen, doris.chen@univie.ac.at.

