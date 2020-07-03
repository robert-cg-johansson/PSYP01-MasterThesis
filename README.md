# PSYP01-MasterThesis
Analysis Code for Masters Thesis

GP_nopauses.txt
Contains the response time and response force data from participant GP.

PF_nopauses.txt
Contains the response time and response force data from participant PF.

RF_RT_Relations.R
Examines relations between response time and response force, and produces a series of plots.

RT_RF_ByCondition.R
This computes bayes factors and 95% highest density intervals for the response time and response force data, and produces a series of plots.
Also subjects the response time and response force data to a joint analysis examining the amount of force transmitted to the response button per unit time ('power').
Requires you to have "Utilities.R" in your working directory.

Redundancy_Models.R
Formally tests the data for the redundant signals effect and motor coactivation.

SFT_Analysis.R
Performs capacity analysis, the survivor interaction contrast, and plots the CDFs of the response time distributions (race model inequality).

Utilities.R
Processes the raw data from .txt-format and contains algorithms for analysis.
