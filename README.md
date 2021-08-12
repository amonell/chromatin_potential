# chromatin_potential
calculating chromatin potential 

Order of running scripts:

Read_in_initial_data.R to get all the files needed from the RDS objects.

Createannotaions.R to annotate the peaks. (Can combine with read_in_initial_data.R if desired).

Start running cistopic.R with the necessary initial data.

Readinannotations python script and then findingspearman python script to find DORCS.

Smoothing_scores python to get the smoothed matrices.

