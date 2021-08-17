# chromatin_potential
calculating chromatin potential 

Order of running scripts:

Prerequisites: The RNA and DNA data should be ordered and named with identical cell barcodes.  RNA expression should be log normalized.  Peak counts shoud be the raw numbers.

Read_in_initial_data.R to get all the files needed from the RDS objects.
  
  -This is set up to read in the rna_exression matrix and peak_expression matrix of the foxA1 data as SummarizedExperiment and RangedSummarizedExperiment Objects.  This file will have to be modified based on the input format and names to create the necessary files for downstream use.
  
  
  
    (In order of how they are written out)
    
    There is a luminalvsbasal file that will be used downstream to plot the celltypes in the final UMAP.
    
    There is a peaks bed file output that contains the list of peak names.
    
    There is a peaks celltype file that contains the cell barcodes given in the ATACseq data.
    
    There is a raw peak counts mtx file that is a sparse matrix of the raw fragment counts data.
    
    There is a normalized peak counts mtx file that is a sparse matrix of the log-normalized fragment counts data.
    
    There is an rna features file that contains the RNA features in order.
    
    There is an rna counts mtx file which is a sparse matrix of the normalized RNA expression.
    
    There is an rna celltypes file which gives the RNA cell barcodes in order.
    
    

Createannotaions.R to annotate the peaks. (Can combine with read_in_initial_data.R if desired).

  
  
    Just feed in the peak names bed file and get the output annotations file.
    
    Can adjust the genome being used if necessary.
   
  
  

Start running cistopic.R with the necessary initial data.  Run the next step (Readinannotations.ipynb and FindingSpearman.ipynb) at the same time.

  
  
    Input the entire peak set (counts, peaks bed names, and peak cell barcodes)
    
    Output: The cisTopic topic-cell probabilities and the UMAP coordinates produced by cisTopic.
  
  
  
Readinannotations.ipynb and then FindingSpearman.ipynb to find DORCS.


  
    In ReadInAnnotations, input the annotation file produced in Createannotations.R and get a peak table as output.
  
    In FindingSpearman read in this peak table. the peak counts matrix, the cell barcodes for both peaks and rna, the rna expression matrix, and the RNA features.
    
      - Next, output the peaks_used_in_barcodes file and read it into generating_background_peaks.R along with the peaks bed file.
      - Take the ouputted background peaks, read them into findingspearman, and continue.  The correlations will take a long time to run, a 5000 cell dataset will take 0.5 days.
      - The significantly correlated peaks will be outputted.
      



Smoothing_scores python to get the smoothed matrices.


    
    Read in the significantly correlated peaks that were output from findingspearman, the normalized peak counts, the rna counts, and the cistopic topic-cell probability matrix.
    
    Ouptut is the smoothed DORC scores, smoothed RNA exression, and the cellbarcodes for the DORC scores.
    



Put into chr_pot.R to see chromatin potential.



    Input the smoothed tables, barcodes, and celltypes (luminal vs. basal) along with the umap coordinates (mapcoords).  
    Make sure peak names are in same format between cellindex and rownames(mapcoords). This means they must be both length 300 with the same spacing characters. 
    If not using to look at luminal vs basal, change the umap.smooth mutate cases.
    
    In ggplot2, change the scale_mag maxsize parameter to adjust the arrow length.  Change the skip parameter in geom_arrow to adjust how sparse the arrow layout is.

