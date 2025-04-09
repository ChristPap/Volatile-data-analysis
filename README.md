# Statistical Analysis of GC/MS data
Pipeline for statistical analysis of MZmine output of GC/MS data used in Papazlatani et al., 2025 (10.1016/j.crmicr.2025.100385).

In this study, the Gram negative bacterium _Burkholderia_ sp. AD24 was stimulated to emit volatiles that suppress the growth of plant pathogenic fungi. For this purpose _Burkholderia_ sp. AD24 was cultivated on 0.1 TSBA and 0.1 TSBA supplemented with amino acids (in mixture or individually). Non-inoculated controls of each growth medium were also prepared to distinguish between bacterial and ambient volatiles.

The "Volatile statistical analysis.R" script firstly normalizes the peak intensity data using the "pretreat command V2.1.R". Then the bacterial volatiles are distinguished from the ambient volatiles by performing pairwise comaprison between the inoculated and non-inoculated treatment for each growth medium.
The bacterial volatiles are further subjected to normalization using the "pretreat command V2.1.R" and differences in the volatiles' intensities when _Burkholderia_ sp. AD24 was cultivated on 0.1 TSBA and 0.1 TSBA supplemented with amino acids were uncovered by pairwise comparison.
Lastly an overview of the data was obtained by performing Partial Least Squares - Discriminant Analysis (PLSDA) analysis on the normalized peak intensity table.

## Steps to take before the analysis
1. Extract the peak intensity table from MZmine in the metaboanalyst form
2. Manually curate the peak intensity table to remove
     - Silica Compounds
     - Analysis artifacts / False compounds
3. Replace missing vlaues with 0
4. Transpose the table
5. Create the "Volatilomic analysis" folder in which prepare the following folders:
     - "0_Data_preparation"
     - "1_Data_norm"
     - "2_Diff.Int"
     - "3_Burk_Comp_Norm"
     - "4_Diff_Int_Media"
6. Place the curated peak intensity table in the "0_Data_preparation" folder
