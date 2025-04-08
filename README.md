# Statistical Analysis of GC/MS data
Pipeline for statistical analysis of MZmine output of GC/MS data used in Papazlatani et al., 2025 (10.1016/j.crmicr.2025.100385).

## Steps to take before the analysis
1. Extract the peak intensity table from MZmine in the metaboanalyst form
2. Manually curate the peak intensity table to remove
     i. Silica Compounds
     ii. Analysis artifacts / False compounds
3. Replace missing vlaues with 0
4. Transpose the table
5.  In the analysis folder prepare the following folders:
    i. "0_Data_preparation"
    ii. "1_Data_norm"
    iii. "2_Diff.Int"
    iv. "3_Burk_Comp_Norm"
    v. "4_Diff_Int_Media"
6. Place the curated peak intensity table in the "0_Data_preparation" folder
