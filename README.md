Re-analysis based on Mathieson et al. (2024), “Cancer-associated fibroblasts expressing fibroblast activation protein and podoplanin in non-small cell lung cancer predict poor clinical outcome.” 

Goal: To isolate and characterize CAF-S1 and CAF-S5 in late NSCLC. 

Wu et al. data (late NSCLC) were taken from GSE148071. 
Code files include 
    bash_script.sh (to load samples)
    read_samples.R (to read samples into seurat objects)
    late_nsclc.R (QC, clustering, and fibroblast isolation)
    DE_analysis.R (DE analysis and pathway analysis)
    
