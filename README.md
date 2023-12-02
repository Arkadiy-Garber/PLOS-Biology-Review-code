## This repository contains the code and data associated with PLOS Biology special issue review:
[LINK TO ARTICLE]

### Identification of cell envelope and translation-related genes in a 85 bacterial genomes (Figure 2)
 - The main script to run is **yfgenie_post_processing.py**

 - All files referenced in that script are present within this GitHub repository, except for the raw output from HMMER analysis, which you can find here: [10.6084/m9.figshare.24716529](10.6084/m9.figshare.24716529)

 - The output of this script is a CSV file called **final_summary.aminoacyl_updated2.ribo_updated.csv**, which was used to make Figure 2.

### PCA analysis of amino acid frequencies in >300,000 prokaryotioc genomes (Figure 3)
 - The **gc.py** script was used to process the roughly 330,000 annotated prokarytic genomes listed in the ncbi_assembly_info.tsv files.

 - All the output files from this analysis were concatanated to make the **combined.aa_gc.txt.zip** file, which was used to generate Figure 3
