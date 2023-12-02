## This repository contains the code and data associated with PLOS Biology special issue review:
[LINK TO ARTICLE]

### Identification of cell envelope and translation-related genes in a 85 bacterial genomes (Figure 2)
 - The main script to run is **yfgenie_post_processing.py**

 - HMMs were queried against a select set of 85 genomes from NCBI's RefSeq and GenBank databases using YfGenie software from the MagicLamp toolkit: [https://github.com/Arkadiy-Garber/MagicLamp](https://github.com/Arkadiy-Garber/MagicLamp)

        GCA_000008885.1
        GCA_000021505.1
        GCA_000043285.1
        GCA_000147015.1
        GCA_000331975.1
        GCA_000412755.1
        GCA_000699475.1
        GCA_000953435.1
        GCA_001029715.1
        GCA_002074035.1
        GCA_002290365.1
        GCA_003391315.1
        GCA_005697395.1
        GCA_015134435.1
        GCA_015139955.1
        GCA_015256915.3
        GCA_015257795.3
        GCA_015697925.1
        GCA_016699915.1
        GCA_017134395.1
        GCA_023898665.1
        GCA_024818515.1
        GCA_025999335.1
        GCA_025999395.1
        GCA_026016225.1
        GCA_027925225.1
        GCA_030060845.1
        GCA_900016775.1
        GCA_900019315.1
        GCA_944407095.1
        GCA_946893605.1
        GCF_000008605.1
        GCF_000011745.1
        GCF_000012345.1
        GCF_000013165.1
        GCF_000021865.1
        GCF_000027325.1
        GCF_000093065.1
        GCF_000195735.1
        GCF_000196075.1
        GCF_000219175.1
        GCF_000247565.1
        GCF_000287295.1
        GCF_000304455.1
        GCF_000709555.1
        GCF_000754265.1
        GCF_000828515.1
        GCF_000829155.1
        GCF_000931575.1
        GCF_000974825.1
        GCF_001242845.1
        GCF_001269425.1
        GCF_002119845.1
        GCF_002776555.1
        GCF_002777075.1
        GCF_002849995.1
        GCF_004296475.1
        GCF_005080705.1
        GCF_009883795.1
        GCF_012562765.1
        GCF_012571345.1
        GCF_013030075.1
        GCF_013463375.1
        GCF_013463555.1
        GCF_014251515.1
        GCF_014251975.1
        GCF_016889585.1
        GCF_017656055.1
        GCF_019669085.1
        GCF_019923565.1
        GCF_020541245.1
        GCF_021222645.1
        GCF_024205245.1
        GCF_025021925.1
        GCF_026644375.1
        GCF_028856405.1
        GCF_900039485.1
        GCF_900044015.1
        GCF_900048035.1
        GCF_900048045.1
        GCF_900090215.1
        GCF_900343015.1
        GCF_900638335.1
        GCF_902713755.1
        GCF_902860225.1

- All files referenced in that script are present within this GitHub repository, except for the raw output from HMMER analysis, which you can find here: [10.6084/m9.figshare.24716529](10.6084/m9.figshare.24716529)

 - The output of this script is a CSV file called **final_summary.aminoacyl_updated2.ribo_updated.csv**, which was used to make Figure 2.

### PCA analysis of amino acid frequencies in >300,000 prokaryotioc genomes (Figure 3)
 - The **gc.py** script was used to process the roughly 330,000 annotated prokarytic genomes listed in the ncbi_assembly_info.tsv files.

 - All the output files from this analysis were concatanated to make the **combined.aa_gc.txt.zip** file, which was used to generate Figure 3
