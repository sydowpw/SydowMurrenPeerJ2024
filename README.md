Root Variation in Response to Auxin Across Genotypes; From Seedling to Maturity
Sydow & Murren 2023

# Questions for analysis

question #1

Do mutant lines of key auxin pathway genes vary in their response to exogenous 
treatments of auxin in the form of IAA? We predict that lines with insertion 
mutations in auxin pathway genes will vary in their response to exogenous auxin 
treatments showing differences in below ground and above ground phenotypes 
paired with differences in plant fitness quantified as fruit production. 

question #2

Do differences in root traits of mutant lines and natural accessions  of auxin 
genes change across developmental stages? We hypothesized that root traits in
mutant lines of auxin genesall genotypes would differ such that root traits 
would be greater in later developmental stages and vary significantly across 
mutant lines. Additionally, we hypothesized that auxin treatments would not 
affect the timing of life stage transitions in mutant lines in comparison to 
wildtype genotype.

question #3

Do the root phenotypic responses of mutant lines of auxin genes correspond with 
the number of alternative polyadenylation sites? Is this pattern the same in 
seedlings and reproductive adult plants? We predicted that genes with fewer 
poly(A) sites would function similarly across plant development assessed via 
root trait phenotypes when compared to mutants with inserts in genes with many 
poly(A) sites that may confer functional variation across the lifecycle of the 
plant. Specifically, we expect mutants of genes with fewer poly(A) sites to 
exhibit both fewer lateral roots and less overall root length.   

question #4

Is the phenotypic variation across auxin related gene mutant lines comparable to
that produced by natural variation across multiple natural accessions? We 
predict that natural accessions that exhibit variation across the genome will 
produce a wider range of phenotypes than insert mutants which share genetic 
backgrounds, and only differ at one locus. 

# The structure of the code-base

- data (contains all .csv files)
  
  - combined_data_clean.csv (contains genotypic and phenotypic data from 
                                summer 2021 experiment)
  - filtered_PAC_df.csv (contains PlantAPAdb PAC data)
  - T-DNA.SALK.NG.csv (contains SALK insert data for mutants)
  - APA_files (folder contains a .csv for every mutant locus of Araport11 
                genome annotation data)
                
- figs
  
  - Containing .pdfs of diagnostic plots of ANOVA tests analyzing models created 
    that were used to determine which traits/models to use for analysis.
                  
- source

  - question_1.R (Analysis of question #1)
  - question_2.R (Analysis of question #2)
  - question_3.R (Analysis of question #3)
  - ms_figures.R - Script for creating key figures and Levene's Test analysis of
                    question #4
  - gene_annotation.R - R script that creates visualizes the position and number
                        of APA sites and SALK inserts
  - insert_PAC_dist.R - R script that creates the variables "PAC_count" and
                        "mean_PAC_insert_dist" for all, "neg" and "pos" strands 
                        PACs in combined_data_clean.csv df
