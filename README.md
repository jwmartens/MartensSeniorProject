# MartensSeniorProject
Program calculates an immune induction score for the bacteria specific to multiple sclerosis (MS). This is accomplished within the script through multiple methods, these being main(), bacterialprevalence(), upregulated(), downregulated(), and geneexpression(). The batcerialprevalence, upregulated, and downregulated methods work by gathering information stored in pyhton dictionaries specific to each bacteria. This information was compiled through comprehensive literature review, pathway analysis, and protein homolog research in order to determine prevalence in MS and metabolite profiles for each bacteria. The geneexpression method takes the input of a csv file of differentially expressed genes with headings (blank), logFC, AveExpr, t, P.Value, adj.P.Val, B, and symbol. For this study, gene expression data was calculated using R (Affy, limma) using GSE21942 as a data source. Calculated values are found first by finding the sum of the inflammatory and anti-inflammatory pathways that are induced (+1 for inflammatory, -1 for anti-inflammatory) and the sum of the inflammatory and anti-inflammatory pathways that are repressed (+1 for anti-inflammatory, -1 for inflammatory) for each metabolite-derived pathway effect. This number is then weighted according to the average host gene expression of the pathways induced by the bacteria (either inflammatory or anti-inflammatory). Lastly, these values are weighted additionally by multiplication of the prevalence of the bacteria in MS. This culmination of data gives the immune induction scores for each bacterial species.

CONTACT
  If any issues or methodology questions arise in the running of this program, contact is jwmartens@unomaha.edu.

RUNNING THE PROGRAM
Sample input :
  python immuneinduction.py -i Martens_MS_DE_result.csv

Requires -i option corresponding to the gene expression data that the user wishes to use for the calculations. Must be in CSV format.

Sample output:
  Immune Induction Score
  Methanobrevibacter: 1.351
  Proteobacteria-Caulobacteraceae: 1.110
  Proteobacteria-Mycoplana: 1.110
  Proteobacteria-Pseudomonas: 1.003
  Acinetobacter calcoaceticus: 0.767
  Akkermansia muciniphila: 0.747
  Sutterella wadsworthensis 2_1_59BFAA: 0.000
  Desulfotomaculum sp. CYP1: 0.000
  Anaerostipes hadrus: -0.077
  Streptococcus thermophiles/salivarius: -0.094
  butyrate-producing bacterium SL7/1: -0.107
  Roseburia sp.1120: -0.177
  Eubacterium rectale ATCC 33656: -0.183
  Clostridium sp.: -0.196
  butyrate-producing bacterium A2-175: -0.236
  Eggerthella lenta: -0.288
  Faecalibacterium prausnitzii: -0.294
  Clostridium sp. ID5: -0.392
  Actinobacteria; Adlercreutzia: -0.486
  Clostridium sp. RT8: -0.490
  Bacteroides coprocola: -0.524
  Lactobacillus rogosae: -0.524
  Bacteroides stercoris: -0.562
  Bacteroidetes-Pedobacter: -0.584
  Bacteroides coprophilus: -0.912
  Prevotella copri DSM 18205: -0.914
  Megamonas funiformis YIT 11815: -1.131
  Bacteroidetes-Flavobacterium: -1.406

 
SYSTEM REQUIREMENTS
  Python 2.7.6
  Packages/Libraries required: import sys, operator, itemgetter, argparse


