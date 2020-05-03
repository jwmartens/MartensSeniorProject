#Analysis of bacterial prevalence in disease
import sys
from operator import itemgetter
import argparse

def main():
	parser = argparse.ArgumentParser(description='Take commands for program')
        parser.add_argument('-i', required=True, type=str, help='CSV File with gene expression for host response analysis.')
        args = parser.parse_args()


        expressioninput = args.i
          
	IIS = []
	expressionWeight = 0
	MSarray = []
	expression = []
	counter = 0
	metabolitescore = 0


	bacterialprevalence(MSarray)
	geneexpression(expression, expressioninput)
	#print(expression)

	for bacteria in MSarray:
		metabolitescore =float(upregulated(bacteria[0])) + float(downregulated(bacteria[0]))
		if (float(bacteria[1]) < 0 and metabolitescore > 0):
			bacteria[1] = float(bacteria[1]) * metabolitescore * expression[0]
		elif (float(bacteria[1]) < 0 and metabolitescore < 0):
			bacteria[1] = abs(float(bacteria[1])) * metabolitescore * abs(expression[1])
		elif (metabolitescore == 0):
                        bacteria[1] = metabolitescore
		else:
			if (metabolitescore > 0):
				bacteria[1] = float(bacteria[1]) * metabolitescore * expression[0]
			elif (metabolitescore < 0):
				bacteria[1] = float(bacteria[1]) * metabolitescore * abs(expression[1])

	print("\nImmune Induction Score")
	
	output = sorted(MSarray, key=itemgetter(1), reverse=True)
	for sortedentry in output:
		print(sortedentry[0] + ": " + str('%.3f'%(sortedentry[1])))

	print("")

#Returns library of MS prevalence values for each bacterial species documented by populating the input array MSarray.
def bacterialprevalence(MSarray):
	MSdictionary = {
		"Bacteroidetes-Pedobacter":"2.74",
		"Methanobrevibacter":"2.64",
		"Bacteroidetes-Flavobacterium":"2.20",
		"Proteobacteria-Mycoplana":"2.17",
		"Proteobacteria-Caulobacteraceae":"2.17",
		"Akkermansia muciniphila":"1.46",
		"Proteobacteria-Pseudomonas":"0.98",
		"Acinetobacter calcoaceticus":"0.75",
		"Eggerthella lenta":"0.45",
		"Streptococcus thermophiles/salivarius":"0.44",
		"Faecalibacterium prausnitzii":"-0.23",
		"Anaerostipes hadrus":"-0.36",
		"Eubacterium rectale ATCC 33656":"-0.43",
		"Clostridium sp.":"-0.46",
		"butyrate-producing bacterium SL7/1":"-0.50",
		"Bacteroides stercoris":"-0.66",
		"Actinobacteria; Adlercreutzia":"-0.76",
		"Bacteroides coprocola":"-0.82",
		"Lactobacillus rogosae":"-0.82",
		"Roseburia sp.1120":"-0.83",
		"Sutterella wadsworthensis 2_1_59BFAA":"-0.86",
		"Clostridium sp. ID5":"-0.92",
		"Bacteroides coprophilus":"-1.07",
		"butyrate-producing bacterium A2-175":"-1.11",
		"Desulfotomaculum sp. CYP1":"-1.11",
		"Clostridium sp. RT8":"-1.15",
		"Prevotella copri DSM 18205":"-1.43",
		"Megamonas funiformis YIT 11815":"-1.77"
	}
	row = []
	for entry in MSdictionary:
		row =[]
		row.append(entry)
		row.append(MSdictionary[entry])
		MSarray.append(row)

	return MSarray

#Returns upregulated pathway numbers for each bacterial species as described here. Values are returned according to the bacterial species input 'bacteria' and consist of the name of the bacteria and the corresponding value. Add one for every literature-backed pathway influence on inflammatory response, subtract for anti-inflammatory.
def upregulated(bacteria):
	MetaboliteDictionary = {
      
		"Bacteroidetes-Pedobacter":"-1",
                "Methanobrevibacter":"1",
                "Bacteroidetes-Flavobacterium":"-3",
                "Proteobacteria-Mycoplana":"1",
                "Proteobacteria-Caulobacteraceae":"1",
                "Akkermansia muciniphila":"0",
                "Proteobacteria-Pseudomonas":"2",
                "Acinetobacter calcoaceticus":"1",
                "Eggerthella lenta":"-1",
                "Streptococcus thermophiles/salivarius":"1",
                "Faecalibacterium prausnitzii":"-3",
                "Anaerostipes hadrus":"-1",
                "Eubacterium rectale ATCC 33656":"-1",
                "Clostridium sp.":"-2",
                "butyrate-producing bacterium SL7/1":"-1",
                "Bacteroides stercoris":"-1",
                "Actinobacteria; Adlercreutzia":"0",
                "Bacteroides coprocola":"-3",
                "Lactobacillus rogosae":"-2",
                "Roseburia sp.1120":"-1",
                "Sutterella wadsworthensis 2_1_59BFAA":"0",
                "Clostridium sp. ID5":"-2",
                "Bacteroides coprophilus":"-1",
                "butyrate-producing bacterium A2-175":"-1",
                "Desulfotomaculum sp. CYP1":"0",
                "Clostridium sp. RT8":"-2",
                "Prevotella copri DSM 18205":"-2",
                "Megamonas funiformis YIT 11815":"-1"

	}	
        for entry in MetaboliteDictionary:
               	if entry == bacteria:
			return(MetaboliteDictionary[entry])


#Returns upregulated pathway numbers for each bacterial species as described here. Values are returned according to the bacterial species input 'bacteria' and consist of the name of the bacteria and the corresponding value. Add one for every literature-backed pathway influence on anti-inflammatory response, subtract for inflammatory.
def downregulated(bacteria):
        MetaboliteDictionary = {

                "Bacteroidetes-Pedobacter":"0",
                "Methanobrevibacter":"0",
                "Bacteroidetes-Flavobacterium":"0",
                "Proteobacteria-Mycoplana":"0",
                "Proteobacteria-Caulobacteraceae":"0",
                "Akkermansia muciniphila":"1",
                "Proteobacteria-Pseudomonas":"0",
                "Acinetobacter calcoaceticus":"1",
                "Eggerthella lenta":"-2",
                "Streptococcus thermophiles/salivarius":"-2",
                "Faecalibacterium prausnitzii":"-3",
                "Anaerostipes hadrus":"0",
                "Eubacterium rectale ATCC 33656":"-1",
                "Clostridium sp.":"0",
                "butyrate-producing bacterium SL7/1":"0",
                "Bacteroides stercoris":"-3",
                "Actinobacteria; Adlercreutzia":"-3",
                "Bacteroides coprocola":"0",
                "Lactobacillus rogosae":"-1",
                "Roseburia sp.1120":"0",
                "Sutterella wadsworthensis 2_1_59BFAA":"0",
                "Clostridium sp. ID5":"0",
                "Bacteroides coprophilus":"-3",
                "butyrate-producing bacterium A2-175":"0",
                "Desulfotomaculum sp. CYP1":"0",
                "Clostridium sp. RT8":"0",
                "Prevotella copri DSM 18205":"-1",
                "Megamonas funiformis YIT 11815":"-2"

        }
        for entry in MetaboliteDictionary:
                if entry == bacteria:
                        return(MetaboliteDictionary[entry])              		
#Returns array of the average of inflammatory and anti-inflammatory gene expression according to genes listed in the form of the 'expression' array. The 'expressioninput' input file taken must be in the form of a csv with the corresponding headings (blank), logFC, AveExpr, t, P.Value, adj.P.Val, B, and symbol.

def geneexpression(expression, expressioninput):
	lineArray = []
	Inftotal = 0.0
	Infnumber = 0.0
	Antitotal = 0.0
	Antinumber = 0.0

	InflammatoryGenes = {"IL6", "TRAF4", "NFKBIZ", "IL17C", "NFKB2"}
	AntiinflammatoryGenes = ("IL10RB", "STAT3", "JAK1", "JAK2", "TGFBR1")
	#MSexpression = open("Martens_MS_DE_result.csv", "r" )
	MSexpression = open(expressioninput, "r" )	
	for line in MSexpression:
		line = line[:-2]
		lineArray = line.split(',')
		if lineArray[7] in InflammatoryGenes:
			Inftotal += float(lineArray[1])
			Infnumber = Infnumber + 1

		elif lineArray[7] in AntiinflammatoryGenes:
			Antitotal += float(lineArray[1])
			Antinumber += 1				
	#Average gene expression for each pathway.
	expression.append(float(Inftotal) / float(Infnumber))
	expression.append(float(Antitotal) / float(Antinumber))
	return expression

if __name__ == "__main__":
	main()

