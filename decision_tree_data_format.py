import sys

#argv[1]=txt file of list of genes want to look at for decision tree
#argv[2]=tab delimited file of clinical data for TCGA in TCGA format
#argv[3]=tab delimited file of gene variant data with patientID

def main():

	#a list of genes to analyze in decision tree
	geneOrder=[]
	with open(sys.argv[1]) as input:
		for line in input:
			geneOrder.append(line.rstrip('\n'))
	
	#stores all clincinal patient data
	patientID=[]
	with open(sys.argv[2]) as input:
		for line in input:
			data=line.split('\t')
			patientID.append(data)

	#stores all the data from sequencing
	geneVariants=[]
	with open(sys.argv[3]) as input:
		for line in input:
			data=line.split('\t')
			geneVariants.append(data)

		
	#/////////////////////////////////////////////////////////////////////////////////#
	#Goal: match patientIDs to individual patient outcome 							  #
	#Input:  None; utiziles files read from [patientID]								  #
	#Output: returns outcome, a dictionary of key as patient ID and 				  #
	#		 value=[key][0]=days alive/until death [key][1]:0=alive; 1=deceased       #
	#/////////////////////////////////////////////////////////////////////////////////#
	
	def patientOutcome():

		outcome={}
		for x in range(2, len(patientID)-1):
			#[35] stores days alive; [36] stores days until deceased
			if patientID[x][35]=='[Not Available]':
				outcome[patientID[x][1]]=[patientID[x][36]]+[int(1)]
			else:
				outcome[patientID[x][1]]=[patientID[x][35]]+[int(0)]

		return outcome

	
	#/////////////////////////////////////////////////////////////////////////////////#
	#Goal:  matches patient ID with gene Variants seen in each individual 			  #
	#Input:  None; utilzes file read from [geneVariants]                              #
	#Output:  returns variants, a dictionary of key as patient ID and                 #
	#         value=list of genes that show a variant in at least one allele		  #  
	#/////////////////////////////////////////////////////////////////////////////////#
	
	def patientGenes():

		variants={}
		for x in range(2, len(geneVariants)-1):
			#checks if patient ID already exists
			if geneVariants[x][15][:12] in variants:
				variants[geneVariants[x][15][:12]]=variants[geneVariants[x][15][:12]]+[geneVariants[x][0]]
			#if patient ID does not exits, makes a new key in variants library
			else:
				variants[geneVariants[x][15][:12]]=[geneVariants[x][0]]	

		return variants

	#/////////////////////////////////////////////////////////////////////////////////#
	#Goal: checks if each patient carries on the of the gene variants from [geneOrder]#
	#Input: variants dictionary from def patientGenes and uses [geneOrder] file read  #
	#Output: returns variantCheck a dictionary key=patient ID and 					  #
	#        value=list of +'s and -'s +=has variant; -=does not have varian          #
	#/////////////////////////////////////////////////////////////////////////////////#

	def checkGene(variants):

		variantCheck={}
		for key in variants:
			#checks which genes in geneOrder are seen in patient gene variant profile
			for x in range(0, len(geneOrder)):
				#if the variant is found, checks if key exists in dictionary, if not, makes one
				if geneOrder[x] in variants[key]:
					if key in variantCheck:
						variantCheck[key]=variantCheck[key]+['+']
					else:
						variantCheck[key]=['+']
				#if variant is not found, checks if key exists in dictionary, if not, make one
				else:
					if key in variantCheck:
						variantCheck[key]=variantCheck[key]+['-']
					else:
						variantCheck[key]=['-']

		return variantCheck


	#/////////////////////////////////////////////////////////////////////////////////#
	#Goal: merges variant check with outcome for each patient 						  #
	#Input: variantCheck dictionary and outcome dictionary							  #
	#Output: returns variantsWithOutcomes, a dictionary where key=patientID and       #
	#        value=+/- list for variants in order of [geneOrder], number of days      #
	#        alive/deceased, and an int 1 or 0, 1=deceased; 0=alive                   #
	#/////////////////////////////////////////////////////////////////////////////////#

	def merge(variantCheck, outcome):

		variantsWithOutcomes={}
		for key in variantCheck:
			#only picks out patients that have both sequence and clinical information
			if key in outcome:
				variantsWithOutcomes[key]=variantCheck[key]+outcome[key]

		return variantsWithOutcomes

	#///////////////////////////////////////////////////////////////////////////////////////////////////////#
	#																										#
	#										Calls made to run defs 											#
  	#																										#
  	#///////////////////////////////////////////////////////////////////////////////////////////////////////#
  	
  	outcome=patientOutcome()
  	variants=patientGenes()
  	variantCheck=checkGene(variants)
  	variantsWithOutcomes=merge(variantCheck, outcome)

  	
  	#///////////////////////////////////////////////////////////////////////////////////////////////////////#
  	#																										#
  	#											Output Files                                                #
  	#																										#
  	#///////////////////////////////////////////////////////////////////////////////////////////////////////#

  	#creates output txt file with patient IDs, +/- for gene variant, survival days, and outcome
  	f=open('variantsWithOutcomes.txt', 'w')
  	f.write('patient ID'+'\t')
  	for x in range(0, len(geneOrder)):
  		f.write(str(geneOrder[x])+'\t')
	f.write('days'+'\t'+'outcome'+'\n')
  	for key in variantsWithOutcomes:
  		f.write(str(key)+'\t')
  		for i in range(0, len(geneOrder)):
  			f.write(str(variantsWithOutcomes[key][i])+'\t')
  		f.write(str(variantsWithOutcomes[key][len(geneOrder)])+'\t'+str(variantsWithOutcomes[key][len(geneOrder)+1])+'\n')


if __name__=='__main__':
	main();
