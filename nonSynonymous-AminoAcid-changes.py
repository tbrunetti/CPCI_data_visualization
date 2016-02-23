import matplotlib.pyplot as plt
import os
import sys
import pandas
from collections import Counter 
import seaborn as sns
from collections import OrderedDict
import operator
import plotly
import plotly.plotly as pyplot
 
pyplot.sign_in('tbrunetti', '72r07ya8dr')

from plotly.graph_objs import Scatter
from plotly.graph_objs import Bar
from plotly.graph_objs import Layout
import plotly.graph_objs as go

#sys.argv[1]=path to UC retrospective clinical patient information (Andrew Schneider)
#sys.argv[2]=path to UC retrospective offtargets removed (all pancreatic cancers) text files properly formatted


#The code directly in main() is just data formatting
#No data analysis happen here, all analysis is defined in new functions
def main():
	noPDAC=[]
	#key=patientID
	#value=dataframe of all amino acid non-synonmous variants
	varsInAA={}
	#headers should be the order of the offtargets removed text files
	sequenceHeaders=['chromosome', 'position', 'ref', 'alt', 'ref_count', 'alt_count', 'pct_norm_alt', 'tumor_ref_count', 'tumor_alt_count', 'pct_tumor', 'tn_pct_alt_ratio', 'gene', 'context', 'dbsnp_id', 'effect', 'coding', 'codon_change', 'amino_acid_change', 'amino_acid_length', 'mutect_call', 'on/off-target']
	#a running list of all genes that have an amino acid change
	genesWithAminoAcidChanges=[]
	#stores the patient IDs of patients diagnosed with PDAC
	patientIDPDAC=[]
	#subset of genesWithAminoAcidChanges, only PDAC diagnosed patients
	genesWithAminoAcidChangesPDAConly=[]
	
	#stores Andrew Schnieders clinical data as a dataframe call patientInfo with headers provided in txt file
	temp=[]
	with open(sys.argv[1]) as input:
		for line in input:
			temp.append(line.split('\t'))
	clinicalHeaders=[temp[0][x] for x in range(len(temp[0]))]
	temp.pop(0)
	patientInfo=pandas.DataFrame(temp, columns=clinicalHeaders)

	for patientMutationFiles in os.listdir(sys.argv[2]):
		os.chdir(sys.argv[2])
		with open(patientMutationFiles) as input:
			#gets the patients ID associated with mutation file
			for x in range(0, len(patientMutationFiles)):
				if patientMutationFiles[x]=='-':
					patientID=patientMutationFiles[0:x]
					break;
			
			#stores patientID of those with PDAC
			if patientInfo[patientInfo['study number']==patientID]['Tumor type'].item()=='PDAC':
				patientIDPDAC.append(patientID)
			else:
				noPDAC.append(patientID)
	
			#only stores mutations that are associated with amino acid changes
			temp=[]
			for line in input:
				data=line.split('\t')
				if data[14]=='NON_SYNONYMOUS_CODING':
					temp.append(data)
					#adds the gene name to the list if an amino acid change has occurred
					genesWithAminoAcidChanges.append(data[11])
					#if gene variant results in a non-synonymous mutation check if PDAC
					if patientID in patientIDPDAC:
						genesWithAminoAcidChangesPDAConly.append(data[11])

			varsInAA[patientID]=pandas.DataFrame(temp, columns=sequenceHeaders)	
	
	#outputs some basic count and patient total information
	def basicInformation():
		#print 'The total number of amino acid mutations per gene across all pancreatic cancers seen are: '+str(Counter(genesWithAminoAcidChanges))+'\n'
		
		print 'The total different number of genes that are associated with a non-synonymous amino acid change across all pancreatic cancers are ' + str(len(set(genesWithAminoAcidChanges)))+'\n'

		print 'The top 50 most frequently mutated genes associated with a non-synonymous amino acid change across all pancreatic cancers are: '+ str(Counter(genesWithAminoAcidChanges).most_common(50))+'\n'

		print 'The total number of amino acid mutations per gene across PDAC patients are: '+str(Counter(genesWithAminoAcidChangesPDAConly))+'\n'

		print 'The total different number of genes that are associated with a non-synonymous amino acid change across PDACs are ' + str(len(set(genesWithAminoAcidChangesPDAConly)))+'\n'

		print 'The top 50 most frequently mutated genes associated with a non-synonymous amino acid change PDACs are: '+ str(Counter(genesWithAminoAcidChangesPDAConly).most_common(50))+'\n'

		print 'The total number of PDAC patients with sequences at UC retro are '+str(len(patientIDPDAC))
		print 'The total number of patients that are non-PDAC at UC retro are '+str(len(noPDAC))
		#print patientIDPDAC
		#print noPDAC

		#-------------------------graphical overview of most frequent gene mutations----------------------------------
		#only PDAC patients, each sublist stores the gene in the first index
		#and the frequency (# of hits, not patients) in the 2nd index
		mostCommonPDAC=[]
		#mostCommonPDACgeneOnly=just extracts the gene name, no fequency measure
		mostCommonPDACgeneOnly=[]
		headers=['gene', 'number of mutations']
		
		for gene, freqeuncy in Counter(genesWithAminoAcidChangesPDAConly).most_common(20):
			temp=[]
			temp.append(gene)
			mostCommonPDACgeneOnly.append(gene)
			temp.append(freqeuncy)
			mostCommonPDAC.append(temp)
		#commonDataFrame is essentially mostCommonPDAC but converted into a convenient dataframe
		commonDataFrame=pandas.DataFrame(mostCommonPDAC, columns=headers)
		
		#this graph will show the genes with the most frequent nonsynon hits (not most frequently occurring in patients) 
		#in PDAC patients
		sns.set(font_scale=2)
		sns.barplot(x='number of mutations', y='gene', data=commonDataFrame)
		sns.plt.title('Genes with most non-synonymous mutations', fontsize=25)
		sns.plt.show()

		#mostCommon is a list of lists, each length 2
		#mostCommon stores the most frequenly nonsyn hit genes across all patients, NOT JUST PDAC!
		#also this list is not gene mutations most frequenly seen in patients just most frequently hit/mutated
		mostCommon=[]
		#mostCommonGeneOnly-is just a list of the gene names, no frequencies in the mostCommon list
		mostCommonGeneOnly=[]
		for gene, frequency in Counter(genesWithAminoAcidChanges).most_common(20):
			temp=[]
			temp.append(gene)
			mostCommonGeneOnly.append(gene)
			temp.append(frequency)
			mostCommon.append(temp)

		#acrossAllDataFrame=essentially mostCommon list but converted into convenient dataframe
		acrossAllDataFrame=pandas.DataFrame(mostCommon, columns=headers)
		#makes a bar graph of the most frequently hit/mutated genes across ALL patients
		#note: as stated before, not most frequently seen among patients, just most nonsyn mutations/hits seen
		sns.set(font_scale=2)
		sns.barplot(x='number of mutations', y='gene', data=acrossAllDataFrame)
		sns.plt.title('Genes with most non-synonymous mutations', fontsize=25)
		sns.plt.show()

		
		return mostCommonPDACgeneOnly, mostCommonGeneOnly



	#analyzes particular variant in gene of interest
	def aminoacidAnalysisPDAC(mostCommonPDACgeneOnly):
		#key: patientID, value: pandas dataframe of filtered variant viewer MUTECt nonsyn calls
		noPDACInfo={}
		#removes patients in varsInAA that are not PDAC
		#NOTE! after this loop, varsInAA becomes a PDAC only set of patients
		#To get the patients w/o PDAC, use the dictionary noPDACInfo
		for patients in noPDAC:
			if patients in varsInAA:
				noPDACInfo[patients]=varsInAA[patients]
				del varsInAA[patients]
		
	#---------------Overview of most common PDAC genes-------------------------------------------------
		#stores pateint IDs of most commonly mutated PDAC genes
		#key=gene name (one of the most commonly mutated in PDAC (again, most nonsyn hit/mutated, not fequency among patients))
		#value=list of patient IDs of PDAC patients only! with a nonsyn mutation in given gene (key)
		commonGeneWithPatientsIDs={}

		for x in range(0, len(mostCommonPDACgeneOnly)):
			commonGeneWithPatientsIDs[mostCommonPDACgeneOnly[x]]=[]
			#remember varsInAA now is just PDAC patients, removed nonPDACs from this dictionary above
			for key in varsInAA:
				if len(list(varsInAA[key].loc[varsInAA[key].gene==mostCommonPDACgeneOnly[x]]['gene']))>0:
					commonGeneWithPatientsIDs[mostCommonPDACgeneOnly[x]]=commonGeneWithPatientsIDs[mostCommonPDACgeneOnly[x]]+[key]
		
		#frequencyPatWithMutation=stores percent of PDAC patients with at least one nonsyn mutation (list of lists)
		#each sublist representative of one gene
		#list index[0]: gene name
		#list index[1]: a float, representing percent of PDAC patients that exhibit at least one nonsyn mutation in this gene
		#list index[2]: type of pancancer, in this case, all should be PDAC
		frequencyPatWithMutation=[[key,float(len(commonGeneWithPatientsIDs[key]))/len(varsInAA),'PDAC']for key in commonGeneWithPatientsIDs]
		#a sorted frequencyPatWithMutation by percent value, converted into a data frame
		commonPDACdataframe=pandas.DataFrame(sorted(frequencyPatWithMutation, key=operator.itemgetter(1), reverse=True))
		
		
		#commonGeneWithNOPDACIDs=same as commonGeneWithPatientIDs except patients DO NOT HAVE PDAC!
		#key: gene name
		#value:list of patients IDs, WITHOUT PDAC, that have at least one nonsyn mutation in gene
		commonGeneWithNoPDACIDs={}
		for i in range(0, len(mostCommonPDACgeneOnly)):
			commonGeneWithNoPDACIDs[mostCommonPDACgeneOnly[i]]=[]
			for key in noPDACInfo:
				if len(list(noPDACInfo[key].loc[noPDACInfo[key].gene==mostCommonPDACgeneOnly[i]]['gene']))>0:
						commonGeneWithNoPDACIDs[mostCommonPDACgeneOnly[i]]=commonGeneWithNoPDACIDs[mostCommonPDACgeneOnly[i]]+[key]

		#frequencyNoPDACmutations is the same as frequencyPatWithMutation, except in NON-PDAC patients
		#each sublist representative of one gene
		#list index[0]: gene name
		#list index[1]: a float, representing percent of NON-PDAC patients that exhibit at least one nonsyn mutation in this gene
		#list index[2]: type of pancancer, in this case, all should be No PDAC
		frequencyNoPDACmutations=[[key, float(len(commonGeneWithNoPDACIDs[key]))/len(noPDACInfo), 'No PDAC'] for key in commonGeneWithNoPDACIDs]
		
		#finalList, is sorted by highest percent and is a concatenation of list fequencyPatWithMutation
		#and frequencyNoPDACmutations, which is why list index[2] is critical in order to compare
		#PDAC with non-PDAC nonsyn patient mutations
		finalList=sorted(frequencyPatWithMutation+frequencyNoPDACmutations, key=operator.itemgetter(1), reverse=True)
		
		#converst finalList from above into a convenient to use dataframe
		finalDataframe=pandas.DataFrame(finalList, columns=['gene', 'percent of patients', 'diagnosis'])
		#graphs (bargraph) the percent of patients that have at least one mutation in given gene and compares it to PDAC
		#patients versus non-PDAC patients
		sns.set(font_scale=0.5)
		sns.barplot(x='gene', y='percent of patients', hue='diagnosis', data=finalDataframe, palette=['midnightblue', 'orange'])
		sns.plt.title('Most frequent PDAC non-synonymous mutations', fontsize=20)
		sns.plt.show()

	#---------------end of overview-------------------------------------------------------------------		
		
	#**this function requires the basicInformation() is run prior to calling this function**
	def specificAminoAcidPDACanalysis(mostCommonPDACgeneOnly):
		#variant variable should be updated to variant of interest
		variant='U2AF1'
		#codon_changes is a list of lists, each list represents a patient, and the values in the list
		#are a list of all the nonsyn mutations that patient has of the given gene above
		#NOTE: one cannot derive with particular patient is responsible for the list (see patientIDwithVariant variable)
		codon_changes=[]

		#removes patients in varsInAA that are not PDAC
		#NOTE: varsInAA now only consists of PDAC ONLY patients for duration of function
		for patients in noPDAC:
			if patients in varsInAA:
				del varsInAA[patients]
		
		#a list of patient IDs that have at least one nonsyn mutation for the variant gene listed above
		patientIDwithVariant=[]
		for key in varsInAA:
			codon_changes.append(list(varsInAA[key].loc[varsInAA[key].gene==variant]['amino_acid_change']))
			if len(list(varsInAA[key].loc[varsInAA[key].gene==variant]['amino_acid_change']))>0:
				patientIDwithVariant.append(key)

		#patientsWithMuts removes lists that are empty, each list should have the nonsyn codon changes listed for a particular patientID
		#as a check, the len(patientsWithMuts)==len(patientIDwithVariant))	
		patientsWithMuts=[codon_changes[x] for x in range(len(codon_changes)) if len(codon_changes[x])!=0]
		print patientsWithMuts
		print 'The number of PDAC patients with at least one variant/mutation in '+str(variant)+' is '+ str(len(patientIDwithVariant))
	
		#removes codon_changes from list and combines all into a single list so can turn into dataframe
		amino_acid_changes=[]
		for x in range(0, len(patientsWithMuts)):
			for i in range(0, len(patientsWithMuts[x])):
				amino_acid_changes.append(patientsWithMuts[x][i])
		
		#ordered=Counts and orders amino_acid_changes from most frequent counts to least frequent counts
		#ordered is a dictionary where key=nonsyn mutation, value=number of occurennces of that mutation
		ordered=OrderedDict(sorted(Counter(amino_acid_changes).items(), key=lambda x: -x[1]))

		#extracts just the mutation name from ordered, so it can be passed through sns to
		#illustrate the order of mutations in graph in decreasing frequency
		order_in_graph=[]
		for i in ordered:
			order_in_graph.append(i)
	
		
		#this will plot the number of patients that have a particular non-syn mutation of the
		#given variant variable gene listed at the beginning of function
		#countplot parameters
		headers=['non-synonymous amino acid changes']
		#converts amino_acid_changes into a convenient dataframe for graphing purposes,
		#although amino_acid_changes is not orders, it is ok, because when using countplot, we specify
		#the order in which to plot the mutations
		variant_dataframe=pandas.DataFrame(amino_acid_changes, columns=headers)
		sns.countplot(x='non-synonymous amino acid changes', data=variant_dataframe, order=order_in_graph)
		sns.plt.title('Variants in '+str(variant)+' across '+str(len(patientsWithMuts))+' PDAC patients')
		sns.plt.show()

		return variant, patientIDwithVariant, patientsWithMuts

#----------------end of specificAminoAcidPDACanalysis(mostCommonPDACgeneOnly) -----------------------------	

	#NOTE! Requires output from specificAminoAcidPDACAnalysis(mostCommonPDACgeneOnly)
	def specificNucleotidePDACanalysis(variant, patientIDwithVariant, patientsWithMuts):

		#NOTE! After this loop varsInAA is ONLY PDAC patients
		for patients in noPDAC:
			if patients in varsInAA:
				del varsInAA[patients]

		#for the variant/gene selected, it stores the nucleotide position and allele mutation
		#key=nonsyn mutation
		#value=list of lists, each sublist corresponds to one nonsyn mutation; sublist are ordered as follows:
		#sublist index[0]=position
		#sublist index[1]=ref allele/nt
		#sublist index[2]=alt allele/nt
		nucleotideVariants={}
		
		#patientIDwithVariant are the IDs of the patients that have at least one nonsyn mutation at given gene
		#Therefore the len(patientIDwithVariant) should equal the total number of patients that have a nonsyn mutation
		#at the given variant/gene of interest
		for i in patientIDwithVariant:
			nonsynMutations=list(varsInAA[i].loc[varsInAA[i].gene==variant]['amino_acid_change'])		
			for x in nonsynMutations:
				if x in nucleotideVariants:
					nucleotideVariants[x]=nucleotideVariants[x]+([list(varsInAA[i].loc[varsInAA[i].amino_acid_change==x]['position'])+list(varsInAA[i].loc[varsInAA[i].amino_acid_change==x]['ref'])+list(varsInAA[i].loc[varsInAA[i].amino_acid_change==x]['alt'])])
				else:
					nucleotideVariants[x]=[list(varsInAA[i].loc[varsInAA[i].amino_acid_change==x]['position'])+list(varsInAA[i].loc[varsInAA[i].amino_acid_change==x]['ref'])+list(varsInAA[i].loc[varsInAA[i].amino_acid_change==x]['alt'])]
			
		
		##len(patientIDwithVariant) should be total number of patients with at least one mutation in gene
		keys=[key for key in nucleotideVariants]
		#print keys
		nucChanges=[]
		for key in nucleotideVariants:
			temp=[]
			for nt in range(0, len(nucleotideVariants[key])):
				temp.append(nucleotideVariants[key][nt][1]+nucleotideVariants[key][nt][2])
			nucChanges.append(temp)
		print nucChanges
		print nucleotideVariants

		#each plotly graph needs to modified
		#DOES NOT UPDATE WITH VARIANT, MUST DO MANUAL
		#temp.html file must also be removed when running next iteration
		plotly.offline.plot({
		"data":[{"direction": "clockwise",
		"domain":{"x":[0.25, 0.75], "y":[0, 1]},
		#"hole":0.83,
		"insidetextfont":{"size":15},
		"labels":["G -> A or C -> T", "A -> G or T -> C"],
		"marker":{"colors":['#FF6633', '#66CCCC']},
		"textposition":"inside",
		"type":"pie",
		"values":["7", "1"]
		},
		{
		"direction":"clockwise",
		"domain":{"x":[0.3, 0.7], "y":[0, 1]},
		"insidetextfont":{"size":14},
		"textinfo":"label+percent",
		"labels":['S158L (7)', 'Y114H (1)'],
		"marker":{"colors":['#B2B2B2', '#CCCCCC']},
		"type":"pie",
		"values":["7", "1"]
		}],

		"layout":Layout(
			title='Nucleotide variant distribution over '+ str(len(patientIDwithVariant))+ ' patients with a '+str(variant)+' variant')
		})

		
	#**this function requires the basicInformation() is run prior to calling this function**
	def panAAchangesTopPDACMutations(mostCommonPDACgeneOnly):
		#percentAminoAcidChangesInPDAC only gets made once, it is NEVER cleared out, auumulates throughout
		#duration of the function
		#stores a list of lists.  Each sublist is made of three components.
		#sublist index[0]=nonsyn mutation
		#sublist index[1]=percent of PDAC patients that have this particular nonsyn mutation
		#sublist index[2]=the gene/variant the mutation belongs to
		percentAminoAcidChangesInPDAC=[]
		#removes patients in varsInAA that are not PDAC
		#NOTE! after this loop, varsInAA becomes a PDAC only set of patients for remainder of this function
		for patients in noPDAC:
			if patients in varsInAA:
				del varsInAA[patients]
		
		#cycles through top 20 PDAC genes mutated
		for gene in range(0, len(mostCommonPDACgeneOnly)):
			#codon_changes resets every time a new variant is selected
			#codon_changes is a list of lists, each sublist represents a single patients nonsyn mutations
			#exhibited in the gene variant of interest
			codon_changes=[]
			variant=mostCommonPDACgeneOnly[gene]

			#patientIDwithVariant is reset everytime a new variant is selected
			#patientIDwithVariant is a list of patient IDs of PDAC only that have at least one nonsyn mutation of given variant
			patientIDwithVariant=[]
			for key in varsInAA:
				codon_changes.append(list(varsInAA[key].loc[varsInAA[key].gene==variant]['amino_acid_change']))
				if len(list(varsInAA[key].loc[varsInAA[key].gene==variant]['amino_acid_change']))>0:
					patientIDwithVariant.append(key)

			#removes any sublist in codon_changes that are empty as to not a false number for the total
			#number of patients that have at least one non-syn mutation in the given gene variant of interest
			patientsWithMuts=[codon_changes[x] for x in range(len(codon_changes)) if len(codon_changes[x])!=0]

			#amino_acid_changes is a single list of the patientsWithMuts
			#essentially it removes all sublisting to create one continuous list of nonsyn mutations
			#amino_acid_changes clears out and resets every time a new variant/gene is selected
			amino_acid_changes=[]
			for x in range(0, len(patientsWithMuts)):
				#the i loop is critical because if a patient has mutliple nonsyn mutations in the gene
				#it extracts all of them, not just one
				for i in range(0, len(patientsWithMuts[x])):
					amino_acid_changes.append(patientsWithMuts[x][i])
			

			#Counter function will now count all the different nonsyn mutations with their frequencies
			#and stores like a dictionary
			#key=nonsyn mutation
			#value=number of times mutation is seen
			for key in Counter(amino_acid_changes):
				temp=[]
				temp.append(key)
				temp.append(float(Counter(amino_acid_changes)[key])/len(patientsWithMuts))
				temp.append(variant)
				percentAminoAcidChangesInPDAC.append(temp)

	#----strip/dot plot to show the differnt mutations each patient has given gene---------------------	
		
		#converts the list of lists percentAminoAcidChangesInPDAC into a convenient to use dataframe
		#for more information on percentAminoAcidChangesInPDAC, so to beginning of function to see explanation		
		percentAAdataframe=pandas.DataFrame(percentAminoAcidChangesInPDAC, columns=['mutation', 'frequency', 'gene'])
		
		#uses a strip plot of plot each gene/variant with the percent of PDAC patients that show each
		#nonsyn mutation
		sns.set(font_scale=1.0)
		#sns.set_style("whitegrid")
		sns.factorplot(x='gene', y='frequency', hue='mutation', data=percentAAdataframe, kind='strip', jitter=0.05, edgecolor="gray", legend=False)
		sns.plt.ylim(0, 1.05)
		sns.plt.title("Percent same mutation across PDAC patients", fontsize=30)
		sns.plt.show()

	
#-----------------------------end of panAAchangesTopPDACMutations(mostCommonPDACgeneOnly)--------------------------------





	


	#-----------------------------------calls to functions--------------------------------------------
	
	mostCommonPDACgeneOnly, mostCommonGeneOnly=basicInformation();
	
	#-------------------in order to run the calls below, basicInformation() must be called first-------------------
	#---------------------------------output of basicInformation() is used ------------------------------------------
	#--------------------------continue with funciton calls after basicInformation() is called first---------------
	
	#aminoacidAnalysisPDAC(mostCommonPDACgeneOnly);
	#panAAchangesTopPDACMutations(mostCommonPDACgeneOnly);
	#variant, patientIDwithVariant, patientsWithMuts=specificAminoAcidPDACanalysis(mostCommonPDACgeneOnly);

	#-------------------------------in order to fun specificNucleotidePDACanalysis------------------------------
	#-------------------specificAminoAcidPDACanalysis(mostCommonPDACgeneOnly) must be called first--------------
	#---input of specificNucleotidePDACanalysis is dependent on output of specificAminoAcidPDACanlaysis---------
	
	variant, patientIDwithVariant, patientsWithMuts=specificAminoAcidPDACanalysis(mostCommonPDACgeneOnly);
	specificNucleotidePDACanalysis(variant, patientIDwithVariant, patientsWithMuts);


if __name__=='__main__':
	main();