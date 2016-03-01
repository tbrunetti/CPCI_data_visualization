import sys
import numpy as np
from collections import Counter
import operator
#sys.argv[1]=output file called variantsWithOutcomes.txt which is tab-delimited
def main():
	tree=[]
	#stores information regarding patient outcome, if there is a mutation in gene of interest or not, patientID
	#patients[0]=headings
	patients=[]
	with open(sys.argv[1]) as input:
		for line in input:
			data=line.split('\t')
			patients.append(data)
	print len(patients[0])	

	#/////////////////////////////////////////////////////////////////////////////////////#
	#Goal: calculate overall entropy of root tree node(tree does not yet exist) or 		  #
	#      subtree root node when a tree already exists									  #
	#Input: [tree] uses global [patients] from read and store 							  #
	#Output: returns entropy a float of the entropy of root tree or subtree root node	  #
	#/////////////////////////////////////////////////////////////////////////////////////#

	def overallEntropy(tree, patients):
		if len(tree)==0:
			#list of 1's and 0's of patient outcome
			outcome=[patients[x][len(patients[0])-1].rstrip('\n') for x in range(1, len(patients))]
			print outcome
			percentAlive=float(Counter(outcome)['0'])/len(outcome)
			percentDeceased=float(Counter(outcome)['1'])/len(outcome)
			if percentAlive==0.0:
				entropy=(-1*percentDeceased*np.log2(percentDeceased))
			elif percentDeceased==0.0:
				entropy=(-1*percentAlive*np.log2(percentAlive))
			
			else:
				#note \ is a line continuation symbol
				entropy=(-1*percentAlive*np.log2(percentAlive)) \
					-(percentDeceased*np.log2(percentDeceased))
				#print entropy
				x=0
			
			return entropy, len(outcome), x

		else:
			#most recently added value tuples in tree
			for x in range(1, 3):
				if tree[len(tree)-1][x][0]=='leftLeaf' or tree[len(tree)-1][x][1]=='rightLeaf':
					continue

				else:
					outcome=[tree[len(tree)-1][x][i][len(patients[0])-1].rstrip('\n') for i in range(len(tree[len(tree)-1][x]))]
					percentAlive=float(Counter(outcome)['0'])/len(outcome)
					percentDeceased=float(Counter(outcome)['1'])/len(outcome)
					#note \ is a line continuation symbol
					entropy=(-1*percentAlive*np.log2(percentAlive)) \
					-(percentDeceased*np.log2(percentDeceased))
					print "about to hit break"
					break
			#print entropy
			return entropy, len(outcome), x


	#/////////////////////////////////////////////////////////////////////////////////////#
	#Goal: sorts attributes by if patient has mutation or not and the relative outcome 	  #
	#Input: float(entropy) from overallEntropy()										  #
	#Output: returns dictionary of attributes with their relative counts 				  #
	#/////////////////////////////////////////////////////////////////////////////////////#
	
	def featureCounts(patients):
		#dictionary keyed by attribute
		print 'feature counts' +str(patients[0])
		attributes={}
		#outer loop= # of attibutes + outcome data
		#inner loop= # of patients
		for x in range(1, len(patients[0])):
			#(0, 0)=(deceasedTotals, aliveTotals)
			attributeTotalsPos=(0, 0)
			attributeTotalsNeg=(0, 0)
			#in this loop for each attribute, counts the number of patients that
			#have the mutation versus wt and furthermore, it sub counts the outcome
			#of the patient for each mutation vs wt
			for i in range(1, len(patients)):
				#patients that have mutation x
				if patients[i][x]=='+':
					#patient is alive with mutation
					if patients[i][len(patients[0])-1]=='0\n':
						attributeTotalsPos=(attributeTotalsPos[0], attributeTotalsPos[1]+1)
					#patient is deceased with mutation
					else:
						attributeTotalsPos=(attributeTotalsPos[0]+1, attributeTotalsPos[1])
				#patients that do not have mutation x
				else:
					#patient is alive without mutation
					if patients[i][len(patients[0])-1]=='0\n':	
						attributeTotalsNeg=(attributeTotalsNeg[0], attributeTotalsNeg[1]+1)
					#patient is deceased without mutation
					else:
						attributeTotalsNeg=(attributeTotalsNeg[0]+1, attributeTotalsNeg[1])
			
			#key=mutations
			#value=[(deceased, alive), (deceased, alive)]
			#1st tuple=tests (+) for mutation; 2nd tuple=tests (-) for mutation
			attributes[patients[0][x]]=[attributeTotalsPos]+[attributeTotalsNeg]
		print attributes
		return attributes

	#/////////////////////////////////////////////////////////////////////////////////////#
	#Goal: find the feature/attribute that has the highest info gain
	#Input: float system entropy from overallEntropy() and dictionary of attibutes from   #
	#       featureCounts()
	#Output: 
	#/////////////////////////////////////////////////////////////////////////////////////#

	def gain(entropy, attributes, tree):
		gainValues={}
		for features in attributes:
	#LEFT SIDE OF TREE		
			#all that have the mutation
			firstTuple=attributes[features][0]
			#total sum of patients that have the mutation for the particular feature/attribute
			sumOfFirst=float(firstTuple[0]+firstTuple[1])
			#means nobody falls into this category so entropy will also be 0
			if sumOfFirst==0:
				entropyMut=float(0)
			else:
				deceasedMuts=firstTuple[0]/sumOfFirst
				aliveMuts=firstTuple[1]/sumOfFirst	
				#another check to make sure denominator of one component does not equal zero
				if (deceasedMuts==0 or aliveMuts==0):
					#LEAF REACHED
					if deceasedMuts>0:
						useForCalc=deceasedMuts
					else:
						useForCalc=aliveMuts
					entropyMut=(-1*useForCalc*np.log2(useForCalc))
				else:
					entropyMut=(-1*deceasedMuts*np.log2(deceasedMuts))-(aliveMuts*np.log2(aliveMuts))

	#RIGHT SIDE OF TREE
			#all that do not have mutation
			secondTuple=attributes[features][1]
			#total sum of patients that do not have the mutation for the particular feature/attribute
			sumOfSecond=float(secondTuple[0]+secondTuple[1])
			#means nobody falls into this category so entropy will also be 0
			if sumOfSecond==0:
				entropyNoMut=float(0)
			else:
				deceasedNoMut=secondTuple[0]/sumOfSecond
				aliveNoMut=secondTuple[1]/sumOfSecond
				#another check to make sure denominator of one component does not equal zero
				if (deceasedNoMut==0 or aliveNoMut==0):
					#LEAF REACHED
					if deceasedNoMut>0:
						useForCalc=deceasedNoMut
					else:
						useForCalc=aliveNoMut
					entropyNoMut=(-1*useForCalc*np.log2(useForCalc))
				else:	
					entropyNoMut=(-1*deceasedNoMut*np.log2(deceasedNoMut))-(aliveNoMut*np.log2(aliveNoMut))

			#key=feature/attribute
			#value=gain
			gainValues[features]=entropy-((sumOfFirst/totalPatients)*entropyMut)-((sumOfSecond/totalPatients)*entropyNoMut)
		highestGain=max(gainValues.iteritems(), key=operator.itemgetter(1))[0]
		#prints out the feature/attribute with the highest gain
		print 'The attribute with the highest gain is '+ str(highestGain)
		return highestGain, attributes[highestGain]

	
	def newPatientPopulation(highestGain, newPatientTotals, treeSide):
		if len(tree)!=0:
			tree[len(tree)-1][treeSide]=highestGain
		
		population=patients[0].index(highestGain)
		leftPop=[]
		rightPop=[]
		print newPatientTotals[0]
		for x in range(len(newPatientTotals)):
			if 0 in newPatientTotals[x]:
				if x==0:
					leftPop.append('leftLeaf')
					
					#stores the patient index to be removed
					#have to do this because poping will decrease patients list size
					#and screw up index numbers due to deletion prematurely
					indicesToRemove=[]
					for i in range(1, len(patients)):
						if patients[i][population]=='+':
							indicesToRemove.append(i)
					#deletes the patients that have reached a leaf
					for offset, index in enumerate(indicesToRemove):
						index-=offset
						del patients[index]	

				else:
					rightPop.append('rightLeaf')
					
					#stores the patient index to be removed
					#have to do this because poping will decrease patients list size
					#and screw up index numbers due to deletion prematurely
					indicesToRemove=[]
					for i in range(1, len(patients)):
						if patients[i][population]=='-':
							indicesToRemove.append(i)
					#deletes the patients that have reached a leaf
					for offset, index in enumerate(indicesToRemove):
						index-=offset
						del patients[index]	
			else:
				
				#removeFeature=patients[0].index(highestGain)

				#for i in range(0, len(patients)):
				#	patients[i].pop(removeFeature)
					
				for x in range(1, len(patients)):
					if patients[x][population]=='+':
						leftPop.append(patients[x])
					elif patients[x][population]=='-':
						rightPop.append(patients[x])

		tree.append([highestGain] + [leftPop, rightPop])
		print 'Patients '+str(len(patients))
		print tree[len(tree)-1][0]
		print tree[len(tree)-1][1]
		return tree, patients


	
	for x in range(0, 6):
		entropy, totalPatients, treeSide=overallEntropy(tree, patients);
		if treeSide==1:
			#this means that it has not reached a leaf node on left side of tree
			attributes=featureCounts(patients);
			highestGain, newPatientTotals=gain(entropy, attributes, tree)
			tree, patients=newPatientPopulation(highestGain, newPatientTotals, treeSide)
			entropy, totalPatients, treeSide=overallEntropy(tree, patients);
			print '1 happens'
		attributes=featureCounts(patients);
		highestGain, newPatientTotals=gain(entropy, attributes, tree)
		tree, patients=newPatientPopulation(highestGain, newPatientTotals, treeSide)
	print tree

	

if __name__=='__main__':
	main();