import os

def main():
	minReads=30
	f=open('offtarget-variant-calls-and-low-read-SNVs.txt', 'w')

	for allFiles in os.listdir('/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/txt_format'):
		os.chdir('/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/txt_format')
		
		with open(allFiles) as input:
			if allFiles=='format_UC_ret_DNA.py' or allFiles=='decisionTree_input_UC.txt' or allFiles=='input_into_violin_plot_UC_patientID_numUniqueMuts_gender.txt':
				#skips this file and moves to next file iteration
				continue;

			#creates a new file for each file corresponding to a variant called	
			f2=open(str(allFiles[:len(allFiles)-4]+'-offtargets-low-reads-SNVs-removed.txt'), 'w')
			#the temp array adds data from individual file line by line so it can be easily scanned
			#temp directory is cleared out each time a new file is read
			temp=[]
			for line in input:
				data=line.split('\t')
				temp.append(data)
#-------------------commment out the chunk of code below if no headers are wanted---------
			
			for headers in range(0, len(temp[0])-1):
				f2.write(temp[0][headers]+'\t')
			f2.write(temp[0][len(temp[0])-1])

#-----------------------------------end of headers code------------------------------------
			
			for x in range(1, len(temp)):
				if temp[x][20].rstrip('\n')=='ON' and (int(temp[x][4])+int(temp[x][5])>=minReads) and (int(temp[x][7])+int(temp[x][8])>=minReads):
					print int(temp[x][4])+int(temp[x][5])
					print int(temp[x][7])+int(temp[x][8])
					for w in range(0, len(temp[0])-1):
						f2.write(str(temp[x][w])+'\t')
					f2.write(str(temp[x][len(temp[0])-1]))
					print len(temp[0])
				elif temp[x][20].rstrip('\n')=='OFF' or (int(temp[x][4])+int(temp[x][5])<minReads) or (int(temp[x][7])+int(temp[x][8])<minReads):
					for i in range(0, len(temp[0])-1):
						f.write(str(temp[x][i])+'\t')
					f.write(str(temp[x][len(temp[0])-1]))
				else:
					print 'check formatting at '+ +str(allFiles)+' ' +str(temp[x])
			

if __name__=='__main__':
	main();