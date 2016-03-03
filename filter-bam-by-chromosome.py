import sys
import os
from subprocess import call


#sys.argv[1] is a file with a list of chromosomes to keep in a bam (chr*** format)
#the list needs to be formated with one chromosome per line
#sys.argv[2] is the original bam file that is to be filetered out

def filterBam(pathToSamtools, pathToBamFile):
	#appends results for each chromosome
	output=open('newBam_chromsomeSpecific_'+str(sys.argv[2]), 'a')
	#change directory into samtools working directory
	os.chdir(pathToSamtools)
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.rstrip('\n')
			#call samtools from samtools directory
			call(['samtools', 'view', '-b', str(pathToBamFile)+str(sys.argv[2]), str(line)], stdout=output)

if __name__=='__main__':
	pathToSamtools='/home/ubuntu/TOOLS/samtools-0.1.19'
	pathToBamFile='/mnt/'
	filterBam(pathToSamtools, pathToBamFile);
