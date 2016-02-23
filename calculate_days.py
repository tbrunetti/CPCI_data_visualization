from datetime import datetime
import sys

def main():
#sys.argv[1]=txt file with three columns tab delimited; col1=patID, col2=surgery date, col3=date of death	
	temp=[]
	with open(sys.argv[1]) as input:
		for line in input:
			data=line.split('\t')
			temp.append(data)
	
	print temp[0]
	print temp[1]
	#capital Y for 4-digit year; lowercase y for 2-digit year
	format="%m/%d/%Y"
	f=open('days.txt', 'w')
	#1 because counts header
	#-2 because last two indicies do not have patient data
	for x in range(1, len(temp)):
		#if no patient infomation is available, cell[1] will be left blank
		if temp[x][1]=='':
			f.write(str(temp[x][0])+'\t'+ ''+ '\t'+''+'\n')
		#if patient is alive, cell will have EOF character in cell [2]
		elif temp[x][2]=='\n':
			start=datetime.strptime(temp[x][1], format)
			#means still alive, so use today's date, confirmed by Andrew Schneider
			end=datetime.strptime("10/13/2015", format)
			#.days gets days only, not hours, min, sec, etc...
			f.write(str(temp[x][0])+'\t'+str((end-start).days)+'\t'+'0'+'\n')
		#if enters into else, means patient is deceased
		else:
			removeEOF=temp[x][2].strip('\n')
			start=datetime.strptime(temp[x][1], format)
			end=datetime.strptime(removeEOF, format)
			f.write(str(temp[x][0])+'\t'+str((end-start).days)+'\t'+'1'+'\n')
	
if __name__=='__main__':
	main();