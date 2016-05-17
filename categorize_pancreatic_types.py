#note this will not open password protected files, so temporarily change to no password
import xlrd
import sys
from openpyxl import load_workbook
from openpyxl import Workbook
from openpyxl.cell import get_column_letter, column_index_from_string

#first line of clincal data is header
def categorize(pathToClinicalData):
	wb=load_workbook(pathToClinicalData)
	workingSheet=wb.get_sheet_by_name("Sheet1")
	
	#will get index number for column 'AK' which is the column of diagnosis
	typeIndex=column_index_from_string('AK')

	for i in range(1, 233):
		print (i, workingSheet.cell(row=i, column=typeIndex).value)

if __name__=='__main__':
	pathToClinicalData='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/Retrospective_panData_UC/RetroPancStudy_PathDatabase1_AS_newdeath_noPass_noScheduleCol.xlsx'
	categorize(pathToClinicalData);