

import os


fitPlotList = ["Dibosonlnujj_muon_Stacked.png", "Dibosonlnujj_muon_Pull.png", "Dibosonlnujj_muon_Subtracted.png", "Dibosonlnujj_electron_Stacked.png", "Dibosonlnujj_electron_Pull.png", "Dibosonlnujj_electron_Subtracted.png", "DibosonBtaglnujj_muon_Stacked.png", "DibosonBtaglnujj_muon_Pull.png", "DibosonBtaglnujj_muon_Subtracted.png", "DibosonBtaglnujj_electron_Stacked.png", "DibosonBtaglnujj_electron_Pull.png", "DibosonBtaglnujj_electron_Subtracted.png"]


## strToInsert = "WpJAlt_"
## or 
## strToInsert = "WpJSherpa_"


def insertStringIntoSpecializedFitResultPlotNames(strToInsert):
    for fileName in fitPlotList:
	if fileName.find('muon') >= 0:
		newName=fileName.replace('muon_','muon_'+strToInsert)
	if fileName.find('electron') >=0:
		newName=fileName.replace('electron_','electron_'+strToInsert)
	os.system("mv "+fileName+" "+newName)
    return

