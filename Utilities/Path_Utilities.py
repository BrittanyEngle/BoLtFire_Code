def getRootFolder():
    return "<Root Folder>"

#-------------------ENTLN PATHS-----------------------------------
def getENTLN_Root():
    return f"{getRootFolder()}ENTLN/"

def getENTLN_CSVDataset():
    return f'{getENTLN_Root()}Downloads/'

def getENTLN_Boreal_Filter():
    return f"{getENTLN_Root()}Boreal_Filter/"

def getENTLN_Continent_Split():
    return f"{getENTLN_Root()}Continent_Split/"

def getENTLN_Final_Filter():
    return f"{getENTLN_Root()}Final_Filter/"


#-------------------GROUNDTRUTH PATHS------------------------------
def getGroundtruth_Root():
    return f"{getRootFolder()}Groundtruth/"

def getGroundtruth_Alaska_Points():
    return f"{getGroundtruth_Root()}/Alaska/AlaskaFireHistory_Points_Original.shp"

def getGroundtruth_Alaska_Polygons():
    return f"{getGroundtruth_Root()}/Alaska/AlaskaFireHistory_Polygons_Original.shp"

def getGroundtruth_Canada_Points():
    return f"{getGroundtruth_Root()}/Canada/NFDB_point_20240409.shp"

def getGroundtruth_Canada_Polygons():
    return f"{getGroundtruth_Root()}/Canada/NFDB_poly_20210707.shp"

def getGroundtruth_Processed():
    return f"{getGroundtruth_Root()}Processed_new/"

def getGroundtruth_Year_Split():
    return f"{getGroundtruth_Root()}Year_Split_new/"


##-------------------MODIS PATHS------------------------------
def getMODIS_Landcover():
    return f"{getRootFolder()}MODIS_Landcover/"

def getGWIS_Main_Dataset():
    return f"{getRootFolder()}GlobFire/"

def getGWIS_Processed():
    return f"{getGWIS_Main_Dataset()}Processed_new/"

#-------------------MATCHING PATHS------------------------------
def getMatchingCasePath(groundtruth:bool, day_filter:int, distance_filter:int):
    output = f"{getRootFolder()}Matching_new/"
    if groundtruth:
        output +=f"Groundtruth/{day_filter}days_{distance_filter}km/"
    else:
        output += f"GWIS_new/{day_filter}days_{distance_filter}km/"
    return output

def getMatchingFilterPath(day_filter:int, distance_filter:int,  ha_filter:int, geoshape_filter:int):
    output = getMatchingCasePath(True,day_filter,distance_filter)
    output +=f"{ha_filter}ha/{geoshape_filter}/"
    return output

#-------------------PUBLIC DATASET PATHS------------------------------
def getPublicDatasetPath():
    return f"{getRootFolder()}BoLtFire_Public_Dataset/"
