
import os
import geopandas as gpd
import glob
from Utilities.Path_Utilities import getRootFolder,getPublicDatasetPath
import pandas as pd

def getUSDateRange():
    start_date = '01/01/2012'
    end_date = '12/31/2022'
    return start_date, end_date

def Export_Shapefile(filtered, output_folder, output_filename):
    #check for invalid polygons, buffer them by 0 to remove the issue (works with self inserted polygons)
    filtered["geometry"] = filtered["geometry"]
   # If output_folder not there, create it
    os.makedirs(output_folder, exist_ok=True)
    # Path for the output shapefile
    output_shapefile_path = os.path.join(output_folder, output_filename + ".shp")
    # Export to shapefile
    filtered.to_file(output_shapefile_path, driver='ESRI Shapefile')
    print(f"Filtered file exported to: {output_shapefile_path}")
   

def loadGlanceProjectionZones():
    continents = []
    for continent_shape in glob.glob(f"{getRootFolder()}Glance_Grid/*.shp"):
        continent = gpd.read_file(continent_shape).to_crs(4326)
        continent = continent[["geometry"]]
        name = os.path.basename(continent_shape).split("_")[0]
        continent["name"] = name
        continents.append(continent)
    continents = gpd.GeoDataFrame(pd.concat(continents),crs=4326)
    return continents

def loadBiomes(names=["Boreal Forests/Taiga"]):
    # Ecoregions
    Ecoregions = gpd.read_file(f"{getRootFolder()}Biomes/Ecoregions2017.shp").to_crs(4326)
    # Select only the necessary columns from Ecoregions dataset
    Ecoregions_subset = gpd.GeoDataFrame(Ecoregions[['BIOME_NAME', 'ECO_BIOME_', 'REALM', "ECO_NAME","ECO_ID", "geometry"]])
    # Filter the Ecoregions_subset to include only "Boreal Forests/Taiga" or "Tundra"
    Ecoregions_subset = Ecoregions[Ecoregions['BIOME_NAME'].isin(names)]
    return Ecoregions_subset

def CountryBoundaries():
    countries = gpd.read_file(f"{getRootFolder()}CountryBoundaries/WB_countries_Admin0_10m.shp").to_crs(4326)
    countries["country"] = countries["NAME_EN"]
    for i, row in countries.iterrows():
         if row["country"] == None:
              countries.loc[i,"country"] = row["NAME_2"]
    countries = countries[["country","geometry"]]
    return countries

    
def getGlanceCRS(name):
    if name == "AF":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - AF - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",20],PARAMETER["latitude_of_center",5],UNIT["meter",1.0]]'
    elif name == "AN":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - AN - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",0],PARAMETER["latitude_of_center",-90],UNIT["meter",1.0]]'
    elif name == "AS":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - AS - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",100],PARAMETER["latitude_of_center",45],UNIT["meter",1.0]]'
    elif name == "EU":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - EU - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",20],PARAMETER["latitude_of_center",55],UNIT["meter",1.0]]'
    elif name == "OC":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - OC - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",135],PARAMETER["latitude_of_center",-15],UNIT["meter",1.0]]'
    elif name == "NA":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - NA - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",-100],PARAMETER["latitude_of_center",50],UNIT["meter",1.0]]'
    elif name == "SA":
         return 'PROJCS["BU MEaSUREs Lambert Azimuthal Equal Area - SA - V01",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["longitude_of_center",-60],PARAMETER["latitude_of_center",-15],UNIT["meter",1.0]]'    

def load_shapefile(file_path, expected_crs):
    """
    Loads a shapefile and applies the expected CRS.
    """
    # Determine the CRS code
    if isinstance(expected_crs, str) and len(expected_crs) == 2:  # Continent abbreviation
        crs_code = getGlanceCRS(expected_crs)
    elif isinstance(expected_crs, int):  # Integer EPSG code
        crs_code = f'EPSG:{expected_crs}'
    else:  # Assume it's a valid CRS string
        crs_code = expected_crs
    
    # Load the shapefile
    gdf = gpd.read_file(file_path)
    
    # Apply the CRS
    gdf = gdf.to_crs(crs_code)
    gdf["geometry"] = gdf["geometry"]
    return gdf


def generate_public_dataset(root_folder):
    short_to_long_continent = {
        "NA":"North_America",
        "AS":"Asia",
        "EU":"Europe"
    }
    output_folder = getPublicDatasetPath()
    columns_to_rename = {
    'id': 'FireID',
    'initialdat': 'StartDate',
    'finaldate':"EndDate",
    "fire_year":"FireYear",
    "area_ha":"AreaHa",
    "fire_size":"ClassSize",
    "BIOME_NAME":"BiomeName",
    "ECO_BIOME_":"EcoBiome",
    "ECO_NAME":"EcoName",
    "ECO_ID":"EcoID",
    "REALM":"Realm",
    "LC_DN":"LCDN",
    "LC_Name":"LCName",
    "country":"Country",
    "Continent":"Continent",
    "Holdover_d":"HoldoverD",
    "Holdoverrd":"HoldoverRD",
    "t_lat":"IgnLat",
    "t_long":"IgnLong",
    "t_pol_dist":"DisPol",
    "t_Spatial":"PerCheck"
    }
    filter_column = 'AreaHa' 
    min_value = 200  # The minimum value for filtering
    required_columns = ['HoldoverD']

    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all files in the root folder
    for file_path in glob.glob(root_folder+"*Merged_Results_*_Forest.shp"):
        continent = short_to_long_continent[file_path.split("_")[4]]
        # Load the shapefile
        gdf = gpd.read_file(file_path)
        
        # Rename specified columns
        gdf = gdf.rename(columns=columns_to_rename)
        
        # Keep only renamed columns and the geometry
        columns_to_keep = list(columns_to_rename.values()) + ['geometry']
        gdf = gdf[columns_to_keep]
        
        # Filter rows where the value in 'filter_column' is >= 'min_value'
        gdf = gdf[gdf[filter_column] >= min_value]
        
        # Filter out rows where any of the 'required_columns' have missing or null values
        if required_columns:
            gdf = gdf.dropna(subset=required_columns)
        
        # Save the modified shapefile to the output folder
        filename = f"BoLtFire_{continent}"
        output_path = os.path.join(output_folder, filename)
        gdf.to_file(output_path)