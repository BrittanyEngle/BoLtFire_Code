import geopandas as gpd
import pandas as pd
import datetime
import math
from Utilities.Path_Utilities import getGroundtruth_Alaska_Points,getGroundtruth_Alaska_Polygons, getGroundtruth_Processed
from Utilities.Most_Used_Functions import getUSDateRange, Export_Shapefile, getGlanceCRS


def alaska_filter_polygons():
    # Define the start and end year for filtering
    start_year = 2012
    end_year = 2022
    
    # Load the Alaska Fire History dataset
    Alaska_Dataset_polygons = gpd.read_file(getGroundtruth_Alaska_Polygons()).to_crs(4326)
    print("Alaska_Dataset CRS:", Alaska_Dataset_polygons.crs)
    print("Amount of fire cause available:", len(Alaska_Dataset_polygons))

    # Convert FIREYEAR to datetime, if not already done
    Alaska_Dataset_polygons["FIREYEAR"] = pd.to_datetime(Alaska_Dataset_polygons["FIREYEAR"], format='%Y', errors='coerce')

    # Filter by the defined year range
    Alaska_Dataset_polygons = Alaska_Dataset_polygons[(Alaska_Dataset_polygons["FIREYEAR"].dt.year >= start_year) & 
                                                      (Alaska_Dataset_polygons["FIREYEAR"].dt.year <= end_year)]

    # Convert FPOUTDATE to datetime
    Alaska_Dataset_polygons["FPOUTDATE"] = pd.to_datetime(Alaska_Dataset_polygons["FPOUTDATE"], errors="coerce") 

    # Further filter by start and end date (already filtered above, but keeping this if needed for additional filtering)
    Alaska_Dataset_polygons = Alaska_Dataset_polygons[(Alaska_Dataset_polygons["FPOUTDATE"].notnull())]

    print("Amount of filtered fires (date) available:", len(Alaska_Dataset_polygons))

    # Filter by FALSEALARM column
    Alaska_Dataset_polygons = Alaska_Dataset_polygons[Alaska_Dataset_polygons["PRESCRIBED"] != "Y"]

    # Define the conversion factor for area
    acres_to_hectares = 0.404686

    # Create a new column with the area in hectares
    Alaska_Dataset_polygons['EST_HA'] = Alaska_Dataset_polygons['ACRES'] * acres_to_hectares

    # Filter fires with an area greater than or equal to 100 hectares
    Alaska_Dataset_polygons = Alaska_Dataset_polygons[Alaska_Dataset_polygons['EST_HA'] >= 100]

    # Convert dates to string format for export
    Alaska_Dataset_polygons["FIREYEAR"] = Alaska_Dataset_polygons["FIREYEAR"].dt.year
    Alaska_Dataset_polygons["FPOUTDATE"] = Alaska_Dataset_polygons["FPOUTDATE"].dt.strftime("%Y-%m-%d")

    #--------------------------- Export Alaska 
    output_folder = getGroundtruth_Processed()
    output_filename = "Alaska_Polygons_Groundtruth"
    Export_Shapefile(Alaska_Dataset_polygons, output_folder, output_filename)

def alaska_filter_points():
    # Load date range
    start_date, end_date = getUSDateRange()
    start_date = datetime.datetime.strptime(start_date, "%m/%d/%Y")
    end_date = datetime.datetime.strptime(end_date, "%m/%d/%Y")

    # Load dataset and reproject
    Alaska_Dataset = gpd.read_file(f"{getGroundtruth_Alaska_Points()}").to_crs(4326)
    print("Alaska_Dataset CRS:", Alaska_Dataset.crs)

    # Filter by cause
    keywords = ["lightning", "natural", "wfu", "LIGHTNING", "Natural Out", "Natural","Natural Out"]
    Alaska_Dataset_Lightning = Alaska_Dataset[
        Alaska_Dataset.apply(lambda row: any(
            kw in str(row['FIRECAUSE']).lower() or 
            kw in str(row['GENERALCAU']).lower() or 
            kw in str(row['SPECIFICCA']).lower() for kw in keywords), axis=1)]
    print("Amount of filtered fire cause available:", len(Alaska_Dataset_Lightning))

    # Filter by date
    Alaska_Dataset_Lightning["DISCOVERYD"] = pd.to_datetime(Alaska_Dataset_Lightning["DISCOVERYD"])
    Alaska_Dataset_Lightning["OUTDATE"] = pd.to_datetime(Alaska_Dataset_Lightning["OUTDATE"], errors="coerce")
    Alaska_Dataset_Lightning_Date = Alaska_Dataset_Lightning[
        (Alaska_Dataset_Lightning['DISCOVERYD'] >= start_date) & 
        (Alaska_Dataset_Lightning["OUTDATE"].notnull()) & 
        (Alaska_Dataset_Lightning['OUTDATE'] <= end_date)]
    print("Amount of filtered fires (date) available:", len(Alaska_Dataset_Lightning_Date))

    # Remove false alarms and prescribed fires
    Alaska_Dataset_Lightning_Date = Alaska_Dataset_Lightning_Date[
        (Alaska_Dataset_Lightning_Date["FALSEALARM"] != "Y") & 
        (Alaska_Dataset_Lightning_Date["PRESCRIBED"] != "Y")]

    # Convert fire size from acres to hectares
    acres_to_hectares = 0.404686
    Alaska_Dataset_Lightning_Date['EST_HA'] = Alaska_Dataset_Lightning_Date['ESTIMATEDT'] * acres_to_hectares

    # Filter by fire size
    Alaska_Dataset_Lightning_Date = Alaska_Dataset_Lightning_Date[Alaska_Dataset_Lightning_Date['EST_HA'] >= 100]

    # Convert dates to string format
    Alaska_Dataset_Lightning_Date["DISCOVERYD"] = Alaska_Dataset_Lightning_Date["DISCOVERYD"].dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
    Alaska_Dataset_Lightning_Date["OUTDATE"] = Alaska_Dataset_Lightning_Date["OUTDATE"].dt.strftime("%Y-%m-%d")

    # Export the filtered dataset
    output_folder = getGroundtruth_Processed()
    output_filename = "Alaska_Points_Groundtruth"
    Export_Shapefile(Alaska_Dataset_Lightning_Date, output_folder, output_filename)


# # ******************  Alaska Polygons & Points Filter
def alaska_merge_polygons():
    Alaska_perimeter= gpd.read_file(f"{getGroundtruth_Processed()}Alaska_Polygons_Groundtruth.shp").to_crs(getGlanceCRS("NA"))
    Alaska_points= gpd.read_file(f"{getGroundtruth_Processed()}Alaska_Points_Groundtruth.shp").to_crs(getGlanceCRS("NA"))
    Alaska_perimeter.rename(columns={"FPOUTDATE":"OUTDATE"},inplace=True)
    Alaska_points.rename(columns={"ID":"FIREID"},inplace=True)
    Alaska_perimeter["OUTDATE"] = pd.to_datetime(Alaska_perimeter["OUTDATE"])
    Alaska_points["OUTDATE"] = pd.to_datetime(Alaska_points["OUTDATE"])
    output = []
    for _, polygon_row in Alaska_perimeter.iterrows():
        found_point = False
        new_row = polygon_row.copy()
        matching_points=Alaska_points[Alaska_points["FIREID"]==polygon_row["FIREID"]]
        if len(matching_points) == 1:
            if abs(matching_points.iloc[0]["OUTDATE"]-polygon_row["OUTDATE"]) <= pd.Timedelta(days = 3):
                new_row["LATITUDE"] = matching_points.iloc[0]["LATITUDE"]
                new_row["LONGITUDE"] = matching_points.iloc[0]["LONGITUDE"]
                new_row["Geo_Shape"] = "Single Ignition Point"
                new_row["DISCOVERYD"] = matching_points.iloc[0]["DISCOVERYD"]
                found_point = True
                Alaska_points.drop(Alaska_points[Alaska_points["FIREID"] == matching_points.iloc[0]["FIREID"]].index, inplace=True)
        elif len(matching_points)>1:
            print("Issue!")
        else:
            for _, point_row in Alaska_points.iterrows():
                if polygon_row.geometry.contains(point_row.geometry) and not pd.isna(point_row["OUTDATE"]) and not pd.isna(polygon_row["OUTDATE"]):
                    if abs(polygon_row["OUTDATE"] - point_row["OUTDATE"]) <= pd.Timedelta(days=7):
                        new_row["LATITUDE"] = point_row["LATITUDE"]
                        new_row["LONGITUDE"] = point_row["LONGITUDE"]
                        new_row["Geo_Shape"] = "Ignition Point Wrong Name"
                        new_row["DISCOVERYD"] = matching_points["DISCOVERYD"]
                        Alaska_points.drop(Alaska_points[Alaska_points["FIREID"] == point_row["FIREID"]].index, inplace=True)
                        found_point = True
                        continue
        if not found_point:
            new_row["Geo_Shape"] = "No Ignition Point"
        output.append(new_row)
    for i in range(len(output)):
        if output[i]["Geo_Shape"] == "No Ignition Point":
            for _, point_row in Alaska_points.iterrows():
                if output[i]["geometry"].buffer(2000).contains(point_row.geometry) and not pd.isna(point_row["OUTDATE"]) and not pd.isna(output[i]["OUTDATE"]):
                    if abs(output[i]["OUTDATE"] - point_row["OUTDATE"]) <= pd.Timedelta(days=7):
                        output[i]["LATITUDE"] = point_row["LATITUDE"]
                        output[i]["LONGITUDE"] = point_row["LONGITUDE"]
                        output[i]["Geo_Shape"] = "Ignition Point Wrong Name Buffer"
                        output[i]["DISCOVERYD"] = matching_points["DISCOVERYD"]
                        Alaska_points.drop(Alaska_points[Alaska_points["FIREID"] == point_row["FIREID"]].index, inplace=True)
                        found_point = True
                        continue

    if len(Alaska_points)>0:
        for _, point_row in Alaska_points.iterrows():
            new_row = point_row.copy()
            new_row["Geo_Shape"] = "Unmatched Ignition Point without perimeter"
            for _, polygon_row in Alaska_perimeter.iterrows():
                if polygon_row.geometry.contains(point_row.geometry) and not pd.isna(point_row["OUTDATE"]) and not pd.isna(polygon_row["OUTDATE"]):
                    if abs(polygon_row["OUTDATE"] - point_row["OUTDATE"]) <= pd.Timedelta(days=7):
                        new_row["Geo_Shape"] = "Unmatched Ignition Point with perimeter"    
            new_row["LATITUDE"] = point_row["LATITUDE"]
            new_row["LONGITUDE"] = point_row["LONGITUDE"]
            new_row["geometry"] = point_row["geometry"].buffer(math.sqrt((float(point_row["ESTIMATEDT"])*4046.86)/math.pi)) 
            output.append(new_row)

    output = gpd.GeoDataFrame(output, crs=getGlanceCRS("NA"))
    output["OUTDATE"] = output["OUTDATE"].dt.strftime("%Y-%m-%d")
    output_folder = getGroundtruth_Processed()
    output_filename = "Alaska_Groundtruth_merged_points_perimeters"
    Export_Shapefile(gpd.GeoDataFrame(output).to_crs(getGlanceCRS("NA")), output_folder, output_filename)