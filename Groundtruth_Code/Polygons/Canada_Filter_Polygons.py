import geopandas as gpd
import pandas as pd
import math
from Utilities.Path_Utilities import getGroundtruth_Canada_Points,getGroundtruth_Canada_Polygons,getGroundtruth_Processed 
from Utilities.Most_Used_Functions import Export_Shapefile, getGlanceCRS,load_shapefile

def getCDateRange():
    """
    Returns the start and end dates for filtering the fire data.
    """
    start_date = '2012-01-01'
    end_date = '2020-12-31'
    return start_date, end_date
def canada_filter_polygons():
    """
    Filters the Canadian wildfire points dataset based on specified criteria and exports the filtered data.
    """
    start_date, end_date = getCDateRange()
    POLY_CWFIS_Dataset = getGroundtruth_Canada_Polygons()
    POLY_CWFIS_Dataset = gpd.read_file(POLY_CWFIS_Dataset).to_crs(4326)
    print("CWFIS_NFDB_Dataset CRS:", POLY_CWFIS_Dataset.crs)
    print("Amount of total fires in Canadian dataset unfiltered:", len(POLY_CWFIS_Dataset))

    # Convert dates to datetime format
    POLY_CWFIS_Dataset['REP_DATE'] = pd.to_datetime(POLY_CWFIS_Dataset['REP_DATE'])
    POLY_CWFIS_Dataset['OUT_DATE'] = pd.to_datetime(POLY_CWFIS_Dataset['OUT_DATE'])
    print("pd.to_datetime REP_DATE:", len(POLY_CWFIS_Dataset['REP_DATE']))

    # Filter Fires by ignition source, start date, end date, fire size, and landcover
    Filtered_POLY_CWFIS = POLY_CWFIS_Dataset[(POLY_CWFIS_Dataset['CAUSE'] == 'L') & 
                                             (POLY_CWFIS_Dataset['REP_DATE'] >= start_date) & 
                                             (POLY_CWFIS_Dataset['OUT_DATE'] <= end_date) &
                                             (POLY_CWFIS_Dataset['SIZE_HA'] >= 100)]

    print("Filtered before string:", len(Filtered_POLY_CWFIS))

    # Format dates as strings
    Filtered_POLY_CWFIS["REP_DATE"] = Filtered_POLY_CWFIS["REP_DATE"].dt.strftime("%Y-%m-%d")
    Filtered_POLY_CWFIS["OUT_DATE"] = Filtered_POLY_CWFIS["OUT_DATE"].dt.strftime("%Y-%m-%d")
    print("Amount of filtered fires (date and size) available:", len(Filtered_POLY_CWFIS))
    print("after string:", len(Filtered_POLY_CWFIS))

    # Export filtered_CWFIS
    output_folder = getGroundtruth_Processed()
    output_filename = "Canada_Polygon_Groundtruth"
    Export_Shapefile(Filtered_POLY_CWFIS, output_folder, output_filename)

def canada_filter_points():
    """
    Filters the Canadian wildfire points dataset based on specified criteria and exports the filtered data.
    """
    start_date, end_date = getCDateRange()
    NFDB_CWFIS_Dataset = getGroundtruth_Canada_Points()
    CWFIS_NFDB_Dataset = gpd.read_file(NFDB_CWFIS_Dataset).to_crs(4326)
    print("CWFIS_NFDB_Dataset CRS:", CWFIS_NFDB_Dataset.crs)
    print("Amount of total fires in Canadian dataset unfiltered:", len(CWFIS_NFDB_Dataset))

    # Convert dates to datetime format
    CWFIS_NFDB_Dataset['REP_DATE'] = pd.to_datetime(CWFIS_NFDB_Dataset['REP_DATE'])
    CWFIS_NFDB_Dataset['OUT_DATE'] = pd.to_datetime(CWFIS_NFDB_Dataset['OUT_DATE'])
    print("pd.to_datetime REP_DATE:", len(CWFIS_NFDB_Dataset['REP_DATE']))

    # Filter Fires by ignition source, start date, end date, fire size, and landcover
    Filtered_NFDB_CWFIS = CWFIS_NFDB_Dataset[(CWFIS_NFDB_Dataset['CAUSE'] == 'N') & 
                                             (CWFIS_NFDB_Dataset['REP_DATE'] >= start_date) & 
                                             (CWFIS_NFDB_Dataset['OUT_DATE'] <= end_date) &
                                             (CWFIS_NFDB_Dataset['SIZE_HA'] >= 100)]

    print("Filtered before string:", len(Filtered_NFDB_CWFIS))

    # Format dates as strings
    Filtered_NFDB_CWFIS["REP_DATE"] = Filtered_NFDB_CWFIS["REP_DATE"].dt.strftime("%Y-%m-%d")
    Filtered_NFDB_CWFIS["OUT_DATE"] = Filtered_NFDB_CWFIS["OUT_DATE"].dt.strftime("%Y-%m-%d")
    print("Amount of filtered fires (date and size) available:", len(Filtered_NFDB_CWFIS))
    print("after string:", len(Filtered_NFDB_CWFIS))

    # Export filtered_CWFIS
    output_folder = getGroundtruth_Processed()
    output_filename = "Canada_Points_Groundtruth"
    Export_Shapefile(Filtered_NFDB_CWFIS, output_folder, output_filename)

def canada_merge_polygons():
    CWFIS_points = load_shapefile(f"{getGroundtruth_Processed()}Canada_Points_Groundtruth.shp","NA")
    CWFIS_perimeter = load_shapefile(f"{getGroundtruth_Processed()}Canada_Polygon_Groundtruth.shp","NA")
    CWFIS_perimeter['REP_DATE'] = pd.to_datetime(CWFIS_perimeter['REP_DATE'])
    CWFIS_points['REP_DATE'] = pd.to_datetime(CWFIS_points['REP_DATE'])

    output = []

    for _, polygon_row in CWFIS_perimeter.iterrows():
        found_point = False
        new_row = polygon_row.copy()
        matching_points=CWFIS_points[CWFIS_points["NFDBFIREID"]==polygon_row["CFS_REF_ID"]]
        if len(matching_points) == 1:
            if abs(matching_points.iloc[0]["REP_DATE"]-polygon_row["REP_DATE"]) <= pd.Timedelta(days = 3):
                new_row["LATITUDE"] = matching_points.iloc[0]["LATITUDE"]
                new_row["LONGITUDE"] = matching_points.iloc[0]["LONGITUDE"]
                new_row["Geo_Shape"] = "Single Ignition Point"
                found_point = True
                CWFIS_points.drop(CWFIS_points[CWFIS_points["NFDBFIREID"] == matching_points.iloc[0]["NFDBFIREID"]].index, inplace=True)
        elif len(matching_points)>1:
            print("Issue!")
        else:
            for _, point_row in CWFIS_points.iterrows():
                if polygon_row.geometry.contains(point_row.geometry) and not pd.isna(point_row["REP_DATE"]) and not pd.isna(polygon_row["REP_DATE"]):
                    if abs(polygon_row["REP_DATE"] - point_row["REP_DATE"]) <= pd.Timedelta(days=7):
                        new_row["LATITUDE"] = point_row["LATITUDE"]
                        new_row["LONGITUDE"] = point_row["LONGITUDE"]
                        new_row["Geo_Shape"] = "Ignition Point Wrong Name"
                        CWFIS_points.drop(CWFIS_points[CWFIS_points["NFDBFIREID"] == point_row["NFDBFIREID"]].index, inplace=True)
                        found_point = True
                        continue
        if not found_point:
            new_row["Geo_Shape"] = "No Ignition Point"
        output.append(new_row)
    for i in range(len(output)):
        if output[i]["Geo_Shape"] == "No Ignition Point":
            for _, point_row in CWFIS_points.iterrows():
                if output[i]["geometry"].buffer(2000).contains(point_row.geometry) and not pd.isna(point_row["REP_DATE"]) and not pd.isna(output[i]["REP_DATE"]):
                    if abs(output[i]["REP_DATE"] - point_row["REP_DATE"]) <= pd.Timedelta(days=7):
                        output[i]["LATITUDE"] = point_row["LATITUDE"]
                        output[i]["LONGITUDE"] = point_row["LONGITUDE"]
                        output[i]["Geo_Shape"] = "Ignition Point Wrong Name Buffer"
                        CWFIS_points.drop(CWFIS_points[CWFIS_points["NFDBFIREID"] == point_row["NFDBFIREID"]].index, inplace=True)
                        found_point = True
                        continue
    if len(CWFIS_points)>0:
        for _, point_row in CWFIS_points.iterrows():
            new_row = point_row.copy()
            new_row["Geo_Shape"] = "Unmatched Ignition Point without perimeter"
            for _, polygon_row in CWFIS_perimeter.iterrows():
                if polygon_row.geometry.contains(point_row.geometry) and not pd.isna(point_row["REP_DATE"]) and not pd.isna(polygon_row["REP_DATE"]):
                    if abs(polygon_row["REP_DATE"] - point_row["REP_DATE"]) <= pd.Timedelta(days=7):
                        new_row["Geo_Shape"] = "Unmatched Ignition Point with perimeter"    
            new_row["LATITUDE"] = point_row["LATITUDE"]
            new_row["LONGITUDE"] = point_row["LONGITUDE"]
            new_row["geometry"] = point_row["geometry"].buffer(math.sqrt((float(point_row["SIZE_HA"])*10000)/math.pi))
            output.append(new_row)

    output = gpd.GeoDataFrame(output, crs=getGlanceCRS("NA"))
    output["REP_DATE"] = output["REP_DATE"].dt.strftime("%Y-%m-%d")
    output_folder = getGroundtruth_Processed()
    output_filename = "Canada_Groundtruth_merged_points_perimeters"
    Export_Shapefile(gpd.GeoDataFrame(output), output_folder, output_filename)
        






