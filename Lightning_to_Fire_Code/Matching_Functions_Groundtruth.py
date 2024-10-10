import geopandas as gpd
from Utilities.Path_Utilities import getGroundtruth_Year_Split,getGroundtruth_Processed,getENTLN_Final_Filter,getMatchingCasePath
from datetime import timedelta, datetime
from shapely import Point
from Utilities.Most_Used_Functions import Export_Shapefile, getGlanceCRS
import pandas as pd
from pyproj import CRS,Transformer
import math
pd.options.mode.chained_assignment = None

# Check for polarity of lightning
def check_negative(value):
    if value<0:
        return 'negative'
    elif value>0:
        return 'positive'
    else:
        return "zero"

# Categorize fires based on their area
def categorize_fire_size(area_ha):
    if area_ha < 1000:
        return 'Small'
    elif 1000 <= area_ha < 10000:
        return 'Moderate'
    elif 10000 <= area_ha < 50000:
        return 'Large'
    elif 50000 <= area_ha < 100000:
        return 'Xlarge'
    else:
        return 'Mega'
    
def computeDistance(lat1,long1, lat2,long2):
    transformer = Transformer.from_crs(CRS.from_epsg(4326),CRS.from_wkt(getGlanceCRS("NA")))
    p1 = Point(transformer.transform(lat1,long1))
    p2 = Point(transformer.transform(lat2,long2))

    return p1.distance(p2)

def compute_t_min(row, fire_initial_date):
    t = (fire_initial_date - row["timestamp"]).total_seconds() / 3600
    return abs(t)

def compute_a_max(row, fire_initial_date, fire_center, time_max, distance_max):
    t = abs(fire_initial_date-row["timestamp"]).total_seconds()/3600
    s = row["geometry"].distance(fire_center)
    a_max = (1-t/time_max)*(1-s/distance_max)
    return a_max

def split_into_years():
    # Define the path to the shapefile
    groundtruth_path = f"{getGroundtruth_Processed()}Merged_Groundtruth_landcover_MCD1Q12v61_Type1_NA.shp"
    
    # Load the shapefile into a GeoDataFrame
    groundtruth_gdf = gpd.read_file(groundtruth_path)
    #groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({'No Ignition Point'})]
    output_folder =getGroundtruth_Year_Split()
    # Filter based on Geo_Shape column
    # Reproject the GeoDataFrame to the required CRS
    Fire_Perimeters = groundtruth_gdf.to_crs(getGlanceCRS("NA"))    
    
    # Split by year and export
    year_rows = {}
    for _, row in Fire_Perimeters.iterrows():
        year = int(row["StartDate"][:4])
        if year in year_rows:
            year_rows[year].append(row)
        else:
            year_rows[year] = [row]
    
    for year in year_rows:
        gdf = gpd.GeoDataFrame(year_rows[year], crs=getGlanceCRS("NA"))
        output_filename = f"merged_groundtruth_{year}_NA"
        Export_Shapefile(gdf, output_folder, output_filename)

def matching_function_all(time_delta_days=14, distance_delta_km=10):

    # ---- Minimum t_min
    tmin = time_delta_days*24
    tmin_km = distance_delta_km 
    tmin_buffer = tmin_km * 1000 # convert km to m

    # ---- Maximum Index A 
    amax = time_delta_days*24 
    amax_km = distance_delta_km
    amax_buffer= amax_km*1000 # convert km to m

    # --- Daily min dis
    dmin_km = distance_delta_km 
    dmin_buffer = dmin_km*1000 # convert km to m
    dmin_max_hours = time_delta_days*24

    for year in range(2012,2023):
        start_time = datetime.now()
        groundtruth_path = f"{getGroundtruth_Year_Split()}merged_groundtruth_{year}_NA.shp"
        Fire_Perimeters = gpd.read_file(groundtruth_path).to_crs(getGlanceCRS("NA"))
        Lightning_Flashes = gpd.read_file(f"{getENTLN_Final_Filter()}ENTLN_Boreal_{year}_NA.shp")

        if year > 2012:
            previous_year = gpd.read_file(f"{getENTLN_Final_Filter()}ENTLN_Boreal_{year-1}_NA.shp")
            previous_year["timestamp"] = pd.to_datetime(previous_year["timestamp"])
            previous_year = previous_year[previous_year["timestamp"].dt.month >= 4]
            previous_year["timestamp"] = previous_year["timestamp"].dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
            Lightning_Flashes = gpd.GeoDataFrame(pd.concat([Lightning_Flashes, previous_year]))

        if year < 2022:
            next_year = gpd.read_file(f"{getENTLN_Final_Filter()}ENTLN_Boreal_{year+1}_NA.shp")
            next_year["timestamp"] = pd.to_datetime(next_year["timestamp"])
            next_year = next_year[next_year["timestamp"].dt.month == 1]
            next_year["timestamp"] = next_year["timestamp"].dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
            Lightning_Flashes = gpd.GeoDataFrame(pd.concat([Lightning_Flashes, next_year]))
        LIW = []

        for index, polygon in Fire_Perimeters.iterrows():
            if polygon["Geo_Shape"] == "No Ignition Point":
                LIW.append(polygon)
                continue
            row = perform_tmin(polygon,Lightning_Flashes,tmin,tmin_buffer)
            row = perform_amax(row,Lightning_Flashes,amax,amax_buffer)
            row = perform_dmin(row, Lightning_Flashes, dmin_max_hours,dmin_buffer)
            LIW.append(row)
        if len(LIW)>0:
            output = gpd.GeoDataFrame(LIW, crs=getGlanceCRS("NA"))
            output_folder = getMatchingCasePath(True,time_delta_days,distance_delta_km)
            output_filename = f"LIW_{year}_Groundtruth_all_NA"
            Export_Shapefile(output,output_folder,output_filename)
        else:
            print(f"THIS SHOULD NOT HAPPEN {year}")
        print(f"{year} took {(datetime.now() - start_time).total_seconds()}s to process")

def perform_tmin(row,lightning_flashes ,time_buffer, distance_buffer):
    initial_date = datetime.strptime(row["StartDate"], "%Y-%m-%d")
    buffered_initial_date = initial_date - timedelta(hours=time_buffer)
    initial_date += timedelta(days=1, microseconds=-1)
    row["fail_tmin"] = ""

    # --- 1. Spatial filter (1) Lightning within Fire Perimeter and then (2) Lightning withing 10km outside Perimeter
    lightning_within_geometry = gpd.clip(lightning_flashes, row.geometry)
    lightning_within_buffered_geometry = gpd.clip(lightning_flashes, row.geometry.buffer(distance_buffer))#

    if len(lightning_within_geometry) == 0 and len(lightning_within_buffered_geometry) == 0:
        row["fail_tmin"] = "spatial"
        return row
    
    amax_temporal_max = time_buffer

    initial_spatial_temporal_filter_within = None
    initial_spatial_temporal_filter_buffered = None
    if len(lightning_within_geometry)>0:
        initial_spatial_filter = lightning_within_geometry
        # --- 2. Temporal filter
        initial_spatial_filter["timestamp"] = pd.to_datetime(initial_spatial_filter["timestamp"])
        initial_spatial_temporal_filter_within = initial_spatial_filter[(initial_spatial_filter["timestamp"] >= buffered_initial_date) &
                                                                    (initial_spatial_filter["timestamp"] <= initial_date)]

    if len(lightning_within_buffered_geometry) > 0:
        initial_spatial_filter = lightning_within_buffered_geometry
        # --- 2. Temporal filter
        initial_spatial_filter["timestamp"] = pd.to_datetime(initial_spatial_filter["timestamp"])
        initial_spatial_temporal_filter_buffered = initial_spatial_filter[(initial_spatial_filter["timestamp"] >= buffered_initial_date) &
                                                                    (initial_spatial_filter["timestamp"] <= initial_date)]
    
    if initial_spatial_temporal_filter_within is not None and len(initial_spatial_temporal_filter_within) > 0:
        row["t_Spatial"] = True
        initial_spatial_temporal_filter = initial_spatial_temporal_filter_within
    elif initial_spatial_temporal_filter_buffered is not None and len(initial_spatial_temporal_filter_buffered) > 0:
        row["t_Spatial"] = False
        initial_spatial_temporal_filter = initial_spatial_temporal_filter_buffered
    else:
        row["fail_tmin"] = "temporal"
        return row

    if row["t_Spatial"]:
        initial_spatial_temporal_filter["t_min"] = initial_spatial_temporal_filter.apply(
            lambda row: compute_t_min(row, initial_date), axis=1)
        min_index = initial_spatial_temporal_filter["t_min"].idxmin()
    else:
        initial_spatial_temporal_filter["t_min"] = initial_spatial_temporal_filter.apply(
            lambda r: compute_a_max(r,initial_date,row.geometry,amax_temporal_max,distance_buffer), axis=1)
        min_index = initial_spatial_temporal_filter["t_min"].idxmax()
    Lightning_Candidate = initial_spatial_temporal_filter.loc[min_index]
    output_row = row.copy()
    ###** tmin and holdover are the same
    output_row["t_min"] = Lightning_Candidate["t_min"]
    output_row["t_holdover"] = (initial_date - Lightning_Candidate["timestamp"]).total_seconds()/(24*3600)
    output_row["t_holdrnd"] = int(output_row["t_holdover"])
    output_row["t_lat"] = Lightning_Candidate["latitude"]
    output_row["t_long"] = Lightning_Candidate["longitude"]
    output_row["t_ts"] = Lightning_Candidate["timestamp"].strftime("%Y-%m-%dT%H:%M:%S.%f")
    output_row["t_src_dist"] = computeDistance(output_row["LATITUDE"],output_row["LONGITUDE"],output_row["t_lat"],output_row["t_long"])
    if not row["t_Spatial"]:
        output_row["t_pol_dist"] = row.geometry.distance(Lightning_Candidate["geometry"])
    return output_row

def perform_amax(row,lightning_flashes ,time_buffer, distance_buffer):
    transformer = Transformer.from_crs(CRS.from_epsg(4326),CRS.from_wkt(getGlanceCRS("NA")))
    ignition_point = Point(transformer.transform(row["LATITUDE"],row["LONGITUDE"]))
    buffered_center = ignition_point.buffer(distance_buffer) # Rep point to ensure the point is within the polygon
    initial_date = datetime.strptime(row["StartDate"], "%Y-%m-%d")
    buffered_initial_date = initial_date - timedelta(hours=time_buffer)
    initial_date += timedelta(days=1, microseconds=-1)
    row["fail_amax"] = ""

    # --- 1. Spatial filter of 10km
    initial_spatial_filter = gpd.clip(lightning_flashes,buffered_center)

    if len(initial_spatial_filter) <= 0 :
        row["fail_amax"] = "spatial"
        return row

    # --- 2. Temporal filter
    initial_spatial_filter["timestamp"]=pd.to_datetime(initial_spatial_filter["timestamp"])
    initial_spatial_temporal_filter = initial_spatial_filter[(initial_spatial_filter["timestamp"]>=buffered_initial_date) &
                                                                    (initial_spatial_filter["timestamp"]<=initial_date)]
    if len(initial_spatial_temporal_filter) <= 0 :
        row["fail_amax"] = "temporal"
        return row

    initial_spatial_temporal_filter["A_max"] = initial_spatial_temporal_filter.apply(
        lambda r: compute_a_max(r, initial_date, ignition_point, time_buffer, distance_buffer), axis=1)

    max_index = initial_spatial_temporal_filter["A_max"].idxmax()
    Lightning_Candidate = initial_spatial_temporal_filter.loc[max_index]
    output_row = row.copy()
    output_row["a_Spatial"] = row.geometry.contains(Lightning_Candidate.geometry)
    output_row["a_max"] = Lightning_Candidate["A_max"]
    output_row["a_holdover"] = (initial_date - Lightning_Candidate["timestamp"]).total_seconds()/(24*3600)
    output_row["a_holdrnd"] = int(output_row["a_holdover"])
    output_row["a_lat"] = Lightning_Candidate["latitude"]
    output_row["a_long"] = Lightning_Candidate["longitude"]
    output_row["a_ts"] = Lightning_Candidate["timestamp"].strftime("%Y-%m-%dT%H:%M:%S.%f")
    output_row["a_src_dist"] = computeDistance(output_row["LATITUDE"],output_row["LONGITUDE"],output_row["a_lat"],output_row["a_long"])
    return output_row

def perform_dmin(row,lightning_flashes ,time_buffer, distance_buffer):
    transformer = Transformer.from_crs(CRS.from_epsg(4326),CRS.from_wkt(getGlanceCRS("NA")))
    ignition_point = Point(transformer.transform(row["LATITUDE"],row["LONGITUDE"]))
    buffered_center = ignition_point.buffer(distance_buffer) # Rep point to ensure the point is within the polygon
    initial_date = datetime.strptime(row["StartDate"],"%Y-%m-%d")
    end_buffered_initial_date = initial_date-timedelta(hours=time_buffer)
    initial_date += timedelta(days=1, microseconds=-1)
    plus_buffered_initial_date = initial_date+timedelta(days=1)
    row["fail_dmin"] = ""

    # --- 1. Spatial filter of 10km
    initial_spatial_filter = gpd.clip(lightning_flashes, buffered_center)
    if len(initial_spatial_filter) <= 0 :
        row["fail_dmin"] = "spatial"
        return row
    
    # --- 2. Temporal filter
    initial_spatial_filter["timestamp"]=pd.to_datetime(initial_spatial_filter["timestamp"])
    initial_spatial_temporal_filter = initial_spatial_filter[(initial_spatial_filter["timestamp"]>=end_buffered_initial_date) &
                                                                    (initial_spatial_filter["timestamp"]<=plus_buffered_initial_date)]
    if len(initial_spatial_temporal_filter) <= 0 :
        row["fail_dmin"] = "temporal"
        return row

    Lightning_Candidate = None
    min_dist = float("inf")
    for day_shift in range(0,time_buffer+1,24):
        # Adjust the date based on day_shift
        adjusted_date = initial_date - timedelta(hours=day_shift)
        # Filter lightning strikes dataset for the adjusted date
        lightning_on_date = initial_spatial_temporal_filter[initial_spatial_temporal_filter['timestamp'].dt.day == adjusted_date.day]
        if len(lightning_on_date) == 0:
            continue
        for i, lightning_on_day in lightning_on_date.iterrows():
            lightning_point = lightning_on_day["geometry"]
            distance = lightning_point.distance(ignition_point)
            if distance <= min_dist:
                min_dist = distance
                Lightning_Candidate = lightning_on_day
        break

    if Lightning_Candidate is None:
        day_shift = -1
        adjusted_date = initial_date + timedelta(hours=24)   
        lightning_on_date =initial_spatial_temporal_filter[initial_spatial_temporal_filter["timestamp"].dt.day==adjusted_date.day]
        for _, lightning_on_day in lightning_on_date.iterrows():
            lightning_point = lightning_on_day["geometry"]
            distance = lightning_point.distance(ignition_point)
            if distance <= min_dist:
                min_dist = distance
                Lightning_Candidate = lightning_on_day     
                
    output_row = row.copy()
    output_row["d_Spatial"] = row.geometry.contains(Lightning_Candidate.geometry)      
    output_row["d_holdover"] = (initial_date - Lightning_Candidate["timestamp"]).total_seconds()/(24*3600)
    output_row["d_holdrnd"] = math.floor(output_row["d_holdover"])
    output_row["d_dist"] = min_dist
    output_row["d_disc_day"] = day_shift
    output_row["d_lat"] = Lightning_Candidate["latitude"]
    output_row["d_long"] = Lightning_Candidate["longitude"]
    output_row["d_ts"] = Lightning_Candidate["timestamp"].strftime("%Y-%m-%dT%H:%M:%S.%f")
    output_row["d_src_dist"] = computeDistance(output_row["LATITUDE"],output_row["LONGITUDE"],output_row["d_lat"],output_row["d_long"])
    return output_row