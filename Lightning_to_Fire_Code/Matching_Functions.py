import os
import geopandas as gpd
from Utilities.Path_Utilities import getGWIS_Processed,getENTLN_Final_Filter,getMatchingCasePath
import glob
from datetime import timedelta, datetime
from Utilities.Most_Used_Functions import Export_Shapefile, getGlanceCRS
import pandas as pd
from pyproj import CRS,Transformer
from shapely.geometry import Point
pd.options.mode.chained_assignment = None

def computeDistance(lat1,long1, lat2,long2):
    transformer = Transformer.from_crs(CRS.from_epsg(4326),CRS.from_wkt(getGlanceCRS("NA")))
    p1 = Point(transformer.transform(lat1,long1))
    p2 = Point(transformer.transform(lat2,long2))

    return p1.distance(p2)

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
    
def categorize_distance(dist):
    if 0 <= dist <= 100:         # 0 to 100 meters
        return '100m'
    elif 100 < dist <= 1000:     # 101 meters to 1 km
        return '1km'
    elif 1000 < dist <= 2000:    # 1 km to 2 km
        return '2km'
    elif 2000 < dist <= 3000:    # 2 km to 3 km
        return '3km'
    elif 3000 < dist <= 4000:    # 3 km to 4 km
        return '4km'
    elif 4000 < dist <= 5000:    # 4 km to 5 km
        return '5km'
    elif 5000 < dist <= 6000:    # 5 km to 6 km
        return '6km'
    elif 6000 < dist <= 7000:    # 6 km to 7 km
        return '7km'
    elif 7000 < dist <= 8000:    # 7 km to 8 km
        return '8km'
    elif 8000 < dist <= 9000:    # 8 km to 9 km
        return '9km'
    elif 9000 < dist <= 10000:   # 9 km to 10 km
        return '10km'
    else:
        return 'NA'           # More than 10 km
    
def compute_t_min(row, fire_initial_date):
    t = (fire_initial_date - row["timestamp"]).total_seconds() / 3600
    return abs(t)

def compute_a_max(row, fire_initial_date, fire_center, time_max, distance_max):
    t = abs(fire_initial_date-row["timestamp"]).total_seconds()/3600
    s = row["geometry"].distance(fire_center)
    a_max = (1-t/time_max)*(1-s/distance_max)
    return a_max

def matching_function_polygon_tmin(tmax=14, smax = 10,unique_lightning = False, start_year = 2012,end_year =  2023):
    # ---- Minimum t_min
    tmax_hours = tmax*24
    smax_buffer = smax * 1000 # convert km to m
    for year in range(start_year, end_year):
        start_time = datetime.now()
        continent_year_shp = glob.glob(f"{getGWIS_Processed()}original_globfire_filtered_{year}_*.shp")
        for shape in continent_year_shp:
            continent_name = os.path.basename(shape).split("_")[-1][:2]
            Fire_Perimeters = gpd.read_file(shape).to_crs(getGlanceCRS(continent_name))
            Lightning_Flashes = gpd.read_file(f"{getENTLN_Final_Filter()}ENTLN_Boreal_{year}_{continent_name}.shp")
           
            if year > 2012:
                previous_year = gpd.read_file(f"{getENTLN_Final_Filter()}ENTLN_Boreal_{year-1}_{continent_name}.shp")
                previous_year["timestamp"] = pd.to_datetime(previous_year["timestamp"])
                previous_year = previous_year[previous_year["timestamp"].dt.month >= 4]
                previous_year["timestamp"] = previous_year["timestamp"].dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
                Lightning_Flashes = gpd.GeoDataFrame(pd.concat([Lightning_Flashes, previous_year]))
            
            # Iterate over each polygon
            LIW = []
            for _, polygon in Fire_Perimeters.iterrows():
                initial_date = datetime.strptime(polygon["initialdat"], "%Y-%m-%d")
                buffered_initial_date = initial_date - timedelta(hours=tmax_hours)
                initial_date += timedelta(days=1, microseconds=-1)
                polygon["fail_tmin"] = ""
                
                # --- 1. Spatial filter (1) Lightning within Fire Perimeter and then (2) Lightning withing 10km outside Perimeter
                lightning_within_geometry = gpd.clip(Lightning_Flashes, polygon.geometry)
                lightning_within_buffered_geometry = gpd.clip(Lightning_Flashes, polygon.geometry.buffer(smax_buffer))#

                if len(lightning_within_geometry) == 0 and len(lightning_within_buffered_geometry) == 0:
                    polygon["fail_tmin"] = "spatial"
                    LIW.append(polygon)
                    continue

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
                    polygon["t_Spatial"] = True
                    initial_spatial_temporal_filter = initial_spatial_temporal_filter_within
                elif initial_spatial_temporal_filter_buffered is not None and len(initial_spatial_temporal_filter_buffered) > 0:
                    polygon["t_Spatial"] = False
                    initial_spatial_temporal_filter = initial_spatial_temporal_filter_buffered
                else:
                    polygon["fail_tmin"] = "temporal"
                    LIW.append(polygon)
                    continue
            
                if polygon["t_Spatial"]:
                    initial_spatial_temporal_filter["t_min"] = initial_spatial_temporal_filter.apply(
                        lambda row: compute_t_min(row, initial_date), axis=1)
                    min_index = initial_spatial_temporal_filter["t_min"].idxmin()
                else:
                    initial_spatial_temporal_filter["t_min"]=initial_spatial_temporal_filter.apply(
                        lambda row: compute_a_max(row,initial_date,polygon.geometry,tmax_hours,smax_buffer),axis=1)
                    min_index = initial_spatial_temporal_filter["t_min"].idxmax()
                
                Lightning_Candidate = initial_spatial_temporal_filter.loc[min_index]
                output_row = polygon.copy()
                if unique_lightning:
                    Lightning_Flashes.drop(Lightning_Flashes[Lightning_Flashes["id"] == Lightning_Candidate["id"]].index, inplace=True)
                ###** tmin and holdover are the same
                output_row["t_min"] = Lightning_Candidate["t_min"]
                output_row["Holdover_d"] = (initial_date - Lightning_Candidate["timestamp"]).total_seconds() / 86400
                output_row["Holdoverrd"] = int(output_row["Holdover_d"])
                output_row["t_lat"] = Lightning_Candidate["latitude"]
                output_row["t_long"] = Lightning_Candidate["longitude"]
                output_row["t_ts"] = Lightning_Candidate["timestamp"].strftime("%Y-%m-%dT%H:%M:%S.%f")
                output_row["t_multipli"] = Lightning_Candidate["cgmultipli"]
                output_row["t_peakcurr"] = Lightning_Candidate["peakcurren"]
                output_row["t_duration"] = Lightning_Candidate["duration"]
                output_row["t_polarity"] = check_negative(Lightning_Candidate['peakcurren'])
                output_row["fire_size"] = categorize_fire_size(polygon["area_ha"])
                if not polygon["t_Spatial"]:
                    output_row["t_pol_dist"] = polygon.geometry.distance(Lightning_Candidate["geometry"])
                    output_row["dist_cat"] = categorize_distance(output_row["t_pol_dist"])
                LIW.append(output_row)
            
            if len(LIW) == 0:
                print(f"No matchings found for {os.path.basename(shape)}") 
                continue
            
            # Export the resulting GeoDataFrame
            output = gpd.GeoDataFrame(LIW, crs=getGlanceCRS(continent_name))
            output_folder = getMatchingCasePath(False,tmax,smax)
            output_filename = f"LIW_{year}_t_min_{continent_name}" 
            Export_Shapefile(output, output_folder, output_filename)
        print(f"{year} took {(datetime.now() - start_time).total_seconds()}s to process")

