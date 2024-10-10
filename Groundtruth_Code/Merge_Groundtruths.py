import geopandas as gpd
import pandas as pd
import numpy as np
from tqdm import tqdm
from Utilities.Most_Used_Functions import Export_Shapefile, getGlanceCRS, loadBiomes
from shapely.validation import make_valid
from Utilities.Path_Utilities import getGroundtruth_Processed,getMODIS_Landcover
import rasterio
from shapely.geometry import box

def landcover_number_to_name_type1(df):
    LC_Names = ["Evergreen Needleleaf Forests", "Evergreen Broadleaf Forests", "Deciduous Needleleaf Forests",
                "Deciduous Broadleaf Forests", "Mixed Forests", "Closed Shrublands", "Open Shrublands", "Woody Savannas",
                "Savannas", "Grasslands", "Permanent Wetlands", "Croplands", "Urban and Built-up Lands",
                "Cropland/Natural Vegetation Mosaics","Permanent Snow and Ice", "Barren","Water Bodies"]
    
    df['LC_Name'] = df.apply(lambda row: LC_Names[int(row["LC_DN"])-1], axis=1)
    return df

def landcover_forest_notforest_type1(df):
    forest_list = [1,2,3,4,5,6,7,8,9,10,11,15]
    return "Forest" if int(df) in forest_list else "NonForest"


def load_datasets():
    alaska_path = getGroundtruth_Processed() + "Alaska_Groundtruth_merged_points_perimeters.shp"
    canada_path = getGroundtruth_Processed() + "Canada_Groundtruth_merged_points_perimeters.shp"
    alaska = gpd.read_file(alaska_path).to_crs(4326)
    canada = gpd.read_file(canada_path).to_crs(4326)
    return alaska, canada

def preprocess_datasets(alaska, canada):
    alaska['Country'] = 'Alaska'
    canada['Country'] = 'Canada'
    alaska['Cause'] = 'Natural/Lightning'
    canada['Cause'] = 'Natural/Lightning'

    alaska.rename(columns={ "DISCOVERYD": "StartDate", "OUTDATE": "EndDate", "EST_HA": "SIZE_HA"}, inplace=True)
    alaska.dropna(subset="StartDate",inplace=True)
    canada.rename(columns={"FIRE_ID": "FIREID", "REP_DATE": "StartDate", "OUT_DATE": "EndDate"}, inplace=True)

    included_columns = ["FIREID", "StartDate", "EndDate", "LATITUDE", "LONGITUDE", "SIZE_HA", "Country", "Cause", "geometry", "Geo_Shape"]

    alaska = alaska[included_columns]
    canada = canada[included_columns]
    return alaska, canada

def merge_datasets(alaska, canada):
    merged = gpd.GeoDataFrame(pd.concat([alaska, canada], ignore_index=True))
    return merged.to_crs(4326)


def load_landcover_data_raster():
    landcovers = {}
    landcover_bounds = {}
    for year in range(2012,2023):
        path = f"{getMODIS_Landcover()}MCD1Q12v61_Type1_{year}.tif"
        lc = rasterio.open(path)
        landcovers[str(year)] = lc
        landcover_bounds[str(year)] = box(*lc.bounds)
    return landcovers, landcover_bounds

def clean_geometry(geom):
    if not geom.is_valid:
        geom = make_valid(geom)
    return geom

def assign_landcover(merged_dataset, landcovers, landcovers_sindex):
    merged_dataset['LC_DN'] = None
    for index, row in tqdm(merged_dataset.iterrows(), total=len(merged_dataset)):
        year = row['StartDate'][:4]
        geom = clean_geometry(row.geometry)
        possible_matches_index = list(landcovers_sindex[year].intersection(geom.bounds))
        possible_matches = landcovers[year].iloc[possible_matches_index]

        max_intersection_area = 0
        corresponding_value = None

        for idx, row2 in possible_matches.iterrows():
            intersection = geom.intersection(clean_geometry(row2.geometry))
            if intersection.area > max_intersection_area:
                max_intersection_area = intersection.area
                corresponding_value = row2['DN']

        merged_dataset.at[index, 'LC_DN'] = int(corresponding_value) if corresponding_value is not None else None

def assign_landcover_raster(merged_dataset, landcovers, bounds):
    merged_dataset["LC_DN"] = 0
    to_remove_fires = []
    for index, row in tqdm(merged_dataset.iterrows(),total=len(merged_dataset)):
        year = row["StartDate"][:4]
        geom = row.geometry
        geom_bounds = box(*geom.bounds)

        LandCover = landcovers[year]
        lc_bounds = bounds[year]

        if not geom_bounds.intersects(lc_bounds):
            to_remove_fires.append(index)
            continue

        lcImage, _ = rasterio.mask.mask(LandCover,[geom],crop=True, nodata=255)
        values, counts = np.unique(lcImage[0],return_counts=True)
        mappedUniques = dict(zip(values,counts))
        currMax = 0
        currVal = 0
        for value in mappedUniques:
            if value == 255: continue
            if mappedUniques[value]>currMax:
                currMax= mappedUniques[value]
                currVal = value
        merged_dataset.at[index,"LC_DN"] = currVal
    merged_dataset = merged_dataset.drop(to_remove_fires)
    return merged_dataset

def merge_groundtruths_raster():
    alaska, canada = load_datasets()
    alaska, canada = preprocess_datasets(alaska, canada)
    merged_dataset = merge_datasets(alaska, canada)
    landcovers, bounds= load_landcover_data_raster()
    merged_dataset = assign_landcover_raster(merged_dataset, landcovers, bounds)

    merged_dataset = merged_dataset[merged_dataset["LC_DN"] != 0] # remove water
    merged_dataset = landcover_number_to_name_type1(merged_dataset)
    merged_dataset["Forest_Type"] = merged_dataset["LC_DN"].apply(landcover_forest_notforest_type1)
    merged_dataset = merged_dataset[merged_dataset["Forest_Type"] == "Forest"]

    Ecoregions_subset = loadBiomes()
    print(f"Biomes dataset loaded: {Ecoregions_subset.head()}")

    # Ensure that Ecoregions_subset is not empty and is a GeoDataFrame
    if isinstance(Ecoregions_subset, gpd.GeoDataFrame) and not Ecoregions_subset.empty:
        merged_dataset = gpd.clip(merged_dataset, Ecoregions_subset)
    else:
        print("Ecoregions_subset is either empty or not a valid GeoDataFrame.")
        return

    output_folder = getGroundtruth_Processed()
    output_filename = "Merged_Groundtruth_landcover_MCD1Q12v61_Type1_NA"
    merged_dataset = merged_dataset.to_crs(getGlanceCRS("NA"))
    Export_Shapefile(merged_dataset, output_folder, output_filename)


