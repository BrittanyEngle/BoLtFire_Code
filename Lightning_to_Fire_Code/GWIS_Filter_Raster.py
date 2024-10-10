import os
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
import glob
from datetime import datetime
import rasterio.mask
from Utilities.Most_Used_Functions import Export_Shapefile, loadBiomes, loadGlanceProjectionZones, getGlanceCRS, CountryBoundaries
from Utilities.Path_Utilities import getGWIS_Main_Dataset,getMODIS_Landcover,getGWIS_Processed

def gwis_filter():
    # Load necessary datasets
    Ecoregions_subset = loadBiomes()
    ECOregion_sindex = Ecoregions_subset.sindex

    glance_continents = loadGlanceProjectionZones()

    country_boundary = CountryBoundaries()
    country_boundary_sindex = country_boundary.sindex

    GWIS_files = glob.glob(os.path.join(getGWIS_Main_Dataset(), "original_globfire_filtered_*.shp"))

    # Process each GWIS file
    for gwis in GWIS_files:
        year = os.path.basename(gwis).split(".")[0].split("_")[-1]
        start = datetime.now()

        if os.path.getsize(gwis) <= 0:
            print(f"Empty file {os.path.basename(gwis)}")
            continue

        gdf = gpd.read_file(gwis)
        gdf = gdf.to_crs(4326)
        print("gdf crs:", gdf.crs)
        gdf = gdf[gdf['area_ha'] >= 100]
        gdf = gpd.clip(gdf,Ecoregions_subset)
        if len(gdf) > 0:
            with rasterio.open(f"{getMODIS_Landcover()}MCD1Q12v61_Type1_{year}.tif") as LandCover:
                for index, row in gdf.iterrows():
                    geom = row.geometry

                    #Ecoregion max clipping
                    max_intersection_area_eco = 0
                    corresponding_row = None
                    possible_matches_index_eco = list(ECOregion_sindex.intersection(geom.bounds))
                    possible_matches_eco = Ecoregions_subset.iloc[possible_matches_index_eco]
                    for idx, row2 in possible_matches_eco.iterrows():
                        intersection = geom.intersection(row2.geometry)
                        if intersection.area > max_intersection_area_eco:
                            max_intersection_area_eco = intersection.area
                            corresponding_row = row2
                    gdf.at[index, 'BIOME_NAME'] = corresponding_row["BIOME_NAME"]
                    gdf.at[index, 'ECO_BIOME_'] = corresponding_row["ECO_BIOME_"]
                    gdf.at[index, 'ECO_NAME'] = corresponding_row["ECO_NAME"]
                    gdf.at[index, 'ECO_ID'] = corresponding_row["ECO_ID"]
                    gdf.at[index, 'REALM'] = corresponding_row["REALM"]
                    
                    #Landcover max clipping
                    lcImage, transform = rasterio.mask.mask(LandCover,[geom],crop=True, nodata=255)
                    values, counts = np.unique(lcImage[0],return_counts=True)
                    mappedUniques = dict(zip(values,counts))
                    currMax = 0
                    currVal = 0
                    for value in mappedUniques:
                        if value == 255: continue
                        if mappedUniques[value]>currMax:
                            currMax= mappedUniques[value]
                            currVal = value
                    gdf.at[index,"LC_DN"] = currVal

                    #Country max clipping
                    max_intersection_area_country = 0
                    corresponding_value_country = None
                    possible_matches_index_country = list(country_boundary_sindex.intersection(geom.bounds))
                    possible_matches_country = country_boundary.iloc[possible_matches_index_country]
                    for idx, row2 in possible_matches_country.iterrows():
                        intersection = geom.intersection(row2.geometry)
                        if intersection.area > max_intersection_area_country:
                            max_intersection_area_country = intersection.area
                            corresponding_value_country = row2['country']

                    gdf.at[index, 'country'] = corresponding_value_country

                #update the landcover names and forest type on the whole dataset
                gdf = gdf[np.isfinite(gdf['LC_DN'])]
                gdf["LC_DN"] = gdf["LC_DN"].astype(int)    
                gdf = gdf[gdf["LC_DN"] != 0] # water is 0
                
                # Check unique values in LC_DN
                print("Unique LC_DN values:", gdf["LC_DN"].unique())
                
                gdf = landcover_number_to_name(gdf)
                gdf["ForestType"] = gdf.apply(lambda row: landcover_forest_notforest(row["LC_DN"]), axis=1)
                for i, continent in glance_continents.iterrows():
                    continent_name = continent["name"]
                    continent_geom = continent["geometry"]
                    data_in_continent = gdf[gdf.geometry.within(continent_geom)]
                    if len(data_in_continent) != 0:
                        data_in_continent = data_in_continent.to_crs(getGlanceCRS(continent_name))
                        output_filename = os.path.basename(gwis).split(".")[0] + "_" + continent_name
                        Export_Shapefile(data_in_continent, getGWIS_Processed(), output_filename)
                
            print(f"Elapsed time for {os.path.basename(gwis)}: {(datetime.now() - start).seconds}s total of {len(gdf)} available")


def landcover_number_to_name(df):
    LC_Names = ["Evergreen Needleleaf Forests", "Evergreen Broadleaf Forests", "Deciduous Needleleaf Forests",
                "Deciduous Broadleaf Forests", "Mixed Forests", "Closed Shrublands", "Open Shrublands", "Woody Savannas",
                "Savannas", "Grasslands", "Permanent Wetlands", "Croplands", "Urban and Built-up Lands",
                "Cropland/Natural Vegetation Mosaics","Permanent Snow and Ice", "Barren","Water Bodies"]
    
    df['LC_Name'] = df.apply(lambda row: LC_Names[int(row["LC_DN"])-1], axis=1)
    return df

def landcover_forest_notforest(df):
    forest_list = [1,2,3,4,5,6,7,8,9,10,11,15]
    return "Forest" if int(df) in forest_list else "NonForest"

"""
LC_Type1 Class Table

Value	Color	Description
1	#05450a	Evergreen Needleleaf Forests: dominated by evergreen conifer trees (canopy >2m). Tree cover >60%.
2	#086a10	Evergreen Broadleaf Forests: dominated by evergreen broadleaf and palmate trees (canopy >2m). Tree cover >60%.
3	#54a708	Deciduous Needleleaf Forests: dominated by deciduous needleleaf (larch) trees (canopy >2m). Tree cover >60%.
4	#78d203	Deciduous Broadleaf Forests: dominated by deciduous broadleaf trees (canopy >2m). Tree cover >60%.
5	#009900	Mixed Forests: dominated by neither deciduous nor evergreen (40-60% of each) tree type (canopy >2m). Tree cover >60%.
6	#c6b044	Closed Shrublands: dominated by woody perennials (1-2m height) >60% cover.
7	#dcd159	Open Shrublands: dominated by woody perennials (1-2m height) 10-60% cover.
8	#dade48	Woody Savannas: tree cover 30-60% (canopy >2m).
9	#fbff13	Savannas: tree cover 10-30% (canopy >2m).
10	#b6ff05	Grasslands: dominated by herbaceous annuals (<2m).
11	#27ff87	Permanent Wetlands: permanently inundated lands with 30-60% water cover and >10% vegetated cover.
12	#c24f44	Croplands: at least 60% of area is cultivated cropland.
13	#a5a5a5	Urban and Built-up Lands: at least 30% impervious surface area including building materials, asphalt and vehicles.
14	#ff6d4c	Cropland/Natural Vegetation Mosaics: mosaics of small-scale cultivation 40-60% with natural tree, shrub, or herbaceous vegetation.
15	#69fff8	Permanent Snow and Ice: at least 60% of area is covered by snow and ice for at least 10 months of the year.
16	#f9ffa4	Barren: at least 60% of area is non-vegetated barren (sand, rock, soil) areas with less than 10% vegetation.
17	#1c0dff	Water Bodies: at least 60% of area is covered by permanent water bodies.
"""