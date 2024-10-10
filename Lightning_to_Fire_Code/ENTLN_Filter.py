import pandas as pd
from datetime import datetime
import geopandas as gpd
import numpy as np 
import glob
import os
from Utilities.Most_Used_Functions import Export_Shapefile, loadBiomes, loadGlanceProjectionZones, getGlanceCRS, CountryBoundaries
from Utilities.Path_Utilities import getENTLN_CSVDataset, getENTLN_Boreal_Filter, getENTLN_Continent_Split, getMODIS_Landcover,getENTLN_Final_Filter
import rasterio
 
def entln_filter():
    Ecoregions_subset = loadBiomes()
    for year in range(2012,2023):
        ENTLN_CSV_files = glob.glob(getENTLN_CSVDataset() + str(year)+ "/*.csv")
        # Read each CSV file into a GeoDataFrame and concatenate them into a single GeoDataFrame
        dfs = []
        for csv_file in ENTLN_CSV_files:
            filter_list = [0,40]
            start = datetime.now()
            if os.path.getsize(csv_file) <= 0:
                print(f"Empty file {os.path.basename(csv_file)}")
                continue
            df = pd.read_csv(csv_file)
            df = df[df["type"].isin(filter_list)]
            df = df[df["peakcurrent"] != 0]
            gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["longitude"], df["latitude"]), crs = "EPSG:4326")
            clipped_gdf = gpd.clip(gdf,Ecoregions_subset)
            clipped_gdf = gpd.sjoin(clipped_gdf,Ecoregions_subset, how="left", predicate='intersects')
            print(f"Elapsed time for {os.path.basename(csv_file)}: {(datetime.now()-start).seconds}s total of {len(clipped_gdf)} points available")
            if len(clipped_gdf) > 0:
                dfs.append(clipped_gdf)

        # Concatenate all GeoDataFrames into a single GeoDataFrame
        merged_gdf =pd.concat(dfs, ignore_index=True)
        print(f"Amount of points in merged_gdf {year}:" )
        print(len(merged_gdf))
        print("merged_gdf crs:", merged_gdf.crs)
        # Export
        output_folder = f"{getENTLN_Boreal_Filter()}"
        output_filename = "ENTLN_Boreal_"+ str(year)
        Export_Shapefile(merged_gdf, output_folder, output_filename)

def entln_continent_split():
    root_folder = getENTLN_Boreal_Filter()
    proj_zones = loadGlanceProjectionZones()

    for year in range(2012,2023):
        start = datetime.now()
        shp_files = glob.glob(os.path.join(root_folder,f"ENTLN_Boreal_{year}.shp"))
        if len(shp_files) > 1:
            print(f"More than 1 Shapefile for the year {year}")
            continue
        gdf = gpd.read_file(shp_files[0])
        print(f"loaded {shp_files[0]} in {(datetime.now()-start).total_seconds()} seconds")
        continent_splits = 0
        continents = 0 
        for i, continent in proj_zones.iterrows():
            continent_name = continent["name"]
            continent_geom = continent["geometry"]
            data_in_continent = gdf[gdf.geometry.within(continent_geom)]
            if len(data_in_continent)!=0:
                continent_splits += len(data_in_continent)
                continents+=1
                data_in_continent = data_in_continent.to_crs(getGlanceCRS(continent_name))
                output_filename = os.path.basename(shp_files[0]).split(".")[0]+"_" + continent_name
                Export_Shapefile(data_in_continent, getENTLN_Continent_Split(), output_filename)

        print(f"Processed {year} in {(datetime.now()-start).total_seconds()} seconds.")
        print(f"Assigned {continent_splits} of {len(gdf)} to {continents} continents.")

def entln_landcover_country_split():
    # Load necessary datasets
    glance_continents = loadGlanceProjectionZones()
    country_boundary = CountryBoundaries()

    ENTLN_files = glob.glob(os.path.join(f"{getENTLN_Continent_Split()}", "*.shp"))

    # Process each ENTLN file
    for entln_file in ENTLN_files:
        year = os.path.basename(entln_file).split(".")[0].split("_")[-2]
        start = datetime.now()

        try:
            gdf = gpd.read_file(entln_file).to_crs(4326)
            if gdf.empty:
                print(f"GeoDataFrame is empty for {os.path.basename(entln_file)}")
                continue

            print("gdf crs:", gdf.crs)
            print(f"Processing file {os.path.basename(entln_file)} with {len(gdf)} records")

            with rasterio.open(f"{getMODIS_Landcover()}MCD1Q12v61_Type1_{year}.tif") as LandCover:
                # Extract landcover values for each point
                coords = [(geom.x, geom.y) for geom in gdf.geometry]
                landcover_values = [val[0] for val in LandCover.sample(coords)]
                gdf['LC_DN'] = landcover_values

                # Update the landcover names and forest type on the whole dataset
                gdf = gdf[np.isfinite(gdf['LC_DN'])]
                gdf["LC_DN"] = gdf["LC_DN"].astype(int)
                gdf = gdf[gdf["LC_DN"] != 0]  # water is 0

                # Check unique values in LC_DN
                print("Unique LC_DN values:", gdf["LC_DN"].unique())

                gdf = landcover_number_to_name(gdf)
                gdf["ForestType"] = gdf.apply(lambda row: landcover_forest_notforest(row["LC_DN"]), axis=1)

                # Spatial join with country boundaries to get the country for each point
                gdf = gpd.sjoin(gdf, country_boundary, how="left", op="within")
                gdf.rename(columns={'country': 'country'}, inplace=True)

                for i, continent in glance_continents.iterrows():
                    continent_name = continent["name"]
                    continent_geom = continent["geometry"]
                    data_in_continent = gdf[gdf.geometry.within(continent_geom)]
                    if len(data_in_continent) != 0:
                        data_in_continent = data_in_continent.to_crs(getGlanceCRS(continent_name))
                        output_filename = os.path.basename(entln_file).split(".")[0]
                        Export_Shapefile(data_in_continent, getENTLN_Final_Filter(), output_filename)

            print(f"Elapsed time for {os.path.basename(entln_file)}: {(datetime.now() - start).seconds}s total of {len(gdf)} available")

        except Exception as e:
            print(f"Error processing {os.path.basename(entln_file)}: {e}")

def landcover_number_to_name(df):
    LC_Names = ["Evergreen Needleleaf Forests", "Evergreen Broadleaf Forests", "Deciduous Needleleaf Forests",
                "Deciduous Broadleaf Forests", "Mixed Forests", "Closed Shrublands", "Open Shrublands", "Woody Savannas",
                "Savannas", "Grasslands", "Permanent Wetlands", "Croplands", "Urban and Built-up Lands",
                "Cropland/Natural Vegetation Mosaics", "Permanent Snow and Ice", "Barren", "Water Bodies"]

    df['LC_Name'] = df.apply(lambda row: LC_Names[int(row["LC_DN"]) - 1] if 0 < int(row["LC_DN"]) <= len(LC_Names) else 'Unknown', axis=1)
    return df

def landcover_forest_notforest(df):
    forest_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15]
    return "Forest" if int(df) in forest_list else "NonForest"

def count_forest_and_landcover_rows_in_shapefiles():
    directory = getENTLN_Final_Filter()
    # List to store data for the main dataframe
    main_data = []
    
    # List to store data for the ECO_NAME dataframe
    eco_name_data = []
    
    # List to store data for the landcover dataframe
    landcover_data = []
    
    # List to store data for the continent dataframe
    continent_data = []
    
    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        # Check if the file is a shapefile
        if filename.endswith('.shp'):
            # Determine the continent based on the filename
            if '_NA' in filename:
                continent = 'North America'
            elif '_EU' in filename or '_AS' in filename:
                continent = 'Eurasia'
            else:
                continue  # Skip files that don't match the naming convention
            
            # Construct the full file path
            filepath = os.path.join(directory, filename)
            # Read the shapefile using geopandas
            gdf = gpd.read_file(filepath)
            # Filter rows where the 'ForestType' column is 'Forest'
            forest_gdf = gdf[gdf['ForestType'] == 'Forest']
            # Get the number of rows in the filtered dataframe
            forest_count = len(forest_gdf)
            
            # Count rows by country
            forest_gdf["country"]=forest_gdf["country"].fillna("No Country")
            country_counts = forest_gdf['country'].value_counts().to_dict()
            
            # Count rows by ECO_NAME
            eco_name_counts = forest_gdf['ECO_NAME'].value_counts().to_dict()
            
            # Count rows by landcover
            landcover_counts = forest_gdf['LC_Name'].value_counts().to_dict()
            
            # Add data to the main_data list
            for country, country_count in country_counts.items():
                main_data.append({
                    'filename': filename,
                    'total_forest': forest_count,
                    'country': country,
                    'country_count': country_count
                })
            
            # Add data to the eco_name_data list
            for eco_name, eco_count in eco_name_counts.items():
                eco_name_data.append({
                    'filename': filename,
                    'ECO_NAME': eco_name,
                    'eco_name_count': eco_count
                })
            
            # Add data to the landcover_data list
            for landcover, landcover_count in landcover_counts.items():
                landcover_data.append({
                    'filename': filename,
                    'landcover': landcover,
                    'landcover_count': landcover_count
                })
            
            # Add data to the continent_data list
            continent_data.append({
                'filename': filename,
                'continent': continent,
                'forest_count': forest_count
            })
    
    # Create DataFrames from the data
    main_df = pd.DataFrame(main_data)
    eco_name_df = pd.DataFrame(eco_name_data)
    landcover_df = pd.DataFrame(landcover_data)
    
    # Create DataFrame for continent data
    continent_df = pd.DataFrame(continent_data)
    continent_grouped_df = continent_df.groupby('continent')['forest_count'].sum().reset_index()
    
    # Add a total row
    total_row = pd.DataFrame({'continent': ['Total'], 'forest_count': [continent_grouped_df['forest_count'].sum()]})
    continent_grouped_df = pd.concat([continent_grouped_df, total_row], ignore_index=True)
    
    # Save the dataframes to an Excel file with separate sheets
    output_file = f'{getENTLN_Final_Filter()}entln_stats.xlsx'
    with pd.ExcelWriter(output_file) as writer:
        main_df.to_excel(writer, sheet_name='Country_Stats', index=False)
        eco_name_df.to_excel(writer, sheet_name='ECO_NAME_Stats', index=False)
        landcover_df.to_excel(writer, sheet_name='Landcover_Stats', index=False)
        continent_grouped_df.to_excel(writer, sheet_name='Continent_Stats', index=False)

    print(f"Forest, landcover, and continent row counts saved to {output_file}")