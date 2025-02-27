import geopandas as gpd
import pandas as pd
import numpy as np
from Utilities.Most_Used_Functions import Export_Shapefile, getGlanceCRS, load_shapefile
from Utilities.Path_Utilities import getRootFolder
import datetime
from datetime import timedelta
import glob
import os

def Export_Excel(gdf, output_folder, output_filename):
    """
    Exports a GeoDataFrame to an Excel file.
    """
    gdf.drop(columns='geometry').to_excel(f"{output_folder}{output_filename}.xlsx", index=False)

def merge_results_by_continent(input_folder):
    """
    Merges shapefiles by continent and exports the merged result.
    """
    m = []
    for continent in ["NA", "AS", "EU"]:
        fires = []
        # Search for shapefiles in the directory
        results = glob.glob(f"{input_folder}LIW_*_t_min_{continent}.shp")
        for file in results:
            dataset = gpd.read_file(file)
            dataset = dataset.to_crs(getGlanceCRS(continent))
            fires.append(dataset)
        
        # Concatenate GeoDataFrames ensuring the geometry column is properly handled
        merged_gdf = gpd.GeoDataFrame(pd.concat(fires, ignore_index=True), crs=getGlanceCRS(continent), geometry='geometry')
            # Remove rows with null values in the "country" column
        merged_gdf = merged_gdf.dropna(subset=['country'])
        # Filter the "ForestType" column to only include "Forest"
        merged_gdf = merged_gdf[merged_gdf['ForestType'] == 'Forest']
        # Add "fire_year" column based on the first four digits of the "initialdat" column
        merged_gdf['fire_year'] = merged_gdf['initialdat'].str[:4]
        
        # Add "Continent" column based on the "country" column
        merged_gdf['Continent'] = merged_gdf['country'].apply(lambda x: 'North America' if x in ['Canada', 'United States of America'] else 'Eurasia')

       
        # Remove MODIS fires in Canada 2021 & 2022 since the agency dataset only goes to 2020
        merged_gdf['initialdat'] = pd.to_datetime(merged_gdf['initialdat'])

        # If necessary, format the 'initialdat' column as a string (this might not be necessary depending on the context)
        merged_gdf['initialdat'] = merged_gdf['initialdat'].dt.strftime('%Y-%m-%d')

        # Define output folder and filename
        output_filename = "Merged_Results_2012_2022_" + continent + "_Forest"
        
        # Export the merged shapefile
        Export_Shapefile(merged_gdf, input_folder, output_filename)
        m.append(merged_gdf.to_crs(4326))
    Export_Excel(pd.concat(m),input_folder,"Merged_Results_2012_2022_all_continents")

def filter_fire_ids(main_gdf,exclude_gdf, output_shp_path):

    if len(exclude_gdf)>0:
    # Extract the IDs to exclude
        exclude_ids = exclude_gdf['id'].unique()
    else:
        exclude_ids = []

    # Filter the main GeoDataFrame to exclude these IDs
    filtered_gdf = main_gdf[~main_gdf['id'].isin(exclude_ids)]

    # Write the filtered GeoDataFrame to a new shapefile
    filtered_gdf.to_file(output_shp_path)

    Export_Excel(filtered_gdf, "/".join(output_shp_path.split("/")[:-1])+"/",output_shp_path.split("/")[-1])

def compute_confusion_matrix(input_folder, groundtruth_file, date_column= "StartDate", spatial_threshold=10, temporal_buffer_back=7,temporal_buffer_forward = 7, ha_filter= 100, startYear = 2012):
    """
    Computes confusion matrix by comparing project results with ground truth data.
    """
    results_list = []  # List to store results for each continent
    file_path = input_folder
    groundtruth_geoshape = groundtruth_file.split("/")[-2]
    output_folder = f"{file_path}" \
        f"-{temporal_buffer_back}days_+{temporal_buffer_forward}days_{spatial_threshold}km_{startYear}/{ha_filter}ha_geoshape{groundtruth_geoshape}/"

    # Load project results
    Fire_Perimeters = load_shapefile(file_path + "Merged_Results_2012_2022_NA_Forest.shp","NA")
    Fire_Perimeters = Fire_Perimeters[~Fire_Perimeters["fail_tmin"].isin({'temporal', 'spatial'})]
    Fire_Perimeters = Fire_Perimeters[(Fire_Perimeters['area_ha'] >= ha_filter)]

    # Filter for Canada and the specified date range
    canada_mask = (Fire_Perimeters['country'] == 'Canada') & \
                    (Fire_Perimeters['initialdat'] >= '2021-01-01') & \
                    (Fire_Perimeters['initialdat'] <= '2022-12-31')

    # Update only the entries for Canada within the date range
    Fire_Perimeters = Fire_Perimeters[~canada_mask]

    Export_Shapefile(Fire_Perimeters,output_folder,"GWIS_Filtered")
    Fire_Perimeters[date_column] = pd.to_datetime(Fire_Perimeters["initialdat"])
    Groundtruth = load_shapefile(groundtruth_file, "NA")
    
    Groundtruth["fire_year"] = pd.DatetimeIndex(Groundtruth[date_column]).year
    start_date = datetime.datetime.strptime(f'{startYear}/01/01', "%Y/%m/%d")
    end_date = datetime.datetime.strptime('2022/12/31', "%Y/%m/%d")
    Groundtruth[date_column] = pd.to_datetime(Groundtruth[date_column])
    Groundtruth["EndDate"] = pd.to_datetime(Groundtruth["EndDate"], errors="coerce")
    Groundtruth = Groundtruth[(Groundtruth[date_column] >= start_date) & (Groundtruth["EndDate"].notnull()) & (Groundtruth['EndDate'] <= end_date)]
    print(f"Amount of Groundtruth fires available in NA: {len(Groundtruth)}")
    # Buffer ground truth polygons by the spatial threshold
    Groundtruth['geometry'] = Groundtruth.geometry.buffer(spatial_threshold*1000)

    unmatched = []
    matched = []
    multi_matched=[]
    for index, row in Groundtruth.iterrows():
        row_geometry = row.geometry
        buffered_initial_date = row[date_column] - timedelta(days=temporal_buffer_back)
        forward_initial_date = row[date_column] +timedelta(days=temporal_buffer_forward)
        spatial_filter = Fire_Perimeters[Fire_Perimeters.intersects(row_geometry)]
        if len(spatial_filter)==0:
            row["gwis_fail"] = "spatial"
            unmatched.append(row)
            continue
        spatial_temporal_filter = spatial_filter[(spatial_filter[date_column]>=buffered_initial_date) &
                                                (spatial_filter[date_column]<=forward_initial_date)]
        if  len(spatial_temporal_filter)==0:
            row["gwis_fail"] = "temporal"
            unmatched.append(row)
            continue
        spatial_temporal_filter = spatial_temporal_filter
        spatial_temporal_filter["area"] = spatial_temporal_filter.geometry.area
        spatial_temporal_filter["area_diff"] = (spatial_temporal_filter["area"]-row_geometry.area).abs()

        matching_one = spatial_temporal_filter.loc[spatial_temporal_filter["area_diff"].idxmin()]

        output_row = row.copy()
        output_row["id"] = matching_one["id"]
        output_row["initialdat"] = matching_one["initialdat"]
        output_row["finaldate"] = matching_one["finaldate"]
        output_row["area_ha"] = matching_one["area_ha"]
        output_row["area_diff"] = matching_one["area_diff"]

        matched.append(output_row)
        for i,r in spatial_temporal_filter.iterrows():
            multi_matched.append(r)
    groundtruth_geoshape = groundtruth_file.split("/")[-2]
        
    if len(matched)>0:
        matched = gpd.GeoDataFrame(matched,crs=getGlanceCRS("NA"))
        matched[date_column] = matched[date_column].dt.strftime("%Y-%m-%d")
        matched["EndDate"] = matched["EndDate"].dt.strftime("%Y-%m-%d")
        Export_Shapefile(matched,output_folder, "matched")
        Export_Excel(matched, output_folder, "matched.shp")

    if len(unmatched)>0:
        unmatched = gpd.GeoDataFrame(unmatched,crs=getGlanceCRS("NA"))
        unmatched[date_column] = unmatched[date_column].dt.strftime("%Y-%m-%d")
        unmatched["EndDate"] = unmatched["EndDate"].dt.strftime("%Y-%m-%d")
        Export_Shapefile(unmatched,output_folder, "unmatched")
        Export_Excel(unmatched, output_folder, "unmatched.shp")
    
    if len(multi_matched)>0:
        multi_matched = gpd.GeoDataFrame(multi_matched,crs=getGlanceCRS("NA"))
        multi_matched[date_column] = multi_matched[date_column].dt.strftime("%Y-%m-%d")
        Export_Shapefile(multi_matched,output_folder, "gwis_multi_match")
        Export_Excel(multi_matched, output_folder, "gwis_multi_match.shp")

    Fire_Perimeters[date_column] = Fire_Perimeters[date_column].dt.strftime("%Y-%m-%d")
    filter_fire_ids(Fire_Perimeters, matched, output_folder + "not_used_gwis.shp")
    filter_fire_ids(Fire_Perimeters,multi_matched,output_folder+"not_used_gwis_multi_match_version.shp")
    
    # Construct confusion matrix
    true_positive = len(matched)
    false_positive = len(Fire_Perimeters) - true_positive
    false_negative = len(unmatched)
    true_negative = np.nan  # Not applicable in this context

    # User Accuracy (Precision)
    user_accuracy = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    # Producer Accuracy (Recall)
    producer_accuracy = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0

    # Commission Error
    commission_error = false_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    # Omission Error
    omission_error = false_negative / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0

    results_list.append({
        'true_positive': true_positive,
        'false_positive': false_positive,
        'false_negative': false_negative,
        'true_negative': true_negative,
        "User Accuracy": user_accuracy,
        "Producer Accuracy": producer_accuracy,
        "Commission Error": commission_error,
        "Omission Error": omission_error
    })

        # Construct confusion matrix
    true_positive = len(matched)
    false_positive = len(Fire_Perimeters) - len(multi_matched)
    false_negative = len(unmatched)
    true_negative = np.nan  # Not applicable in this context

    # User Accuracy (Precision)
    user_accuracy = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    # Producer Accuracy (Recall)
    producer_accuracy = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0

    # Commission Error
    commission_error = false_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    # Omission Error
    omission_error = false_negative / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0

    results_list.append({
        'true_positive': true_positive,
        'false_positive': false_positive,
        'false_negative': false_negative,
        'true_negative': true_negative,
        "User Accuracy": user_accuracy,
        "Producer Accuracy": producer_accuracy,
        "Commission Error": commission_error,
        "Omission Error": omission_error
    })

# Convert results to a DataFrame and save as Excel
    results_df = pd.DataFrame(results_list)
    results_df.to_excel(f"{output_folder}Confusion_Matrix_Results_tmin_multi_match.xlsx", index=False)

def mergeConfusionMatrices():
    tmin_base = f"{getRootFolder()}Matching/GwIS_new/"
    initial_depth = tmin_base.count(os.sep)
    target_filename = "Confusion_Matrix_Results_tmin_multi_match.xlsx"

    gt_base = f"{getRootFolder()}Matching/Groundtruth/"

    all_confusion_matrices = []
    for root, dirs, files, in os.walk(tmin_base):
        depth = root.count(os.sep) - initial_depth
        if depth!=2:
            continue

        if target_filename in files:
            all_confusion_matrices.append(os.path.join(root,target_filename))

    all_results = {
        "case":[],
        "Multimatch":[],
        "tmin":[],
        "Year":[],
        "tmin +1day":[],
        "km":[],
        "back days":[],
        "forward days":[],
        "ha size":[],
        "geo_shape":[],
        #"fire_perimeters just t_min":[],
        #"groundtruth just t_min":[],
        "Notes":[],
        "Total gwis Fires":[],
        "true_positive":[],
        "false_positive":[],
        "false_negative":[],
        "true_negative":[],
        "User Accuracy":[],
        "Producer Accuracy":[],
        "Commission Error":[],
        "Omission Error":[],
        "total gt fires":[],
        "a_success %":[],
        "t_success %":[],
        "d_success %":[],
        "a_success total":[],
        "t_success total":[],
        "d_success total":[],
        "same_stroke %":[],
        "dif_stroke":[],
        "same_stroke":[],
        "all no stroke":[],
        "t_src_dist - mean dist (meters)":[],
        "a_src_dist- mean dist (meters)":[],
        "d_src_dist- mean dist (meters)":[]
    }

    for confusion_matrix in all_confusion_matrices:
        df= pd.read_excel(confusion_matrix)
        case_name = os.sep.join(confusion_matrix.split(os.sep)[initial_depth:-1])
        tmin = confusion_matrix.split(os.sep)[initial_depth]
        tmin_plus1day = "no"

        third_component = confusion_matrix.split(os.sep)[initial_depth+1].split("_")
        year = third_component[3]+"-2022"
        km = third_component[2][0:2]
        back_days = third_component[0][1]
        forward_days = third_component[1][1]

        forth_component = confusion_matrix.split(os.sep)[initial_depth+2].split("_")
        ha_size = ">"+forth_component[0][:-2]
        geoshape = forth_component[1][-1]
        fire_perimeter_just_tmin = "yes"
        groundtruth_just_tmin = "no"

        total_gwis_fires = int(df.iloc[0]["true_positive"]) + int(df.iloc[0]["false_positive"])

        gt_df = pd.read_excel(os.path.join(gt_base,tmin,ha_size[1:]+"ha",geoshape,"Groundtruth_Merged_2012_2022_statistics_summary.xlsx"))
        total_gt_fires = int(gt_df[gt_df["Column"]=="Total GT Count"]["Count"].values[0])
        a_success_percent = gt_df[gt_df["Column"]=="a success"]["Value"].values[0]
        t_success_percent = gt_df[gt_df["Column"]=="t success"]["Value"].values[0]
        d_success_percent = gt_df[gt_df["Column"]=="d success"]["Value"].values[0]

        a_success_total = gt_df[gt_df["Column"]=="a success"]["Count"].values[0]
        t_success_total = gt_df[gt_df["Column"]=="t success"]["Count"].values[0]
        d_success_total = gt_df[gt_df["Column"]=="d success"]["Count"].values[0]

        same_stroke_total = int(gt_df[(gt_df["Column"]=="d_a_t_Comp")&(gt_df["Value"]=="same_stroke")]["Count"].values[0])
        dif_stroke_total = int(gt_df[(gt_df["Column"]=="d_a_t_Comp")&(gt_df["Value"]=="dif_stroke")]["Count"].values[0])
        all_no_stroke = int(gt_df[(gt_df["Column"]=="d_a_t_Comp")&(gt_df["Value"]=="all no stroke")]["Count"].values[0])
        same_stroke_percent = same_stroke_total/total_gt_fires*100

        t_src_dist_mean = gt_df[(gt_df["Column"]=="t_src_dist")&(gt_df["Value"]=="mean dist (meters)")]["Count"].values[0]
        a_src_dist_mean = gt_df[(gt_df["Column"]=="a_src_dist")&(gt_df["Value"]=="mean dist (meters)")]["Count"].values[0]
        d_src_dist_mean = gt_df[(gt_df["Column"]=="d_src_dist")&(gt_df["Value"]=="mean dist (meters)")]["Count"].values[0]

        for i in range(2):
            row = df.iloc[i]
            all_results["case"].append(case_name)
            all_results["Multimatch"].append("no" if i==0 else "yes")
            all_results["tmin"].append(tmin)
            all_results["Year"].append(year)
            all_results["tmin +1day"].append(tmin_plus1day)
            all_results["km"].append(km)
            all_results["back days"].append(back_days)
            all_results["forward days"].append(forward_days)
            all_results["ha size"].append(ha_size)
            all_results["geo_shape"].append(geoshape)
            #all_results["fire_perimeters just t_min"].append(fire_perimeter_just_tmin)
            #all_results["groundtruth just t_min"].append(groundtruth_just_tmin)
            all_results["Notes"].append("")
            all_results["Total gwis Fires"].append(total_gwis_fires)
            all_results["true_positive"].append(row["true_positive"])
            all_results["false_positive"].append(row["false_positive"])
            all_results["false_negative"].append(row["false_negative"])
            all_results["true_negative"].append(row["true_negative"])
            all_results["User Accuracy"].append(row["User Accuracy"])
            all_results["Producer Accuracy"].append(row["Producer Accuracy"])
            all_results["Commission Error"].append(row["Commission Error"])
            all_results["Omission Error"].append(row["Omission Error"])
            all_results["total gt fires"].append(total_gt_fires)
            all_results["a_success %"].append(a_success_percent)
            all_results["t_success %"].append(t_success_percent)
            all_results["d_success %"].append(d_success_percent)
            all_results["a_success total"].append(a_success_total)
            all_results["t_success total"].append(t_success_total)
            all_results["d_success total"].append(d_success_total)
            all_results["same_stroke %"].append(same_stroke_percent)
            all_results["dif_stroke"].append(dif_stroke_total)
            all_results["same_stroke"].append(same_stroke_total)
            all_results["all no stroke"].append(all_no_stroke)
            all_results["t_src_dist - mean dist (meters)"].append(t_src_dist_mean)
            all_results["a_src_dist- mean dist (meters)"].append(a_src_dist_mean)
            all_results["d_src_dist- mean dist (meters)"].append(d_src_dist_mean)
    
    all_results_df = pd.DataFrame(all_results)
    all_results_df.to_excel(f"{tmin_base}all_confusion_matrices_new_entln.xlsx",index=False)
