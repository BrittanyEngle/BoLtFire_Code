import geopandas as gpd
import pandas as pd
import os
import numpy as np
from Utilities.Most_Used_Functions import Export_Shapefile, getGlanceCRS
from Utilities.Path_Utilities import getRootFolder
import glob
import matplotlib.pyplot as plt
from shapely.geometry import Point

def calculate_percent_intervals(results, column_name, intervals):
    """Calculate the percentage of total values falling into specified intervals."""
    data = results[column_name].dropna()
    min_val, max_val = data.min(), data.max()
    interval_width = (max_val - min_val) / intervals
    interval_edges = [min_val + i * interval_width for i in range(intervals + 1)]
    
    # Calculate counts in each interval
    bin_labels = [f'{int(i * 100 / intervals)}-{int((i + 1) * 100 / intervals)}%' for i in range(intervals)]
    binned_data = pd.cut(data, bins=interval_edges, include_lowest=True, labels=bin_labels)

    max_values = []
    for label in bin_labels:
        values_in_bin = data[binned_data == label]
        if not values_in_bin.empty:
            max_values.append(values_in_bin.max())
        else:
            max_values.append(None)
    
    return bin_labels, max_values

def calculate_summary_stats(results, columns):
    summary_stats = {'Column': [], 'Value': [], 'Count': []}
    for column in columns:
        counts = results[column].value_counts()
        for value, count in counts.items():
            summary_stats['Column'].append(column)
            summary_stats['Value'].append(value)
            summary_stats['Count'].append(count)
    return summary_stats

def add_summary_statistics(summary_stats, results, column_name, stat_name, stat_func):
    results = results[~(results["Geo_Shape"] =="No Ignition Point")]
    stat_value = stat_func(results[column_name])
    summary_stats['Column'].append(column_name)
    summary_stats['Value'].append(stat_name)
    summary_stats['Count'].append(stat_value)

def excel_to_shapefile(case_name="14days_10km", ha_filter =100, groundtruth_geoshape= 1):
    output_folder = f"{getRootFolder()}Matching\\Groundtruth\\False_unique\\{ha_filter}ha\\{case_name}\\{groundtruth_geoshape}\\"
    shp_file_path = output_folder + "Groundtruth_Merged_2012_2022_statistics.shp"
    output_excel_path = output_folder + "Groundtruth_Merged_2012_2022_converted.xlsx"
    # Read the shapefile using geopandas
    gdf = gpd.read_file(shp_file_path)
    
    # Convert the GeoDataFrame to a regular DataFrame
    df = pd.DataFrame(gdf.drop(columns='geometry'))
    
    # Optionally, you can add the geometry as WKT (Well-Known Text) if needed
    df['geometry'] = gdf['geometry'].apply(lambda geom: geom.wkt if geom else None)
    
    # Write the DataFrame to an Excel file
    df.to_excel(output_excel_path, index=False)

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
    

def merge_groundtruth_results(case_path,groundtruth_geoshape=1,ha_filter=100,):
    """
    Merge shapefiles and export the merged result.
    """
    all_groundtruths_filepaths = glob.glob(f"{case_path}LIW_*_Groundtruth_all_NA.shp")

    fires = []
    for filepath in all_groundtruths_filepaths:
        # Search for shapefiles in the directory
        groundtruth_gdf = gpd.read_file(filepath)
        groundtruth_gdf = groundtruth_gdf.to_crs(getGlanceCRS("NA"))

        #groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({'No Ignition Point'})]
        groundtruth_gdf = groundtruth_gdf[groundtruth_gdf["SIZE_HA"]>=ha_filter]
    # Filter based on Geo_Shape column
        if groundtruth_geoshape == 1:
            groundtruth_gdf = groundtruth_gdf[groundtruth_gdf["Geo_Shape"].isin({'Ignition Point Wrong Name','Single Ignition Point',"Ignition Point Wrong Name Buffer","Unmatched Ignition Point without perimeter","Unmatched Ignition Point with perimeter",'No Ignition Point'})]                                     
        elif groundtruth_geoshape == 2:
            groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({"Unmatched Ignition Point without perimeter","Unmatched Ignition Point with perimeter"})]
        elif groundtruth_geoshape == 3:
            groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({"Unmatched Ignition Point without perimeter"})]
        elif groundtruth_geoshape == 4:
            groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({"Unmatched Ignition Point with perimeter"})]
        elif groundtruth_geoshape==5:
            groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({'Ignition Point Wrong Name'})]
        elif groundtruth_geoshape==6:
            groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({'No Ignition Point'})]
        elif groundtruth_geoshape==7:
            groundtruth_gdf = groundtruth_gdf[~groundtruth_gdf["Geo_Shape"].isin({'Ignition Point Wrong Name Buffer'})]
        fires.append(groundtruth_gdf)
        
    # Concatenate GeoDataFrames ensuring the geometry column is properly handled
    merged_gdf = gpd.GeoDataFrame(pd.concat(fires, ignore_index=True), crs=getGlanceCRS("NA"), geometry='geometry')
    merged_gdf["fire_size"] = merged_gdf["SIZE_HA"].apply(categorize_fire_size)    
    # Define output folder and filename
    output_folder = f"{case_path}{ha_filter}ha/{groundtruth_geoshape}/"
    output_filename = f"Groundtruth_Merged_2012_2022" 
        
    # Export the merged shapefile
    Export_Shapefile(merged_gdf, output_folder, output_filename)

def plot_distribution(results, column_name, xaxis_name, output_folder, output_filename):
    """Plots the distribution of a specified column and saves the plot."""
    results = results[~(results["Geo_Shape"] =="No Ignition Point")]
    plt.figure(figsize=(10, 6))
    plt.hist(results[column_name].dropna(), bins=30, edgecolor='k', alpha=0.7)
    plt.title(f'Distribution of {column_name}')
    plt.xlabel(xaxis_name)
    plt.ylabel('Frequency')
    plt.grid(True)

    plot_path = os.path.join(output_folder, output_filename)
    plt.savefig(plot_path)
    plt.close()

def plot_frequency_distribution(results, col1, col2,output_folder, output_filename):
    """
    Plots the frequency distribution of two columns from a shapefile.
    """
    results = results[(~results["Geo_Shape"].isin(["No Ignition Point"]))]
    # Extract the columns of interest
    data_col1 = results[col1]
    data_col2 = results[col2].dropna()

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot histograms for both columns

    ax.plot(data_col1, label='a_src_dist', linestyle='-', marker='o')
    ax.plot(data_col2, label='a_holdover', linestyle='-', marker='x')

    # Set labels and title
    ax.set_xlabel('Index')
    ax.set_ylabel('Values')
    ax.set_title('Comparison of a_src_dist and a_holdover')

    # Add a legend
    ax.legend()


    plot_path = os.path.join(output_folder, output_filename)
    plt.savefig(plot_path)
    plt.close()

def plot_distribution_subplots(results, output_folder, output_filename):
    """
    Plots histograms for the given distributions and saves the figure to the specified output path.

    Parameters:
    - results (dict): Dictionary containing data to plot.
    - output_folder (str): Folder where the plot will be saved.
    - output_filename (str): Name of the output file.
    """
    results = results[~(results["Geo_Shape"] =="No Ignition Point")]
    # Extract the columns of interest
    t_holdover = results['t_holdrnd']
    t_src_dist = results['t_src_dist']
    a_holdover = results['a_holdrnd']
    a_src_dist = results['a_src_dist']
    d_holdover = results['d_holdrnd']
    d_src_dist = results['d_src_dist']

    # Create a figure with subplots
    fig, axes = plt.subplots(3, 2, figsize=(12, 18))
    
    # Define the subplot titles, labels, and colors
    plots_info = [
        (axes[0, 0], t_holdover, 'b', 't_holdrnd', 'Holdover [Days]', 'TMin Holdover Distribution'),
        (axes[0, 1], t_src_dist, 'r', 't_src_dist', 'Distance [km]', 'TMin Distance from Ignition Point Distribution'),
        (axes[1, 0], a_holdover, 'g', 'a_holdrnd', 'Holdover [Days]', 'MaxA Holdover Distribution'),
        (axes[1, 1], a_src_dist, 'c', 'a_src_dist', 'Distance [km]', 'MaxA Distance from Ignition Point Distribution'),
        (axes[2, 0], d_holdover, 'm', 'd_holdrnd', 'Holdover [Days]', 'DMin Holdover Distribution'),
        (axes[2, 1], d_src_dist, 'y', 'd_src_dist', 'Distance [km]', 'DMin Distance from Ignition Point Distribution')
    ]
    
    # Plot each subplot
    # Plot each subplot
    hbins = 1
    dbins = 1
    for ax, data, color, label, xlabel, title in plots_info:
        if xlabel == "Holdover [Days]": 
            hbins = max(hbins, max(data))
        elif xlabel == "Distance [km]":
            dbins = max(dbins, max(data))
     
    for ax, data, color, label, xlabel, title in plots_info:
        if xlabel == "Holdover [Days]":
            bins = int(hbins)
        else:
            bins = int(dbins/1000)
            data = data/1000

        if label == "d_holdrnd":
            ax.hist(data,bins=np.arange(-1,bins+2), alpha=0.5, color=color,label=label,rwidth=0.8)
            xticks = np.arange(-1,bins+1)
            xtickslabels = [f"+{x}" if x < 0 else str(x) for x in xticks]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtickslabels)
        elif label in ["a_holdrnd","t_holdrnd"]:
            ax.hist(data,bins=np.arange(0,bins+2), alpha=0.5, color=color,label=label,rwidth=0.8)
            xticks = np.arange(0,bins+1)
            xtickslabels = [f"+{x}" if x < 0 else str(x) for x in xticks]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtickslabels)
        else:
            ax.hist(data, bins=bins, alpha=0.5, color=color, label=label, rwidth=0.8, range = (0,bins))
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Frequency')
        ax.set_title(title)
    # Adjust layout
    plt.tight_layout()

    # Save the figure
    output_path = os.path.join(output_folder, output_filename)
    plt.savefig(output_path)
    plt.close(fig)


def add_g_spatial_column(results):
    """
    Loads a shapefile, adds a 'g_Spatial' column to indicate if the point's coordinates are within its geometry, 
    and saves the result to a new shapefile.
    """
    # Load the shapefile

    perimeter_gdf = results

    # Ensure the shapefile is using the same coordinate reference system (CRS)
    perimeter_gdf = perimeter_gdf.to_crs(epsg=4326)
    
    output = []
    for index, row in perimeter_gdf.iterrows():
        if row["Geo_Shape"] == "No Ignition Point":
            output.append(row)
            continue
        # Coordinates of the point
        latitude = row["LATITUDE"]
        longitude = row["LONGITUDE"]

        # Create a GeoDataFrame for the point
        point = Point(longitude, latitude)

        # Add the result to a new column "g_Spatial"
        row['g_Spatial'] = row.geometry.contains(point)
        output.append(row)

    output_gdf = gpd.GeoDataFrame(output, crs='EPSG:4326')
    return output_gdf

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

def groundtruth_comparison(filter_case_path):
    # Load the shapefile
    results = gpd.read_file(f"{filter_case_path}Groundtruth_Merged_2012_2022.shp")

    results["fire_size"] = results["SIZE_HA"].apply(categorize_fire_size)
    no_shape_mask = results["Geo_Shape"] =="No Ignition Point"
    count_results = len(results)
    # tmin vs amax comparison
    fail_condition_t_a = (results['fail_tmin'].isin(['temporal', 'spatial'])) & (results['fail_amax'].isin(['temporal', 'spatial']))# both fail // no stroke found
    tmin_fail_condition = (results['fail_tmin'].isin(['temporal', 'spatial'])) & (~results['fail_amax'].isin(['temporal', 'spatial'])) # no tmin stroke but amax stroke
    amax_fail_condition = (results['fail_amax'].isin(['temporal', 'spatial'])) & (~results['fail_tmin'].isin(['temporal', 'spatial']))# no amax stroke but tmin stroke
    comparison_t_a = (results['t_lat'] == results['a_lat']) & (results['t_long'] == results['a_long']) & (results['t_ts'] == results['a_ts']) # both found same stroke
    results['t_a_Comp'] = comparison_t_a.apply(lambda x: 'same_stroke' if x else 'dif_stroke')
    results.loc[fail_condition_t_a, 't_a_Comp'] = 'both no stroke'
    results.loc[tmin_fail_condition, 't_a_Comp'] = 'tmin no stroke'
    results.loc[amax_fail_condition, 't_a_Comp'] = 'amax no stroke'
    results.loc[no_shape_mask, "t_a_Comp"] = "no ignition point"

    # tmin vs dmin comparison
    fail_condition_t_d = (results['fail_tmin'].isin(['temporal', 'spatial'])) & (results['fail_dmin'].isin(['temporal', 'spatial'])) # both fail // no stroke found
    tmin_dfail_condition = (results['fail_tmin'].isin(['temporal', 'spatial'])) & (~results['fail_dmin'].isin(['temporal', 'spatial']))# no tmin stroke but dmin stroke
    dmin_fail_condition = (results['fail_dmin'].isin(['temporal', 'spatial'])) & (~results['fail_tmin'].isin(['temporal', 'spatial']))# no dmin stroke but tmin stroke
    comparison_t_d = (results['t_lat'] == results['d_lat']) & (results['t_long'] == results['d_long']) & (results['t_ts'] == results['d_ts']) # both found same stroke
    results['t_d_Comp'] = comparison_t_d.apply(lambda x: 'same_stroke' if x else 'dif_stroke')
    results.loc[fail_condition_t_d, 't_d_Comp'] = 'both no stroke'
    results.loc[tmin_dfail_condition, 't_d_Comp'] = 'tmin no stroke'
    results.loc[dmin_fail_condition, 't_d_Comp'] = 'dmin no stroke'
    results.loc[no_shape_mask, "t_d_Comp"] = "no ignition point"

    # amax vs dmin comparison
    fail_condition_a_d = (results['fail_amax'].isin(['temporal', 'spatial'])) & (results['fail_dmin'].isin(['temporal', 'spatial'])) # both fail // no stroke found
    amaxdmin_dfail_condition = (results['fail_amax'].isin(['temporal', 'spatial'])) & (~results['fail_dmin'].isin(['temporal', 'spatial']))# no amax stroke but dmin stroke
    dminamax_fail_condition = (results['fail_dmin'].isin(['temporal', 'spatial'])) & (~results['fail_amax'].isin(['temporal', 'spatial']))# no dmin stroke but amax stroke
    comparison_a_d = (results['a_lat'] == results['d_lat']) & (results['a_long'] == results['d_long']) & (results['a_ts'] == results['d_ts']) # both found same stroke
    results['a_d_Comp'] = comparison_a_d.apply(lambda x: 'same_stroke' if x else 'dif_stroke')
    results.loc[fail_condition_a_d, 'a_d_Comp'] = 'both no stroke'
    results.loc[amaxdmin_dfail_condition, 'a_d_Comp'] = 'amax no stroke'
    results.loc[dminamax_fail_condition, 'a_d_Comp'] = 'dmin no stroke'
    results.loc[no_shape_mask, "a_d_Comp"] = "no ignition point"


    # amax vs dmin comparison
    fail_condition_a_d_t = (results['fail_amax'].isin(['temporal', 'spatial'])) & (results['fail_dmin'].isin(['temporal', 'spatial']) & results['fail_tmin'].isin(['temporal', 'spatial'])) # all fail // no stroke found
    #amaxdmintmin_dfail_condition = (results['fail_amax'].isin(['temporal', 'spatial'])) & (~results['fail_dmin'].isin(['temporal', 'spatial']))& (~results['fail_tmin'].isin(['temporal', 'spatial']))# no amax stroke but dmin stroke
    #dminamaxtim_fail_condition = (results['fail_dmin'].isin(['temporal', 'spatial'])) & (~results['fail_amax'].isin(['temporal', 'spatial']))&(~results['fail_tmin'].isin(['temporal', 'spatial']))# no dmin stroke but amax stroke
    #tminamaxdim_fail_condition = (results['fail_tmin'].isin(['temporal', 'spatial'])) & (~results['fail_amax'].isin(['temporal', 'spatial']))&(~results['fail_dmin'].isin(['temporal', 'spatial']))# no dmin stroke but amax stroke
    comparison_a_d_t = ((results['a_lat'] == results['d_lat']) &  (results['d_lat'] == results['t_lat'])) & ((results['a_long'] == results['d_long']) & (results['d_long'] == results['t_long'])) & ((results['a_ts'] == results['d_ts']) & (results['d_ts'] == results['t_ts'])) # both found same stroke
    results['d_a_t_Comp'] = comparison_a_d_t.apply(lambda x: 'same_stroke' if x else 'dif_stroke')
    results.loc[fail_condition_a_d_t, 'd_a_t_Comp'] = 'all no stroke'
    results.loc[no_shape_mask, "d_a_t_Comp"] = "no ignition point"
    #results.loc[amaxdmintmin_dfail_condition, 'd_a_t_Comp'] = 'amax no stroke'
    #results.loc[dminamaxtim_fail_condition, 'd_a_t_Comp'] = 'dmin no stroke'
    #results.loc[tminamaxdim_fail_condition, 'd_a_t_Comp'] = 'tmin no stroke'

    # Define output folder and filenames
    output_folder = filter_case_path
    output_shapefile = "Groundtruth_Merged_2012_2022_statistics"
    output_excel = "Groundtruth_Merged_2012_2022_statistics_summary.xlsx"

    # Export the merged shapefile
    results = add_g_spatial_column(results)
    Export_Shapefile(results, output_folder, output_shapefile)
    
    # Calculate summary statistics
    summary_stats = calculate_summary_stats(results, ['t_a_Comp', 't_d_Comp', 'a_d_Comp', "d_a_t_Comp"])

    t_a_val_count = results["t_a_Comp"].value_counts()
    t_d_val_count = results["t_d_Comp"].value_counts()

    a_success = t_a_val_count.get("tmin no stroke",0) + t_a_val_count.get("same_stroke",0) + t_a_val_count.get("dif_stroke",0)
    a_success_percent = a_success/count_results*100
    t_success = t_a_val_count.get("amax no stroke",0) + t_a_val_count.get("same_stroke",0) + t_a_val_count.get("dif_stroke",0) 
    t_success_percent = t_success/count_results*100
    d_success = t_d_val_count.get("tmin no stroke",0) + t_d_val_count.get("same_stroke",0) + t_d_val_count.get("dif_stroke",0) 
    d_success_percent = d_success/count_results*100

    summary_stats["Column"].append("a success")
    summary_stats["Value"].append(a_success_percent)
    summary_stats["Count"].append(a_success)
    summary_stats["Column"].append("t success")
    summary_stats["Value"].append(t_success_percent)
    summary_stats["Count"].append(t_success)
    summary_stats["Column"].append("d success")
    summary_stats["Value"].append(d_success_percent)
    summary_stats["Count"].append(d_success)
    
    # Add total count to summary statistics
    total_count = len(results)
    summary_stats['Column'].append('Total GT Count')
    summary_stats['Value'].append('Total')
    summary_stats['Count'].append(total_count)
    
    # Add mean distances to summary statistics
    distance_columns = ['t_src_dist', 'a_src_dist', 'd_src_dist']
    for column in distance_columns:
        add_summary_statistics(summary_stats, results, column, 'mean dist (meters)', pd.Series.mean)
    
    # Add mean holdover times to summary statistics
    holdover_columns = ['t_holdover', 'a_holdover', 'd_holdover']
    for column in holdover_columns:
        add_summary_statistics(summary_stats, results, column, 'mean holdover rnd to day', pd.Series.mean)

    # Add mean distances to summary statistics
    distance_columns = ['t_src_dist', 'a_src_dist', 'd_src_dist']
    for column in distance_columns:
        add_summary_statistics(summary_stats, results, column, 'median dist (meters)', pd.Series.median)
    
    # Add mean holdover times to summary statistics
    holdover_columns = ['t_holdover', 'a_holdover', 'd_holdover']
    for column in holdover_columns:
        add_summary_statistics(summary_stats, results, column, 'median holdover rnd to day', pd.Series.median)

    # Create a DataFrame for summary statistics
    summary_df = pd.DataFrame(summary_stats)

    # Export summary statistics to Excel
    summary_excel_path = os.path.join(output_folder, output_excel)
    with pd.ExcelWriter(summary_excel_path) as writer:
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
    #plt.rcParams["font.family"] = 'Times New Roman'
    # Plot the distribution 
    plot_distribution(results,'t_holdrnd', 'TMin Holdover Distribution', output_folder, 't_holdover_distribution.png')
    plot_distribution(results, 't_src_dist', "TMin Distance from Ignition Point Distribution", output_folder, 't_distance_distribution.png')
    plot_distribution(results,'a_holdrnd', 'MaxA Holdover Distribution', output_folder, 'a_holdover_distribution.png')
    plot_distribution(results, 'a_src_dist', "MaxA Distance from Ignition Point Distribution", output_folder, 'a_distance_distribution.png')
    plot_distribution(results,'d_holdrnd', 'DMin Holdover Distribution', output_folder, 'd_holdover_distribution.png')
    plot_distribution(results, 'd_src_dist', "DMin Distance from Ignition Point Distribution", output_folder, 'd_distance_distribution.png')
    plot_distribution_subplots(results, output_folder, "subplot_distance_holdover_all")
