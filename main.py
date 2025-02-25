from Groundtruth_Code.Polygons import Alaska_Polygons,Canada_Filter_Polygons
from Groundtruth_Code import Groundtruth_Statistics, Merge_Groundtruths
from Lightning_to_Fire_Code import ENTLN_Filter,Matching_Functions,GWIS_Filter_Raster, Matching_Functions_Groundtruth
from main import Validation
from Utilities.Path_Utilities import getMatchingCasePath,getMatchingFilterPath
from Utilities.Most_Used_Functions import generate_public_dataset

#------------- ENTLN
ENTLN_Filter.entln_filter()
ENTLN_Filter.entln_continent_split()
ENTLN_Filter.entln_landcover_country_split()
ENTLN_Filter.count_forest_and_landcover_rows_in_shapefiles()

#------------- GROUNDTRUTH FILTERING
Alaska_Polygons.alaska_filter_points()
Alaska_Polygons.alaska_filter_polygons()
Alaska_Polygons.alaska_merge_polygons()
Canada_Filter_Polygons.canada_filter_points()
Canada_Filter_Polygons.canada_filter_polygons()
Canada_Filter_Polygons.canada_merge_polygons()
Canada_Filter_Polygons.canada_convert_to_noon()
Merge_Groundtruths.merge_groundtruths_raster()
Matching_Functions_Groundtruth.split_into_years()


#------------- GROUNDTRUTH MATCHIND AND STATS
day_filter = 14 #days
distance_fitler = 10 #km
Matching_Functions_Groundtruth.matching_function_all(time_delta_days=day_filter,distance_delta_km=distance_fitler)

ha_filter_analysis = [200]
geoshape_filter = [6]
case_path = getMatchingCasePath(True,day_filter,distance_fitler)
for ha in ha_filter_analysis:
    for geoshape in geoshape_filter:
        Groundtruth_Statistics.merge_groundtruth_results(case_path=case_path, ha_filter=ha, groundtruth_geoshape= geoshape)
        filter_case_path = getMatchingFilterPath(day_filter,distance_fitler,ha,geoshape)
        Groundtruth_Statistics.groundtruth_comparison(filter_case_path)

#------------- GWIS FILTER
GWIS_Filter_Raster.gwis_filter()
#------------- GWIS MATCHING
Matching_Functions.matching_function_polygon_tmin(tmax=day_filter,smax=distance_fitler)

#------------- GWIS VALIDATION
Validation.merge_results_by_continent(getMatchingCasePath(False,day_filter,distance_fitler))
for ha_filter in ha_filter_analysis:
    for geoshape in geoshape_filter:
        groundtruth_path = getMatchingFilterPath(day_filter,distance_fitler,ha_filter,geoshape)+"Groundtruth_Merged_2012_2022_statistics.shp"
        Validation.compute_confusion_matrix(getMatchingCasePath(False,day_filter,distance_fitler),groundtruth_file=groundtruth_path,ha_filter=ha_filter)

generate_public_dataset(getMatchingCasePath(False,day_filter,distance_fitler))
Validation.mergeConfusionMatrices()