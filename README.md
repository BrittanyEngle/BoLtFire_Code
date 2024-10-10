# BoLtFire
The frequency and severity of fire weather have increased under climate change, particularly in high-latitude boreal forests. Lightning, a key ignition source globally, is also expected to become more frequent with climate change and could significantly increase burn area. Current research on lightning-ignited wildfires (LIW) has a long history in boreal ecosystems but has typically focused on North America due to better data availability, while the lack of publicly available data for Eurasia has hindered our comprehensive understanding of important characteristics of LIW, such as holdover time, lightning-ignition efficiency, frequency, and spatial distribution of lightning-ignited wildfires in boreal forests. This study introduces the Temporal Minimum Distance (TMin) method, a novel approach to matching lightning strikes with wildfires without requiring ignition location, that outperformed current methodologies. As a result, we developed a comprehensive dataset of lightning-ignited wildfires across the entire boreal forest from 2012 to 2022, encompassing 6,228 fires — 4,186 in Eurasia and 2,042 in North America — each over 200 hectares in size. This dataset provides new opportunities to model ignition and spread dynamics of boreal wildfires and offers deeper insights into lightning-driven fire activity globally. 

The code to create this dataset can be found in this account. 

The end product will result in a vector based dataset:
- **FireID:**	Unique fire identification number
- **StartDate:** Start date of the fire
- **EndDat:**	End date of the fire
- **FireYear:**	Year fire was discovered
- **AreaHa:**	Total burned area of the fire in hectares
- **ClassSize:**	Fire class size; Small <1,000 ha, Moderate 1,000 < 10,000 ha, Large 10,000 < 50,000 ha, Extremely Large 50,000 < 100,000 ha, Mega Fires ≥ 100,000 ha
- **BiomeName:**	Biome name based on Olson et al. (2001)
- **EcoBiome:**	Ecoregion Biome number based on Olson et al. (2001)
- **EcoName:**	Ecoregion name based on Olson et al. (2001)
- **EcoID:** Ecoregion ID based on Olson et al. (2001)
- **Realm:**	Realm based on Olson et al. (2001)
- **LCDN:**	Land cover type number based on Friedl and Sulla-Menashe (2022)
- **LCName:**	Land cover type name based on Friedl and Sulla-Menashe (2022)
- **Country:**	Country where fire is located based on World Bank (2020)
- **Continent:** Continent on which the fire was located
- **HoldoverD:**	Holdover in days
- **HoldoverRD:** Holdover in days, rounded to the day
- **IgnLat:**	Latitude location of the candidate lightning


