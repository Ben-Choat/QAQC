'''
BChoat 2025/04/02

Used to extract outer boundary of huc-8 boundaries.

This script may better belong with the precip-frequency code, but
that is mostly in R, and I wanted to use geopandas due to comfort.

Using genenv311 conda environment
'''

# %% load libs define vars and etc.
##################################
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import Polygon

# file_in = "E:/Projects/OWRD_PrecipFreq/GIS/CascadesWestOregon_HUC8_100AR.shp"
# file_out = "E:/Projects/OWRD_PrecipFreq/GIS/CascadesWestOregon_100AR.shp"

# file_in = "E:/Projects/OWRD_PrecipFreq/GIS/Oregon_NAD83_HUC8.shp"
# file_out = "E:/Projects/OWRD_PrecipFreq/GIS/Oregon_HUC8_Outer.shp"

file_in = "E:/Projects/OWRD_PrecipFreq/GIS/CascadesEastOregon_OverlapWest.shp"
file_out = "E:/Projects/OWRD_PrecipFreq/GIS/CascadesEastOregon_OverlapWest_Outer.shp"



# %% read in data
##################################
gdf = gpd.read_file(file_in)

gdf['geometry'] = gdf['geometry'].make_valid()

gdf_merged = gdf.dissolve()


def get_outer_boundary(geom):
    '''Extracts only the outermost polygon boundary.'''
    if geom.geom_type == "MultiPolygon":
        # Keep only the largest polygon (outermost one)
        largest = max(geom.geoms, key=lambda g: g.area)
        return Polygon(largest.exterior)
    elif geom.geom_type == "Polygon":
        return Polygon(geom.exterior)
    return geom  # Return as-is if not a polygon

# Apply function to extract the outer boundary
gdf_outer = gdf_merged.copy()
gdf_outer["geometry"] = gdf_merged["geometry"].apply(get_outer_boundary)

gdf_outer.boundary.plot()


gdf_outer.to_file(file_out)

# for i in range(gdf.shape[0]):
#     gdf.iloc[i, :].plot()

# Merge all polygons into a single geometry
# merged_polygon = unary_union(gdf['geometry'])
# merged_polygon.boundary
# # Create a new GeoDataFrame with the merged geometry
# merged_gdf = gpd.GeoDataFrame({'geometry': [merged_polygon]}, crs=gdf.crs)
# merged_gdf.boundary.plot()

# gdf = gdf.to_crs(epsg="2992")
# gdf.plot()
# gdf = gdf.make_valid()

# gdf = gpd.GeoDataFrame(geometry=gdf.geometry)
# gdf_diss = gdf.dissolve()
# gdf_diss = gdf_diss.make_valid()

# gdf_diss.boundary.plot()


# %% Intermediate step on test of transposition that does
# not shift sea-level pressure