# %%
## add pwd to path
# import os
# import sys
# sys.path.append(os.getcwd())

from ee_extract import *


# %%
import argparse
import ee
import pandas as pd

# Dongdong Kong
# js版本已经测试通过
# 这是GPT3.5自动翻译的代码，可能会存在bug。
# https://code.earthengine.google.com/027ef28c2f7f304be6c7ebed95f5e1cc?noload=1

# Initialize the Earth Engine library
ee.Initialize()

# %%
# Define the function to extract points from an image (internal)

def shp2df(fc, outfile=None):
  data = fc.getInfo().get('features')
  props = [x["properties"] for x in data]
  df = pd.DataFrame(props)
  if None != outfile:
    df.to_csv(outfile, index=False)
  return df

# %%
# points = ee.FeatureCollection("users/kongdd/shp/flux-212")
points = ee.FeatureCollection("projects/gee-hydro/assets/wuhan").select(["site"])

# ee.ImageCollection("MODIS/061/MCD12Q1")

# img = "OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02"
img = "users/kongdd/BEPS/CI_240X_1Y_V1"
img = ee.Image(img)

# proj = img.projection().getInfo()
r = image_extract_points(img, points)
r.getInfo()
df = shp2df(r, "flux212_ClampingIndex.csv")
df

# %%
# Define the function to export an image to Drive

def main():
    parser = argparse.ArgumentParser(description="Google Earth Engine Command Line Tool")
    
    # Add command line arguments
    parser.add_argument('--action', required=True, choices=['img_extract_points', 'map_key', 'export_image', 'col_extract_points'],
                        help='Specify the action to perform.')
    parser.add_argument('--image', help='Path to the image for img_extract_points or export_image actions.')
    parser.add_argument('--task', help='Task name for img_extract_points or export_image actions.')
    parser.add_argument('--points', help='Path to the points for img_extract_points or col_extract_points actions.')
    
    args = parser.parse_args()

    if args.action == 'img_extract_points':
        img = ee.Image(args.image)
        task_options = {'points': ee.FeatureCollection(args.points)}
        image_extract_points(img, args.task, task_options)
    elif args.action == 'map_key':
        # Define your 'images' dictionary and 'options' here
        # map_key(images, image_extract_points, options)
    elif args.action == 'export_image':
        img = ee.Image(args.image)
        task_options = {'points': ee.FeatureCollection(args.points)}
        export_image(img, args.task, task_options)
    elif args.action == 'col_extract_points':
        col = ee.ImageCollection(args.image)
        col_extract_points(col, args.task, ee.FeatureCollection(args.points))

if __name__ == "__main__":
    main()
