import argparse
import ee
import pandas as pd

# Dongdong Kong
# js版本已经测试通过
# 这是GPT3.5自动翻译的代码，可能会存在bug。
# https://code.earthengine.google.com/027ef28c2f7f304be6c7ebed95f5e1cc?noload=1

# Initialize the Earth Engine library
ee.Initialize()


def image_extract_points(img, points, proj=None, set_date=False):
  if proj == None:
    proj = img.projection().getInfo()
  options = {
      'collection': points,
      'reducer': 'first',
      'crs': proj['crs'],
      'crsTransform': proj['transform'],
      'tileScale': 16
  }
  date = ee.Date(img.get('system:time_start')).format('yyyy-MM-dd')

  def tidy_props(f):
    f = ee.Feature(None).copyProperties(f)
    if set_date:
      f = f.set('date', date)
    return f
  return img.reduceRegions(**options).map(tidy_props)


def col_extract_points(col, task, points):
  proj = col.first().select(0).projection().getInfo()
  res = col.map(lambda img: image_extract_points(
      img, points, proj, set_date=True)).flatten()

  ee.batch.Export.table.toDrive(
      collection=res,
      description=task,
      folder="gee",
      fileNamePrefix=task,
      fileFormat="GeoJSON"
  ).start()


def export_image(img, task, prefix=""):
  task = prefix + task
  proj = img.projection().getInfo()

  # Export the image to Drive as a GeoTIFF file
  ee.batch.Export.image.toDrive(
      image=img,
      description=task,
      folder="gee",
      crs=proj['crs'],
      crsTransform=proj['transform'],
      maxPixels=1e13,
      fileNamePrefix=task,
      fileFormat='GeoTIFF'
  ).start()


def shp2df(fc, outfile=None):
  data = fc.getInfo().get('features')
  props = [x["properties"] for x in data]
  df = pd.DataFrame(props)
  if None != outfile:
    df.to_csv(outfile, index=False)
  return df

def hello():
  print("Hello, World!")
