import ee
import os
import requests
import geemap

from ee import batch
from functions_dsa import initialise_credentials

def export_image_to_tif_with_region(image, file_path, scale, lat, lon, band):
    """
    Export an image to a TIF file with a specified scale and region of interest using the geemap.ee_export_image() function.

    Args:
    image (ee.Image): The Earth Engine image to be exported.
    file_path (str): The file path where the TIF file will be saved.
    scale (int): The scale at which to export the image.
    lat (float): Latitude coordinate of the region of interest.
    lon (float): Longitude coordinate of the region of interest.
    band (str): The band of the image to export.

    Returns:
    str: A message indicating the success or failure of the export process.
    """
    try:
        roi = ee.Geometry.Point(lon, lat).buffer(200)  # Adjust buffer size as needed
        image = image.select(band)
        geemap.ee_export_image(image, filename=file_path, scale=scale, crs='EPSG:4326', region=roi)
        return f"Image exported successfully to {file_path}"
    except Exception as e:
        return f"An error occurred: {e}"
    

def main_dataset2image():
    # Initialize Earth Engine
    initialise_credentials()

    # Load a Landsat 8 image collection and filter by date and location
    collection = ee.ImageCollection('LANDSAT/LC08/C01/T1').filterDate('2022-01-01', '2022-01-31')

    # Select the first image in the filtered collection
    image = collection.first()

    # Specify the file path, scale, and region for export
    file_path = 'Data_Set_Accuracy/images/image.tif'
    scale = 30
    lat = -33.674  # Example latitude
    lon = 150.9151  # Example longitude
    band = ['B4', 'B3', 'B2'] # Example band name

    # Export the image
    export_result = export_image_to_tif_with_region(image, file_path, scale, lat, lon, band)
    print(export_result)

    return

if __name__ == "__main__":
   main_dataset2image()