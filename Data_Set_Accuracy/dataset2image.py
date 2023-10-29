import ee
import os
import requests

from ee import batch
from functions_dsa import initialise_credentials

def main_dataset2image():
    # Initialize Earth Engine
    initialise_credentials()

    # Define the region of interest (AOI) as a geometry
    aoi = ee.Geometry.Polygon(
        [[[-122.5, 37.5],
        [-122.0, 37.5],
        [-122.0, 37.75],
        [-122.5, 37.75]]])

    # Define the date range
    start_date = '2021-01-01'
    end_date = '2021-12-31'

    # Load the Landsat 8 dataset for the specified time period and location
    dataset = ee.ImageCollection('LANDSAT/LC08/C01/T1') \
        .filterBounds(aoi) \
        .filterDate(ee.Date(start_date), ee.Date(end_date)) \
        .sort('CLOUD_COVER') \
        .first()  # Get the image with the least cloud cover

    # Select the bands you want to export (e.g., Red, Green, Blue)
    bands_to_export = ['B4', 'B3', 'B2']  # RGB bands

    # Define export parameters
    export_params = {
        'image': dataset.select(bands_to_export),
        'description': 'Landsat8_Image',
        'scale': 30,  # Spatial resolution in meters per pixel
        'region': aoi.getInfo()['coordinates'],
        'crs': 'EPSG:4326',
        'fileFormat': 'GeoTIFF'  # Export format (e.g., GeoTIFF)
    }

    # Initiate the export task
    task = batch.Export.image.toDrive(**export_params)

    task.start()

    # Monitor the task status
    print('Exporting...', task)

    print(task.status())

    # Optionally, wait for the task to complete
    while task.active():
        pass

    # Check the task status
    print('Task Status:', task.status())

    # Get the download URL for the exported image
    # After the task is completed, get the download URL from the assets

    asset_id = task.status()['destination_uris'][0]
    download_url = asset_id
    
    # Define the local output path for the downloaded image
    output_path = 'Data_Set_Accuracy/images/Landsat8_Image.tif'

    # Download the exported image from the URL
    response = requests.get(download_url)
    with open(output_path, 'wb') as f:
        f.write(response.content)

    # Check if the image was successfully downloaded
    if os.path.exists(output_path):
        print('Image downloaded successfully.')
    else:
        print('Image download failed.')

    return

if __name__ == "__main__":
   main_dataset2image()