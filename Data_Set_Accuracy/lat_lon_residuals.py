import rasterio
import ee
import csv
import functions_dsa

from rasterio.plot import show
from matplotlib import pyplot as plt

def main_lat_lon():
    
    # Connect to GEE
    functions_dsa.initialise_credentials()
    
    # Defines region of interest as a polygon (San francisco, switch)
    roi = ee.Geometry.Polygon(
      [[[-123.1205, 37.6173],
        [-123.0760, 37.6173],
        [-123.0760, 37.7125],
        [-123.1205, 37.7125]]])
    
    # Choose dataset
    dataset = ee.ImageCollection('LANDSAT/LC08/C01/T1')

    # Define the date range
    start_date = '2022-01-01'
    end_date = '2022-01-31'

    # Filter the dataset by date and region
    filtered_dataset = dataset.filterDate(start_date, end_date).filterBounds(roi)

    # Merge the images in the collection into a single image
    merged_image = filtered_dataset.mosaic()

    # Extract lat/lon data using reduceRegion
    lat_lon_data = merged_image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=roi,
        scale=30  # Resolution in meters
    )

    # Convert the lat/lon data to a Python dictionary
    lat_lon_dict = lat_lon_data.getInfo()

    # Extract lat/lon values
    latitude = lat_lon_dict['latitude']
    longitude = lat_lon_dict['longitude']

    # Export lat/lon data to a CSV file
    with open('lat_lon_data.csv', 'w', newline='') as csvfile:
        fieldnames = ['Latitude', 'Longitude']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({'Latitude': latitude, 'Longitude': longitude})

    return

if __name__ == "__main__":
   main_lat_lon()