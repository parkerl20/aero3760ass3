import ee
import geemap
import folium
import webbrowser
import os
import functions

# Initialise credentials
functions.initialise_credentials()

# Set some points of interest
# Mackay, Qld
lat = -21.1434
lon = 149.1868

# Get URL for our image
Data_Set_Clearing = ee.ImageCollection("WRI/GFW/FORMA/raw_output_firms")

# Create an interactive map
# Map = folium.Map(center=[lat, lon], zoom_start=7)
Map = geemap.Map(center=[lat, lon], zoom=7)

# Define the dataset
dataset = ee.ImageCollection('WRI/GFW/FORMA/raw_output_firms')
dataset.filter(ee.Filter.date('2018-08-01', '2018-08-15'))

# Select the 'nday' band
percentageOfClearing = dataset.select('nday')

# Visualization parameters
visParams = {
'min': 0.0,
'max': 0.01
}

# Add layer to map
Map.addLayer(percentageOfClearing, visParams, 'Percentage of clearing')

# Show the map
functions.show_map(Map, "initial_map")