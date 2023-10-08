import ee
import os
import webbrowser

def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "5-GEE/key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def show_map(map_object, map_name):
    # Directory
    directory = os.path.join("5-GEE", "maps")
    file_extension = ".html"
    
    # Build the full file path
    file_path = os.path.join(directory, map_name + file_extension)

    # Save and open the map
    map_object.save(file_path)
    webbrowser.open('file:///' + os.path.realpath(file_path))