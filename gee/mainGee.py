from . import functions
import datetime as dt
import ee
import matplotlib.pyplot as plt

def extract_roi(Map, lon_lat):
    # Assuming Map is a data structure containing the map data and lon_lat is the region of interest
    # Implement the extraction logic here based on the structure of your Map data
    # Extract the region of interest using the provided lon_lat coordinates

    # For example, if Map is a 2D array or an image, and lon_lat contains the coordinates of the region of interest
    # Extract the region based on the lon_lat coordinates
    # You may need to use some image processing techniques or slicing operations depending on the structure of Map

    # Example of extracting a region from a 2D array using slicing
    # Assuming lon_lat contains the coordinates of the region of interest in the format [x_start, x_end, y_start, y_end]
    x_start, x_end, y_start, y_end = lon_lat
    extracted_roi = Map[x_start:x_end, y_start:y_end]

    # Return the extracted region
    return extracted_roi


def mainGee(results):
    """Currently running this code will take about 300 seconds, but will make a html map of nsw in the NDVI index,
       Select a few locations and write their specific NDVI values to a csv file and then display the html map
       in your browser.
    """
    # Load in data
    r_obv = results[0]['r'][:]
    t_obv = results[0]['t'][:]
    epoch = dt.datetime(2023, 1, 1)

    # Connect to GEE
    functions.initialise_credentials()

    # Latitude and longitude conversion
    # lon_lat = functions.eci_to_llh(r_obv, t_obv, epoch)
    lon_lat = functions.eci_to_llh_nsw(r_obv, t_obv, epoch, num_points=100)

    print("Length of lon_lat:", len(lon_lat))
    print("lon_lat:", lon_lat)
    print("lon_lat[0], lon_lat[50]:", lon_lat[6], lon_lat[50])

    # Sentinel-2A satellite
    # Map = functions.S2A("2019-12-01", "2020-01-31")
    # Map = functions.fires()
    # Map = functions.S2A_NDVI("2019-12-01", "2020-01-31")
    Map = functions.S2A_coverage("2019-12-01", "2020-01-31", lon_lat, circle_radius=345088) # Radius corresponding to a 6.2 degree swathe width
    # Map = functions.plot_red_points(lon_lat, circle_radius=10)

    # Create map
    # functions.create_map(Map, "NSW Infrared")
    # functions.create_map(Map, "Fires")
    # functions.create_map(Map, "NDVI")
    functions.create_map(Map, "S2A coverage")
    # functions.create_map(Map, "Red")

    # Show map
    # functions.show_map("NSW Infrared")
    # functions.show_map("Fires")
    # functions.show_map("NDVI")
    functions.show_map("S2A coverage")
    # functions.show_map("Red")



if __name__ == "__main__":
    mainGee()