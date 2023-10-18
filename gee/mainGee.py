from . import functions

def mainGee():
    """Currently running this code will take about 300 seconds, but will make a html map of nsw in the NDVI index,
       Select a few locations and write their specific NDVI values to a csv file and then display the html map
       in your browser.
    """
    # Connect to GEE
    functions.initialise_credentials()

    # Sentinel-2A satellite
    # Map = functions.S2A("2019-12-01", "2020-01-31")
    # Map = functions.surface_temperature()
    Map = functions.fires()
    # Map = functions.S2A_NDVI("2019-12-01", "2020-01-31")
    # Create map
    # functions.create_map(Map, "NSW Infrared")
    # functions.create_map(Map, "Surface temperature")
    functions.create_map(Map, "Fires")
    # functions.create_map(Map, "NDVI")

    # Show map
    # functions.show_map("NSW Infrared")
    # functions.show_map("Surface temperature")
    functions.show_map("Fires")
    # functions.show_map("NDVI")


if __name__ == "__main__":
    mainGee()