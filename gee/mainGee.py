from . import functions

def mainGee(r_obv, t_obv, epoch):
    """Currently running this code will take about 300 seconds, but will make a html map of nsw in the NDVI index,
       Select a few locations and write their specific NDVI values to a csv file and then display the html map
       in your browser.
    """
    # Connect to GEE
    functions.initialise_credentials()

    # Latitude and longitude conversion
    lon_lat = functions.eci_to_llh(r_obv, t_obv, epoch)

    # Sentinel-2A satellite
    # Map = functions.S2A("2019-12-01", "2020-01-31")
    # Map = functions.surface_temperature()
    # Map = functions.fires()
    # Map = functions.S2A_NDVI("2019-12-01", "2020-01-31")
    Map = functions.S2A_coverage("2019-12-01", "2020-01-31", lon_lat)
    # Map = functions.heat_map(lon_lat)

    # Create map
    # functions.create_map(Map, "NSW Infrared")
    # functions.create_map(Map, "Surface temperature")
    # functions.create_map(Map, "Fires")
    # functions.create_map(Map, "NDVI")
    functions.create_map(Map, "S2A coverage")
    # functions.create_map(Map, "Heatmap")

    # Show map
    # functions.show_map("NSW Infrared")
    # functions.show_map("Surface temperature")
    # functions.show_map("Fires")
    # functions.show_map("NDVI")
    functions.show_map("S2A coverage")
    # functions.show_map("Heatmap")


if __name__ == "__main__":
    mainGee()