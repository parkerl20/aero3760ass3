import functions

def main():
    # Connect to GEE
    functions.initialise_credentials()

    # Sentinel-2A satellite
    # Map = functions.elevation_5m()
    # Map = functions.S2A("2019-12-01", "2020-01-31")
    Map = functions.surface_temperature()
    
    # Create map
    # functions.create_map(Map, "Elevation 5m")
    # functions.create_map(Map, "NSW Infrared")
    functions.create_map(Map, "Surface temperature")

    # Show map
    # functions.show_map("Elevation 5m")
    # functions.show_map("NSW Infrared")
    functions.show_map("Surface temperature")


if __name__ == "__main__":
    main()