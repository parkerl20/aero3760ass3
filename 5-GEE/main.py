import functions

def main():
    # Connect to GEE
    functions.initialise_credentials()

    # Sentinel-2A satellite
    Map = functions.S2A()
    
    # Create map
    functions.create_map(Map, "S2A-RGB-Winter")

    # Show map
    functions.show_map("S2A-RGB-Winter")


if __name__ == "__main__":
    main()