import functions

def main():
    # Connect to GEE
    functions.initialise_credentials()

    # Sentinel-2A satellite
    Map = functions.S2A()
    
    # Show map
    functions.show_map(Map, "S2A-RGB")


if __name__ == "__main__":
    main()