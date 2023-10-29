import ee
import geemap
import rasterio
import matplotlib.pyplot as plt


def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def export_image_to_tif_with_region(image, file_path, scale, lat, lon, band):
    """
    Export an image to a TIF file with a specified scale and region of interest using the geemap.ee_export_image() function.

    Args:
    image (ee.Image): The Earth Engine image to be exported.
    file_path (str): The file path where the TIF file will be saved.
    scale (int): The scale at which to export the image.
    lat (float): Latitude coordinate of the region of interest.
    lon (float): Longitude coordinate of the region of interest.
    band (str): The band of the image to export.

    Returns:
    str: A message indicating the success or failure of the export process.
    """
    try:
        roi = ee.Geometry.Point(lon, lat).buffer(10000)  # Adjust buffer size as needed
        image = image.select(band)
        geemap.ee_export_image(image, filename=file_path, scale=scale, region=roi)
        return f"Image exported successfully to {file_path}"
    except Exception as e:
        return f"An error occurred: {e}"
    


def save_to_png(file_path):
    with rasterio.open(file_path) as src:
        # Read the image as a numpy array
        img = src.read(1)  # Replace '1' with the band you want to visualize

        # Get the spatial transform (geotransform)
        transform = src.transform

        # Get the spatial extent (bounding box)
        extent = [transform[2], transform[2] + transform[0] * src.width,
                transform[5] + transform[4] * src.height, transform[5]]

        # Create a figure and plot the image
        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(img, extent=extent)

        # Add a colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Value')

        # Set axis labels
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        # Show the plot
        print("Plot saved!")
        plt.show()
        plt.savefig('figures/sydney.png')



def main():

    initialise_credentials()

    # Example usage
    # Load an image (you should replace this with your own Earth Engine image)
    image = ee.Image('USGS/SRTMGL1_003')

    # Specify the file path, scale, and region for export
    file_path = 'tifs/image.tif'
    scale = 30
    lat = -33.7749  # Example latitude
    lon = 151.4194  # Example longitude
    band = 'elevation'  # Example band name

    # Export the image
    export_result = export_image_to_tif_with_region(image, file_path, scale, lat, lon, band)
    print(export_result)

    save_to_png(file_path)

    return 0


if __name__ == "__main__":
    main()
