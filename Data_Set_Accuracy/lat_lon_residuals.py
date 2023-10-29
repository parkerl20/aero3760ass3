import rasterio
import ee
import csv
import functions_dsa

from rasterio.plot import show
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

def main_lat_lon():
    with rasterio.open('/home/tpgrep/UNI/2023/SEM_2/AERO3760_SPACE_ENG/Assignments/Assignment_3/aero3760ass3/Data_Set_Accuracy/Images/Opera_house_test.tif') as src:
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
        plt.show()
        plt.savefig('Test.png')

if __name__ == "__main__":
   main_lat_lon()