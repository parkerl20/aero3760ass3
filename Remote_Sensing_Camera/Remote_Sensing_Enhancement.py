import subprocess
import os
import shutil

from PIL import Image

def main_remote_sensing() -> None:
    # Real Time image enhancement, moves captured image into folder for image enhancement
    # Pre
    source_folder_pre = r"D:\UNi\USYD_units\3rd year\SEM 2\AERO3760\Assignment_3\aero3760ass3\Remote_Sensing_Camera\Images\Captured"
    destination_folder_pre = r"D:\UNi\USYD_units\3rd year\SEM 2\AERO3760\Assignment_3\aero3760ass3\Remote_Sensing_Camera\AI\ESRGAN\LR"

    # Post
    source_folder_post = r"D:\UNi\USYD_units\3rd year\SEM 2\AERO3760\Assignment_3\aero3760ass3\Remote_Sensing_Camera\AI\ESRGAN\results"
    destination_folder_post = r"D:\UNi\USYD_units\3rd year\SEM 2\AERO3760\Assignment_3\aero3760ass3\Remote_Sensing_Camera\Images\Enhanced"
        
    file_list_pre = [f for f in os.listdir(source_folder_pre)]
    file_list_post = [f for f in os.listdir(source_folder_post)]

    for image_name in file_list_pre:
        source_image = os.path.join(source_folder_pre, image_name)
        destination_image = os.path.join(destination_folder_pre, image_name)
        shutil.move(source_image, destination_image)
            
        # Call the test.py function to perform actual image enhancing
        subprocess.Popen(["python", r"Remote_Sensing_Camera\AI\ESRGAN\test_rsc.py"]).wait()
        
        for image_name in file_list_post:
            source_image_post = os.path.join(source_folder_post, image_name)
            destination_image_post = os.path.join(destination_folder_post, image_name)
            shutil.move(source_image_post, destination_image_post)

    return

if __name__ == "__main__":
    main_remote_sensing()


