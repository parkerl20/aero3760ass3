import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

from spacesim import satellite as sat
from spacesim import celestial_body as cb
from spacesim import constants as const
from spacesim import estimation as est
from spacesim import sensor

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

def main() -> None:
    pass

if __name__ == "__main__":
    main()