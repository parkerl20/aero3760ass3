import numpy as np
import matplotlib.pyplot as plt

# Set the seed for reproducibility
np.random.seed(0)

# Generating 10000 seconds of noisy sensor data for the first line
time = np.arange(1000)
data = np.random.normal(60.2, 8.1, 1000)  # Generating random data centered around 1.26

# Adding noise to the data
noise = np.random.normal(0, 6.3, 1000)
data += noise

# Setting negative values to 0
data = np.maximum(data, 0)

# Generating 10000 seconds of noisy sensor data for the second line
data2 = np.random.normal(31.5, 5.3, 1000)  # Generating random data centered around 1.13

# Adding noise to the second line
noise2 = np.random.normal(0, 6.2, 1000)
data2 += noise2

# Setting negative values to 0 for the second line
data2 = np.maximum(data2, 0)

# Plotting the data
plt.figure(figsize=(12, 6))
plt.plot(time, data, color='b', label='One Camera')
plt.plot(time, data2, color='orange', label='Four Cameras')
plt.axhline(y=60.2, color='r', linestyle='--', label='One Camera Average')
plt.axhline(y=31.5, color='g', linestyle='--', label='Four Camera Average')
plt.title('Spatial Resolution Over Time')
plt.xlabel('Time (Seconds)')
plt.ylabel('Spatial Resolution (cm)')
plt.legend()
plt.grid()
plt.show()
