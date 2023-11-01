import numpy as np
import matplotlib.pyplot as plt

# Set the seed for reproducibility
np.random.seed(0)

# Generating 10000 seconds of noisy sensor data
time = np.arange(1000)
data = np.random.normal(1.26, 0.5, 1000)  # Generating random data centered around 1.26

# Adding noise to the data
noise = np.random.normal(0, 0.2, 1000)
data += noise

# Plotting the data
plt.figure(figsize=(12, 6))
plt.plot(time, data, color='b', label='Sensor Data')
plt.axhline(y=1.26, color='r', linestyle='--', label='Center Line')
plt.title('Mapping Accuracy Over Time')
plt.xlabel('Time (Seconds)')
plt.ylabel('Mapping Accuracy (m)')
plt.legend()
plt.grid()
plt.show()
