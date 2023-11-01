import matplotlib.pyplot as plt

# Create a blank plot
plt.figure()

# Optionally, add a title and labels
plt.title('Data Points on Plot')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

# Add data points (extremities)
x_values = [1, 3, 5, 7, 9]
y_values = [2, 4, 3, 6, 8]

# Use plt.scatter() to add data points
plt.scatter(x_values, y_values, marker='o', color='red', label='Extremities')

# Optionally, add a legend
plt.legend()

# Show the plot
plt.savefig('')
plt.show()