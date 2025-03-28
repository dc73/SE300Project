import numpy as np

print ("hello SE300")
Name = "SE300"
print(Name) 

import numpy as np
import matplotlib.pyplot as plt

# Given data for turbine expansion process
T1 = 1100  # Initial temperature (K)
T2 = 550   # Final temperature (K)
gamma = 1.25

# Generate temperature range for visualization
T_values = np.linspace(T1, T2, 100)

# Compute corresponding pressure ratio using isentropic relation
P_ratio = (T_values / T1) ** (gamma / (gamma - 1))

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(T_values, P_ratio, label="Expansion Ratio", color="b")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure Ratio (P1/P2)")
plt.title("Turbine Expansion Process")
plt.grid(True)
plt.legend()
plt.show()