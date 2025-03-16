import numpy as np
import matplotlib.pyplot as plt


# Change the value to Energy_level to get wave funtion at different Energy
Energy_level = 1


L = 1e-10  # Width of the potential well in meters

# Number of points and step size
N = 1000
dx = L / N

# Create the spatial grid
x = np.linspace(0, L, N)


# Define the potential energy function
def V(x):
    return 0  # Infinite potential well


# Construct the Hamiltonian matrix
def construct_hamiltonian(N, dx, V):
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = -2 / (dx**2) + V(i * dx)
        if i > 0:
            H[i, i - 1] = -1 / (dx**2)
        if i < N - 1:
            H[i, i + 1] = -1 / (dx**2)
    return H


# Solve for the eigenvalues and eigenvectors
H = construct_hamiltonian(N, dx, V)
eigenvalues, eigenvectors = np.linalg.eigh(H)

# Plot the ground state wavefunction
ground_state_wavefunction = eigenvectors[:, Energy_level]
ground_state_wavefunction /= np.sqrt(dx)  # Normalize the wavefunction

plt.plot(x, ground_state_wavefunction)
plt.xlabel("Position (m)")
plt.ylabel("Wavefunction Amplitude")

## If you change the Energy level values ,  don't forget to change it in the plot title


plt.title("Ground State Wavefunction")
plt.show()


# To plot Probability density

plt.plot(x, ground_state_wavefunction**2)
plt.xlabel("Position (m)")
plt.ylabel("Wavefunction Amplitude")

## If you change the Energy level values don't forget to chnge it in the plot title


plt.title("Ground State Probability Density")
plt.show()
