import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

OUTPUT_FILE = 'output/energies.dat'

TIME_LABEL = 'Time Step'
POTENTIAL_LABEL = 'Potential'
KINETIC_LABEL = 'Kinetic'
TOTAL_LABEL = 'Total'
TEMPERATURE_LABEL = 'Temperature'
PRESSURE_LABEL = 'Pressure'
ENERGY_DRIFT_LABEL = 'Energy Drift'


def plot_energies():
    """Plot energy evolution over time."""
    data = np.loadtxt(OUTPUT_FILE, skiprows=1)
    steps = data[:, 0]
    potential = data[:, 1]
    kinetic = data[:, 2]
    total = data[:, 3]
    temperature = data[:, 4]
    pressure = data[:, 5]

    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Molecular Dynamics Simulation Results')

    # Plot energies
    ax1.plot(steps, potential, label=POTENTIAL_LABEL)
    ax1.plot(steps, kinetic, label=KINETIC_LABEL)
    ax1.plot(steps, total, label=TOTAL_LABEL)
    ax1.set_xlabel(TIME_LABEL)
    ax1.set_ylabel('Energy')
    ax1.set_title('Energy Evolution')
    ax1.legend()
    ax1.grid(True)

    # Plot temperature
    ax2.plot(steps, temperature, 'r-')
    ax2.set_xlabel(TIME_LABEL)
    ax2.set_ylabel(TEMPERATURE_LABEL)
    ax2.set_title(TEMPERATURE_LABEL + ' Evolution')
    ax2.grid(True)

    # Plot pressure
    ax3.plot(steps, pressure, 'g-')
    ax3.set_xlabel(TIME_LABEL)
    ax3.set_ylabel(PRESSURE_LABEL)
    ax3.set_title(PRESSURE_LABEL + ' Evolution')
    ax3.grid(True)

    # Plot energy conservation
    energy_drift = abs(total - total[0]) / abs(total[0])
    ax4.plot(steps, energy_drift, 'k-')
    ax4.set_xlabel(TIME_LABEL)
    ax4.set_ylabel(ENERGY_DRIFT_LABEL)
    ax4.set_title(ENERGY_DRIFT_LABEL)
    ax4.set_yscale('log')
    ax4.grid(True)

    plt.tight_layout()

def plot_rdf():
    """Plot radial distribution function evolution."""
    rdf_files = sorted(glob('output/rdf_*.dat'))
    if not rdf_files:
        return

    # Plot the first and last RDF for comparison
    first_rdf = np.loadtxt(rdf_files[0], skiprows=1)
    last_rdf = np.loadtxt(rdf_files[-1], skiprows=1)

    plt.figure(figsize=(8, 6))
    plt.plot(first_rdf[:, 0], first_rdf[:, 1], 'b-', label='Initial')
    plt.plot(last_rdf[:, 0], last_rdf[:, 1], 'r-', label='Final')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.title('Radial Distribution Function')
    plt.legend()
    plt.grid(True)

def print_statistics():
    """Print statistical analysis of the simulation."""
    data = np.loadtxt(OUTPUT_FILE, skiprows=1)
    total_energy = data[:, 3]
    temperature = data[:, 4]
    pressure = data[:, 5]

    print("=== Statistical Analysis ===")
    print("\nEnergy Statistics:")
    print(f"Mean total energy: {np.mean(total_energy):.6f}")
    print(f"Energy standard deviation: {np.std(total_energy):.6f}")
    print(f"Max energy drift: {np.max(abs(total_energy - total_energy[0])/abs(total_energy[0])):.2e}")

    print("\nTemperature Statistics:")
    print(f"Mean temperature: {np.mean(temperature):.2f}")
    print(f"Temperature standard deviation: {np.std(temperature):.2f}")

    print("\nPressure Statistics:")
    print(f"Mean pressure: {np.mean(pressure):.2f}")
    print(f"Pressure standard deviation: {np.std(pressure):.2f}")

def main():
    """Main analysis routine."""
    # Check if required files exist
    if not os.path.exists(OUTPUT_FILE):
        print("Error: energies.dat not found in output directory!")
        return

    # Print statistics
    print_statistics()

    # Create and display plots
    plot_energies()
    if glob('output/rdf_*.dat'):
        plot_rdf()
    
    # Show all plots
    plt.show()

if __name__ == "__main__":
    main() 