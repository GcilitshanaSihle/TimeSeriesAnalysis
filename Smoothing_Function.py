import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

def plot_flux_time_error(flux, timestamp, error):
    # Assuming 'flux', 'time', and 'errors' are your data arrays
    
    # Create the plot
    plt.figure(figsize=(25, 6))
    
    # Plot the flux with error bars
    plt.errorbar(timestamp, flux*1e+6, yerr=error*1e+6, fmt='o', markersize=2, color='b', ecolor='r', capsize=0, elinewidth=0.5)
    
    # Add labels and title
    plt.xlabel('Time')
    plt.ylabel('Flux')
    plt.title('Flux vs. Time')
    # Show the plot
    plt.show()


def plot_flux_time_error_smoothed(flux_list, timestamp, labels, window_length=301, polyorder=3):
    # Assuming 'flux_list', 'timestamp', 'error_list', and 'labels' are your data arrays
    
    # Create the plot
    plt.figure(figsize=(30, 6))
    
    for i, (flux, error, label) in enumerate(zip(flux_list, error_list, labels)):
        # Apply Savitzky-Golay filter to each set of data
        smoothed_flux = savgol_filter(flux, window_length, polyorder)
        
        # Plot the smoothed flux for each light curve with error bars
        plt.errorbar(timestamp, smoothed_flux*1e+6, fmt='o', markersize=2, capsize=0, elinewidth=0.5, label=f'Light Curve {label}')
    
    # Add labels and title
    plt.xlabel('Time')
    plt.ylabel('Smoothed Flux')
    plt.title('Smoothed Flux vs. Time')
    
    # Add legend
    plt.legend()
    
    # Show the plot
    plt.show()
