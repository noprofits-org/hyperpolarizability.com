import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd

# Set up figure parameters for publication-quality output
plt.rcParams['font.family'] = 'DejaVu Serif'  # Use a font that's likely installed
plt.rcParams['font.size'] = 12
plt.rcParams['figure.figsize'] = (10, 12)

# Virial coefficients (cm³/mol) for different gases at different temperatures
# Data based on literature values
virial_data = {
    'N2': {
        '200K': {'B': -160.0, 'C': 4800},
        '273K': {'B': -10.5, 'C': 1200},
        '300K': {'B': 4.2, 'C': 1050},
        '400K': {'B': 21.3, 'C': 800},
        '600K': {'B': 35.1, 'C': 520}
    },
    'Ar': {
        '200K': {'B': -187.0, 'C': 5600},
        '273K': {'B': -21.7, 'C': 1200},
        '300K': {'B': -5.3, 'C': 950},
        '400K': {'B': 19.4, 'C': 760},
        '600K': {'B': 32.6, 'C': 490}
    },
    'CO2': {
        '273K': {'B': -142.0, 'C': 3400},
        '300K': {'B': -114.8, 'C': 3100},
        '400K': {'B': -49.2, 'C': 2300},
        '500K': {'B': -13.3, 'C': 1800},
        '600K': {'B': 11.9, 'C': 1400}
    },
    'CH4': {
        '200K': {'B': -148.2, 'C': 4200},
        '273K': {'B': -53.6, 'C': 2300},
        '300K': {'B': -42.3, 'C': 2100},
        '400K': {'B': -9.1, 'C': 1400},
        '600K': {'B': 28.4, 'C': 780}
    }
}

# Gas constants
R = 0.08314  # L·bar/mol·K (convenient for pressure in bar)

# Function to calculate Z using the virial equation with B and C
def calculate_Z(pressure, temperature, B, C):
    # First calculate molar volume of ideal gas at these conditions
    Vm_ideal = R * temperature / pressure
    
    # Calculate Z using the virial equation truncated after C term
    Z = 1 + (B * 1e-3) / Vm_ideal + (C * 1e-6) / (Vm_ideal**2)
    
    return Z, Vm_ideal

# Pressure range for calculation (in bar)
pressures = np.linspace(1, 200, 100)

# Create the figure with subplots
fig = plt.figure()
gs = GridSpec(2, 2, figure=fig)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), 
        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

# Colors for different temperatures
colors = {
    '200K': 'blue',
    '273K': 'green',
    '300K': 'red',
    '400K': 'purple',
    '500K': 'orange',
    '600K': 'brown'
}

# Store compression factor values at 100 bar for the table
table_data = {}

# Plot Z vs pressure for each gas
for idx, (gas, temp_data) in enumerate(virial_data.items()):
    ax = axes[idx]
    table_data[gas] = {}
    
    for temp, coef in temp_data.items():
        B = coef['B']
        C = coef['C']
        Z_values = []
        temperature = float(temp[:-1])  # Remove 'K' and convert to float
        
        # Calculate Z for each pressure
        for p in pressures:
            Z, Vm = calculate_Z(p, temperature, B, C)
            Z_values.append(Z)
            
            # Store Z at 100 bar for the table
            if abs(p - 100) < 1.0:  # Close to 100 bar
                table_data[gas][temp] = Z
        
        # Plot Z vs pressure
        ax.plot(pressures, Z_values, label=f"{temp} (B={B})", color=colors[temp], 
                linestyle='-' if B > 0 else '--')
    
    # Add horizontal line at Z=1 (ideal gas behavior)
    ax.axhline(y=1, color='black', linestyle=':', alpha=0.7)
    
    # Customize the plot
    ax.set_title(f"{gas}")
    ax.set_xlabel("Pressure (bar)")
    ax.set_ylabel("Compression Factor (Z)")
    ax.set_ylim(0.4, 1.6)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

plt.tight_layout()

# Save the figure
plt.savefig('compression_factor_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('compression_factor_comparison.pdf', bbox_inches='tight')

# Create a DataFrame for the tabular data
# First, get all temperatures used
all_temps = sorted(list(set([temp for gas_data in virial_data.values() for temp in gas_data.keys()])))

# Create a DataFrame
df_data = []
for gas in virial_data.keys():
    row = [gas]
    for temp in all_temps:
        row.append(table_data[gas].get(temp, float('nan')))
    df_data.append(row)

columns = ['Gas'] + [f"Z at {temp}, 100 bar" for temp in all_temps]
df = pd.DataFrame(df_data, columns=columns)

# Save to CSV and create markdown table
df.to_csv('compression_factor_data.csv', index=False, float_format='%.4f')

# Create a markdown table
with open('compression_factor_table.md', 'w') as f:
    f.write(df.to_markdown(index=False))

# Print the table
print("\nCompression Factor (Z) Values at 100 bar:\n")
print(df.to_string(index=False, na_rep='-', float_format='%.4f'))

print("\nFigures and data saved successfully.")
