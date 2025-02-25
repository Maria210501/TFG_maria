import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Any
import logging
import numpy as np

def load_json(file_path: str) -> Dict[str, Any]:
    """Load JSON data from a file."""
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        raise
    except json.JSONDecodeError:
        logging.error(f"Invalid JSON in file: {file_path}")
        raise

def create_dataframe(data: Dict[str, Any]) -> pd.DataFrame:
    """Create a DataFrame from JSON data."""
    return pd.DataFrame(data)

def convert_to_centimeters(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    """Convert specified columns from millimeters to centimeters."""
    for col in columns:
        if col in df.columns:
            df[col] = df[col] / 10
    return df

def convert_to_liters(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    """Convert specified columns from milliliters to liters."""
    for col in columns:
        if col in df.columns:
            df[col] = df[col] / 1000
    return df

def plot_correlation_matrix(df: pd.DataFrame, title: str):
    """Create an enhanced correlation matrix plot optimized for PDF output."""
    # Compute the correlation matrix
    corr_matrix = df.corr()

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Create a mask to hide the upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    
    # Draw the heatmap
    sns.heatmap(corr_matrix, mask=mask, cmap=cmap, vmax=1, vmin=-1, center=0,
                annot=True, fmt='.2f', square=True, linewidths=.5, cbar_kws={"shrink": .8},
                ax=ax)

    # Improve the readability of the plot
    plt.title(title, fontsize=16, pad=20)
    plt.tight_layout()
    
    # Rotate the x-axis labels
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.2, left=0.15)
    
    # Save the figure
    plt.savefig(f"{title.lower().replace(' ', '_')}.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def main():
    # Load data
    control_data = load_json('src/data/datos_control.json')
    rupture_data = load_json('src/data/datos_rupture.json')

    # Create DataFrames
    df_control = create_dataframe(control_data)
    df_rupture = create_dataframe(rupture_data)

    # Columns to convert
    length_columns = ['HYPOTHETIC AORTA LENGTH', 'ANEURYSM LENGTH', 'DAMAX', 'D_EXT_MAX', 'CENTERLINE LENGTH']
    volume_columns = ['ANEURYSM VOLUME', 'THROMBUS VOLUME', 'TOTAL VOLUME']

    # Convert units
    df_control = convert_to_centimeters(df_control, length_columns)
    df_control = convert_to_liters(df_control, volume_columns)
    df_rupture = convert_to_centimeters(df_rupture, length_columns)
    df_rupture = convert_to_liters(df_rupture, volume_columns)

    # Set index
    df_control.set_index('id', inplace=True)
    df_rupture.set_index('id', inplace=True)

    # Plot correlation matrices
    plot_correlation_matrix(df_control, 'Control Group Correlation Matrix')
    plot_correlation_matrix(df_rupture, 'Rupture Group Correlation Matrix')

    logging.info("Analysis completed successfully. Enhanced plots have been saved as PDFs and displayed.")

if __name__ == "__main__":
    main()