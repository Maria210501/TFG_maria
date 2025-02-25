import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk
from tkinter import font as tkfont

# Load the data
control_df = pd.read_json('src/data/datos_control.json')
rupture_df = pd.read_json('src/data/datos_rupture.json')

control_df['status'] = 'control'
rupture_df['status'] = 'rupture'
combined_df = pd.concat([control_df, rupture_df])

# Convert to centimeters
combined_df['HYPOTHETIC AORTA LENGTH'] = combined_df['HYPOTHETIC AORTA LENGTH']/10
combined_df['ANEURYSM LENGTH'] = combined_df['ANEURYSM LENGTH']/10
combined_df['DAMAX'] = combined_df['DAMAX']/10
combined_df['D_EXT_MAX'] = combined_df['D_EXT_MAX']/10
combined_df['CENTERLINE LENGTH'] = combined_df['CENTERLINE LENGTH']/10
combined_df['ANEURYSM VOLUME'] = combined_df['ANEURYSM VOLUME']/1000
combined_df['THROMBUS VOLUME'] = combined_df['THROMBUS VOLUME']/1000
combined_df['TOTAL VOLUME'] = combined_df['TOTAL VOLUME']/1000

if 'id' in combined_df.columns:
    combined_df = combined_df.drop('id', axis=1)

variables = [col for col in combined_df.columns if col != 'status']
color_palette = {'control': '#3498db', 'rupture': '#e74c3c'}  # Updated colors

class DataVisualizationApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Data Visualization App")
        self.master.state('zoomed')  # Start maximized
        self.master.configure(bg='#f0f0f0')  # Light gray background

        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure('TFrame', background='#f0f0f0')
        self.style.configure('TLabel', background='#f0f0f0', font=('Segoe UI', 10))
        self.style.configure('TButton', font=('Segoe UI', 10))
        self.style.configure('TCombobox', font=('Segoe UI', 10))

        self.title_font = tkfont.Font(family='Segoe UI', size=16, weight='bold')

        self.create_widgets()

    def create_widgets(self):
        main_frame = ttk.Frame(self.master, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        title_label = ttk.Label(main_frame, text="Data Visualization Dashboard", font=self.title_font)
        title_label.pack(pady=(0, 20))

        # Create a frame for controls
        control_frame = ttk.Frame(main_frame, padding="10")
        control_frame.pack(fill=tk.X, pady=(0, 20))

        # Dropdown for selecting visualization type
        viz_frame = ttk.Frame(control_frame)
        viz_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(viz_frame, text="Visualization Type:").pack(anchor='w')
        self.viz_type = tk.StringVar(value="boxplot")
        viz_dropdown = ttk.Combobox(viz_frame, textvariable=self.viz_type, 
                                    values=["boxplot", "violin", "histogram", "scatter"], width=15)
        viz_dropdown.pack(anchor='w')
        viz_dropdown.bind("<<ComboboxSelected>>", self.update_controls)

        # Dropdown for selecting variable
        var_frame = ttk.Frame(control_frame)
        var_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(var_frame, text="Variable:").pack(anchor='w')
        self.variable = tk.StringVar(value=variables[0])
        self.var_dropdown = ttk.Combobox(var_frame, textvariable=self.variable, values=variables, width=15)
        self.var_dropdown.pack(anchor='w')
        self.var_dropdown.bind("<<ComboboxSelected>>", self.update_plot)

        # Dropdown for selecting second variable (initially hidden)
        self.var2_frame = ttk.Frame(control_frame)
        ttk.Label(self.var2_frame, text="Second Variable:").pack(anchor='w')
        self.variable2 = tk.StringVar(value=variables[1] if len(variables) > 1 else variables[0])
        self.var2_dropdown = ttk.Combobox(self.var2_frame, textvariable=self.variable2, values=variables, width=15)
        self.var2_dropdown.pack(anchor='w')
        self.var2_dropdown.bind("<<ComboboxSelected>>", self.update_plot)

        # Entry for number of bins (initially hidden)
        self.bins_frame = ttk.Frame(control_frame)
        ttk.Label(self.bins_frame, text="Number of Bins:").pack(side=tk.LEFT)
        self.bins_var = tk.StringVar(value="10")
        self.bins_entry = ttk.Entry(self.bins_frame, textvariable=self.bins_var, width=5)
        self.bins_entry.pack(side=tk.LEFT)
        self.bins_entry.bind("<Return>", self.update_plot)

        # Create a frame for the plot
        self.plot_frame = ttk.Frame(main_frame, padding="10")
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

        # Add an Exit button
        exit_button = ttk.Button(main_frame, text="Exit", command=self.master.quit)
        exit_button.pack(pady=(10, 0))

        self.create_plot()

    def create_plot(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        fig, ax = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor('#f0f0f0')
        ax.set_facecolor('#ffffff')
        self.plot_visualization(ax)

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, self.plot_frame)
        toolbar.update()
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)

    def update_controls(self, event=None):
        if self.viz_type.get() == "scatter":
            self.var2_frame.pack(side=tk.LEFT)
            self.bins_frame.pack_forget()
        elif self.viz_type.get() == "histogram":
            self.var2_frame.pack_forget()
            self.bins_frame.pack(side=tk.LEFT)
        else:
            self.var2_frame.pack_forget()
            self.bins_frame.pack_forget()
        self.update_plot()

    def update_plot(self, event=None):
        self.create_plot()

    def plot_visualization(self, ax):
        viz_type = self.viz_type.get()
        var = self.variable.get()
        var2 = self.variable2.get()

        if viz_type == "boxplot":
            sns.boxplot(x='status', y=var, data=combined_df, ax=ax, hue='status', palette=color_palette)
            ax.set_title(f'Boxplot of {var}', fontsize=14, fontweight='bold')
        elif viz_type == "violin":
            sns.violinplot(x='status', y=var, data=combined_df, ax=ax, hue='status', palette=color_palette, split=True)
            ax.set_title(f'Violin Plot of {var}', fontsize=14, fontweight='bold')
        elif viz_type == "histogram":
            bins = int(self.bins_var.get())
            sns.histplot(data=combined_df, x=var, hue='status', palette=color_palette, kde=True, ax=ax, bins=bins)
            ax.set_title(f'Histogram of {var}', fontsize=14, fontweight='bold')
        elif viz_type == "scatter":
            sns.scatterplot(data=combined_df, x=var, y=var2, hue='status', palette=color_palette, ax=ax)
            ax.set_title(f'Scatter Plot of {var} vs {var2}', fontsize=14, fontweight='bold')

        ax.set_xlabel(var, fontsize=12)
        if viz_type in ["boxplot", "violin"]:
            ax.set_xlabel('Status', fontsize=12)
        elif viz_type == "scatter":
            ax.set_ylabel(var2, fontsize=12)
        else:
            ax.set_ylabel('Count', fontsize=12)

        ax.tick_params(axis='both', which='major', labelsize=10)
        sns.set_style("whitegrid")
        plt.tight_layout()

if __name__ == "__main__":
    root = tk.Tk()
    app = DataVisualizationApp(root)
    root.mainloop()