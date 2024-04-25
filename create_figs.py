import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import globals
from globals import feature_dict

def craete_subfigs(seq_length, indel_rate, result_path, modules, number_of_simulations, fig_name):
    feature_df = pd.DataFrame()
    files = []
    for mod in modules:
        file = os.path.join(result_path, mod.upper(), "features.csv")
        df = pd.read_csv(file, index_col=0)
        df["Module"] = [f"{mod}"] * number_of_simulations
        feature_df = pd.concat([feature_df, df], ignore_index=True)
    feature_df.to_csv(os.path.join(result_path, "all_modules_features.csv"), index=False)
    annotation_str = f"Indel Rate:{indel_rate}\nSequence Length:{seq_length}"
    rows = 4
    cols = 4
    fig, axs = plt.subplots(rows,cols, figsize=(15, 15))  # Create a 4x4 grid for 16 subplots
    axs = axs.flatten()  # Flatten the subplot array to access each subplot easily

    # Iterate over the first 13 features and plot on each subplot
    for i, (feature_number, feature_name) in enumerate(list(feature_dict.items())[:13]):
            hist = sns.histplot(data=feature_df, x=f'{feature_number}', hue='Module', alpha=0.7, ax=axs[i], legend=True)  # Plot the histogram
            axs[i].set_title(f"{feature_name}")
            axs[i].spines['top'].set_visible(True)  # Show top spine
            axs[i].spines['right'].set_visible(True)  # Show right spine
            axs[i].set_ylabel('')  # Hide the y-axis label
            axs[i].set_xlabel('')  # Hide the x-axis label

    # Remove the last two subplots (index 13 and 14) since they are not needed
    j = len(feature_dict)
    while j<rows*cols:
        fig.delaxes(axs[j])
        j+=1

    plt.subplots_adjust(hspace=0.5)  # Adjust the vertical space between subplots
    plt.tight_layout()
    plt.savefig(f"{result_path}/{fig_name}_figures.png")
    plt.show()

    
