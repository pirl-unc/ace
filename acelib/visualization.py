#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to visualization.

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan

Last updated date: Aug 18, 2022
"""


import networkx as nx
import functools
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from networkx.drawing.nx_pydot import graphviz_layout
from matplotlib.colors import LinearSegmentedColormap


class Pool:

    def __init__(self, name, prev_pool_ids):
        self.__name = name
        self.__prev_pool_ids = prev_pool_ids

    def get_name(self):
        return self.__name


def __compare_id(item1, item2):
    item1_id = int(str(item1).split('_')[1])
    item2_id = int(str(item2).split('_')[1])
    return item1_id - item2_id


def plot_configuration_graph(df_configuration):
    """
    Plots configuration.

    Args
    ----
    df_configuration        :   DataFrame with the following columns:
                                'coverage_id', 'pool_id', 'peptide_id'

    """
    # Step 1. Generate a graph
    G = nx.Graph()
    for curr_peptide_id in df_configuration['peptide_id'].unique():
        df_matched = df_configuration.loc[df_configuration['peptide_id'] == curr_peptide_id,:]
        curr_pool_ids = df_matched['pool_id'].values.tolist()
        curr_pool_ids = sorted(curr_pool_ids, key=functools.cmp_to_key(__compare_id))
        curr_pool_ids = [i.replace('pool_', 'g') for i in curr_pool_ids] # rename pools
        curr_peptide_id = curr_peptide_id.replace('peptide_', 'p')

        for j in range(1, len(curr_pool_ids)):
            curr_pool_id = curr_pool_ids[j]
            prev_pool_ids = []
            for k in range(0, j):
                prev_pool_ids.append(curr_pool_ids[k])
            node_a = '_'.join(prev_pool_ids)
            node_b = '_'.join(prev_pool_ids) + '_' + curr_pool_id
            G.add_edge(node_a, node_b)

            if j == len(curr_pool_ids) - 1:
                G.add_edge(node_b, curr_peptide_id)

    # Step 2. Relabel nodes
    labels_dict = {}
    for curr_node in G.nodes():
        if '_' in curr_node:
            curr_node_name = str(curr_node).split('_')[-1]
            labels_dict[str(curr_node)] = curr_node_name
        else:
            labels_dict[str(curr_node)] = str(curr_node)

    # Step 3. Map colors
    palette = sns.color_palette("bright", len(df_configuration['pool_id'].unique()))
    palette_idx = 0
    node_colors_dict = {}
    for curr_node in G.nodes():
        if 'p' != str(curr_node)[0]:
            curr_node_id = str(curr_node).split('_')[-1]
            if curr_node_id not in node_colors_dict.keys():
                node_colors_dict[curr_node_id] = palette[palette_idx]
                palette_idx += 1

    node_colors = []
    for curr_node in G.nodes():
        if 'p' != str(curr_node)[0]:
            curr_node_id = str(curr_node).split('_')[-1]
            node_colors.append(node_colors_dict[curr_node_id])
        else:
            node_colors.append('#E32633')

    # Step 3. Visualize graph
    pos = graphviz_layout(G, prog="dot")
    nx.draw(G, pos,
            node_size=120,
            font_size=4,
            font_color='white',
            with_labels=True,
            labels=labels_dict,
            node_color=node_colors)
    plt.show()



def plot_configuration_table(df_configuration, output_pdf_file, save_figure=False):
    """
    Plots configuration.

    Args
    ----
    df_configuration        :   DataFrame with the following columns:
                                'coverage_id', 'pool_id', 'peptide_id'
    output_pdf_file         :   Output PDF file.
    save_figure             :   Save figure as a PDF file if True (default: False).
    """
    # Step 1. Generate matrix
    peptide_ids = sorted(df_configuration['peptide_id'].unique(), key=functools.cmp_to_key(__compare_id))
    pool_ids = sorted(df_configuration['pool_id'].unique(), key=functools.cmp_to_key(__compare_id))
    df_matrix = pd.DataFrame(columns=pool_ids, index=peptide_ids)
    df_matrix = df_matrix.fillna(0)
    for curr_peptide_id in peptide_ids:
        df_matched = df_configuration.loc[df_configuration['peptide_id'] == curr_peptide_id,:]
        for curr_pool_id in df_matched['pool_id'].values.tolist():
            df_matrix.loc[curr_peptide_id, curr_pool_id] = 1

    # Step 2. Plot
    figure_height = 8.27 + (0.5 * (len(peptide_ids) - 25))
    figure_width = 11.7 + (1.0 * (len(pool_ids) - 15))
    sns.set(font_scale=0.75)
    sns.set(rc={'figure.figsize': (figure_width, figure_height)})
    ax = sns.heatmap(data=df_matrix,
                     annot=False,
                     fmt='.1f',
                     linewidths=1,
                     linecolor='black',
                     cbar=False,
                     cmap=LinearSegmentedColormap.from_list('Custom', ['#FFFFFF', '#d23b68'], 2))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.yticks(rotation=0)
    plt.tight_layout()
    if save_figure:
        plt.savefig(output_pdf_file)
    else:
        plt.show()

