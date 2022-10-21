# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to
implement functions related to visualization.
"""


from __future__ import print_function, division, absolute_import


import functools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def __compare_id(item1: str, item2: str) -> bool:
    item1_id = int(str(item1).split('_')[1])
    item2_id = int(str(item2).split('_')[1])
    return item1_id - item2_id


def plot_configuration_table(df_configuration: pd.DataFrame) -> plt:
    """
    Plots an ELIspot configuration.

    Parameters
    ----------
    df_configuration    :   DataFrame with the following columns:
                            'coverage_id'
                            'pool_id'
                            'peptide_id'
    Returns
    -------
    plt                 :   matplotlib.pyplot object.
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
    return plt

