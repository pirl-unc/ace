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


import tkinter
import tkinter.filedialog
import eel
import pandas as pd
import multiprocessing
from acelib.main import run_ace_generate, run_ace_identify
from acelib import utils
from acelib import third_party_api
from acelib import constants
from acelib.sequence_features import *


eel.init('views')


# root.withdraw()
# root.wm_attributes('-topmost', True)


@eel.expose
def download_xls_bttn():
    root = tkinter.Tk()
    root.wm_withdraw()
    root.wm_attributes('-topmost', True)
    child_window = tkinter.Toplevel(root)
    child_window.title("Export File")
    child_window.focus()
    child_window.update()
    child_window.destroy()
    file_path = tkinter.filedialog.asksaveasfilename(
        title="Save ACE Configuration File",
        defaultextension='.xlsx'
    )
    root.update()
    root.destroy()
    return file_path


@eel.expose
def upload_csv_bttn():
    root = tkinter.Tk()
    root.wm_withdraw()
    root.attributes('-topmost', True)
    child_window = tkinter.Toplevel(root)
    child_window.title("Select File")
    child_window.focus()
    child_window.update()
    child_window.destroy()
    file_path = tkinter.filedialog.askopenfilename(
        title="Select Peptide List (*CSV)",
        filetypes=(("CSV Files","*.csv"),)
    )
    root.update()
    root.destroy()
    return file_path


@eel.expose
def upload_xls_bttn():
    root = tkinter.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', True)
    child_window = tkinter.Toplevel(root)
    child_window.title("Export File")
    child_window.focus()
    child_window.update()
    child_window.destroy()
    file_path = tkinter.filedialog.askopenfilename(
        title="Select Peptide List (*XLS)",
        filetypes=(("Excel Files","*.xlsx"),
                   ("Excel Files","*.xls"))
    )
    root.update()
    root.destroy()
    return file_path


@eel.expose
def write_elispot_configuration_excel(configuration_dict, file_path):
    # Step 1. ELIspot configuration
    df = pd.DataFrame(configuration_dict)

    # Step 2. ELIspot configuration for bench-side
    data = {
        'peptide_id': [],
        'sequence': [],
        'plate_unique_ids': []
    }
    for curr_peptide_id in df['peptide_id'].unique():
        df_matched = df.loc[df['peptide_id'] == curr_peptide_id,:]
        curr_sequence = df_matched['sequence'].values.tolist()[0]
        curr_plate_unique_ids = ','.join(df_matched['plate_unique_id'].unique())
        data['peptide_id'].append(curr_peptide_id)
        data['sequence'].append(curr_sequence)
        data['plate_unique_ids'].append(curr_plate_unique_ids)
    df_bench_output = pd.DataFrame(data)

    # Step 3. Write to file
    df_bench_output["order"] = df_bench_output["peptide_id"].str.split("_", 1).str[1].astype(int)
    df_bench_output = df_bench_output.sort_values(by=['order'], ascending=True)
    df_bench_output = df_bench_output.loc[:,['peptide_id', 'sequence', 'plate_unique_ids']]
    with pd.ExcelWriter(file_path) as writer:
        df.to_excel(writer, sheet_name="ACE_config", index=False)
        df_bench_output.to_excel(writer, sheet_name="Bench_ready", index=False)

    return True


@eel.expose
def read_peptide_sequences_csv_file(file_path):
    df = pd.read_csv(file_path)
    peptide_sequences = df['Peptide'].values.tolist()
    return peptide_sequences


@eel.expose
def download_probert_data() -> bool:
    try:
        third_party_api.get_bert_model(model=constants.ROSTLAB_PROT_BERT)
    except:
        return False
    return True


@eel.expose
def generate_configuration(num_peptides,
                           num_peptides_per_pool,
                           num_coverage,
                           num_cores,
                           peptide_sequences):
    num_peptides = int(num_peptides)
    num_peptides_per_pool = int(num_peptides_per_pool)
    num_coverage = int(num_coverage)
    num_cores = int(num_cores)
    disallowed_peps = get_disallowed_peptides(peptide_sequences, dummy_embed, ratio_similarity, 0.8)
    df_configuration = run_ace_generate(
        n_peptides=num_peptides,
        n_peptides_per_pool=num_peptides_per_pool,
        n_coverage=num_coverage,
        num_processes=num_cores,
        disallowed_peptide_pairs=disallowed_peps
    )
    df_configuration = utils.assign_96_well_plate_physical_ids(
        df_configuration=df_configuration,
        peptide_sequences=peptide_sequences
    )
    return df_configuration.to_dict()


def reformat_plate_reader(plate_readout):
    spot_counts = plate_readout.to_numpy().flatten()
    pool_ids = ["pool_" + str(i) for i in np.arange(1, len(spot_counts)+1)]
    d = {'pool_id': pool_ids, 'spot_count': spot_counts}
    return pd.DataFrame(d)


@eel.expose
def ace_identify_helper(
        dict_readout,
        dict_config,
        min_spot_count
    ):
    df_readout = pd.DataFrame(dict_readout)
    df_config = pd.DataFrame(dict_config)
    df_hits = run_ace_identify(df_readout, df_config, min_spot_count)
    return df_hits.to_dict()


@eel.expose
def identify_positives(plate_readout_path,
                       config_path):
    plate_readout = pd.read_excel(
                                plate_readout_path,
                                skiprows=2,
                                index_col="Unnamed: 0",
                                nrows=8
    )
    config_df = pd.read_excel(config_path)
    df_hits = run_ace_identify(
                                df_readout=reformat_plate_reader(plate_readout),
                                df_configuration=config_df
    )

    return (df_hits.to_dict(), config_df.to_dict(), reformat_plate_reader(plate_readout).to_dict())


if __name__ == '__main__':
    multiprocessing.freeze_support()
    eel.start('landing.html', size=(1920, 1080), port=utils.get_open_port())
