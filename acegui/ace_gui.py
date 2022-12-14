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
from acelib.main import run_ace_generate, run_ace_identify
from acelib import utils
from acelib import third_party_api
from acelib import constants
from acelib.sequence_features import *


eel.init('views')


def get_tk_root():
    root = tkinter.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)


@eel.expose
def upload_csv_bttn():
    get_tk_root()
    return tkinter.filedialog.askopenfilename(title="Select Peptide List (*CSV)",
                                              filetypes=(("CSV Files","*.csv"),))

@eel.expose
def upload_xls_bttn():
    get_tk_root()
    print(get_tk_root())
    return tkinter.filedialog.askopenfilename(title="Select Peptide List (*XLS)",
                                              filetypes=(("Excel Files","*.xlsx"),
                                                         ("Excel Files","*.xls")))

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
    df_configuration = utils.assign_96_well_plate_physical_ids(df_configuration=df_configuration)
    return df_configuration.to_dict()

def reformat_plate_reader(plate_readout):
    spot_counts = plate_readout.to_numpy().flatten()
    pool_ids = ["pool_" + str(i) for i in np.arange(1, len(spot_counts)+1)]
    d = {'pool_id': pool_ids, 'spot_count': spot_counts}
    return pd.DataFrame(d)
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
    #df_hits = utils.assign_96_well_plate_physical_ids(df_configuration=df_hits)

    return df_hits.to_dict()


eel.start('index.html', size=(1920, 1080), port=utils.get_open_port())
