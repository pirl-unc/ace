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


import eel
from acelib.main import run_ace_generate
from acelib import devtools


eel.init('')


@eel.expose
def generate_configuration(num_peptides,
                           num_peptides_per_pool,
                           num_coverage,
                           num_cores):
    num_peptides = int(num_peptides)
    num_peptides_per_pool = int(num_peptides_per_pool)
    num_coverage = int(num_coverage)
    num_cores = int(num_cores)
    disallowed_peps = [],
    df_configuration = run_ace_generate(
        n_peptides=num_peptides,
        n_peptides_per_pool=num_peptides_per_pool,
        n_coverage=num_coverage,
        num_processes=num_cores,
        disallowed_peptide_pairs=disallowed_peps
    )
    return df_configuration.to_dict()


eel.start('res/index.html', size=(1200, 1000), port=devtools.get_open_port())

