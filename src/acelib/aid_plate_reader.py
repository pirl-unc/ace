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
The purpose of this python3 script is to implement the AidPlateReader dataclass.
"""


import pandas as pd
from dataclasses import dataclass, field
from .logger import get_logger


logger = get_logger(__name__)


@dataclass(frozen=True)
class AidPlateReader:

    @staticmethod
    def load_readout_file(excel_file, plate_id) -> pd.DataFrame:
        """
        Loads an AID plate reader Excel file.

        Parameters
        ----------
        excel_file  :   Excel file.
        plate_id    :   Plate ID.

        Returns
        -------
        df_counts   :   DataFrame with the following columns:
                        'plate_id'
                        'well_id'
                        'spot_count'
        """
        df = pd.read_excel(excel_file, sheet_name='PlateData', header=None)
        if df.iloc[1,0] == ' Plate Data - Spots Number ':
            df = df.loc[3:10,1:12]
            data = {
                'plate_id': [],
                'well_id': [],
                'spot_count': []
            }
            df.index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            for index, value in df.iterrows():
                for i in range(1, len(value) + 1):
                    data['plate_id'].append(plate_id)
                    data['well_id'].append('%s%i' % (index, i))
                    data['spot_count'].append(df.loc[index,i])

            df_counts = pd.DataFrame(data)
            return df_counts
        else:
            logger.error('Could not locate the spots data.')
            exit(1)