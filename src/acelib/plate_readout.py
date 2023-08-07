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
The purpose of this python3 script is to implement the SpotCount dataclass.
"""


import pandas as pd
from dataclasses import dataclass, field
from typing import List, Type
from .logger import get_logger
from .types import *


logger = get_logger(__name__)


@dataclass(frozen=True)
class PlateReadout:
    observations: PlateReadoutObservations = field(default_factory=list)

    def add_observation(
            self,
            plate_id: PlateId,
            well_id: WellId,
            spot_count: SpotCount
    ):
        """
        Adds a readout observation.

        Parameters
        ----------
        plate_id        :   Plate ID.
        well_id         :   Well ID.
        spot_count      :   Spot count.
        """
        self.observations.append((plate_id, well_id, spot_count))

    def to_dataframe(self) -> pd.DataFrame:
        data = {
            'plate_id': [],
            'well_id': [],
            'spot_count': []
        }
        for observation in self.observations:
            data['plate_id'].append(observation[0])
            data['well_id'].append(observation[1])
            data['spot_count'].append(observation[2])
        return pd.DataFrame(data)

    @staticmethod
    def merge(plate_readouts: List[Type['PlateReadout']]) -> Type['PlateReadout']:
        """
        Merges

        Parameters
        ----------
        plate_readouts  :   List of PlateReadout objects.

        Returns
        -------
        plate_readout   :   PlateReadout object.
        """
        plate_readout = PlateReadout()
        for plate_readout_ in plate_readouts:
            for readout in plate_readout_.observations:
                plate_readout.add_observation(
                    plate_id=readout[0],
                    well_id=readout[1],
                    spot_count=readout[2]
                )
        return plate_readout

    @staticmethod
    def read_aid_plate_reader_file(
            excel_file: str,
            plate_id: int
    ) -> 'PlateReadout':
        """
        Loads an AID plate reader Excel file.

        Parameters
        ----------
        excel_file      :   Excel file.
        plate_id        :   Plate ID.

        Returns
        -------
        plate_readout   :   PlateReadout object.
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

            # Convert to a PlateReadout object
            plate_readout = PlateReadout()
            for index, value in df_counts.iterrows():
                plate_readout.add_observation(
                    plate_id=value['plate_id'],
                    well_id=value['well_id'],
                    spot_count=value['spot_count']
                )
            return plate_readout
        else:
            logger.error('Could not locate the spots data.')
            exit(1)

