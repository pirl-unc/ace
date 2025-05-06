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
The purpose of this python3 script is to implement the PlateReadout dataclass.
"""


import os
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Type
from .block_assignment import BlockAssignment
from .logger import get_logger
from .plate_readout_entry import PlateReadoutEntry


logger = get_logger(__name__)


@dataclass(frozen=True)
class PlateReadout:
    entries: List[PlateReadoutEntry] = field(default_factory=list)
    pool_ids: Dict[str,int] = field(default_factory=dict)

    def add_entry(
            self,
            plate_id: int,
            well_id: str,
            spot_count: int
    ):
        """
        Adds a plate readout entry.

        Parameters:
            plate_id        :   Plate ID.
            well_id         :   Well ID.
            spot_count      :   Spot count.
        """
        entry = PlateReadoutEntry(plate_id=plate_id, well_id=well_id, spot_count=spot_count)
        self.entries.append(entry)

    def assign_pool_ids(self, block_assignment: BlockAssignment):
        df_assignment = block_assignment.to_dataframe()
        for entry in self.entries:
            df_matched = df_assignment[
                (df_assignment['plate_id'] == entry.plate_id) &
                (df_assignment['well_id'] == entry.well_id)
            ]
            if len(df_matched) == 0:
                logger.info('Could not find matching pool ID for %i-%s (plate ID - well ID)' % (entry.plate_id, entry.well_id))
            else:
                self.pool_ids['%i-%s' % (entry.plate_id, entry.well_id)] = df_matched['pool_id'].values[0]

    def to_dataframe(self) -> pd.DataFrame:
        assert len(self.pool_ids.keys()) > 0
        data = {
            'plate_id': [],
            'well_id': [],
            'pool_id': [],
            'spot_count': []
        }
        for entry in self.entries:
            if '%i-%s' % (entry.plate_id, entry.well_id) in self.pool_ids:
                data['plate_id'].append(entry.plate_id)
                data['well_id'].append(entry.well_id)
                data['pool_id'].append(self.pool_ids['%i-%s' % (entry.plate_id,entry.well_id)])
                data['spot_count'].append(entry.spot_count)
        return pd.DataFrame(data)

    @staticmethod
    def merge(plate_readouts: List['PlateReadout']) -> 'PlateReadout':
        """
        Merge PlateReadout objects into one PlateReadout object.

        Parameters:
            plate_readouts  :   List of PlateReadout objects.

        Returns:
            plate_readout   :   PlateReadout object.
        """
        merged_plate_readout = PlateReadout()
        for plate_readout in plate_readouts:
            for entry in plate_readout.entries:
                merged_plate_readout.add_entry(
                    plate_id=entry.plate_id,
                    well_id=entry.well_id,
                    spot_count=entry.spot_count
                )
        return merged_plate_readout

    @staticmethod
    def read_aid_plate_reader_file(
            excel_file: str,
            plate_id: int
    ) -> 'PlateReadout':
        """
        Load an AID plate reader Excel file.

        Parameters:
            excel_file      :   Excel file.
            plate_id        :   Plate ID.

        Returns:
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
                plate_readout.add_entry(
                    plate_id=int(value['plate_id']),
                    well_id=str(value['well_id']),
                    spot_count=int(value['spot_count'])
                )
            return plate_readout
        else:
            raise Exception("Could not locate the spot count data.")

    @staticmethod
    def read_pool_id_file(
            file: str
    ) -> 'PlateReadout':
        """
        Load a pool ID Excel file.

        Parameters:
            excel_file      :   Excel or CSV file.

        Returns:
            plate_readout   :   PlateReadout object.
        """
        file_name, file_extension = os.path.splitext(file)
        if file_extension == '.xlsx':
            df = pd.read_excel(file)
        elif file_extension == '.csv':
            df = pd.read_csv(file)
        else:
            raise Exception('Unexpected file extension: %s' % file_extension)
        plate_readout = PlateReadout()
        for _, row in df.iterrows():
            plate_readout.add_entry(
                plate_id=int(row['plate_id']),
                well_id=str(row['well_id']),
                spot_count=int(row['spot_count'])
            )

        return plate_readout
