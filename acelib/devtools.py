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
serve as a developer toolkit to be ignored by the user
"""

import psutil
import random


def find_port_addr(connection):
    try:
        return connection.raddr.port
    except:
        return


def get_open_port():
    connections = psutil.net_connections()
    existing_ports = [find_port_addr(conx) for conx in connections]
    proposed_port = 1111
    while proposed_port in existing_ports:
        proposed_port = random.randint(1111, 9999)
    return proposed_port

