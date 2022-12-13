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


import random
import socket
from acelib.logger import get_logger


logger = get_logger(__name__)


def find_port_addr(connection):
    try:
        return connection.raddr.port
    except:
        return


def get_open_port() -> int:
    port = 1
    for i in range(10000):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.bind(("127.0.0.1", i))
        except socket.error as e:
             continue
        s.close()
        port = i
        break
    return port
