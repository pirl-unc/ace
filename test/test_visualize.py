from .data import get_data_path
from acelib.main import *


def test_visualize():
    tsv_file = get_data_path(name='elispot_configuration_sample.tsv')
    df_configuration = pd.read_csv(tsv_file, sep='\t')
    plt = run_ace_visualize(
        df_configuration=df_configuration
    )
