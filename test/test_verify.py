import pandas as pd
from .data import get_data_path
from acelib.main import run_ace_verify


def test_verify():
    configuration_file = get_data_path(name='25peptides_5perpool_3x.csv')
    df_configuration = pd.read_csv(configuration_file)
    is_valid = run_ace_verify(
        df_configuration=df_configuration,
        num_coverage=3
    )

    assert is_valid, \
        "'25peptides_5perpool_3x.csv' is a valid ELIspot configuration."
