import pandas as pd
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign


def test_verify():
    excel_file = get_data_path(name='25peptides_5perpool_3x_configuration.xlsx')
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=excel_file
    )
    block_design = BlockDesign.read_excel_file(
        excel_file=excel_file
    )
    is_optimal = block_assignment.is_optimal(
        num_coverage=block_design.num_coverage,
        num_peptides_per_pool=block_design.num_peptides_per_pool
    )
    assert is_optimal, "'25peptides_5perpool_3x_configuration.xlsx' " \
                       "is an optimal ELISpot assignment."
