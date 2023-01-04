from enum import Enum

class EndType(Enum):
    paired_end = "paired-end"
    single_end = "single-end"
    single_cell = "single-cell"
