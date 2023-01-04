from enum import Enum

class EndType(Enum):
    paired_end = "PE"
    single_end = "SE"
    single_cell = "SLC"
