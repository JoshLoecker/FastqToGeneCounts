import os
import csv
import sys
import snakemake
from pathlib import Path
from typing import Literal
from . import perform
from .constants import EndType


def from_master_config(config: dict, attribute: Literal["SRR", "tissue", "tag", "PE_SE"]) -> list[str]:
    valid_attributes = ["SRR", "tissue", "tag", "PE_SE"]
    sub_attribute = ["tissue", "tag"]
    if attribute not in valid_attributes:
        sys.exit(f"\nInvalid attribute input. '{attribute}' is not one of: {valid_attributes}\n")
    else:
        collect_attributes = []
        index_value = valid_attributes.index(attribute)

        # We have to subtract one because "tissue" and "tag" are in the same index, thus the index value in valid_inputs is increased by one
        if index_value >= 2:
            index_value -= 1

        control_lines = open(config["MASTER_CONTROL"], "r").readlines()
        reader = csv.reader(control_lines)

        for line in reader:

            # Get the column from master_control we are interested in
            column_value = line[index_value]
            PE_SE_value = line[2]  # get PE or SE

            # test if we are looking for "tissue" or "tag", as these two values are located at master_control index 1
            if attribute in sub_attribute:
                sub_index = sub_attribute.index(attribute)
                tissue_and_tag: list[str] = str(line[index_value]).split("_")  # Get the tissue and tag value

                # We must append the target attribute twice if it is paired end, once if it is single end
                if PE_SE_value in [EndType.paired_end.value, EndType.single_cell.value]:
                    target_attribute = [tissue_and_tag[sub_index], tissue_and_tag[sub_index]]
                elif PE_SE_value == EndType.single_end.value:
                    target_attribute = [tissue_and_tag[sub_index]]

            elif attribute == "PE_SE":
                # We must append the target attribute twice if it is paired end, once if it is single end
                if column_value in [EndType.paired_end.value, EndType.single_cell.value]:
                    target_attribute = ["1", "2"]
                elif column_value == EndType.single_end.value:
                    target_attribute = ["S"]

            else:
                if PE_SE_value in [EndType.paired_end.value, EndType.single_cell.value]:
                    target_attribute = [line[index_value], line[index_value]]
                elif PE_SE_value == EndType.single_end.value:
                    target_attribute = [line[index_value]]

            collect_attributes += target_attribute

        return collect_attributes


def srr_code(config: dict) -> list[str]:
    """
    Only should be getting SRR values if we are performing prefetch
    """
    if perform.trim(config=config):
        return from_master_config(config=config, attribute="SRR")


def tissue_name(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        return from_master_config(config=config, attribute="tissue")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tissue_name


def tags(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        return from_master_config(config=config, attribute="tag")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tag


def PE_SE(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        return from_master_config(config=config, attribute="PE_SE")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.PE_SE

def end_type(config: dict, tissue_name: str, tag: str) -> EndType:
    """
    This function is responsible for determining if the input to star align is paired-end, single-end, or single-cell data
    """
    with open(config["MASTER_CONTROL"], "r") as i_stream:
        for line in i_stream:
            sample: str = f"{tissue_name}_{tag}"
            if sample in line:
                if str(EndType.paired_end.value) in line:
                    return EndType.paired_end
                elif str(EndType.single_end.value) in line:
                    return EndType.single_end
                elif str(EndType.single_cell.value) in line:
                    return EndType.single_cell
        else:
            raise ValueError(f"Tissue name of '{tissue_name}' and tag of '{tag}' could not be found in the control file. Please double check your control file")


def tag_from_filename(file_path: str | Path) -> str:
    file_name = os.path.basename(file_path)
    purge_extension = file_name.split(".")[0]
    tag = purge_extension.split("_")[-1]
    return str(tag)

def direction_from_name(file: str):
    file_name = os.path.basename(file)
    purge_extension = file_name.split(".")[0]
    direction = purge_extension.split("_")[-1]
    return direction

def sample(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        tag = from_master_config(config=config, attribute="tag")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        tag = fastq_input.tag

    sample = []
    for t in tag:
        sample.append(t.split("R")[0])
    return sample
