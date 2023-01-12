import os
import csv
import sys
import snakemake
from pathlib import Path
from typing import Literal
from . import perform
from .constants import EndType
from . import single_cell

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
            # These values will be used if we are collecting tissue/tag data
            if attribute in sub_attribute:
                sub_index = sub_attribute.index(attribute)
                tissue_and_tag: list[str] = str(line[index_value]).split("_")  # Get the tissue and tag value
            
            # If we are collecting single cell data, get it from NCBI database
            if PE_SE_value == EndType.single_cell.value:
                srr_code = line[0]
                srr_data: single_cell.SRR = single_cell.collect(srr_code)
            
            target_attribute: list[str] = []
            if PE_SE_value == EndType.paired_end.value:
                # Append everything twice
                if attribute in sub_attribute:
                    target_attribute.extend([tissue_and_tag[sub_index], tissue_and_tag[sub_index]])
                elif attribute == "PE_SE":
                    target_attribute.extend(["1", "2"])
                else:
                    target_attribute.extend([line[index_value], line[index_value]])
            
            elif PE_SE_value == EndType.single_end.value:
                # Append everything once
                if attribute in sub_attribute:
                    target_attribute.append(tissue_and_tag[sub_index])
                elif attribute == "PE_SE":
                    target_attribute.append("S")
                else:
                    target_attribute.append(line[index_value])
            
            elif PE_SE_value == EndType.single_cell.value:
                # Append based on number of reads collected
                if attribute in sub_attribute:
                    for _ in range(int(srr_data.num_reads)):
                        target_attribute.append(tissue_and_tag[sub_index])
                elif attribute == "PE_SE":
                    if srr_data.num_reads == "variable":
                        target_attribute.extend(["1", "2"])
                    else:
                        for i in range(1, int(srr_data.num_reads) + 1):
                            target_attribute.append(str(i))
                else:
                    for _ in range(int(srr_data.num_reads)):
                        target_attribute.append(line[index_value])

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
                end_type = line.split(",")[2]
                if end_type == "PE":
                    return EndType.paired_end
                elif end_type == "SE":
                    return EndType.single_end
                elif end_type == "SLC":
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
