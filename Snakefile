"""
This is an attempt to convert HCC_Data/align_liver_PE.bash into a snakefile
"""

import glob
import os
import subprocess
import csv
configfile: "snakemake_config.yaml"

def get_tissue_name():
    """
    Looking to return the base filename from the controls/init_files
    Example:
    controls/init_files/naiveB_S1R1.csv
    controls/init_files/naiveB_S1R2.csv

    We would return: ["naiveB", "naiveB"]
    :return:
    """
    tissue_data = []
    with open(config["MASTER_CONTROL"], "r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            id = line[1].split("_")[0]  # naiveB_S1R1 -> naiveB
            tissue_data.append(id)

    return tissue_data

def get_tag_data():
    """
    Return tag from cell ID
    Example:
        input: naiveB_S1R1
        output: S1R1
    :return:
    """
    tag_data = []
    with open(config["MASTER_CONTROL"], "r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            tag = line[1].split("_")[-1]
            tag_data.append(tag)  # example: S1R1

    return tag_data

def get_srr_data():
    """
    Get the SRR information from the master init_file
    Example:
        input:
            SRR14231328,naiveB_S1R1,PE
            SRR14231329,naiveB_S1R2,PE
        output: ["SRR14231329", "SRR14231328"]

    :return:
    """
    srr_data = []
    with open(config["MASTER_CONTROL"], "r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            srr = line[0]
            srr_data.append(srr)
    return srr_data

def PE_SE_Data():
    """
    This function will read from the config[MASTER_CONTROL] file and return the paired_end or single_end variable
    Example:
        input:
            SRR14231328,naiveB_S1R1,PE
            SRR14231329,naiveB_S1R2,SE
        output:
            ["PE", "SE"]
    :return:
    """
    pe_se_data = []
    with open(config["MASTER_CONTROL"], "r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            pe_se = line[2]
            pe_se_data.append(pe_se)
    return pe_se_data

rule all:
    input:
        # download fastq
        # do not request distribute_init_files or prefetch_fastq output as these are marked as temp()
        expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw_temp", "Bulk_{PE_SE}"), tissue_name=get_tissue_name(), PE_SE=PE_SE_Data()),

        # Rename srr
        # expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","raw","{tissue_name}_{tag}.fastq.gz"), tissue_name=get_tissue_name(), tag=get_tag_data())

        # STAR output
        # config["STAR"]["GENOME_DIR"]

"""
STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir genome/star \
--genomeFastaFiles genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile genome/Homo_sapiens.GRCh38.104.gtf \
--sjdbOverhang 99
"""
rule generate_genome:
    input:
        genome_fasta_file = config["STAR"]["GENERATE_GENOME"]["GENOME_FASTA_FILE"],
        gtf_file = config["STAR"]["GENERATE_GENOME"]["GTF_FILE"]
    output:
        genome_dir = directory(config["STAR"]["GENERATE_GENOME"]["GENOME_DIR"]),
        rule_complete = touch(os.path.join("temp", "rule_complete", "generate_genome.complete"))
    threads: workflow.cores * 0.35
    params:
        run_mode = config["STAR"]["GENERATE_GENOME"]["RUN_MODE"],
        overhang = config["STAR"]["GENERATE_GENOME"]["OVERHANG"]
    shell:
        """
        module load star/2.7
        
        STAR \
        --runThreadN {threads} \
        --runMode {params.run_mode} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang {params.overhang}
        """

rule distribute_init_files:
    input: config["MASTER_CONTROL"]
    output: temp(os.path.join(config["ROOTDIR"], "controls", "init_files", "{tissue_name}_{tag}.csv"))
    params:
        id = "{tissue_name}_{tag}"
    run:
        # Get lines in master control file
        # Open output for writing
        lines = open(str(input), "r").readlines()
        wfile = open(str(output), "w")
        for line in lines:
            # Only write line if the output file has the current tissue-name_tag (naiveB_S1R1) in the file name
            if params.id in line:
                wfile.write(line)
        wfile.close()

checkpoint prefetch_fastq:
    input: rules.distribute_init_files.output
    output:
        data = temp(directory(os.path.join("temp", "prefetch", "{tissue_name}_{tag}/"))),
        rule_complete = touch(os.path.join("temp", "rule_complete", "prefetch_{tissue_name}_{tag}.complete"))
    params:
        id = "{tissue_name}_{tag}"
    shell:
        """
        module load SRAtoolkit
        IFS=","
        while read srr name endtype; do
            prefetch $srr --output-directory {output.data}
        done < {input}
        """

rule dump_fastq:
    input:
        prefetch_data = expand(rules.prefetch_fastq.output.data, tag=get_tag_data(), allow_missing=True),
        prefetch_data_complete = expand(rules.prefetch_fastq.output.rule_complete, tag=get_tag_data(), allow_missing=True)
    output:
        data = directory(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw_temp", "Bulk_{PE_SE}")),
        rule_complete = touch(os.path.join("temp", "rule_complete", "dump_fastq_{tissue_name}_{PE_SE}.complete"))
    threads: workflow.cores * 0.9  # max threads
    params:
        srr_ids = get_srr_data()
    shell:
        """
        module load parallel-fastq-dump
        parallel-fastq-dump \
        --split-files \
        --gzip \
        --sra-id {params.srr_ids} \
        --threads {threads} \
        --outdir {output}
        """

"""
dump_fastq output: raw_temp/[SRR14231328_1.fastq.gz, SRR14231328_2.fastq.gz]
"""
rule rename_srr:
    input: rules.dump_fastq.output.data  # wildcards: tissue_tag, PE_SE
    output: os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}.fastq.gz")
    params:
        tags = get_tag_data()

    shell: """
        for file in {input}/*; do
            mv "$file" "{output}"
        done
    """


"""
Trimming Plan
Have trim pull each file from rename_srr
Perform trimming on each file
"""
if config["PERFORM_TRIM"]:
    rule trim:
        input: expand(rules.dump_fastq.output.data, tissue_name=get_tissue_name(), PE_SE=PE_SE_Data())
        output: directory(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads"))
        shell: """
            module load gnu-parallel # DOI https://doi.org/10.5281/zenodo.1146014
            module load trim_galore  # trimGalore 10.5281/zenodo.5127898 
                                     # Cutadapt DOI:10.14806/ej.17.1.200  
            parallel --xapply trim_galore --illumina --paired --fastqc -o trim_galore/ ::: *_1.fastq.gz ::: *_2.fastq.gz
            
            """

def collect_star_align_input(wildcards):
    if config["PERFORM_TRIM"]:
        return rules.trim.output
    else:
        return rules.dump_fastq.output.data
rule star_align:
    input: collect_star_align_input
    output: directory(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads"))
    threads: workflow.cores * 0.90
    shell:
        """
        STAR --runThreadN {threads} \
		--readFilesCommand {config[STAR][ALIGN_READS][READ_COMMAND]} \
		--readFilesIn $file1 $file2 \
		--genomeDir {config[STAR][GENERATE_GENOME][GENOME_DIR]} \
		--outFileNamePrefix {output} \
		--outSAMtype {config[STAR][ALIGN_READS][OUT_SAM_TYPE]} \
		--outSAMunmapped {config[STAR][ALIGN_READS][OUT_SAM_UNMAPPED} \
		--outSAMattributes {config[STAR][ALIGN_READS][OUT_SAM_ATTRIBUTES]} \
		--quantMode {config[STAR][ALIGN_READS][QUANT_MODE]}
        """
