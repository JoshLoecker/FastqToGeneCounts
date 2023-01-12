import os
import csv
import warnings
import sys
import pandas as pd
from typing import Literal
from utils import get, perform
from utils.constants import EndType
from utils import single_cell

configfile: "snakemake_config.yaml"
# samples: pd.DataFrame = pd.read_csv(config["MASTER_CONTROL"])

# Ensure the results directory is made
os.makedirs(config["ROOTDIR"], exist_ok=True)


# Validate users are using conda. This is important for temporary conda environments defined in the workflow
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass the '--use-conda' flag to snakemake.\nExample: snakemake --cores 10 --use-conda\n\n")
    sys.exit(1)

def fastqc_trimmed_reads(wildcards):
    """
    If we are going to trim, return output for rule fastqc_trim
    """
    if perform.trim(config=config):
        return expand(
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        return []


def perform_screen_rule(wildcards):
    """
    If screening for contamination, return fastq_screen output
    """
    if perform.screen(config=config):
        return snakemake.io.expand(
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fq_screen", "{tissue_name}_{tag}_{PE_SE}_screen.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        return []


def perform_get_insert_size_rule(wildcards):
    """
    If getting insert sizes with picard, return GetinsertSizeMetrics output
    """
    if perform.get_insert_size(config=config):
        return snakemake.io.expand(
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "picard", "insert", "{tissue_name}_{tag}_insert_size.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        )
    else:
        return []


def perform_get_fragment_size_rule(wildcards):
    """
    If getting fragment sizes with deeptools, return RNA_fragment_size.py output
    """
    if perform.get_fragment_size(config=config):
        return snakemake.io.expand(
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fragmentSizes", "{tissue_name}_{tag}_fragment_length.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        )
    else:
        return []


def perform_trim_rule(wildcards):
    """
    If we are performing trimming, return trim's output
    """
    if perform.trim(config=config):
        return snakemake.io.expand(
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        return []


def perform_dump_fastq(wildcards):
    if perform.prefetch(config=config):
        return expand(
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        return []


rule_all = [
    # Preroundup
    expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "layouts", "{tissue_name}_{tag}_layout.txt"), tissue_name=get.tissue_name(config=config), tag=get.tags(config=config)),
    expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "prepMethods", "{tissue_name}_{tag}_prep_method.txt"), tissue_name=get.tissue_name(config=config), tag=get.tags(config=config)),

    # Screen Genome data
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Adapters"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Arabidopsis"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Drosophila"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "E_coli"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Human"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Lambda"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Mitochondria"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Mouse"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "PhiX"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Rat"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Vectors"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Worm"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Yeast"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "rRNA"),
    os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "fastq_screen.conf"),

    # Generate Genome
    config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],
    perform_dump_fastq,  # dump_fastq
    perform_screen_rule,  # fastq_screen
    perform_trim_rule,  # trim reads

    # FastQC
    expand(
        os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        PE_SE=get.PE_SE(config=config)
    ),

    fastqc_trimmed_reads,

    # STAR aligner
    expand(
        os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.tab"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    ),

    # FastQ aligned reads
    expand(
        os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.bam.bai"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    ),

    # copy .tab
    expand(
        os.path.join("MADRID_input", "{tissue_name}", "geneCounts", "{sample}", "{tissue_name}_{tag}.tab"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        sample=get.sample(config=config)
    ),

    # get rnaseq metrics
    expand(
        os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "picard", "rnaseq", "{tissue_name}_{tag}_rnaseq.txt"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    ),

    # copy strandedness
    expand(
        os.path.join("MADRID_input", "{tissue_name}", "strandedness", "{sample}", "{tissue_name}_{tag}_strandedness.txt"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        sample=get.sample(config=config)
    ),

    # MultiQC
    expand(
        os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "multiqc", "{tissue_name}_multiqc_report.html"),
        tissue_name=get.tissue_name(config=config)
    ),
]

if perform.get_insert_size(config=config):
    rule_all.extend([
        # Get Insert sizes
        perform_get_insert_size_rule,

        # copy insert
        expand(os.path.join("MADRID_input", "{tissue_name}", "insertSizeMetrics", "{sample}", "{tissue_name}_{tag}_insert_size.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            sample=get.sample(config=config)
        )
    ])

if perform.get_fragment_size(config=config):
    rule_all.extend([
        # get fragment lengths
        perform_get_fragment_size_rule,

        # copy fragment
        expand(os.path.join("MADRID_input", "{tissue_name}", "fragmentSizes", "{sample}", "{tissue_name}_{tag}_fragment_size.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            sample=get.sample(config=config)
        )
    ])

rule all:
    input: rule_all

rule preroundup:
    input: config["MASTER_CONTROL"]
    output:
        layout=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "layouts", "{tissue_name}_{tag}_layout.txt"),
        preparation=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "prepMethods", "{tissue_name}_{tag}_prep_method.txt"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 200 * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt
    run:
        # SRR12873784,effectorcd8_S1R1,PE,total
        with open(str(input),"r") as i_stream:
            reader = csv.reader(i_stream)
            for line in reader:
                # Collect the required data
                srr_code: str = line[0]                 # SRR123
                name: str = line[1]                     # naiveB_S1R1
                tissue_name: str = name.split("_")[0]   # naiveB
                tag: str = name.split("_")[1]           # S1R1
                end_type: str = line[2]                  # PE
                prep: str = line[3].lower()             # total

                # Set the study
                # If the tag is S1R1, extract the "S1" component
                study: str = ""
                for char in tag:
                    if char == "R":  # Once we reach the replicate, we can exit the loop
                        break
                    elif char == "S":
                        study += char
                    elif char.isdigit():
                        study += char

                # Make the required directories
                os.makedirs(os.path.join("MADRID_input", tissue_name, "geneCounts"), exist_ok=True)
                os.makedirs(os.path.join("MADRID_input", tissue_name, "insertSizeMetrics"), exist_ok=True)
                os.makedirs(os.path.join("MADRID_input", tissue_name, "layouts"), exist_ok=True)
                os.makedirs(os.path.join("MADRID_input", tissue_name, "layouts", study), exist_ok=True)
                os.makedirs(os.path.join("MADRID_input", tissue_name, "fragmentSizes"), exist_ok=True)
                os.makedirs(os.path.join("MADRID_input", tissue_name, "prepMethods"), exist_ok=True)
                os.makedirs(os.path.join("MADRID_input", tissue_name, "prepMethods", study), exist_ok=True)
                os.makedirs(os.path.join(config["ROOTDIR"], "data", tissue_name, "layouts"), exist_ok=True)
                os.makedirs(os.path.join(config["ROOTDIR"], "data", tissue_name, "prepMethods"), exist_ok=True)

                # Write single/paired end to the appropriate location
                end_type_write_root = open(os.path.join(config["ROOTDIR"], "data", tissue_name, "layouts", f"{name}_layout.txt"), "w")
                end_type_write_madrid = open(os.path.join("MADRID_input",tissue_name,"layouts",study,f"{name}_layout.txt"), "w")
                match end_type:
                    case EndType.single_end.value:
                        end_type_write_root.write("single-end")
                        end_type_write_madrid.write("single-end")
                    case EndType.paired_end.value:
                        end_type_write_root.write("paired-end")
                        end_type_write_madrid.write("paired-end")
                    case EndType.single_cell.value:
                        end_type_write_root.write("single-cell")
                        end_type_write_madrid.write("single-cell")
                    case _:
                        print(f"Rule preroundup: Invalid endtype. Must be one of: {[i.name for i in EndType.__members__.values()]}")
                        raise ValueError(f"Rule preroundup: Invalid endtype: {end_type}")
                end_type_write_root.close()
                end_type_write_madrid.close()

                # Write mrna/total to the appropriate location
                prep_method_root = open(os.path.join(config["ROOTDIR"],"data",tissue_name,"prepMethods",f"{name}_prep_method.txt"), "w")
                prep_method_madrid = open(os.path.join("MADRID_input",tissue_name,"prepMethods",study,f"{name}_prep_method.txt"), "w")
                if prep == "mrna":
                    prep_method_root.write("mrna")
                    prep_method_madrid.write("mrna")
                elif prep == "total":
                    prep_method_root.write("total")
                    prep_method_madrid.write("total")
                else:
                    raise ValueError(f"Rule preroundup: Invalid library preparation: {prep}. Must be `total` or `mrna`")
                prep_method_root.close()
                prep_method_madrid.close()


rule generate_genome:
    input:
        genome_fasta_file=config["GENERATE_GENOME"]["GENOME_FASTA_FILE"],
        gtf_file=config["GENERATE_GENOME"]["GTF_FILE"]
    output:
        genome_dir=directory(os.path.join(config["GENERATE_GENOME"]["GENOME_SAVE_DIR"])),
        rule_complete=touch(os.path.join(config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],"generate_genome.complete"))
    threads: 40
    resources:
        mem_mb=50000,# 50 GB
        runtime=lambda wildcards, attempt: 120 * attempt
    conda: "envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang 99
        """

if perform.screen(config=config):
    rule get_screen_genomes:
        """
        Download genomes to screen against
        """
        output:
            genomes=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes")),
            Adapters=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Adapters")),
            Arabidopsis=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Arabidopsis")),
            Drosophila=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Drosophila")),
            E_coli=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "E_coli")),
            Human=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Human")),
            Lambda=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Lambda")),
            Mitochondria=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Mitochondria")),
            Mouse=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Mouse")),
            PhiX=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "PhiX")),
            Rat=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Rat")),
            Vectors=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Vectors")),
            Worm=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Worm")),
            Yeast=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "Yeast")),
            rRNA=directory(os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "rRNA")),
            config=os.path.join(config["ROOTDIR"], "FastQ_Screen_Genomes", "fastq_screen.conf")
        threads: 15
        params:
            output_dir = config["ROOTDIR"],
            download_paths = [
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Adapters/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Arabidopsis/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Drosophila/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/E_coli/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Human/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Lambda/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Mitochondria/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Mouse/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/PhiX/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Rat/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Vectors/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Worm/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Yeast/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/rRNA/",
                "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/fastq_screen.conf"
            ]
        resources:
            mem_mb=lambda wildcards, attempt: 1500 * attempt, # 1.5 GB * attempt
            runtime=lambda wildcards, attempt: 240 * attempt  # 240 minutes * attempt (4 hours)
        shell:
            """
            for path in {params.download_paths}; do
                # Get the species from the URL (i.e., E_coli, Human, Rat, etc.)
                species=$(echo $path | cut -d "/" -f 6)
                
                # Only download directories that are not currently existing
                # Use `--reject` to ignore any index files
                # If it does exist, touch the file so snakeamke doesn't complain about output files not existing
                if [[ ! -d "{output.genomes}/$species" ]]; then
                    # Test if working on fastq_screen.conf, we don't want to make this a directory because it is a file
                    if [[ ! $path =~ .*fastq_screen\.conf.* ]]; then
                        mkdir -p /scratch/$species
                    fi
                    wget --quiet --recursive --no-parent --no-host-directories --cut-dirs=2 --reject="index.html*" -P /scratch $path && echo "Finished $species" &
                else
                    find "{output.genomes}/$species" -t file -exec touch {{}} \; &
                fi
            done
            
            # Wait for all downloads to be done
            wait
            
            # Move scratch downloads into the results directory
            mv /scratch/* {output.genomes}
            
            # Replace "[FastQ_Screen_Genomes_Path]" with the sed_dir
            sed -i 's.\[FastQ_Screen_Genomes_Path\].{output.genomes}.g' {output.config}
            """


if perform.prefetch(config=config):
    rule distribute_init_files:
        input: ancient(config["MASTER_CONTROL"])  # Always run this rule to update its output
        output: report(os.path.join(config["ROOTDIR"],"temp","init_files","{tissue_name}", "{tissue_name}_{tag}.csv"))
        params: id="{tissue_name}_{tag}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 300 * attempt,
            runtime=lambda wildcards, attempt: 5 * attempt
        benchmark: repeat(os.path.join( "benchmarks", "{tissue_name}", "distribute_init_files", "{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
        run:
            # Get lines in master control file
            # Open output for writing
            lines = open(str(input),"r").readlines()
            for line in lines:

                # Only write line if the output file has the current tissue-name_tag (naiveB_S1R1) in the file name
                if params.id in line:
                    with open(str(output), "w") as o_stream:
                        o_stream.write(line)


    rule prefetch:
        input: rules.distribute_init_files.output
        output: os.path.join(config["ROOTDIR"],"temp","prefetch","{tissue_name}","{tissue_name}_{tag}","{tissue_name}_{tag}.sra")
        conda: "envs/SRAtools.yaml"
        threads: 1
        params:
            temp_directory="/scratch",
            temp_file="/scratch/{tissue_name}_{tag}.sra",
            output_directory=os.path.join(config["ROOTDIR"],"temp","prefetch","{tissue_name}_{tag}"),
        resources:
            mem_mb=lambda wildcards, attempt: 10000 * attempt,
            runtime=lambda wildcards, attempt: 30 * attempt,
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","prefetch","{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            # If the SRA file lock exists, remove it
            lock_file={output}.lock
            if [ -f "$lock_file" ]; then
                rm $lock_file
            fi
                    
            IFS=", "
            curr_dir=$(pwd)
            while read srr name endtype prep; do
                # Change into the "scratch" directory so temp files do not populate in the working directory
                cd {params.temp_directory}
                
                # set unlimited max size for prefetch
                prefetch --max-size u --progress --resume yes --output-file {params.temp_file} $srr
            done < {input}
            
            # Change back to the working directory before moving files
            cd $curr_dir
            mv {params.temp_file} {output}
            
            # Move dependencies into the output directory, checking if files exist in "/scratch"
            if [ -n "$(find {params.temp_directory} -prune -empty)" ]; then
                mv {params.temp_directory}/* {params.output_directory}
            fi
                
            
            """

    checkpoint fasterq_dump:
        input:
            prefetch=rules.prefetch.output,
            srr=rules.distribute_init_files.output
        output: fastq=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        threads: 40
        conda: "envs/SRAtools.yaml"
        params:
            # split_command=lambda wildcards: "--split-files" if wildcards.PE_SE in ["1", "2"] else "--concatenate-reads",
            temp_dir="/scratch",
            temp_filename=lambda wildcards: f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq" if wildcards.PE_SE in ["1", "2"]
                                            else f"{wildcards.tissue_name}_{wildcards.tag}.fastq",
            gzip_file=lambda wildcards: f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq.gz" if wildcards.PE_SE in ["1", "2"]
                                        else f"{wildcards.tissue_name}_{wildcards.tag}.fastq.gz",
            split_files=lambda wildcards: True if wildcards.PE_SE in ["1", "2"] else False,
            end_type= lambda wildcards: get.end_type(config=config,tissue_name=wildcards.tissue_name,tag=wildcards.tag)
        resources:
            mem_mb=lambda wildcards, attempt: 25600 * attempt,  # 25 GB
            time_min=lambda wildcards, attempt: 45 * attempt,
            disk_mb=lambda wildcards, attempt: 256000 * attempt,  # 250 GB
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","fasterq_dump","{tissue_name}_{tag}_{PE_SE}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            command='fasterq-dump --force --progress --threads {threads} --temp {params.temp_dir} --outdir {params.temp_dir}'
            
            # Set the split/concatenate based on paired end or single end data
            if [[ "{params.split_files}" == "True" ]]; then
                command+=' --split-files'
            else
                command+=' --concatenate-reads'
            fi
            
            if [[ "{params.end_type}" == "single-cell" ]]; then
                command+=" --include-technical"
            fi
            
            # Add the SRA file path to the command
            command+=' {input.prefetch}'
            
            echo $command
            eval $command
        
            # gzip the output
            echo "scratch"
            ls /scratch/
            pigz --synchronous --processes {threads} {params.temp_dir}/{params.temp_filename}
        
            mv {params.temp_dir}/{params.gzip_file} {output}
            """


    def dump_fastq_input(wildcards):
        output_files = expand(
            rules.prefetch.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            srr_code=get.srr_code(config=config)
        )
        for file in output_files:
            if wildcards.tissue_name in file and wildcards.tag in file:
                return file

    def get_dump_fastq_srr_code(wildcards, input):
        """Get SRR codes corresponding to dump_fastq output"""
        file_name = os.path.basename(str(input))
        srr_code = file_name.split(".")[0]
        return srr_code


def fastqc_dump_fastq_input(wildcards):
    """
    This function will return the input for fastqc_dump_fastq
    It is going to return forward read AND reverse read if the input file is the forward read
    If input is the reverse read, it will only return the reverse read
    If input is a single end read, it will only return the single end read
    """
    if perform.prefetch(config=config):
        if str(wildcards.PE_SE) == "1":
            forward_read: str = checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="1").output[0]
            reverse_read: str = checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="2").output[0]
            return [forward_read, reverse_read]
        else:
            return checkpoints.fasterq_dump.get(**wildcards).output
    else:
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                if (wildcards.tissue_name in file) and (wildcards.tag in file) and (f"_{wildcards.PE_SE}" in file):
                    file_one = os.path.join(path,file)
                    if str(wildcards.PE_SE) == "1":
                        forward_read: str = file_one
                        reverse_read: str = forward_read.replace("_1.fastq.gz","_2.fastq.gz")
                        return [forward_read, reverse_read]
                    else:
                        return file_one

rule fastqc_dump_fastq:
    input:
        fastq = fastqc_dump_fastq_input,
    output: os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    params:
        file_one_zip=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        file_one_html=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "{tissue_name}_{tag}_{PE_SE}_fastqc.html"),
        file_two_zip=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "{tissue_name}_{tag}_2_fastqc.zip"),
        file_two_html=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "{tissue_name}_{tag}_2_fastqc.html"),

        file_one_zip_rename=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        file_one_html_rename=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.html"),
        file_two_zip_rename=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_2_fastqc.zip"),
        file_two_html_rename=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_2_fastqc.html")
    threads: 8
    conda: "envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 15000 * attempt, # 15 GB * attempt number
        runtime=lambda wildcards, attempt: 150 * attempt  # 150 minutes * attempt number
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","fastqc_dump_fastq","{tissue_name}_{tag}_{PE_SE}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
        output_directory="$(dirname {output})"
        mkdir -p "$output_directory"

        if [ "{wildcards.PE_SE}" == "1" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

            mv "{params.file_one_zip}" "{params.file_one_zip_rename}"
            mv "{params.file_one_html}" "{params.file_one_html_rename}"
            mv "{params.file_two_zip}" "{params.file_two_zip_rename}"
            mv "{params.file_two_html}" "{params.file_two_html_rename}"

        # Touch the output file because Snakemake will complain about missing files otherwise
        # This file will be created when PE_SE == "1"
        elif [ "{wildcards.PE_SE}" == "2" ]; then
            touch {output}

        elif [ "{wildcards.PE_SE}" == "S" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

            mv "{params.file_one_zip}" "{params.file_one_zip_rename}"
            mv "{params.file_one_html}" "{params.file_one_html_rename}"
        fi
        """


if perform.screen(config=config):
    def get_screen_input(wildcards):
        """
        aggregate filesnames of all fastqs
        """
        if str(wildcards.PE_SE) == "1":
            forward_read: str = str(checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="1").output)
            reverse_read: str = str(checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="2").output)
            # reverse_read: str = forward_read.replace("_1.fastq.gz", "_2.fastq.gz")
            return [forward_read, reverse_read]
        else:
            return checkpoints.fasterq_dump.get(**wildcards).output


    rule contaminant_screen:
        input:
            files=get_screen_input,
            genomes=rules.get_screen_genomes.output.genomes,
        output: os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fq_screen", "{tissue_name}_{tag}_{PE_SE}_screen.txt")
        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            PE_SE="{PE_SE}",
            genomes_config=rules.get_screen_genomes.output.config,
            output_directory=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fq_screen")
        conda: "envs/screen.yaml"
        threads: lambda wildcards: 20 if str(wildcards.PE_SE) in ["1", "S"] else 1
        resources:
            # TODO: Limit ram if on reverse read
            mem_mb=lambda wildcards, attempt: 20000 * attempt if str(wildcards.PE_SE) in ["1", "S"] else 200, # 20 GB
            runtime=lambda wildcards, attempt: 30 * attempt
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","contaminant_screen","{tissue_name}_{tag}_{PE_SE}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            # Run fastq screen if PE_SE is 1 or S
            if [[ "{params.PE_SE}" == "1" ]]; then
                fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {params.genomes_config} --outdir {params.output_directory} {input.files}
                    # Run on single strand

            elif [[ "{params.PE_SE}" == "S" ]]; then
                fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {params.genomes_config} --outdir {params.output_directory} {input.files}

            # Only touch reverse read, will be created by forward read
            elif [[ "{params.PE_SE}" == "2" ]]; then
                touch {output}
            fi
            """

if perform.trim(config=config):
    def get_trim_input(wildcards):
        output_files = checkpoints.fasterq_dump.get(**wildcards).output
        if str(wildcards.PE_SE) in ["1", "2"]:
            forward_read: str = str(checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="1").output)
            reverse_read: str = str(checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="2").output)
            return [forward_read, reverse_read]
        else:
            return output_files.fastq

    checkpoint trim:
        input: get_trim_input,
        output: os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        params:
            temp_dir="/scratch",
        # Trim galore call uses 4 threads for forward/single reads. Request more because Trim can use UP TO this many
        threads: 16
        conda: "envs/trim.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: 10000 * attempt,# 10 GB
            runtime=lambda wildcards, attempt: 120 * attempt
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","trim","{tissue_name}_{tag}_{PE_SE}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """            
            output_directory="$(dirname {output})"

            if [[ "{wildcards.PE_SE}" == "1" ]]; then
                file_out_1="{params.temp_dir}/{wildcards.tissue_name}_{wildcards.tag}_1_val_1.fq.gz"    # final output paired end, forward read
                trim_galore --paired --cores 4 -o {params.temp_dir} {input}
                mv "$file_out_1" "{output}"

            # Skip over reverse-reads. Create the output file so snakemake does not complain about the rule not generating output
            elif [[ "{wildcards.PE_SE}" == "2" ]]; then
                file_out_2="{params.temp_dir}/{wildcards.tissue_name}_{wildcards.tag}_2_val_2.fq.gz"    # final output paired end, reverse read
                trim_galore --paired --cores 4 -o {params.temp_dir} {input}
                mv "$file_out_2" "{output}"

            # Work on single-end reads
            elif [[ "{wildcards.PE_SE}" == "S" ]]; then
                file_out="$output_directory/{wildcards.tissue_name}_{wildcards.tag}_S_trimmed.fq.gz"   # final output single end
                trim_galore --cores 4 -o "$output_directory" {input}
                mv "$file_out" "{output}"
            fi
            """

    def get_fastqc_trim_input(wildcards):
        if wildcards.PE_SE == "1":
            forward: str = str(checkpoints.trim.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="1").output)
            reverse: str = str(checkpoints.trim.get(tissue_name=wildcards.tissue_name, tag=wildcards.tag, PE_SE="2").output)
            return [forward, reverse]
        else:
            return checkpoints.trim.get(**wildcards).output
    rule fastqc_trim:
        input: get_fastqc_trim_input  # Original: rules.trim.output
        output: os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "fastqc", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
        params:
            file_two_input=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_2.fastq.gz"),
            file_two_out=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "trimmed_reads", "trimmed_{tissue_name}_{tag}_2_fastqc.zip")
        threads: 8
        conda: "envs/fastqc.yaml"
        resources:
            # Allocate 250MB per thread, plus extra to be safe
            # threads * 250 * 2 ~= 500 to 1000 GB
            mem_mb=lambda wildcards, attempt, threads: attempt * threads * 1000,
            runtime=lambda wildcards, attempt: 150 * attempt  # 2.5 hours * attempt
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","fastqc_trim","{tissue_name}_{tag}_{PE_SE}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            output_directory="$(dirname {output})"
            mkdir -p "$output_directory"

            if [ "{wildcards.PE_SE}" == "1" ]; then
                # send fastqc commands to background so we can run both at the same time 
                fastqc {input} --threads {threads} -o "$output_directory"

            # Skip reverse reads, but create the output file so Snakemake does not complain about missing files
            # This file will be created when wildcards.PE_SE == "1"
            elif [ "{wildcards.PE_SE}" == "2" ]; then
                touch {output}

            elif [ "{wildcards.PE_SE}" == "S" ]; then
                fastqc {input} --threads {threads} -o "$output_directory"
            fi
            """


def collect_star_align_input(wildcards):
    if perform.trim(config=config):
        # Have not expanded output from rule trim, need to expand it here
        in_files = sorted(
            expand(
                rules.trim.output,
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config)
            )
        )
    else:
        # already expanding output from dump_fastq, no need to expand it here
        in_files = sorted(
            expand(
                rules.fasterq_dump.output,
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config)
            )
        )

    grouped_reads = []
    for i, in_file in enumerate(in_files):
        direction = get.direction_from_name(in_file)
        try:
            next_file = in_files[i + 1]
            next_direction = get.direction_from_name(next_file)
        except:
            if direction == "S":
                grouped_reads.append(in_file)
                continue
            elif direction == "2":
                continue
            else:
                warnings.warn(f"{in_file} expects additional paired-end read! Skipping....")
                continue

        if direction == "S":
            grouped_reads.append(in_file)
        elif direction == "1" and next_direction == "2":
            if in_file[:-10] == next_file[:-10]:  # remove _1.fastq.gz to make sure they are same replicate
                both_reads = " ".join([in_file, next_file])
                grouped_reads.append(both_reads)
            else:
                warnings.warn(f"{in_file} and {next_file} are incorrectly called together, either the file order is getting scrambled or one end of {in_file} and one end of {next_file} failed to download")

        elif direction == "1" and not next_direction == "2":
            warnings.warn(f"{in_file} expects additional paired-end read! Skipping....")
        elif direction == "2":
            continue
        else:
            warnings.warn(f"{in_file} not handled, unknown reason!")

    """
    We need to return a string, or list of strings. If we return "grouped_reads" directly, some values within are not actually valid files, such as:
        ["results/data/naiveB/naiveB_S1R1_1.fastq.gz results/data/naiveB/naiveB_S1R1_2.fastq.gz", "results/data/naiveB/naiveB_S1R2_S.fastq.gz"]
    Index 0 is taken literally, as a string to a file location. Thus, it does not exist
    Because of this, we are going to filter through each input file and return it if it matches our desired tissue_name and tag
    This is much like what was done in the function get_dump_fastq_output, located above rule all
    """
    for read in grouped_reads:
        if wildcards.tissue_name in read and wildcards.tag in read:
            return read.split(" ")



def new_star_input(wildcards):
    # Open the control file to determine which samples are paired end or not
    is_paired_end: bool = False
    is_single_cell: bool = False
    is_single_end: bool = False

    end_type: EndType = get.end_type(config=config, tissue_name=wildcards.tissue_name, tag=wildcards.tag)
    match end_type:
        case EndType.paired_end:
            is_paired_end = True
        case EndType.single_end:
            is_single_cell = True
        case EndType.single_cell:
            is_single_cell = True

    # Get the output files, using our determined paired/single ends
    if is_paired_end:
        forward = checkpoints.trim.get(**wildcards, PE_SE="1").output
        reverse = checkpoints.trim.get(**wildcards, PE_SE="2").output
        returnal = forward + reverse
    elif is_single_end:
        returnal = checkpoints.trim.get(**wildcards, PE_SE="S").output
    elif is_single_cell:
        srr_code: list[str] = get.srr_code(config)
        srr_data: single_cell.SRR = single_cell.collect(srr_code[0])

        forward = checkpoints.trim.get(**wildcards, PE_SE="1").output
        reverse = checkpoints.trim.get(**wildcards, PE_SE="2").output
        returnal = forward + reverse
        if srr_data.num_reads is int and int(srr_data.num_reads) == 3:
            returnal += checkpoints.trim.get(**wildcards, PE_SE="3").output

    return returnal

def get_star_align_runtime(wildcards, input, attempt):
    """
    This function will return the length of time required for star_align to complete X number of reads
    Using 40 threads, it takes ~9 minutes per input file
    Round this value to 30 minutes (to be on the very safe side)
    We are also going to multiply by the attempt that the workflow is on.
    If on the second/third/etc. attempt, double/triple/etc. time is requested
    Return an integer of: len(input) * 30 minutes * attempt number = total runtime
    """
    # Convert input to a list, it is sent as a string
    input_list = str(input.reads).split(" ")

    # Max time is 7 days (10,080 minutes). Do not let this function return more than this time
    return min(len(input_list) * 20 * attempt, 10079)


rule star_align:
    input:
        # reads=collect_star_align_input,
        reads=new_star_input,
        genome_dir=rules.generate_genome.output.genome_dir,
        generate_genome_complete=rules.generate_genome.output.rule_complete
    output:
        gene_table=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.tab"),
        bam_file=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.bam")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        end_type=lambda wildcards: get.end_type(config=config, tissue_name=wildcards.tissue_name, tag=wildcards.tag),
        gene_table_output=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}_ReadsPerGene.out.tab"),
        bam_output=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}_Aligned.sortedByCoord.out.bam"),
        prefix=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}_"),
    threads: 40
    conda: "envs/star.yaml"
    resources:
        mem_mb=51200,  # 50 GB
        runtime=get_star_align_runtime
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","star_align","{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
        
        if [[ "{params.end_type.name}" == "single-cell" ]]; then
            command="STAR --outFilterType BySJout --outFilterMultimapNmax 20"
        else
            command="STAR"
        fi
        
        $command \
        --runThreadN {threads} \
		--readFilesCommand "zcat" \
		--readFilesIn {input.reads} \
		--genomeDir {input.genome_dir} \
		--outFileNamePrefix {params.prefix} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard \
		--quantMode GeneCounts
		
		mv {params.gene_table_output} {output.gene_table}
		mv {params.bam_output} {output.bam_file}
        """

rule copy_geneCounts:
    input: rules.star_align.output.gene_table
    output: os.path.join("MADRID_input", "{tissue_name}", "geneCounts", "{sample}", "{tissue_name}_{tag}.tab")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        sample=os.path.join("MADRID_input", "{tissue_name}", "geneCounts", "{sample}")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,# 0.5 GB
        runtime=5
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","copy_geneCounts","{sample}_{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
            mkdir -p {params.sample}
            cp {input} {output}
        """


rule index_bam_file:
    input: rules.star_align.output.bam_file
    output: os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.bam.bai")
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt, # 2 GB per attempt
        runtime=lambda wildcards, attempt: 5 * attempt  # 5 minutes per attempt
    conda: "envs/samtools.yaml"
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","index_bam_file","{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
        samtools index -@ {threads} {input} {output}
        """


rule get_rnaseq_metrics:
    input:
        bam=rules.star_align.output.bam_file,
        tab=rules.star_align.output.gene_table
    output:
        metrics=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "picard", "rnaseq", "{tissue_name}_{tag}_rnaseq.txt"),
        strand=os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "strand", "{tissue_name}_{tag}_strand.txt")
    params:
        ref_flat=config["REF_FLAT_FILE"],
        ribo_int_list=config["RRNA_INTERVAL_LIST"]
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 2500 * 5 * attempt,# 5 GB / attempt
        runtime=lambda wildcards, attempt: 60 * attempt
    conda: "envs/picard.yaml"
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","get_rnaseq_metrics","{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
        # Create the parent output directories
        #mkdir -p $(dirname -- "{output.metrics}")
        #mkdir -p $(dirname -- "{output.strand}")
        
        # Get the column sums and store them in unst, forw, and rev, respectively
        # We are interested in columns 2, 3, and 4, which correspond to the number of reads in the unstranded, forward, and reverse strand, respectively
        # Column 1: Gene ID 
        # Column 2: Counts for unstranded RNA-seq
        # Column 3: Counts for  1st read strand aligned with RNA
        # Column 4:Counts for 2nd read strand aligned with RNA
        colsums=$(grep -v "N_" {input.tab} | awk '{{unstranded+=$2;forward+=$3;reverse+=$4}}END{{print unstranded,forward,reverse}}') || colsums="0 1 2"
        
        # Split colsums based on space (create array of three items)
        IFS=" "
        read -ra arr <<< "$colsums"
        
        # Declare these variables as integers
        declare -i unstranded=${{arr[0]}}
        declare -i forward=${{arr[1]}}
        declare -i reverse=${{arr[2]}}
        
        # Increment the denominator by 1 to prevent "divide by 0"
        if [[ $(( reverse / (forward+1) )) -gt 2 ]]; then
            strand_spec="SECOND_READ_TRANSCRIPTION_STRAND"
        elif [[ $(( forward / (reverse+1) )) -gt 2 ]]; then
            strand_spec="FIRST_READ_TRANSCRIPTION_STRAND"
        else
            strand_spec="NONE"
        fi
        
        echo $strand_spec > {output.strand}
        
        picard CollectRnaSeqMetrics I={input.bam} O={output.metrics} REF_FLAT={config[REF_FLAT_FILE]} STRAND_SPECIFICITY=$strand_spec RIBOSOMAL_INTERVALS={config[RRNA_INTERVAL_LIST]}
        """

rule copy_strandedness:
    input: rules.get_rnaseq_metrics.output.strand
    output: os.path.join("MADRID_input", "{tissue_name}", "strandedness", "{sample}", "{tissue_name}_{tag}_strandedness.txt")
    params:
        sample=os.path.join("MADRID_input", "{tissue_name}", "strandedness", "{sample}")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 200 * attempt,  # 200 MB * attempt
        runtime=5
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","copy_strandedness","{sample}_{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
        mkdir -p "{params.sample}"
        cp {input} {output}
        """


if perform.get_insert_size(config=config):
    def insert_size_get_star_data(wildcards):
        return_files = []
        for file in expand(rules.star_align.output.bam_file,zip,tissue_name=get_tissue_name(),tag=get_tags()):
            if wildcards.tissue_name in file:
                return_files.append(file)
        return return_files


    rule get_insert_size:
        input:
            bam=rules.star_align.output.bam_file,
            preround=rules.preroundup.output.layout
        output:
            txt=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "picard", "insert", "{tissue_name}_{tag}_insert_size.txt"),
            pdf=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "picard", "hist", "{tissue_name}_{tag}_insert_size_histo.pdf"),
        # params: layout=os.path.join(config["ROOTDIR"],"MADRID_input", "{tissue_name}", "layouts", "{tissue_name}_{tag}_layout.txt")
        threads: 4
        resources:
            mem_mb=lambda wildcards, attempt: 1000 * 5 * attempt,# 5 GB / attempt
            runtime=lambda wildcards, attempt: 60 * attempt
        conda: "envs/picard.yaml"
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","get_insert_size","{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            lay=$(cat {input.preround})
            if [ $lay == "paired-end"]; then
                picard CollectinsertSizeMetrics \
                I={input.bam} \
                O={output.txt} \
                H={output.pdf} \
                M=0.05
            else
                echo "cannot collect metrics for single-end data" > {output.txt}
                touch {output.pdf}
            fi
            """

    rule copy_insert_size:
        input: rules.get_insert_size.output.txt
        output: os.path.join("MADRID_input", "{tissue_name}", "insertSizeMetrics", "{sample}", "{tissue_name}_{tag}_insert_size.txt")
        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            sample=os.path.join("MADRID_input", "{tissue_name}", "insertSizeMetrics", "{sample}")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 500 * attempt,# 0.5 GB
            runtime=5
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","copy_insert_size","{sample}_{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            mkdir -p {params.sample}
            cp {input} {output}
            """

if perform.get_fragment_size(config=config):
    rule get_fragment_size:
        input:
            bam=rules.star_align.output.bam_file,
            bai=rules.index_bam_file.output
        output:
            os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fragmentSizes", "{tissue_name}_{tag}_fragment_length.txt")
        params:
            layout=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "layouts", "{tissue_name}_{tag}_layout.txt")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 8192 * attempt,  # 8 GB / attempt
            runtime=lambda wildcards, attempt: 90 * attempt  # 90 minutes, should never take this long
        conda: "envs/rseqc.yaml"
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","get_fragment_size","{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            # get matches of script file ( should only be one, but just to be safe run it anyway )
            files=(.snakemake/conda/*/bin/RNA_fragment_size.py)
            
            # run first match      
            python3 ${{files[0]}} -r {config[BED_FILE]} -i {input.bam} > {output}
            """

    rule copy_fragment_size:
        input: rules.get_fragment_size.output
        output: os.path.join("MADRID_input", "{tissue_name}", "fragmentSizes", "{sample}", "{tissue_name}_{tag}_fragment_size.txt")

        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            sample=os.path.join("MADRID_input", "{tissue_name}", "fragmentSizes", "{sample}")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 500 * attempt,# 0.5 GB
            runtime=5
        benchmark: repeat(os.path.join("benchmarks","{tissue_name}","resources","{sample}_{tissue_name}_{tag}.benchmark"), config["BENCHMARK_TIMES"])
        shell:
            """
            mkdir -p {params.sample}
            cp {input} {output}
            """


def multiqc_get_dump_fastq_data(wildcards):
    if perform.prefetch(config=config):
        output = expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "raw",
                "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
            ),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        output = []
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                output.append(os.path.join(path,file))
    return_files = []
    for file in output:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_fastqc_data(wildcards):
    if perform.trim(config=config):
        output_files = expand(
            rules.fastqc_trim.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        output_files = expand(
            rules.fastqc_dump_fastq.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    return_files = []
    for file in output_files:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_star_data(wildcards):
    return_files = []
    for file in expand(
            rules.star_align.output.gene_table,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_screen_data(wildcards):
    if perform.screen(config=config):
        output_files = expand(
            rules.contaminant_screen.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        output_files = []
    return_files = []
    for file in output_files:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_insertsize_data(wildcards):
    return_files = []
    if perform.get_insert_size(config=config):
        for file in expand(
                rules.get_insert_size.output.txt,
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config)
        ):
            if wildcards.tissue_name in file:
                return_files.append(file)
    return return_files


def multiqc_get_fragmentsize_data(wildcards):
    return_files = []
    if perform.get_fragment_size(config=config):
        for file in expand(
                rules.get_fragment_size.output,
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config)
        ):
            if wildcards.tissue_name in file:
                return_files.append(file)
    return return_files


def multiqc_get_rnaseq_data(wildcards):
    return_files = []
    for file in expand(
            rules.get_rnaseq_metrics.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


rule multiqc:
    input:
        fastqc_data=multiqc_get_fastqc_data,
        star_data=multiqc_get_star_data,
        dump_fastq_data=multiqc_get_dump_fastq_data,
        screen_data=multiqc_get_screen_data,
        insertsize_data=multiqc_get_insertsize_data,
        rnaseq_data=multiqc_get_rnaseq_data,
        fragment_size_data=multiqc_get_fragmentsize_data
    output:
        output_file=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "multiqc", "{tissue_name}_multiqc_report.html"),
        output_directory=directory(os.path.join(config["ROOTDIR"],"data", "{tissue_name}", "multiqc"))
    params:
        input_directory=os.path.join(config["ROOTDIR"],"data", "{tissue_name}")
    threads: 1
    conda: "envs/multiqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 10240 * attempt,# 10 GB * attempt
        runtime=lambda wildcards, attempt: int(30 * (attempt * 0.75))  # 30 minutes, don't need much more time than this if it fails
    benchmark: repeat(os.path.join("benchmarks","{tissue_name}","multiqc","{tissue_name}.benchmark"), config["BENCHMARK_TIMES"])
    shell:
        """
        mkdir -p "{output.output_directory}"
        multiqc --force --title "{wildcards.tissue_name}" --filename {wildcards.tissue_name}_multiqc_report.html --outdir {output.output_directory} "{params.input_directory}"
        """
