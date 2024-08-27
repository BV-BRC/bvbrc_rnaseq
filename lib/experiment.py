#!/usr/bin/python3

import sys
import os
import subprocess
import json
import io
import bvbrc_api as bvb


class Genome:
    valid_data_types = [
        "fasta",
        "annotation",
        "bed",
        "hisat_index",
        "hisat_prefix",
        "bowtie_prefix",
        "counts_table",
        "merged_gtf",
        "sample_metadata_file",
    ]
    valid_genome_types = ["bacteria", "host"]
    # if genome_type is bacteria: HTSeq quant
    # if genome_type is host: Stringtie quant
    # if cufflinks==true in job_data, run old pipeline

    genome_id = None
    genome_ref_id = None
    genome_name = None  # TODO: incorporate into paths
    genome_type = None  # bacteria or host
    genome_dir = None
    genome_data = (
        None  # Dictionary that holds fasta, annotation, hisat index, etc
    )
    sample_path_dict = None

    def __init__(self, gi, gt, session, genome_query=True):
        self.genome_id = gi
        self.genome_ref_id = gi
        self.genome_type = gt
        if self.genome_type not in self.valid_genome_types:
            print(
                "{0} is not a valid genome type:\n{1}".format(
                    self.genome_type, ",".join(self.valid_genome_types)
                )
            )
            return None
        self.sample_path_dict = {}
        self.genome_data = {}
        # get genome id prefix
        if genome_query:
            base = "https://www.bv-brc.org/api/genome/?"
            query = f"eq(genome_id,{gi})"
            headers = {
                "accept": "application/json",
                "content-type": "application/rqlquery+x-www-form-urlencoded",
                "Authorization": session.headers["Authorization"],
            }
            print(
                "genome_query:\nurl = {0}\n{1}\n".format(base + query, headers)
            )
            genome_res = bvb.getQueryDataText(base, query, headers)
            # res['response']['docs'][0]['common_name']
            response_json = json.load(io.StringIO(genome_res))
            if len(response_json) == 0:
                self.genome_name = self.genome_id
            else:
                if "common_name" in response_json[0]:
                    self.genome_name = response_json[0]["common_name"]
                else:
                    self.genome_name = self.genome_id
        else:
            self.genome_name = self.genome_id

    def add_genome_data(self, key, data):
        # TODO: (maybe not)check if key is valid
        # if key in self.valid_data_types:
        self.genome_data[key] = data
        # else:
        #    print("{0} not a valid genome data type:\n{1}".format(key,",".join(self.valid_data_types)))

    def get_genome_data(self, key):
        if key not in self.genome_data:
            print("{0} not in genome {1}".format(key, self.genome_id))
            return None
        else:
            return self.genome_data[key]

    def set_genome_dir(self, gd):
        self.genome_dir = gd

    def get_genome_dir(self):
        return self.genome_dir

    def get_ref_id(self):
        return self.genome_ref_id

    def get_id(self):
        return self.genome_id

    def get_genome_name(self):
        if self.genome_name:
            return self.genome_name
        else:
            return self.genome_id

    def get_path_dict(self):
        return self.sample_path_dict

    def get_sample_path(self, sample_id):
        if sample_id in self.sample_path_dict:
            return self.sample_path_dict[sample_id]
        else:
            return None

    def get_genome_type(self):
        return self.genome_type

    def create_path_entry(self, sample_id, path):
        self.sample_path_dict[sample_id] = path

    def setup_genome_database(self):
        if self.genome_type == "bacteria":
            # setup bowtie index
            index_prefix = os.path.join(
                self.get_genome_dir(), self.get_id()
            )  # TODO: check to make sure it's correct
            bowtie_build_cmd = [
                "bowtie2-build",
                self.get_genome_data("fasta"),
                index_prefix,
            ]
            try:
                print(" ".join(bowtie_build_cmd))
                subprocess.check_call(bowtie_build_cmd)
                self.add_genome_data("bowtie_prefix", index_prefix)
                # return True
            except Exception as e:
                sys.stderr.write(
                    f"bowtie build failed for genome {self.get_id()}\n"
                )
                print("error: {0}".format(e))
                return False
        elif self.genome_type == "host":
            print(self.genome_data)
            index_prefix = os.path.join(
                self.get_genome_dir(),
                os.path.basename(
                    self.get_genome_data("hisat_index").replace(".ht2.tar", "")
                ),
            )
            tar_cmd = [
                "tar",
                "-xvf",
                self.get_genome_data("hisat_index"),
                "-C",
                os.path.dirname(self.get_genome_data("hisat_index")),
            ]
            try:
                print("unpacking hisat2 index:\n{0}".format(" ".join(tar_cmd)))
                subprocess.check_call(tar_cmd)
                print("finished unpcking hisat2 index")
            except Exception as e:
                sys.stderr.write("Error unpacking hisat2 index: exiting")
                print("error: {0}".format(e))
                sys.exit(-1)
            self.add_genome_data("hisat_prefix", index_prefix)
        # setup bed file
        bed_file = self.get_genome_data("annotation").replace(".gff", ".bed")
        # TODO: add --do-not-sort after gff2bed?
        bed_cmd = "gff2bed < " + self.get_genome_data("annotation")
        try:
            print(bed_cmd)
            with open(bed_file, "w") as bf:
                subprocess.check_call(bed_cmd, stdout=bf, shell=True)
            self.add_genome_data("bed", bed_file)
        except Exception as e:
            sys.stderr.write(
                "bed file generation failed for genome {0}\n".format(
                    self.get_id()
                )
            )
            print("error: {0}".format(e))
            return False

        # setup gtf file
        # TPMCalculator -g ~/Genomes/9606/GCF_000001405.39_GRCh38.p13_genomic.mod.gtf -b Test_Htseq/9606/avirulent/replicate1/SRR10307420.bam -e
        gff_to_gtf_cmd = [
            "gffread",
            self.get_genome_data("annotation"),
            "-T",
            "-o-",
        ]
        gtf_output = self.get_genome_data("annotation").replace(".gff", ".gtf")
        try:
            print(" ".join(gff_to_gtf_cmd))
            gtg_proc = subprocess.Popen(gff_to_gtf_cmd, stdout=subprocess.PIPE)
            gtf_text_lines = gtg_proc.stdout.read().decode().split('\n')
            gtf_text = []
            for line in gtf_text_lines:
                line = line.strip() 
                if 'gene_id' in line:
                    gtf_text.append(line)
            gtf_text = '\n'.join(gtf_text).replace("\ttranscript\t", "\texon\t")
            with open(gtf_output, "w") as o:
                o.write(gtf_text)
                o.write('\n')
            self.add_genome_data("gtf", gtf_output)
        except Exception as e:
            sys.stderr.write("Error converting gff to gtf:\n{0}\n".format(e))
            return False

    def get_genome_database_prefix(self):
        return None

    # TODO: make more robust
    def __repr__(self):
        return "Test print genome repr: {0}".format(self.genome_id)

    # TODO: make more robust
    def __str__(self):
        return "Test print genome str: {0}".format(self.genome_id)


class Sample:
    # Class Enums
    valid_sample_types = ["paired", "single", "sra"]

    sample_id = None
    sample_type = None  # paired, single, sra
    sample_path = None
    reads_list = None  # list (paired or single) or none (sra)
    accession = (
        None  # TODO: maybe remove?, string (sra) or none (paired or single)
    )
    condition = None
    sample_data = None
    # TODO: sample_output_folder
    # TODO: think about how to store pipeline commands
    command_dict = None
    command_status_dict = (
        None  # TODO: values "running", "finished" or an execution error
    )
    # flags for different parts in the pipeline
    alignment_success = False
    alignment_check = False

    # Initialize Sample class object
    # Parameters:
    #   - si: sample id
    #   - st: sample type
    #   - fp: file path
    #   - ac: accession
    #   - c: condition
    def __init__(self, si, st, rl, ac, c):
        self.sample_id = si
        self.sample_type = st
        # TODO: change sample type based on SRA single or paired
        if self.sample_type not in self.valid_sample_types:
            print(
                "{0} not a valid sample type: {1}".format(
                    self.sample_type, ",".join(self.valid_sample_types)
                )
            )
            return None
        self.reads_list = rl
        self.accession = ac
        self.condition = c
        self.command_dict = {}
        self.command_status_dict = {}
        self.sample_data = {}

    def add_command(self, key, command, status):
        if key in self.command_dict:
            print(
                "{0} already exists in command_dict, status {1}:\n{2}".format(
                    key, self.command_status_dict[key], self.command_dict[key]
                )
            )
        else:
            self.command_dict[key] = command
            self.command_status_dict[key] = status

    def set_command_status(self, key, status):
        if key not in self.command_dict:
            print("{0} does not exist in command dict".format(key))
        else:
            self.command_status_dict[key] = status

    def add_sample_data(self, key, data):
        # TODO: check if key is valid
        self.sample_data[key] = data

    def get_sample_data(self, key):
        if key not in self.sample_data:
            print("{0} not in sample {1}".format(key, self.sample_id))
            return None
        else:
            return self.sample_data[key]

    def check_key(self, key):
        if key in self.sample_data:
            return True
        else:
            return False

    def set_path(self, path):
        self.sample_path = path

    def get_path(self):
        return self.sample_path

    def get_id(self):
        return self.sample_id

    def get_type(self):
        return self.sample_type

    def get_reads_as_list(self):
        return self.reads_list

    def set_reads_list(self, rl):
        self.reads_list = rl

    def get_condition(self):
        return self.condition

    def get_alignment_status(self):
        return self.alignment_success

    def set_alignment_status(self, status):
        self.alignment_success = status

    def set_alignment_check(self, status):
        self.alignment_check = status

    def get_alignment_check(self):
        return self.alignment_check

    # TODO: make more robust
    def __repr__(self):
        return "Test print sample repr: {0}".format(self.sample_id)

    # TODO: make more robust
    def __str__(self):
        return "Test print sample str: {0}".format(self.sample_id)


class Condition:
    condition = None
    cond_path = None
    sample_list = None

    def __init__(self, c):
        self.condition = c
        self.sample_list = []

    def add_sample(self, sample):
        self.sample_list.append(sample)

    def set_path(self, path):
        self.cond_path = path

    def get_path(self):
        return self.cond_path

    def get_sample_list(self):
        return self.sample_list

    def get_condition(self):
        return self.condition

    def __str__(self):
        return "Test print condition str: {0}\n{1}\n".format(
            self.condition, self.sample_list
        )

    def __repr__(self):
        return "Test print condition repr: {0}\n{1}\n".format(
            self.condition, self.sample_list
        )


class Comparison:
    contrast_list = None

    def __init__(self):
        self.contrast_list = []

    def add_contrast(self, condition1, condition2):
        self.contrast_list.append([condition1, condition2])

    def get_contrast_list(self):
        ret_list = []
        for pair in self.contrast_list:
            ret_list.append(":".join(pair))
        return ret_list

    def check_diffexp(self):
        return len(self.contrast_list) > 0
