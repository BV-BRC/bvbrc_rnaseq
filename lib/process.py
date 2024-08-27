#!/usr/bin/python3

# library modules
import sys
import os
import subprocess
import shutil
import concurrent.futures
import glob
import json
import tempfile
import gzip
from threading import Lock
from Bio import SeqIO
from math import log

import pandas as pd
import numpy as np

from bvbrc_api import getQueryDataText


class DifferentialExpression:
    comparisons = None
    genome = None
    recipe = None

    def __init__(self, c):
        print("Creating DifferentialExpression manager")
        self.comparisons = c

    # used to set correct counts matrix delimiter and
    # to run gene+transcript diffexp if host (just gene for bacteria)
    def set_genome(self, g):
        self.genome = g

    def set_recipe(self, r):
        self.recipe = r

    def run_differential_expression(self, output_dir, sample_list):
        if self.recipe == "HTSeq-DESeq" or self.recipe == "Host":
            return self.run_deseq2(output_dir)
        elif self.recipe == "cufflinks":
            return self.run_cuffdiff(output_dir, sample_list)
        else:
            sys.stderr.write(
                "Invalid recipe for differential expression: ", str(self.recipe)
            )
            return False

    def run_cuffdiff(self, output_dir, sample_list):
        threads = 8
        merged_gtf = self.genome.get_genome_data("merge_gtf")
        cuffdiff_cmd = ["cuffdiff", "-p", str(threads), merged_gtf]
        # TODO: replace with contrast or condition class?
        sam_dict = {}
        for sample in sample_list:
            sample_condition = sample.get_condition()
            if sample_condition not in sam_dict:
                sam_dict[sample_condition] = []
            sam_dict[sample_condition].append(sample.get_sample_data("bam"))
        for sample_condition in sam_dict:
            cuffdiff_cmd += [",".join(sam_dict[sample_condition])]
        print("Running Command:\n{0}".format(" ".join(cuffdiff_cmd)))
        try:
            subprocess.check_call(cuffdiff_cmd)
            gmx_file = os.path.join(output_dir, "gene_exp.gmx")
            diff_file = os.path.join(output_dir, "gene_exp.diff")
            self.create_gmx_file([diff_file], gmx_file)
            self.genome.add_genome_data("gmx", gmx_file)
        except Exception as e:
            sys.stderr.write("Error running cuffdiff:\n{0}".format(e))
            return -1

    def create_metadata_file(self, sample_list, output_dir):
        meta_file = os.path.join(output_dir, "sample_metadata.tsv")
        output_text_list = ["Sample\tCondition"]
        for sample in sample_list:
            text = sample.get_id() + "\t" + sample.get_condition()
            output_text_list.append(text)
        output_text = "\n".join(output_text_list) + "\n"
        with open(meta_file, "w") as o:
            o.write(output_text)
        return meta_file

    # TODO: don't let colons exist in the condition names on UI
    def run_deseq2(self, output_dir):
        contrast_list = self.comparisons.get_contrast_list()
        gene_counts = self.genome.get_genome_data(
            f"{self.genome.get_id()}_gene_counts"
        )
        genome_type = self.genome.get_genome_type()
        meta_file = self.genome.get_genome_data("sample_metadata_file")

        if output_dir[-1] != "/":
            output_dir = output_dir + "/"

        deseq_cmd = [
            "run_deseq2_bvbrc",
            gene_counts,
            meta_file,
            output_dir,
            self.genome.get_genome_data("report_img_path"),
            genome_type,
        ]
        # deseq_cmd = ['run_deseq2_bvbrc',gene_counts,meta_file,output_dir,self.genome.get_genome_data('report_img_path'),genome_type]
        # vp_figure = os.path.join(self.genome.get_genome_data('report_img_path'),output_prefix+'volcano_plot.svg')
        vp_figure = os.path.join(
            self.genome.get_genome_data("report_img_path"), "volcano_plot.png"
        )
        ev_cmd = [
            "rnaseq_volcano_plots",
            self.genome.get_genome_data("report_img_path"),
        ]
        contrast_file_list = []
        for contrast in contrast_list:
            deseq_cmd = deseq_cmd + [contrast]
            cond1 = contrast.split(":")[0]
            cond2 = contrast.split(":")[1]
            diffexp_file = os.path.join(
                output_dir, cond1 + "_vs_" + cond2 + ".deseq2.tsv"
            )
            ev_cmd = ev_cmd + [diffexp_file, contrast.replace(":", "_vs_")]
            contrast_file_list.append(diffexp_file)
        try:
            print("Running Command:\n{0}".format(" ".join(deseq_cmd)))
            subprocess.check_call(deseq_cmd)
            self.genome.add_genome_data(
                "contrast_file_list", contrast_file_list
            )
            # invoke volcano plots script
            # rnaseq_volcano_plots.R <output_prefix> <deseq2_file1> <contrast_name1> <deseq2_file2> <contrast_name2>...
            print("Running Command:\n{0}".format(" ".join(ev_cmd)))
            subprocess.check_call(ev_cmd)
            self.genome.add_genome_data("rnaseq_volcano_plots", vp_figure)
        except Exception as e:
            sys.stderr.write("Error running run_deseq2:\n{0}".format(e))
            return -1

        # TODO: run transcript counts
        if self.genome.get_genome_type() == "host":
            print("implement")
            return self.genome.get_genome_data(
                self.genome.get_id() + "_transcript_counts"
            )

    def create_gmx_file(self, init_args, output_file):
        output_handle = open(output_file, "w")
        if len(init_args) < 1:
            print("Usage cuffdiff_to_genematrix.py  <cuffdiff files>")
            exit(0)
        master_list_genes = set()
        master_list_comparisons = set()
        log_lookup = {}
        for input_file in init_args:
            input_handle = open(input_file, "r")
            lines = input_handle.readlines()
            for line in lines[1:]:
                parts = line.strip().split("\t")
                status = "NOT OK"
                if len(parts) == 14:
                    try:
                        # test_id gene_id gene    locus   sample_1        sample_2        status  value_1 value_2 log2(fold_change)       test_stat       p_value q_value significant
                        (
                            gene_col,
                            sample1,
                            sample2,
                            status,
                            value1,
                            value2,
                            log_change,
                        ) = (
                            parts[2],
                            parts[4],
                            parts[5],
                            parts[6],
                            float(parts[7]),
                            float(parts[8]),
                            float(parts[9]),
                        )
                    except ValueError:
                        sys.stderr.write(
                            "One of the input files does not match the formatting of a CuffDiff gene differential expression testing file\n"
                        )
                        sys.exit()
                else:
                    sys.stderr.write(
                        "One of the input files does not match the formatting of a CuffDiff gene differential expression testing file\n"
                    )
                    sys.exit()
                if status != "OK":
                    continue
                gene_ids = []
                if "," in gene_col:
                    gene_ids = gene_col.split(",")
                else:
                    gene_ids = [gene_col]
                changed = False
                if value1 == 0:
                    value1 = 0.01
                    changed = True
                if value2 == 0:
                    value2 = 0.01
                    changed = True
                if changed:
                    log_change = log(value2 / value1) / log(2)
                # pandas would be better for this
                for gene_id in gene_ids:
                    master_list_genes.add(gene_id)
                    comp_id = sample1 + " vs " + sample2
                    master_list_comparisons.add(comp_id)
                    if comp_id not in log_lookup:
                        log_lookup[comp_id] = {}
                    log_lookup[comp_id][gene_id] = log_change
            input_handle.close()
        comparisons = list(master_list_comparisons)
        comparisons.sort()
        headers = ["Gene ID"] + comparisons
        genes = list(master_list_genes)
        genes.sort()
        output_handle.write("\t".join(headers) + "\n")
        for g in genes:
            value_list = []
            for c in comparisons:
                try:
                    current_val = str(log_lookup[c][g])
                except:
                    current_val = "NaN"
                value_list.append(current_val)
            output_handle.write("\t".join([g] + value_list) + "\n")
        output_handle.close()


class GenomeData:
    genome = None
    recipe = None

    def __init__(self):
        print("Creating GenomeData manager")

    def set_genome(self, g):
        self.genome = g

    def set_recipe(self, r):
        self.recipe = r

    def run_queries(self, output_dir, session):
        self.run_subsystems(output_dir, session)
        self.run_pathway(output_dir, session)

    def create_system_figures(self, output_dir):
        if self.recipe == "HTSeq-DESeq":
            self.create_tpm_figures(output_dir)
        elif self.recipe == "cufflinks":
            self.create_fpkm_figures(output_dir)
        else:
            sys.stderr.write("Invalid recipe: exiting create system figures\n")
            return -1

    def create_tpm_figures(self, output_dir):
        superclass_mapping = self.genome.get_genome_data("superclass_mapping")
        pathway_mapping = self.genome.get_genome_data("pathway_mapping")
        genome_counts = self.genome.get_genome_data("tpm")
        if genome_counts is None:
            sys.stderr.write(
                "No tpm's matrix in genome data: exiting create_tpm_figures\n"
            )
            return False
        metadata = self.genome.get_genome_data("sample_metadata_file")

        try:
            superclass_figure = os.path.join(
                self.genome.get_genome_data("report_img_path"),
                "Superclass_Distribution",
            )
            superclass_cmd = [
                "rnaseq_grid_violin_plots",
                superclass_mapping,
                genome_counts,
                metadata,
                superclass_figure,
            ]
            if os.path.exists(superclass_mapping):
                print("Running command:\n{0}".format(" ".join(superclass_cmd)))
                # TODO: ENABLE
                subprocess.check_call(superclass_cmd)
                # self.genome.add_genome_data('superclass_figure',superclass_figure+'.svg')
                self.genome.add_genome_data(
                    "superclass_figure", superclass_figure + ".png"
                )
        except Exception as e:
            sys.stderr.write(
                "Error creating superclass violin plots:\n{0}\n".format(e)
            )

        try:
            pathway_figure = os.path.join(
                self.genome.get_genome_data("report_img_path"),
                "PathwayClass_Distribution",
            )
            pathway_cmd = [
                "rnaseq_grid_violin_plots",
                pathway_mapping,
                genome_counts,
                metadata,
                pathway_figure,
            ]
            if os.path.exists(pathway_mapping):
                print("Running command:\n{0}".format(" ".join(pathway_cmd)))
                # TODO: ENABLE
                subprocess.check_call(pathway_cmd)
                # self.genome.add_genome_data('pathway_figure',pathway_figure+'.svg')
                self.genome.add_genome_data(
                    "pathway_figure", pathway_figure + ".png"
                )
        except Exception as e:
            sys.stderr.write(
                "Error creating pathway violin plots:\n{0}\n".format(e)
            )

    def create_fpkm_figures(self, output_dir):
        superclass_mapping = self.genome.get_genome_data("superclass_mapping")
        pathway_mapping = self.genome.get_genome_data("pathway_mapping")
        genome_counts = self.genome.get_genome_data("fpkm")
        if genome_counts is None:
            sys.stderr.write(
                "No fpkm's matrix in genome data: ",
                "exiting create_fpkm_figures\n",
            )
            return False
        metadata = self.genome.get_genome_data("sample_metadata_file")
        try:
            superclass_figure = os.path.join(
                self.genome.get_genome_data("report_img_path"),
                "Superclass_Distribution",
            )
            superclass_cmd = [
                "rnaseq_grid_violin_plots_cufflinks",
                superclass_mapping,
                genome_counts,
                metadata,
                superclass_figure,
            ]
            if os.path.exists(superclass_mapping):
                print("Running command:\n{0}".format(" ".join(superclass_cmd)))
                # TODO: ENABLE
                subprocess.check_call(superclass_cmd)
                # self.genome.add_genome_data('superclass_figure',superclass_figure+'.svg')
                self.genome.add_genome_data(
                    "superclass_figure", superclass_figure + ".png"
                )
        except Exception as e:
            sys.stderr.write(
                "Error creating superclass violin plots:\n{0}\n".format(e)
            )
        try:
            pathway_figure = os.path.join(
                self.genome.get_genome_data("report_img_path"),
                "PathwayClass_Distribution",
            )
            pathway_cmd = [
                "rnaseq_grid_violin_plots_cufflinks",
                pathway_mapping,
                genome_counts,
                metadata,
                pathway_figure,
            ]
            if os.path.exists(pathway_mapping):
                print("Running command:\n{0}".format(" ".join(pathway_cmd)))
                # TODO: ENABLE
                subprocess.check_call(pathway_cmd)
                # self.genome.add_genome_data('pathway_figure',pathway_figure+'.svg')
                self.genome.add_genome_data(
                    "pathway_figure", pathway_figure + ".png"
                )
        except Exception as e:
            sys.stderr.write(
                "Error creating pathway violin plots:\n{0}\n".format(e)
            )

    def run_pathway(self, output_dir, session):
        # pathway_df = getPathwayDataFrame([self.genome.get_id()], session)
        base = "https://www.bv-brc.org/api/pathway/?http_download=true"
        query = f"eq(genome_id,{self.genome.get_id()})&sort(+id)&limit(2500000)"
        headers = {
            "accept": "application/json",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }
        pathway_df = pd.DataFrame(
            json.loads(getQueryDataText(base, query, headers))
        )
        if pathway_df is not None:
            mapping_table = pathway_df[["patric_id", "pathway_class"]]
            mapping_output = os.path.join(output_dir, "pathway_mapping.tsv")
            mapping_table.to_csv(mapping_output, sep="\t", index=False)
            self.genome.add_genome_data("pathway_mapping", mapping_output)
        else:
            sys.stderr.write("Error, pathway_df is None")
            return -1

    # subsystem_df is a pandas dataframe
    def run_subsystems(self, output_dir, session):
        # subsystem_df = getSubsystemsDataFrame([self.genome.get_id()],session)
        base = "https://www.bv-brc.org/api/subsystem/?http_download=true"
        query = f"eq(genome_id,{self.genome.get_id()})&sort(+id)&limit(2500000)"
        headers = {
            "accept": "application/json",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }
        try:
            subsystem_df = pd.DataFrame(
                json.loads(getQueryDataText(base, query, headers))
            )
        except Exception as e:
            sys.stderr.write(f"Error retrieving subsystems data:\n{e}\n")
            return -1
        if subsystem_df is not None:
            mapping_table = subsystem_df[["patric_id", "superclass"]]
            mapping_output = os.path.join(output_dir, "superclass_mapping.tsv")
            mapping_table.to_csv(mapping_output, sep="\t", index=False)
            self.genome.add_genome_data("superclass_mapping", mapping_output)
        else:
            sys.stderr.write("Error, subsystem_df is None")
            return -1


class Quantify:
    genome = None
    recipe = None

    def __init__(self):
        print("Creating Quantify manager")

    def set_genome(self, g):
        self.genome = g

    def set_recipe(self, r):
        self.recipe = r

    def run_quantification(self, sample_list, threads, output_dir):
        if self.recipe is None:
            sys.stderr.write("Recipe is None: set recipe with set_recipe()")
            return False
        if self.recipe == "HTSeq-DESeq":
            htseq_ret = self.run_htseq(sample_list, threads, output_dir)
            if htseq_ret != 0:
                return htseq_ret
            return self.run_tpmcalc(sample_list, threads)
        elif self.recipe == "Host":
            self.run_stringtie(sample_list, threads)
        elif self.recipe == "cufflinks":
            self.run_cufflinks(sample_list, threads)
        else:
            sys.stderr.write("Invalid recipe: {0}".format(self.recipe))
            return -1

    def run_tpmcalc(self, sample_list, threads):
        tpm_calc_list = []
        sample_details_list = []
        genome_gtf = self.genome.get_genome_data("gtf")
        for sample in sample_list:
            bam_file = sample.get_sample_data("bam")
            tpm_cmd = ["TPMCalculator", "-g", genome_gtf, "-b", bam_file]
            if self.genome.get_genome_type() == "host":
                tpm_cmd += ["-e"]
            if self.genome.get_genome_type() == "bacteria":
                tpm_cmd += ["-k", "gene_id"]
            tpm_output = os.path.join(
                self.genome.get_sample_path(sample.get_id()), "tpm_counts.tsv"
            )
            sample_details_list.append([tpm_output, sample])
            tpm_calc_list.append(tpm_cmd)
        tpm_args_list = list(zip(tpm_calc_list, sample_details_list))
        future_returns = []
        # TPMCalculator outputs to current directory
        if not os.path.exists("TPMCalculator"):
            os.mkdir("TPMCalculator")
        os.chdir("TPMCalculator")
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            future_returns = list(pool.map(self.run_tpmcalc_job, tpm_args_list))
        os.chdir("../")
        for f in future_returns:
            if f != 0:
                sys.stderr.write("Error in HTSeq-count: check logs\n")
                sys.stderr.write("{0}\n".format(future_returns))
                return -1
        return 0

    def run_tpmcalc_job(self, cmd_details):
        cmd = cmd_details[0]
        sample_details = cmd_details[1]
        sample = sample_details[1]
        sample.add_command("tpmcalc_" + self.genome.get_id(), cmd, "running")
        print("Running command:\n{0}\n".format(" ".join(cmd)))
        try:
            subprocess.check_call(cmd)
            sample.set_command_status(
                "tpmcalc" + "_" + self.genome.get_id(), "finished"
            )
            output_file = os.path.abspath(sample.get_id() + "_genes.out")
            if not os.path.exists(output_file):
                sys.stderr.write(
                    "TPMCalculator output file does not exist:\n{0}\n".format(
                        output_file
                    )
                )
            sample.add_sample_data(
                f"{self.genome.get_id()}_tpm_out", output_file
            )
        except Exception as e:
            sys.stderr.write("Error running concurrent tpm job:{0}".format(e))
            sample.set_command_status("tpm" + "_" + self.genome.get_id(), e)
            return -1
        return 0

    def run_htseq(self, sample_list, threads, output_dir):
        # TODO: add strandedness parameter: -s
        # featurey_type: CDS or Gene
        quant_cmd_list = []
        annotation_file = self.genome.get_genome_data("annotation")
        sample_details_list = []
        for sample in sample_list:
            bam_file = sample.get_sample_data("bam")
            quant_cmd = [
                "htseq-count",
                "-n",
                str(threads),
                "-t",
                "gene",
                "-f",
                "bam",
                "-r",
                "pos",
                "-i",
                "ID",
                bam_file,
                annotation_file,
            ]
            quant_cmd_list.append(quant_cmd)
            sample_dir = self.genome.get_sample_path(sample.get_id())
            sample_output_file = os.path.join(
                sample_dir, sample.get_id() + ".counts"
            )
            sample_err_file = os.path.join(
                sample_dir, sample.get_id() + ".htseq_err"
            )
            sample_details_list.append(
                [sample_output_file, sample, sample_err_file]
            )
        quant_args_list = list(zip(quant_cmd_list, sample_details_list))
        future_returns = []
        # redirect stderr
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            future_returns = list(pool.map(self.run_htseq_job, quant_args_list))
        for f in future_returns:
            if f != 0:
                sys.stderr.write("Error in HTSeq-count: check logs\n")
                sys.stderr.write("{0}\n".format(future_returns))
                return -1
        return 0

    def run_htseq_job(self, cmd_details):
        cmd = cmd_details[0]
        sample_details = cmd_details[1]
        output_file = sample_details[0]
        sample = sample_details[1]
        err_file = sample_details[2]
        sample.add_command("htseq" + "_" + self.genome.get_id(), cmd, "running")
        print("Running command:\n{0}\n".format(" ".join(cmd)))
        try:
            with open(output_file, "w") as o, open(err_file, "w") as e:
                subprocess.check_call(cmd, stdout=o, stderr=e)
            sample.set_command_status(
                "htseq" + "_" + self.genome.get_id(), "finished"
            )
            sample.add_sample_data(
                self.genome.get_id() + "_gene_counts", output_file
            )
        except Exception as e:
            sys.stderr.write("Error running concurrent htseq job:{0}".format(e))
            sample.set_command_status("htseq" + "_" + self.genome.get_id(), e)
            return -1
        return 0

    def create_genome_counts_table(self, output_dir, sample_list):
        if self.recipe == "HTSeq-DESeq":
            return self.create_genome_counts_table_htseq(
                output_dir, sample_list
            )
        elif self.recipe == "Host":
            return self.create_genome_counts_table_stringtie(
                output_dir, sample_list
            )
        elif self.recipe == "cufflinks":
            return self.create_genome_counts_table_cufflinks(
                output_dir, sample_list, 8
            )
        else:
            sys.stderr.write(
                "No counts table method found for recipe {0}\n".format(
                    self.recipe
                )
            )
            return None

    def create_genome_counts_table_cufflinks(
        self, output_dir, sample_list, threads
    ):
        gtf_list = []
        for sample in sample_list:
            gtf_list.append(sample.get_sample_data("cuff_gtf"))
        # run cuffmerge to produce a top-level transcripts.gtf file
        cuffmerge_cmd = ["cuffmerge", "-p", str(threads)]
        ref_gtf = self.genome.get_genome_data("gtf")
        if ref_gtf is not None:
            cuffmerge_cmd += ["-g", ref_gtf]
        ref_fasta = self.genome.get_genome_data("fasta")
        if ref_fasta is not None:
            cuffmerge_cmd += ["-s", ref_fasta]
        with open("assembly_list.txt", "w") as o:
            o.write("{0}\n".format("\n".join(gtf_list)))
        cuffmerge_cmd.append("assembly_list.txt")
        print("Running command:\n{0}".format(" ".join(cuffmerge_cmd)))
        merged_gtf = os.path.join(output_dir, "merged_asm/merged.gtf")
        try:
            subprocess.check_call(cuffmerge_cmd)
            self.genome.add_genome_data("merge_gtf", merged_gtf)
        except Exception as e:
            sys.stderr.write("Error running cuffmerge:\n{0}\n".format(e))
            return False

        # cuffquant [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]
        cuffquant_cmd = [
            "cuffquant",
            "-o",
            output_dir,
            "-p",
            str(threads),
            self.genome.get_genome_data("annotation"),
        ]
        # TODO: add strandedness to command
        for sample in sample_list:
            cuffquant_cmd += [sample.get_sample_data("bam")]
        abundance_file = os.path.join(output_dir, "abundances.cxb")
        print("Running command:\n{0}\n".format(" ".join(cuffquant_cmd)))
        try:
            subprocess.check_call(cuffquant_cmd)
            self.genome.add_genome_data("cxb", abundance_file)
        except Exception as e:
            sys.stderr.write("Error running cuffquant\n")
            return -1
        # TODO: rename output file?
        return 0

    # Call rnaseqPrepDE.py script which formats data in the gtf files as a table
    # https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    def create_genome_counts_table_stringtie(self, output_dir, sample_list):
        # create input file for -i option
        output_text_list = []
        avg_len_list = []
        for sample in sample_list:
            gtf_path = sample.get_sample_data(
                self.genome.get_id() + "_merged_transcripts"
            )
            text = sample.get_id() + " " + gtf_path
            output_text_list.append(text)
            avg_len_list.append(int(sample.get_sample_data("avg_read_length")))
        output_text = "\n".join(output_text_list) + "\n"
        path_file = os.path.join(output_dir, "sample_transcript_paths.txt")
        with open(path_file, "w") as o:
            o.write(output_text)
        # run rnaseqPrepDE.py
        gene_matrix_file = os.path.join(output_dir, "gene_counts_matrix.csv")
        transcript_matrix_file = os.path.join(
            output_dir, "transcript_counts_matrix.csv"
        )
        avg_read_length = int(np.average(avg_len_list))
        if avg_read_length == 0:
            sys.stderr.write(
                "Error creating gene and transcript counts table: average read length is 0.\n"
            )
            return None
        # TODO: set path or import or something
        prepde_cmd = [
            "rnaseqPrepDE",
            "-i",
            path_file,
            "-g",
            gene_matrix_file,
            "-t",
            transcript_matrix_file,
            "l",
            str(avg_read_length),
        ]
        try:
            print(" ".join(prepde_cmd))
            subprocess.check_call(prepde_cmd)
            self.genome.add_genome_data(
                self.genome.get_id() + "_gene_counts", gene_matrix_file
            )
            self.genome.add_genome_data(
                self.genome.get_id() + "_transcript_counts",
                transcript_matrix_file,
            )
        except Exception as e:
            sys.stderr.write(
                "Error in rnaseqPrepDE.py: ",
                "cannot generate genome counts or transcript counts file",
                e,
            )
            sys.exit(-1)

    def create_genome_counts_table_htseq(self, output_dir, sample_list):
        genome_df = None
        for sample in sample_list:
            sample_df = pd.read_csv(
                sample.get_sample_data(self.genome.get_id() + "_gene_counts"),
                delim_whitespace=True,
                index_col=0,
                header=None,
                names=[sample.get_id()],
            )
            if genome_df is None:
                genome_df = sample_df
            else:
                genome_df = genome_df.join(sample_df)
        output_file = os.path.join(output_dir, "gene_counts_matrix.tsv")
        genome_df.to_csv(output_file, sep="\t")
        self.genome.add_genome_data(
            self.genome.get_id() + "_gene_counts", output_file
        )
        return output_file

    def create_genome_quant_table(self, output_dir, sample_list):
        if self.recipe == "HTSeq-DESeq":
            return self.create_tpm_table_tpmcalculator(output_dir, sample_list)
        elif self.recipe == "Host":
            return self.create_tpm_table_stringtie(output_dir, sample_list)
        elif self.recipe == "cufflinks":
            # TODO: testing cuffnorm
            return self.create_fpkm_table_cufflinks(output_dir, sample_list)
            # return None

    # outputs to directory 'cuffnorm_output'
    def create_fpkm_table_cufflinks(self, output_dir, sample_list):
        threads = 8
        merged_gtf = self.genome.get_genome_data("merge_gtf")
        cuffnorm_outdir = "cuffnorm_output"
        if not os.path.exists(cuffnorm_outdir):
            os.mkdir(cuffnorm_outdir)
        # TODO: add library-type
        cuffnorm_cmd = [
            "cuffnorm",
            "-p",
            "1",
            "-o",
            cuffnorm_outdir,
            "-library-norm-method",
            "classic-fpkm",
            merged_gtf,
        ]
        sam_dict = {}
        condition_list = []
        sample_id_list = []
        for sample in sample_list:
            sample_condition = sample.get_condition()
            sample_id_list.append(sample.get_id())
            if sample_condition not in sam_dict:
                sam_dict[sample_condition] = []
                condition_list.append(sample_condition)
            sam_dict[sample_condition].append(sample.get_sample_data("bam"))
        for sample_condition in sam_dict:
            cuffnorm_cmd += [",".join(sam_dict[sample_condition])]
        try:
            cuffnorm_err = os.path.join(output_dir, "cuffnorm_output.err")
            print("Running command:\n{0}\n".format(" ".join(cuffnorm_cmd)))
            with open(cuffnorm_err, "w") as err:
                subprocess.check_call(cuffnorm_cmd, stderr=err)
        except Exception as e:
            sys.stderr.write("Error running cuffnorm:\n{0}\n".format(e))
            return -1
        # create fpkm matrix
        try:
            fpkm_file = os.path.join(
                output_dir, "cuffnorm_output/genes.fpkm_table"
            )
            if not os.path.exists(fpkm_file):
                sys.stderr.write(
                    f"Error running cuffnorm: {fpkm_file} does not exist\n"
                )
                return -1
            attr_file = os.path.join(
                output_dir, "cuffnorm_output/genes.attr_table"
            )
            if not os.path.exists(attr_file):
                sys.stderr.write(
                    f"Error running cuffnorm: {attr_file} oes not exist\n"
                )
                return -1
            # print(f'fpkm_file = {fpkm_file}')
            attr_dict = {}
            with open(attr_file, "r") as af:
                af_data = af.readlines()
                af_headers = af_data[0].strip().split("\t")
                for idx, af_line in enumerate(af_data):
                    if idx == 0:
                        continue
                    af_line_parts = af_line.strip().split("\t")
                    attr_dict[af_line_parts[0]] = af_line_parts[4]
            fpkm_output_list = []
            fpkm_header = ["Gene_ID"] + sample_id_list
            fpkm_output_list.append("\t".join(fpkm_header))
            with open(fpkm_file, "r") as ff:
                ff_data = ff.readlines()
                for idx, ff_line in enumerate(ff_data):
                    if idx == 0:
                        continue
                    ff_line_parts = ff_line.strip().split()
                    new_line = (
                        attr_dict[ff_line_parts[0]]
                        + "\t"
                        + "\t".join(ff_line_parts[1:])
                    )
                    fpkm_output_list.append(new_line)
            fpkm_output = os.path.join(output_dir, "fpkm_counts_matrix.tsv")
            fpkm_output_data = "\n".join(fpkm_output_list)
            with open(fpkm_output, "w") as o:
                o.write(fpkm_output_data)
            self.genome.add_genome_data("fpkm", fpkm_output)
        except Exception as e:
            sys.stderr.write("Error parsing fpkm table:\n{0}\n".format(e))
            return -1

    def create_tpm_table_tpmcalculator(self, output_dir, sample_list):
        genome_df = None
        for sample in sample_list:
            sample_df = pd.read_csv(
                sample.get_sample_data(self.genome.get_id() + "_tpm_out"),
                delim_whitespace=True,
            )
            sample_df = sample_df[["Gene_Id", "TPM"]]
            sample_df.rename(
                columns={"Gene_Id": "Gene_Id", "TPM": sample.get_id()},
                inplace=True,
            )
            if genome_df is None:
                genome_df = sample_df
            else:
                genome_df = genome_df.merge(
                    sample_df, how="outer", on="Gene_Id"
                )
        genome_df = genome_df.fillna(0)
        genome_df.set_index("Gene_Id", inplace=True)
        output_file = os.path.join(output_dir, "tpm_counts_matrix.tsv")
        genome_df.to_csv(output_file, sep="\t")
        self.genome.add_genome_data("tpm", output_file)

    def create_tpm_table_stringtie(self, output_dir, sample_list):
        genome_df = None
        for sample in sample_list:
            sample_df = pd.read_csv(
                sample.get_sample_data(
                    f"{self.genome.get_id()}_merged_gene_counts"
                ),
                sep="\t",
            )
            sample_df = sample_df.loc[sample_df["Gene ID"] != "."]
            sample_df = sample_df[["Gene ID", "TPM"]]
            sample_df.rename(
                columns={"Gene ID": "Gene_Id", "TPM": sample.get_id()},
                inplace=True,
            )
            if genome_df is None:
                genome_df = sample_df
            else:
                genome_df = genome_df.merge(
                    sample_df, how="outer", on="Gene_Id"
                )
            genome_df = genome_df.fillna(0)
            genome_df.set_index("Gene_Id", inplace=True)
            output_file = os.path.join(output_dir, "tpm_counts_matrix.tsv")
            genome_df.to_csv(output_file, sep="\t")
            self.genome.add_genome_data("tpm", output_file)

    # TODO: for all sample stuff replace sample.get_id() with self.genome.get_id()
    def run_stringtie(self, sample_list, threads):
        annotation_file = self.genome.get_genome_data("annotation")
        gtf_list = []
        print("sample_list = [{0}],".format(sample_list))
        for sample in sample_list:
            bam_file = sample.get_sample_data("bam")
            sample_dir = self.genome.get_sample_path(sample.get_id())

            gene_output = os.path.join(sample_dir, "gene_abund.tab")
            gtf_output = os.path.join(sample_dir, "transcripts.gtf")
            gtf_list.append(gtf_output)
            # not including -e (included below), including restricts novel feature prediction
            quant_cmd = [
                "stringtie",
                bam_file,
                "-p",
                str(threads),
                "-A",
                gene_output,
                "-G",
                annotation_file,
                "-o",
                gtf_output,
            ]
            if not os.path.exists(gtf_output) or True:
                sample.add_command(
                    "stringtie_" + self.genome.get_id(), quant_cmd, "running"
                )
                print("Running command:\n{0}\n".format(" ".join(quant_cmd)))
                try:
                    subprocess.check_call(quant_cmd)
                    sample.set_command_status(
                        "stringtie_" + self.genome.get_id(), "finished"
                    )
                    sample.add_sample_data(
                        self.genome.get_id() + "_gene_counts", gene_output
                    )
                    sample.add_sample_data(
                        self.genome.get_id() + "_transcripts", gtf_output
                    )
                except Exception as e:
                    sys.stderr.write(
                        "Error running stringtie:\n{0}\n".format(e)
                    )
                    sample.set_command_status(
                        f"stringtie_{self.genome.get_id()}", e
                    )
                    return -1
            else:
                sys.stderr.write(
                    "{0} already exists: skipping stringtie".format(gtf_output)
                )
        # merge reconstructed transcriptomes
        merge_file = os.path.join(self.genome.get_genome_dir(), "merged.gtf")
        if not os.path.exists(merge_file) or True:
            merge_cmd = [
                "stringtie",
                "--merge",
                "-G",
                annotation_file,
                "-o",
                merge_file,
            ] + gtf_list
            try:
                print("Running command:\n{0}\n".format(" ".join(merge_cmd)))
                subprocess.check_call(merge_cmd)
            except Exception as e:
                sys.stderr.write(
                    "ERROR running stringtie-merge:\n{0}".format(e)
                )
                return -1
        self.genome.add_genome_data("merged_gtf", merge_file)

        for sample in sample_list:
            bam_file = sample.get_sample_data("bam")
            sample_dir = self.genome.get_sample_path(sample.get_id())

            gene_output = os.path.join(sample_dir, "merged_gene_abund.tab")
            gtf_output = os.path.join(sample_dir, "merged_transcripts.gtf")
            # -e requires using trancripts found in -G
            # (turns off novel feature prediction)
            quant_cmd = [
                "stringtie",
                "-e",
                bam_file,
                "-p",
                str(threads),
                "-A",
                gene_output,
                "-G",
                merge_file,
                "-o",
                gtf_output,
            ]
            if not os.path.exists(gtf_output) or True:
                sample.add_command(
                    f"stringtie_merged_{self.genome.get_id()}",
                    quant_cmd,
                    "running",
                )
                print("Running command:\n{0}\n".format(" ".join(quant_cmd)))
                try:
                    subprocess.check_call(quant_cmd)
                    sample.set_command_status(
                        "stringtie_merged_" + self.genome.get_id(), "finished"
                    )
                    sample.add_sample_data(
                        self.genome.get_id() + "_merged_transcripts", gtf_output
                    )
                    sample.add_sample_data(
                        self.genome.get_id() + "_merged_gene_counts",
                        gene_output,
                    )
                except Exception as e:
                    sys.stderr.write(
                        "Error running stringtie-merged:\n{0}\n".format(e)
                    )
                    sample.set_command_status(
                        "stringtie_merged_" + self.genome.get_id(), e
                    )
                    return -1
            else:
                sys.stderr.write(
                    "{0} already exists: skipping stringtie merged annotation".format(
                        gtf_output
                    )
                )

    def run_cufflinks(self, sample_list, threads):
        reference = self.genome.get_genome_data("fasta")
        annotation = self.genome.get_genome_data("annotation")
        threads = 8
        for sample in sample_list:
            sample_bam = sample.get_sample_data("bam")
            cufflinks_cmd = [
                "cufflinks",
                "--quiet",
                "-p",
                str(threads),
                "-G",
                annotation,
                "-b",
                reference,
                "-I",
                "50",
                "-o",
                self.genome.get_sample_path(sample.get_id()),
            ]
            # Review: Memory mapped location on each system
            # Attempt to copy to /dev/shm. cufflinks seeks a lot in the file.
            # If that fails, try tmp.
            #
            bam_to_use = None

            try:
                tmpfd, bam_tmp = tempfile.mkstemp(
                    prefix="CUFFL.", dir="/dev/shm"
                )
                os.close(tmpfd)
                shutil.copy(sample_bam, bam_tmp)
                print("Copy succeeded to %s" % (bam_tmp))
                bam_to_use = bam_tmp
            except IOError as err:
                os.unlink(bam_tmp)
                bam_to_use = None
                bam_tmp = None
                sys.stderr.write(
                    "Can't copy %s to %s: %s\n" % (sample_bam, bam_tmp, err)
                )

            if bam_to_use is None:
                try:
                    tmpfd, bam_tmp = tempfile.mkstemp(prefix="CUFFL.", dir=".")
                    os.close(tmpfd)
                    shutil.copy(sample_bam, bam_tmp)
                    bam_to_use = bam_tmp

                except IOError as err:
                    os.unlink(bam_tmp)
                    bam_to_use = None
                    bam_tmp = None
                    sys.stderr.write(
                        "Can't copy %s to %s: %s\n" % (sample_bam, bam_tmp, err)
                    )

            if bam_to_use is None:
                sys.stderr.write("Can't copy %s to tmp space\n" % (sample_bam))
                bam_to_use = sample_bam
                bam_tmp = None

            cufflinks_cmd += [bam_to_use]
            cuff_gtf = os.path.join(
                self.genome.get_sample_path(sample.get_id()), "transcripts.gtf"
            )

            print("Running command:\n{0}".format(" ".join(cufflinks_cmd)))
            try:
                subprocess.check_call(cufflinks_cmd)
                sample.add_sample_data("cuff_gtf", cuff_gtf)
            except Exception as e:
                sys.stderr.write(
                    "Error running cufflinks on sample {0}:\n{1}\n".format(
                        sample.get_id(), e
                    )
                )
                return False

            if bam_tmp is not None:
                sys.stderr.write("remove temp %s\n" % (bam_tmp))
                os.unlink(bam_tmp)
        return True


class Alignment:
    # TODO: check for genome before running functions
    genome = None
    sampled_reads = 2000
    strand_reads = int(sampled_reads / 2)

    def __init__(self):
        print("Creating Alignment manager")

    def set_genome(self, g):
        self.genome = g

    def run_alignment(self, sample, threads):
        # TODO: if exists, add strand param to alignment
        sample_dir = self.genome.get_sample_path(sample.get_id())
        print("sample_dir = {0}".format(sample_dir))
        reads_list = sample.get_reads_as_list()
        sam_file = os.path.join(sample_dir, sample.get_id() + ".sam")
        if self.genome.get_genome_type() == "bacteria":
            align_cmd = [
                "bowtie2",
                "-x",
                self.genome.get_genome_data("bowtie_prefix"),
            ]
        else:
            align_cmd = [
                "hisat2",
                "-x",
                self.genome.get_genome_data("hisat_prefix"),
                "--mp",
                "1,0",
                "--pen-noncansplice",
                "20",
            ]
        if sample.sample_type == "paired":
            align_cmd += ["-1", reads_list[0], "-2", reads_list[1]]
        elif sample.sample_type == "single":
            align_cmd += ["-U", reads_list[0]]
        align_cmd += ["-S", sam_file, "-p", str(threads)]
        sample.add_command(
            f"align_{self.genome.get_id()}", align_cmd, "running"
        )
        align_output_file = sam_file.replace(".sam", ".align_stdout")
        print("Running command:\n{0}".format(" ".join(align_cmd)))
        try:
            # capture stdout for stats later
            # TODO: ENABLE
            with open(align_output_file, "w") as o:
                subprocess.check_call(align_cmd, stderr=o)
            # print captured stdout
            with open(align_output_file, "r") as aof:
                print(aof.read())
            sample.add_sample_data("sam", sam_file)
            sample.add_sample_data(
                self.genome.get_id() + "_align_stats", align_output_file
            )
            sample.set_alignment_status(True)
            sample.set_command_status(
                "align_" + self.genome.get_id(), "finished"
            )
        except Exception as e:
            sys.stderr.write(
                "Sample-alignment encountered an error in Sample ",
                f"{sample.get_id()}:\ncheck error log file\n",
            )
            sample.set_command_status("align" + "_" + self.genome.get_id(), e)
            return False

        bam_file = self.convert_sam_to_bam(sam_file, threads)
        if bam_file:
            sample.add_sample_data("bam", bam_file)
        else:
            sys.stderr.write(
                "Bam file entry does not exist for Sample ",
                f"{sample.get_id()}:\ncheck error log file\n",
            )
            return False
        if not os.path.exists(bam_file):
            sys.stderr.write(
                "Bam file does not exist for Sample ",
                str(sample.get_id()),
                ":\ncheck error log file\n",
            )
            return False
        # remove sam file
        if os.path.exists(sam_file):
            os.remove(sam_file)
        return True

    def check_alignment(self, sample):
        # threshold for number of unique read alignments
        # required by each sample to continue the pipeline
        counts_threshold = 1000
        sample_align_file = sample.get_sample_data(
            self.genome.get_id() + "_align_stats"
        )
        if os.path.exists(sample_align_file):
            try:
                with open(sample_align_file) as saf:
                    align_data = saf.readlines()
                for line in align_data:
                    if "aligned concordantly exactly 1 time" in line:
                        unique_counts = line.split()[0]
                        if int(unique_counts) < counts_threshold:
                            sample.set_alignment_check(False)
                            return False
                sample.set_alignment_check(True)
                return True
            # 407 (0.04%) aligned exactly 1 time
            except Exception as err:
                sys.stderr.write(
                    f"Error checking alignment stats file {sample.get_id()}:\n{err}\n"
                )
                sys.stderr.write("Skipping assessment and hope it works\n")
                sample.set_alignment_check(False)
                return True
        else:
            print(
                f"alignmnt stats output file does not exist for sample {sample.get_id()}, skipping assessment and hope it works"
            )
            sample.set_alignment_check(False)
            return True

    def run_alignment_stats(self, sample, threads):
        sample_dir = self.genome.get_sample_path(sample.get_id())
        sample_bam = sample.get_sample_data("bam")
        # samtools stats
        stats_cmd = ["samtools", "stats", "--threads", str(threads), sample_bam]
        stats_output = os.path.join(
            sample_dir, sample.get_id() + ".samtools_stats"
        )
        sample.add_command(
            "samtools_stats_" + self.genome.get_id(), stats_cmd, "running"
        )
        try:
            # TODO: ENABLE
            print("Running command:\n{0}".format(" ".join(stats_cmd)))
            with open(stats_output, "w") as so:
                subprocess.check_call(stats_cmd, stdout=so)
            sample.set_command_status(
                "samtools_stats_" + self.genome.get_id(), "finished"
            )
        except Exception as e:
            sys.stderr.write(
                "Samtools stats encountered an error in Sample {0}:\ncheck error log file\n".format(
                    sample.get_id()
                )
            )
            sample.set_command_status(
                "samtools_stats_" + self.genome.get_id(), e
            )

        avg_len = self.get_average_read_length_per_file(stats_output)
        sample.add_sample_data("avg_read_length", avg_len)

        # samstat
        samstat_cmd = ["samstat", sample_bam]
        sample.add_command(
            "samstat_" + self.genome.get_id(), samstat_cmd, "running"
        )
        try:
            print("Running command:\n{0}".format(" ".join(samstat_cmd)))
            # TODO: ENABLE
            subprocess.check_call(samstat_cmd)
            sample.set_command_status(
                "samstat_" + self.genome.get_id(), "finished"
            )
        except Exception as e:
            sys.stderr.write(
                "Samstat encountered an error in Sample {0}:\ncheck error log file".format(
                    sample.get_id()
                )
            )
            sample.set_command_status("samstat_" + self.genome.get_id(), e)

    # Reads the output from samtools stat and
    # grabs the average read length value
    def get_average_read_length_per_file(self, stats_file):
        avg_len = 0
        with open(stats_file, "r") as sf:
            for line in sf:
                if "average length:" in line:
                    line = line.strip().split()
                    avg_len = line[-1]
                    break
        return avg_len

    # TODO: incorporate checking the sample alignment results
    def run_sample_alignment(self, sample, threads):
        # TODO: here change sample_dir to the genome directory?
        #   - sample.get_path() to genome.get_sample_path(sample.get_path())
        # sample reads
        sample_dir = sample.get_path()
        reads = sample.get_reads_as_list()
        sampled_reads_list = []
        readNum = 1
        for r in reads:
            sample_file = r.split(".")  # TODO: change where this is output?
            sample_file[len(sample_file) - 1] = "sampled"
            sample_file.append("fq")
            sample_file = ".".join(sample_file)
            sampled_reads_list.append(sample_file)
            sample_cmd = [
                "seqtk",
                "sample",
                "-s",
                "42",
                r,
                str(self.sampled_reads),
            ]
            sample.add_command("sample" + str(readNum), sample_cmd, "running")
            print("Running command:\n{0}".format(" ".join(sample_cmd)))
            try:
                # TODO: ENABLE
                with open(sample_file, "w") as so:
                    subprocess.check_call(sample_cmd, stdout=so)
                sample.set_command_status("sample" + str(readNum), "finished")
            except Exception as e:
                sys.stderr.write(
                    "Sampling encountered an error in Sample {0}:\ncheck error log file".format(
                        sample.get_id()
                    )
                )
                sample.set_command_status("sample" + str(readNum), e)
                return False
            readNum = readNum + 1
        print("sampled reads list:\n{0}".format(sampled_reads_list))
        # align sampled reads
        # TODO: enable for host
        sampled_sam = os.path.join(sample_dir, sample.get_id() + "_sample.sam")
        if self.genome.get_genome_type() == "bacteria":
            sample_align_cmd = [
                "bowtie2",
                "-x",
                self.genome.get_genome_data("bowtie_prefix"),
            ]
        else:
            sample_align_cmd = [
                "hisat2",
                "-x",
                self.genome.get_genome_data("hisat_prefix"),
                "--mp",
                "1,0",
                "--pen-noncansplice",
                "20",
            ]
        if sample.get_type() == "paired":
            sample_align_cmd += [
                "-1",
                sampled_reads_list[0],
                "-2",
                sampled_reads_list[1],
            ]
        elif sample.sample_type == "single":
            sample_align_cmd += ["-U", sampled_reads_list[0]]
        sample_align_cmd += ["-S", sampled_sam, "-p", str(threads)]
        sample.add_command("sample_align", sample_align_cmd, "running")
        print("Running command:\n{0}".format(" ".join(sample_align_cmd)))
        try:
            # TODO: ENABLE
            subprocess.check_call(sample_align_cmd)
            sample.set_command_status("sample_align", "finished")
        except Exception as e:
            sys.stderr.write(
                "Sample-alignment encountered an error ",
                f"in Sample {sample.get_id()}:",
                "\ncheck error log file",
            )
            sample.set_command_status("sample_align", e)
            return False

        # TODO: condition for assigning strandedness?
        infer_cmd = [
            "infer_experiment.py",
            "-i",
            sampled_sam,
            "-r",
            self.genome.get_genome_data("bed"),
            "-s",
            str(self.strand_reads),
        ]
        infer_file = os.path.join(sample_dir, sample.get_id() + "_strand.infer")
        sample.add_command("infer_strand", infer_cmd, "running")
        print("Running command:\n{0}".format(" ".join(infer_cmd)))
        try:
            # TODO: ENABLE
            with open(infer_file, "w") as o:
                subprocess.check_call(infer_cmd, stdout=o)
            sample.set_command_status("infer_strand", "finished")
            sample.add_sample_data("infer_strand_file", infer_file)
            strand = self.infer_strand_from_file(
                sample.get_sample_data("infer_strand_file")
            )
            sample.add_sample_data("strand", strand)
        except Exception as e:
            sys.stderr.write(
                "Infer strand encountered an error in ",
                f"Sample {sample.get_id()}:\ncheck error log file",
            )
            sample.set_command_status("infer_strand", e)
            return False

        for sampled_reads_file in sampled_reads_list:
            if os.path.exists(sampled_reads_file):
                os.remove(sampled_reads_file)
        if os.path.exists(sampled_sam):
            os.remove(sampled_sam)

    def convert_sam_to_bam(self, sam_file, threads):
        bam_file = sam_file.replace(".sam", ".bam")
        print("bam_file = {0}".format(bam_file))
        sam_to_bam_cmd = (
            "samtools view -Su "
            + sam_file
            + " | samtools sort -o - - -@ "
            + str(threads)
            + " > "
            + bam_file
        )
        try:
            print("Running command:\n{0}".format(sam_to_bam_cmd))
            # TODO: ENABLE
            subprocess.check_call(sam_to_bam_cmd, shell=True)
        except Exception as e:
            sys.stderr.write(
                "Error in converting sam to bam file:\n{0}\n".format(e)
            )
            return None
        index_cmd = "samtools index " + bam_file
        try:
            print("Running command:\n{0}".format(index_cmd))
            # TODO: ENABLE
            subprocess.check_call(index_cmd, shell=True)
        except Exception as e:
            sys.stderr.write("Error indexing bam file:\n{0}\n".format(e))
            return None
        return bam_file

    def infer_strand_from_file(self, infer_file):
        strand = None
        with open(infer_file, "r") as handle:
            for x in range(0, 3):
                next(handle)
            undetermined = float(next(handle).split(":")[1].strip())
            rf_stranded = float(next(handle).split(":")[1].strip())
            fr_stranded = float(next(handle).split(":")[1].strip())
            if undetermined > fr_stranded and undetermined > rf_stranded:
                strand = "undetermined"
            elif (
                abs(fr_stranded - rf_stranded) <= 0.1
            ):  # subjective threshold: should work in most cases
                strand = "undetermined"
            else:
                strand = "FR" if fr_stranded > rf_stranded else "RF"
        return strand


class DiffExpImport:
    genome = None
    recipe = None

    def __init__(self):
        print("Creating differential expression import manager")

    def set_genome(self, g):
        self.genome = g

    def set_recipe(self, r):
        self.recipe = r

    def write_gmx_file(self, output_dir):
        print("writing gmx file")
        contrast_file_list = self.genome.get_genome_data("contrast_file_list")
        contrast_list = []
        gene_count_dict = {}
        gene_set = set()
        if contrast_file_list is None:
            print("Skipping contrast file export; contrast_file_list is None")
        else:
            for contrast_file in contrast_file_list:
                contrast_name = os.path.basename(contrast_file).replace(
                    ".tsv", ""
                )
                contrast_list.append(contrast_name)
                gene_count_dict[contrast_name] = {}
                with open(contrast_file, "r") as cf:
                    next(cf)
                    for line in cf:
                        (
                            gene,
                            baseMean,
                            log2FC,
                            lfcSE,
                            stat,
                            pvalue,
                            padj,
                        ) = line.strip().split("\t")
                        # strip 'gene-' from identifiers for host
                        gene_set.add(gene.replace("gene-", ""))
                        gene_count_dict[contrast_name][gene] = log2FC
        gmx_output = os.path.join(output_dir, "gene_exp.gmx")
        self.genome.add_genome_data("gmx", gmx_output)
        # TODO: rewrite this?
        with open(gmx_output, "w") as o:
            o.write("Gene_ID\t%s\n" % "\t".join(contrast_list))
            for gene in gene_set:
                o.write(gene)
                for contrast in contrast_list:
                    if gene in gene_count_dict[contrast]:
                        o.write("\t%s" % gene_count_dict[contrast][gene])
                    else:
                        o.write("\t0")
                o.write("\n")

    def run_diff_exp_import(self, output_dir, map_args):
        if (
            self.recipe == "HTSeq-DESeq" or self.recipe == "Host"
        ):  # create gmx file from DESeq2 results
            self.write_gmx_file(output_dir)
        # elst cufflinks, file should already exist
        gmx_file = self.genome.get_genome_data("gmx")
        transform_script = "expression_transform_bvbrc"
        if gmx_file is None:
            sys.stderr.write(
                "gmx_file is null, exiting differential expression import\n"
            )
            return False
        if os.path.exists(gmx_file):
            experiment_path = os.path.join(output_dir, map_args.d)
            print(f"experiment_path={experiment_path}")
            subprocess.call(["mkdir", "-p", experiment_path])
            transform_params = {
                "output_path": experiment_path,
                "xfile": gmx_file,
                "xformat": "tsv",
                "xsetup": "gene_matrix",
                "source_id_type": "patric_id",
                "data_type": "Transcriptomics",
                "experiment_title": "RNA-Seq",
                "experiment_description": "RNA-Seq",
                "organism": self.genome.get_id(),
            }
            diffexp_json = self.setup_diffexp_json()
            params_file = os.path.join(output_dir, "diff_exp_params.json")
            # TODO: incorporate check for sstring
            with open(params_file, "w") as params_handle:
                params_handle.write(json.dumps(transform_params))
            convert_cmd = [
                transform_script,
                "--ufile",
                params_file,
                "--sstring",
                map_args.sstring,
                "--output_path",
                experiment_path,
                "--xfile",
                gmx_file,
            ]
            print(" ".join(convert_cmd))
            try:
                subprocess.check_call(convert_cmd)
            except subprocess.CalledProcessError:
                sys.stderr.write(
                    "Running differential expression import failed.\n"
                )
                # subprocess.call(["rm","-rf",experiment_path])
                return
            diffexp_obj_file = os.path.join(
                output_dir, os.path.basename(map_args.d).lstrip(".")
            )
            with open(diffexp_obj_file, "w") as diffexp_job:
                diffexp_job.write(json.dumps(diffexp_json))
            return True
        else:
            sys.stderr.write(
                "GMX file does not exist, ",
                "exiting differential expression import",
            )
            return False

    def setup_diffexp_json(self):
        # job template for differential expression object
        diffexp_json = json.loads(
            """                    {
                        "app": {
                            "description": "Parses and transforms users differential expression data",
                            "id": "DifferentialExpression",
                            "label": "Transform expression data",
                            "parameters": [
                                {
                                    "default": null,
                                    "desc": "Comparison values between samples",
                                    "id": "xfile",
                                    "label": "Experiment Data File",
                                    "required": 1,
                                    "type": "wstype",
                                    "wstype": "ExpList"
                                },
                                {
                                    "default": null,
                                    "desc": "Metadata template filled out by the user",
                                    "id": "mfile",
                                    "label": "Metadata File",
                                    "required": 0,
                                    "type": "wstype",
                                    "wstype": "ExpMetadata"
                                },
                                {
                                    "default": null,
                                    "desc": "User information (JSON string)",
                                    "id": "ustring",
                                    "label": "User string",
                                    "required": 1,
                                    "type": "string"
                                },
                                {
                                    "default": null,
                                    "desc": "Path to which the output will be written. Defaults to the directory containing the input data. ",
                                    "id": "output_path",
                                    "label": "Output Folder",
                                    "required": 0,
                                    "type": "folder"
                                },
                                {
                                    "default": null,
                                    "desc": "Basename for the generated output files. Defaults to the basename of the input data.",
                                    "id": "output_file",
                                    "label": "File Basename",
                                    "required": 0,
                                    "type": "wsid"
                                }
                            ],
                            "script": "App-DifferentialExpression"
                        },
                        "elapsed_time": null,
                        "end_time": null,
                        "hostname": "",
                        "id": "",
                        "is_folder": 0,
                        "job_output": "",
                        "output_files": [
                        ],
                        "parameters": {
                            "mfile": "",
                            "output_file": "",
                            "output_path": "",
                            "ustring": "",
                            "xfile": ""
                        },
                        "start_time": ""
                    }
        """
        )
        return diffexp_json


# Preprocessing: fastqc, trimming, sampled alignment(alignment?),
# strandedness(alignment?)
class Preprocess:
    def __init__(self):
        print("Creating Preprocess manager")

    def check_reads_worker(self,sample,reads_errors):
        reads = sample.get_reads_as_list()
        minReads = 2000
        all_good = True
        curr_errors = []
        if len(reads) == 2: # paired
            if reads[0].endswith('.gz'):
                with gzip.open(reads[0],'rt') as r1, gzip.open(reads[1],'rt') as r2:        
                    r1_ids = {record.id.split()[0] for record in SeqIO.parse(r1,'fastq')}
                    r2_ids = {record.id.split()[0] for record in SeqIO.parse(r2,'fastq')}
            else:
                with open(reads[0],'r') as r1, open(reads[1],'r') as r2:        
                    r1_ids = {record.id.split()[0] for record in SeqIO.parse(r1,'fastq')}
                    r2_ids = {record.id.split()[0] for record in SeqIO.parse(r2,'fastq')}

            if r1_ids != r2_ids:
                all_good = False
                curr_errors.append(f'sample {sample.get_id()} paired reads file is not paired correctly') 

            if len(r1_ids) < minReads:
                all_good = False
                curr_errors.append(f'too few reads in file {reads[0]} ')
            if len(r2_ids) < minReads:
                all_good = False
                curr_errors.append(f'too few reads in file {reads[1]}')

            unpaired_r1 = r1_ids - r2_ids
            unpaired_r2 = r2_ids - r1_ids
            if unpaired_r1 or unpaired_r2:
                print(f'unpaired reads found in sample {sample.get_id()}')
                if unpaired_r1:
                    curr_errors.append(f'unpaired reads found in {reads[0]}')
                if unpaired_r2:
                    curr_errors.append(f'unpaired reads found in {reads[1]}')
                all_good = False
        else: # single
            with open(reads[0],'r') as r:
                read_ids = {record.id.split()[0] for record in SeqIO.parse(r,'fastq')}
                if len(read_ids) < minReads:
                    all_good = False
                    curr_errors.append(f'too few reads in file {reads[0]}')
        if len(curr_errors) > 0:
            with thread_lock:
                reads_list += curr_errors
        return all_good

    # reads_errors is a list passed in from main
    def check_reads(self, sample_list, reads_errors, threads):
        thread_lock = Lock()
        reads_results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            future_returns = [pool.submit(self.check_reads_worker, sample, reads_errors) for sample in sample_list] 
            for result in concurrent.futures.as_completed(future_returns):
                reads_results.append(result)
        return all(reads_results)

    def run_fastqc(self, sample):
        sample_dir = sample.get_path()
        reads = sample.get_reads_as_list()
        fastqc_cmd = ["fastqc", "--outdir", sample_dir]
        fastqc_cmd += reads
        sample.add_command("fastqc", fastqc_cmd, "running")
        print("Running command:\n {0}".format(" ".join(fastqc_cmd)))
        try:
            # TODO: ENABLE
            subprocess.check_call(fastqc_cmd)
            sample.set_command_status("fastqc", "finished")
        except Exception as e:
            sys.stderr.write(
                "FastQC encountered an error in Sample {0}:\ncheck error log file".format(
                    sample.get_id()
                )
            )
            sample.set_command_status("fastqc", e)
            return False

    # runs trim_galore on samples
    # TODO: option for rerunning fastqc?
    # TODO: option for adapter specification
    def run_trimming(self, sample, threads):
        sample_dir = sample.get_path()
        reads = sample.get_reads_as_list()
        trimmed_reads = []
        trim_cmd = [
            "trim_galore",
            "--output_dir",
            sample_dir,
            "--cores",
            str(threads),
        ]
        if sample.get_type() == "paired":
            # trimmed_reads.append(os.path.join(sample_dir,os.path.basename(reads[0]).split(".")[0]+"_val_1.fq.gz"))
            # trimmed_reads.append(os.path.join(sample_dir,os.path.basename(reads[1]).split(".")[0]+"_val_2.fq.gz"))
            trim_cmd += ["--paired"]
        # if sample.get_type() == "single":
        #    trimmed_reads.append(os.path.join(sample_dir,os.path.basename(reads[0]).split('.')[0]+"_trimmed.fq.gz"))
        trim_cmd += reads
        sample.add_command("trim", trim_cmd, "running")
        trim_galore_stderr = os.path.join(sample_dir, "trim_galore.stderr")
        try:
            print("Running command:\n{0}".format(" ".join(trim_cmd)))
            # TODO: ENABLE
            with open(trim_galore_stderr, "w") as err:
                subprocess.check_call(trim_cmd, stderr=err)
            sample.set_command_status("trim", "finished")
            if sample.get_type() == "paired":
                read_parts1 = os.path.basename(reads[0]).split(".")
                read_parts2 = os.path.basename(reads[1]).split(".")
                new_r1 = glob.glob(os.path.join(sample_dir, "*_val_1.*"))[0]
                new_r2 = glob.glob(os.path.join(sample_dir, "*_val_2.*"))[0]
                trimmed_reads.append(new_r1)
                trimmed_reads.append(new_r2)
            else:
                read_parts = os.path.basename(reads[0]).split(".")
                new_r = glob.glob(os.path.join(sample_dir, "*_trimmed.*"))[0]
                trimmed_reads.append(new_r)
            sample.set_reads_list(trimmed_reads)
            for cutadapt_file in glob.glob("./*cutadapt.log"):
                os.remove(cutadapt_file)
        except Exception as e:
            sys.stderr.write(
                "Trimming encountered an error in Sample {0}:\ncheck error log file".format(
                    sample.get_id()
                )
            )
            sample.set_command_status("trim", e)
            return False
        # TODO: remove cutadapt log files: _cutadapt.log
        for r in sample.get_reads_as_list():
            if not os.path.exists(r):
                print("{0} does not exist".format(r))
