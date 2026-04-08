#!/usr/bin/python3

# library modules
import argparse, sys, os, json
import requests
import subprocess

# pipeline modules
import experiment
import process
import report

# import process_statistics
from bvbrc_api import authenticateByEnv

# https://www.nature.com/articles/s41598-020-76881-x

# valid recipes
valid_recipes = ["HTSeq-DESeq", "cufflinks", "Host"]


def validate_job_data(job_data):
    """Validate required fields in job JSON before pipeline starts."""
    errors = []
    for field in ["reference_genome_id", "genome_type", "recipe", "output_path"]:
        if field not in job_data:
            errors.append("Missing required field: {0}".format(field))
    if "genome_type" in job_data and job_data["genome_type"] not in ["bacteria", "host"]:
        errors.append("Invalid genome_type: {0}".format(job_data["genome_type"]))
    if "recipe" in job_data and job_data["recipe"] not in valid_recipes:
        errors.append("Invalid recipe: {0}".format(job_data["recipe"]))
    # Validate contrasts structure if present
    contrasts = job_data.get("contrasts")
    if contrasts is not None:
        for i, con in enumerate(contrasts):
            if not isinstance(con, list) or len(con) < 2:
                errors.append("Contrast {0} must be a list of 2 condition names".format(i))
    if errors:
        sys.stderr.write("Job validation errors:\n")
        for e in errors:
            sys.stderr.write("  - {0}\n".format(e))
        sys.exit(1)


def main(
    genome, experiment_dict, tool_params, output_dir, comparisons, session, map_args
):
    # setup folder structure and genome databases
    setup(output_dir, experiment_dict, genome)
    diffexp_flag = comparisons.check_diffexp()

    # get all samples as a list for certain functions 
    sample_list = []
    for condition in experiment_dict:
        samples = experiment_dict[condition].get_sample_list()
        sample_list = sample_list + samples
    print("sample_list = {0}".format(sample_list))

    # process data independently of genomes
    # assess reads and terminate gracefully if reads are not set up correctly
    # Fastqc
    preprocess = process.Preprocess()
    reads_errors = []
    reads_all_good = preprocess.check_reads(sample_list, reads_errors, 8)
    for condition in experiment_dict:
        for sample in experiment_dict[condition].get_sample_list():
            preprocess.run_fastqc(sample)
    if not reads_all_good:
        # report manager and terminate
        report_manager = report.ReportManager()
        report_stats = {}
        sample_count = 0
        condition_count = 0
        for condition in experiment_dict:
            sample_count += len(experiment_dict[condition].get_sample_list())
            if condition != "no_condition":
                condition_count += 1
        report_stats["num_samples"] = sample_count
        report_stats["num_conditions"] = condition_count
        report_stats["recipe"] = map_args.recipe
        report_stats["reads_failure"] = reads_errors
        report_manager.run_multiqc(output_dir)
        report_manager.create_report(
            genome,
            output_dir,
            experiment_dict,
            report_stats,
            map_args.workspace_dir,
            diffexp_flag,
        )
        sys.exit(0)

    # Trimming (non-fatal: use original reads as fallback on failure)
    if map_args.trimming:
        for condition in experiment_dict:
            for sample in experiment_dict[condition].get_sample_list():
                trim_result = preprocess.run_trimming(sample, 8)
                if trim_result is False:
                    sys.stderr.write(
                        f"Warning: trimming failed for sample {sample.get_id()}, "
                        "continuing with original reads\n"
                    )

    # Sampled align against genome
    # TODO: assess strandedness with one genome?
    alignment = process.Alignment()
    alignment.set_genome(genome)
    for condition in experiment_dict:
        for sample in experiment_dict[condition].get_sample_list():
            try:
                alignment.run_sample_alignment(sample, 8)
            except Exception as e:
                sys.stderr.write(f'Error running sample alignment on sample {sample.get_id()}:\n{e}\n')

    # Align against genome
    alignment.set_genome(genome)
    alignment_all_good = True
    for condition in experiment_dict:
        for sample in experiment_dict[condition].get_sample_list():
            try:
                align_complete = alignment.run_alignment(sample, 8)
                align_result = alignment.check_alignment(sample)
                if align_complete and align_result:
                    alignment.run_alignment_stats(sample, 8)
                else:
                    alignment_all_good = False
            except Exception as e:
                sys.stderr.write(f'Error running alignment on sample {sample.get_id()}:\n{e}\n')
                alignment_all_good = False

    if not alignment_all_good:
        report_manager = report.ReportManager()
        # number of samples and conditions
        report_stats = {}
        sample_count = 0
        condition_count = 0
        for condition in experiment_dict:
            sample_count += len(experiment_dict[condition].get_sample_list())
            if condition != "no_condition":
                condition_count += 1
        report_stats["num_samples"] = sample_count
        report_stats["num_conditions"] = condition_count
        # get recipe
        report_stats["recipe"] = map_args.recipe
        report_manager.run_multiqc(output_dir)
        report_manager.create_report(
            genome,
            output_dir,
            experiment_dict,
            report_stats,
            map_args.workspace_dir,
            diffexp_flag,
        )
        sys.exit(0)

    # HTSeq(bacteria), Stringtie(host)
    quantifier = process.Quantify()
    quantifier.set_genome(genome)
    quantifier.set_recipe(map_args.recipe)
    quant_result = quantifier.run_quantification(sample_list, 8, output_dir)
    if quant_result is False:
        sys.stderr.write("Quantification failed: check logs\n")
        sys.exit(1)
    genome_quant_file = quantifier.create_genome_counts_table(output_dir, sample_list)
    if genome_quant_file is None:
        sys.stderr.write("Failed to create genome counts table: check logs\n")
        sys.exit(1)
    genome.add_genome_data("counts_table", genome_quant_file)
    genome_quant_file = quantifier.create_genome_quant_table(output_dir, sample_list)

    # sample_list used in function below
    sample_list = []
    for condition in experiment_dict:
        samples = experiment_dict[condition].get_sample_list()
        sample_list = sample_list + samples

    # Differential expression
    diff_exp = process.DifferentialExpression(comparisons)
    diff_exp.set_recipe(map_args.recipe)
    meta_file = diff_exp.create_metadata_file(sample_list, output_dir)
    genome.add_genome_data("sample_metadata_file", meta_file)
    if diffexp_flag:
        diffexp_import = process.DiffExpImport()
        diffexp_import.set_recipe(map_args.recipe)
        diff_exp.set_genome(genome)
        diffexp_result = diff_exp.run_differential_expression(output_dir, sample_list)
        if diffexp_result is False:
            sys.stderr.write("Warning: differential expression analysis failed, continuing pipeline\n")
        elif genome.get_genome_type() == "bacteria":
            diffexp_import.set_genome(genome)
            import_result = diffexp_import.run_diff_exp_import(output_dir, map_args)
            if import_result is False:
                sys.stderr.write("Warning: differential expression import failed, continuing pipeline\n")

    # Queries: subsystems, kegg
    # output files are used in creating figures
    genome_data = process.GenomeData()
    genome_data.set_recipe(map_args.recipe)
    if genome.get_genome_type() == "bacteria":
        genome_data.set_genome(genome)
        genome_data.run_queries(output_dir, session)
        genome_data.create_system_figures(output_dir)

    # call multiqc without any adjustments
    if not map_args.disable_reports:
        report_manager = report.ReportManager()
        # number of samples and conditions
        report_stats = {}
        sample_count = 0
        condition_count = 0
        for condition in experiment_dict:
            sample_count += len(experiment_dict[condition].get_sample_list())
            if condition != "no_condition":
                condition_count += 1
        report_stats["num_samples"] = sample_count
        report_stats["num_conditions"] = condition_count
        # get recipe
        report_stats["recipe"] = map_args.recipe
        report_manager.run_multiqc(output_dir)
        report_manager.create_report(
            genome,
            output_dir,
            experiment_dict,
            report_stats,
            map_args.workspace_dir,
            diffexp_flag,
        )

    # TODO: Add command output and status
    # File cleanup in perl
    print("done")
    return True


# sets up initial condition, sample, genome folder structure
# folder stucture is: output_dir/condition/sample/genome
# TODO: diffexp object folder
def setup(output_dir, experiment_dict, genome):
    # create output directory
    if not os.path.exists(output_dir):
        print("Creating output directory: {0}".format(output_dir))
        os.mkdir(output_dir)
    else:
        print("Output directory alread exists: {0}".format(output_dir))

    # create subfolder list and add paths to relevant class objects
    subfolder_list = []
    for condition in experiment_dict:
        cond_path = os.path.abspath(os.path.join(output_dir, condition))
        experiment_dict[condition].set_path(cond_path)
        subfolder_list.append(cond_path)
        for sample in experiment_dict[condition].get_sample_list():
            sample_path = os.path.abspath(
                os.path.join(experiment_dict[condition].get_path(), sample.get_id())
            )
            sample.set_path(sample_path)
            subfolder_list.append(sample.get_path())
            genome_path = sample.get_path()
            genome.create_path_entry(sample.get_id(), genome_path)
            subfolder_list.append(genome.get_sample_path(sample.get_id()))

    # create subfolders
    print("Creating {0} subfolders:".format(len(subfolder_list)))
    for folder in subfolder_list:
        if not os.path.exists(folder):
            print("Creating subfolder: {0}".format(folder))
            os.mkdir(folder)
        else:
            print("Subfolder already exists: {0}".format(folder))

    # create report images folder
    report_img_folder = os.path.join(output_dir, "report_images/")
    if not os.path.exists(report_img_folder):
        os.mkdir(report_img_folder)
    genome.add_genome_data("report_img_path", report_img_folder)

    # TODO: genome_data checking, delete genome if not enough data??? Throw errors???
    # setup genome index objects
    """
    genome_links_dir = os.path.join(output_dir,"genome_links")
    if not os.path.exists(genome_links_dir):
        print("Creating {0}".format(genome_links_dir))
        os.mkdir(genome_links_dir)
    else:
        print("{0} already exists".format(genome_links_dir))
    for genome in genome_list:
        genome_data_dir = os.path.join(genome_links_dir,genome.get_id())
        if not os.path.exists(genome_data_dir):
            print("Creating {0}".format(genome_data_dir))
            os.mkdir(genome_data_dir)
        else:
            print("{0} already exists".format(genome_data_dir))
        genome.setup_genome_database(genome_data_dir)
    """
    if not genome.setup_genome_database():
        sys.stderr.write(
            "Genome database setup failed for {0}: exiting\n".format(genome.get_id())
        )
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # TODO: update help
    parser.add_argument(
        "--jfile",
        help='json file for job {"reference_genome_id": "1310806.3", "experimental_conditions":\
                    ["c1_control", "c2_treatment"], "output_file": "rnaseq_baumanii_1505311", \
                    "recipe": "RNA-Rocket", "output_path": "/anwarren@patricbrc.org/home/test",\
                    "paired_end_libs": [{"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R2.fq.gz", "condition": 1},\
                    {"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R2.fq.gz", "condition": 2}], "contrasts": [[1, 2]]}',
        required=True,
    )
    parser.add_argument(
        "-o",
        help="output directory. defaults to current directory.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-g",
        help="csv list of directories each containing all genome data: fna, gff, and hisat indices",
        required=True,
    )
    parser.add_argument(
        "--sstring",
        help='json server string specifying api {"data_api":"url"}',
        required=False,
        default='{"data_api":"https://p3.theseed.org/services/data_api/genome_feature"}',
    )
    parser.add_argument(
        "-p", help="tool parameters", required=False, type=str, default="{}"
    )
    parser.add_argument(
        "-d", help="differential expression folder", required=False, default=".diff_exp"
    )
    # TODO:
    # link to genome files (gff, fa)???
    # link to hisat2 index
    # parameters as input json file OR inline
    # unit testing parameters

    map_args = parser.parse_args()
    print("map_args:\n{0}".format(map_args))

    # set path to use correct version of samtools
    os.environ["PATH"] = (
        "/disks/patric-common/runtime/samtools-1.9/bin:" + os.environ["PATH"]
    )
    os.environ["R_LIBS"] = (
        "/disks/patric-common/runtime/lib/R/library:" + os.environ.get("R_LIBS", "")
    )
    print("R_LIBS path = {0}".format(os.environ["R_LIBS"]))

    # output directory
    output_dir = map_args.o
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # load library
    job_data = None
    try:
        with open(map_args.jfile, "r") as job_handle:
            job_data = json.load(job_handle)
    except Exception as e:
        sys.stderr.write("Error in opening job json file:\n{0}\n".format(e))
        sys.exit(1)
    if not job_data:
        sys.stderr.write("job_data is null\n")
        sys.exit(1)

    # Validate job data
    validate_job_data(job_data)

    # Setup session
    s = requests.Session()
    try:
        authenticateByEnv(s)
        print("authentication success")
    except Exception as e:
        sys.stderr.write("Error during authentication, exiting:\n{0}\n".format(e))
        sys.exit(1)

    # load genome ids
    # genome_list = []
    # for i in range(0,len(job_data["reference_genome_id"])):
    # genome_list.append(experiment.Genome(job_data["reference_genome_id"][i],job_data["genome_type"][i],s))
    # genome_list.append(experiment.Genome(job_data["reference_genome_id"],job_data["genome_type"],s))
    genome = experiment.Genome(
        job_data["reference_genome_id"], job_data["genome_type"], s
    )

    # Load genome data
    genome_dir = map_args.g.strip()
    if genome_dir.endswith("/"):
        genome_key = os.path.basename(os.path.dirname(genome_dir))
    else:
        genome_key = os.path.basename(genome_dir)
    if genome_key == genome.get_id():
        print(os.listdir(genome_dir))
        for f in os.listdir(genome_dir):
            data_key = None
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                data_key = "fasta"
            elif f.endswith(".gff"):  # TODO: other annotation types?
                data_key = "annotation"
            elif f.endswith(".ht2.tar"):
                data_key = "hisat_index"
            else:
                continue
            genome.set_genome_dir(genome_dir)
            genome.add_genome_data(
                data_key, os.path.abspath(os.path.join(genome_dir, f))
            )

    # sample_list = [] # maybe don't store this, access samples by condition like in original
    experiment_dict = {}
    condition_list = []
    for cond_str in job_data.get("experimental_conditions", []):
        cond_str = cond_str.replace(" ", "_")
        condition_list.append(cond_str)
        new_condition = experiment.Condition(cond_str)
        experiment_dict[cond_str] = new_condition

    # paired_end_libs
    if "paired_end_libs" in job_data:
        for paired_sample in job_data["paired_end_libs"]:
            if "condition" in paired_sample:
                condition = paired_sample["condition"]
            else:
                if "no_condition" not in experiment_dict:
                    no_condition = experiment.Condition("no_condition")
                    experiment_dict["no_condition"] = no_condition
                condition = "no_condition"
            condition = condition.replace(" ", "_")
            sample_reads = [paired_sample["read1"]]
            if "read2" in paired_sample:
                sample_reads.append(paired_sample["read2"])
                sample_type = "paired"
            else:
                sample_type = "single"
            sample_id = paired_sample["sample_id"].replace(" ", "_")
            new_sample = experiment.Sample(
                sample_id, sample_type, sample_reads, None, condition
            )
            if condition:
                experiment_dict[condition].add_sample(new_sample)

    # single_end_libs
    if "single_end_libs" in job_data:
        for single_sample in job_data["single_end_libs"]:
            if "condition" in single_sample:
                condition = single_sample["condition"]
            else:
                if "no_condition" not in experiment_dict:
                    no_condition = experiment.Condition("no_condition")
                    experiment_dict["no_condition"] = no_condition
                condition = "no_condition"
            condition = condition.replace(" ", "_")
            sample_read = [single_sample["read"]]
            sample_id = single_sample["sample_id"].replace(" ", "_")
            new_sample = experiment.Sample(
                sample_id, "single", sample_read, None, condition
            )
            if condition:
                experiment_dict[condition].add_sample(new_sample)

    # put sra-fastq files in <output_dir>/SRA_Fastq/
    # put sra-metadata files in <output_dir>/SRA_Meta/
    # set type to paired or single here
    if "srr_libs" in job_data:
        sra_fastq_dir = os.path.join(output_dir, "SRA_Fastq")
        sra_meta_dir = os.path.join(output_dir, "SRA_Meta")
        if not os.path.exists(sra_fastq_dir):
            os.mkdir(sra_fastq_dir)
        if not os.path.exists(sra_meta_dir):
            os.mkdir(sra_meta_dir)
        for sra_sample in job_data["srr_libs"]:
            if "condition" in sra_sample:
                condition = sra_sample["condition"]
            else:
                if "no_condition" not in experiment_dict:
                    no_condition = experiment.Condition("no_condition")
                    experiment_dict["no_condition"] = no_condition
                condition = "no_condition"
            condition = condition.replace(" ", "_")
            srr_id = sra_sample["srr_accession"]
            meta_file = os.path.join(sra_meta_dir, srr_id + "_meta.txt")
            reads_dir = {}
            try:
                subprocess.check_call(
                    [
                        "p3-sra",
                        "--out",
                        sra_fastq_dir,
                        "--metadata-file",
                        meta_file,
                        "--id",
                        srr_id,
                    ]
                )
                with open(meta_file, "r") as meta_handle:
                    job_meta = json.load(meta_handle)
                    files = job_meta[0].get("files", [])
                    for i, f in enumerate(files):
                        if f.endswith("_2.fastq"):
                            read2 = os.path.join(sra_fastq_dir, f)
                            reads_dir["read2"] = read2
                        elif f.endswith("_1.fastq"):
                            read1 = os.path.join(sra_fastq_dir, f)
                            reads_dir["read1"] = read1
                        elif f.endswith(".fastq"):
                            read = os.path.join(sra_fastq_dir, f)
                            reads_dir["read"] = read
                missing_read = False
                for read_key in reads_dir:
                    if not os.path.exists(reads_dir[read_key]):
                        sys.stderr.write(
                            "Error: fastq file doesn't exist:\n{0}\n".format(
                                reads_dir[read_key]
                            )
                        )
                        missing_read = True
                        continue
                # TODO: exit or something?
                if missing_read:
                    continue
                if "read2" in reads_dir:
                    reads_list = [reads_dir["read1"], reads_dir["read2"]]
                    new_sample = experiment.Sample(
                        srr_id, "paired", reads_list, None, condition
                    )
                    if condition:
                        experiment_dict[condition].add_sample(new_sample)
                else:  # single end
                    reads_list = [reads_dir["read"]]
                    new_sample = experiment.Sample(
                        srr_id, "single", reads_list, None, condition
                    )
                    if condition:
                        experiment_dict[condition].add_sample(new_sample)

            except Exception as e:
                sys.stderr.write(
                    "Error in downloading SRA ID {0}:\n{1}\n".format(srr_id, e)
                )
                # TODO: continue or exit??
                sys.exit(-1)

    tool_params = json.loads(map_args.p)
    print("tool_params = {0}".format(tool_params))

    comparisons = experiment.Comparison()
    for con in job_data.get("contrasts", []):
        con = [x.replace(" ", "_") for x in con]
        comparisons.add_contrast(con[0], con[1])

    # change into output directory
    os.chdir(output_dir)

    # set recipe in map_args
    map_args.recipe = job_data["recipe"]

    # set trimming in map_args
    if "trimming" in job_data:
        map_args.trimming = job_data["trimming"]
    else:
        map_args.trimming = True

    # set report variable
    if "disable_reports" in job_data:
        map_args.disable_reports = job_data["disable_reports"]
    else:
        map_args.disable_reports = False

    # workspace dir for links
    map_args.workspace_dir = job_data["output_path"]

    # If not cufflinks, run pipeline
    main(genome, experiment_dict, tool_params, output_dir, comparisons, s, map_args)
