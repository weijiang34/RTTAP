import os

import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
import time
import glob

import envs
from .api_qc import run_fastp
from .api_host_removal import remove_host
from .api_read_recruitment import run_kraken2
from .api_read_extraction import split_reads
from .api_taxonomy_bacteria import run_metaphlan4, generate_mpa_report
from .api_taxonomy_bowtie2 import (
    run_mapping_and_extract, 
    append_taxid, 
    convert_taxid_to_path, 
    run_lca, 
    merge_and_cleanup, 
    generate_bt2_report
)
from .api_downstream import run_rgi, run_virstrain, generate_overall_report


# def timing_step(step_name: str):
#     def decorator(func):
#         def wrapper(self, *args, **kwargs):
#             logging.info(f"Starting {step_name} ...")
#             start = time.time()
#             try:
#                 result = func(self, *args, **kwargs)
#             except Exception as e:
#                 logging.error(f"Error in {step_name}: {e}")
#                 raise
#             elapsed = time.time() - start
#             logging.info(f"{step_name} finished in {elapsed:.2f} seconds.")
#             return result
#         return wrapper
#     return decorator

class Step:
    def __init__(
        self, 
        step_name, 
        func, 
        check_files_dependencies: list = [],
        check_files_completed: list = [],
        force: bool = True,
        *args, 
        **kwargs
    ):
        self.step_name = step_name
        self.func = func
        self.check_files_dependencies = check_files_dependencies
        self.check_files_completed = check_files_completed
        self.force = force
        self.completed = False

        self.args = args
        self.kwargs = kwargs

    def check_step_dependencies(self):
        missing = False
        for f in self.check_files_dependencies:
            if not os.path.exists(f):
                logging.error(f"\tMissing dependency file: {f}")
                missing = True
        if missing:
            raise FileNotFoundError(f"Missing dependency files found.")
        return missing

    def check_step_completed(self, log=False):
        files_completed = {}
        for f in self.check_files_completed:
            if not os.path.exists(f):
                if log is True:
                    logging.warning(f"\tMissing file: {f}")
                files_completed[f] = False
            else:
                files_completed[f] = True
        self.completed = all(files_completed.values())
        return self.completed

    def timing_step(self, func, *args, **kwargs):
        info = f"\tStarting {self.step_name.lower()} ..."
        logging.info(info)
        start = time.time()
        try:
            result = func(*args, **kwargs)
        except Exception as e:
            logging.error(f"\tError in {self.step_name.lower()}: {e}")
            raise
        elapsed = time.time() - start
        logging.info(f"\t{self.step_name.capitalize()} finished in {elapsed:.2f} seconds.")
        return result

    def run(self):
        # pre-run - check completed
        self.check_step_completed(log=False)
        # main process
        if not self.force and self.completed:
            logging.info(f"\t{self.step_name.lower()} already completed. Skipping...")
            return
        if self.force:
            logging.info(f"\tForce re-running {self.step_name.lower()}.")
        self.timing_step(self.func, *self.args, **self.kwargs)
        # post-run check completed
        self.check_step_completed(log=True)

class Workflow:
    def __init__(
        self, 
        wf_name: str, 
        steps: list = [],
        force=True, 
        completed=False
    ):
        self.wf_name = wf_name
        self.steps = steps
        self.force = force
        self.completed = completed

    def check_completed(self):
        steps_completed = [step.check_step_completed() for step in self.steps]
        self.completed = all(steps_completed)
        return self.completed

    def timing_workflow(self, func, *args, **kwargs):
        info = f"Starting {self.wf_name.lower()} ..."
        logging.info(info)
        start = time.time()
        try:
            result = func(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error in {self.wf_name.lower()}: {e}")
            raise
        elapsed = time.time() - start
        logging.info(f"{self.wf_name.capitalize()} finished in {elapsed:.2f} seconds.")
        return result
    
    def run_steps(self):
        for step in self.steps:
            step.run()

    def run(self):
        # pre-run check
        self.check_completed()
        # main process
        # if skip
        if not self.force and self.completed:
            logging.info(f"{self.wf_name.capitalize()} already completed. Skipping...")
            return
        # if force
        if self.force:
            logging.info(f"Force re-running {self.wf_name.lower()}.")
        self.timing_workflow(self.run_steps)


class QualityControl(Workflow):
    """
    Quality control workflow using fastp, with stepwise timing and safety checks.
    """
    def __init__(
        self, 
        input_1, 
        out_dir, 
        fileHeader, 
        tool_path, 
        conda_dir, 
        env_name,
        threads=1, 
        min_len=15, 
        input_2=None, 
        cmd='',
        force=True
    ):
        # Parameter validation
        assert isinstance(input_1, str) and os.path.exists(input_1), f"Input file {input_1} not found."
        if input_2 is not None:
            assert isinstance(input_2, str) and os.path.exists(input_2), f"Input file {input_2} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(fileHeader, str), "fileHeader must be str"

        assert isinstance(tool_path, str) and os.path.exists(tool_path), f"fastp path {tool_path} not found."
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"

        assert isinstance(threads, int) and threads > 0, "threads must be positive int"
        assert isinstance(min_len, int) and min_len > 0, "min_len must be positive int"
        assert isinstance(force, bool), "force must be bool"

        assert isinstance(cmd, str), "cmd must be str"

        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"

        super().__init__(
            wf_name="Quality control by Fastp",
            force=force, 
            completed=False
        )

        # FUNC PARAMETERS START
        self.input_1 = input_1
        self.input_2 = input_2
        self.out_dir = out_dir
        self.fileHeader = fileHeader

        self.tool_path = tool_path
        self.conda_dir = conda_dir
        self.env_name = env_name

        self.threads = threads
        self.min_len = min_len
        self.cmd = cmd
        # FUNC PARAMETERS END

        # check files - completed 
        check_files_completed = [os.path.join(out_dir, f"{fileHeader}.json"), os.path.join(out_dir, f"{fileHeader}.html")]
        if input_2 is None:
            check_files_completed.append(os.path.join(out_dir, f"{fileHeader}.clean.fq.gz"))
        else:
            check_files_completed.append(os.path.join(out_dir, f"{fileHeader}.clean.R1.fq.gz"))
            check_files_completed.append(os.path.join(out_dir, f"{fileHeader}.clean.R2.fq.gz"))

        self.steps = [
            Step(
                step_name="Fastp quality control", 
                func=run_fastp, 
                force=self.force,
                # FUNC PARAMETERS START
                input_1=self.input_1, 
                input_2=self.input_2, 
                out_dir=self.out_dir, 
                fileHeader=self.fileHeader, 

                tool_path=self.tool_path, 
                conda_dir=self.conda_dir,
                env_name=self.env_name,

                threads=self.threads, 
                min_len=self.min_len, 
                cmd=self.cmd, 
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=check_files_completed
            )
        ]

        os.makedirs(out_dir, exist_ok=True)


class HostRemoval(Workflow):
    def __init__(
        self, 
        input_1,
        out_dir,
        tool_path,
        conda_dir,
        env_name,
        db_dir,
        fileHeader,
        threads=1,
        input_2=None,
        cmd='',
        force=True
    ):
        # Parameter validation
        assert isinstance(input_1, str) and os.path.exists(input_1), f"Input file {input_1} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(tool_path, str) and os.path.exists(tool_path), f"Tool path {tool_path} not found."
        # check db_dir contains .bt2 or .bt2l files
        bt2_files = glob.glob(os.path.join(db_dir, "*.bt2"))
        bt2l_files = glob.glob(os.path.join(db_dir, "*.bt2l"))
        assert isinstance(db_dir, str) and (bt2_files, bt2l_files), f"Bowtie2 db path {db_dir} not found."
        assert isinstance(fileHeader, str), "fileHeader must be str"
        assert isinstance(threads, int) and threads > 0, "threads must be positive int"
        if input_2 is not None:
            assert isinstance(input_2, str) and os.path.exists(input_2), f"Input file {input_2} not found."
        assert isinstance(cmd, str), "cmd must be str"
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"

        super().__init__(
            wf_name="Host read removal", 
            force=force, 
            completed=False
        )

        # FUNC PARAMETERS START
        self.input_1 = input_1
        self.input_2 = input_2
        self.out_dir = out_dir
        self.fileHeader = fileHeader

        self.tool_path = tool_path
        self.conda_dir = conda_dir
        self.env_name = env_name
        self.db_dir = db_dir
        self.threads = threads
        
        self.cmd = cmd
        # FUNC PARAMETERS END

        os.makedirs(out_dir, exist_ok=True)

        # check files - completed
        check_files_completed = [os.path.join(out_dir, f"{fileHeader}.nohost.fq.gz")]
        self.steps = [
            Step(
                step_name="Host removal by Bowtie2", 
                func=remove_host, 
                force=self.force,
                # FUNC PARAMETERS START
                input_1=self.input_1, 
                input_2=self.input_2, 
                out_dir=self.out_dir, 
                fileHeader=self.fileHeader, 

                tool_path=self.tool_path, 
                db_dir=self.db_dir, 
                conda_dir=self.conda_dir,
                env_name=self.env_name,
                
                threads=self.threads, 
                
                cmd=self.cmd,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=check_files_completed
            )
        ]

        

class ReadRecruitment(Workflow):
    """
    Recruit reads workflow using Kraken2, with stepwise timing and safety checks.
    """
    def __init__(
        self, 
        input: str, 
        out_dir: str,
        fileHeader: str,  

        tool_path: str, 
        conda_dir: str, 
        env_name: str, 
        db_path: str, 

        threads: int = 1, 
        cmd: str = '',
        force: bool = True
    ):
        # Parameter validation
        assert isinstance(input, str) and os.path.exists(input), f"Input file {input} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(tool_path, str) and os.path.exists(tool_path), f"Kraken2 path {tool_path} not found."
        assert isinstance(fileHeader, str), "fileHeader must be str"
        assert isinstance(db_path, str) and os.path.exists(db_path), f"Kraken2 db path {db_path} not found."
        assert isinstance(threads, int) and threads > 0, "threads must be positive int"
        assert isinstance(cmd, str), "cmd must be str"
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"

        super().__init__(
            wf_name="Read recruitment by Kraken2", 
            force=force, 
            completed=False
        )
        # FUNC PARAMETERS START
        self.input = input
        self.out_dir = out_dir
        self.tool_path = tool_path
        self.conda_dir = conda_dir
        self.env_name = env_name
        self.fileHeader = fileHeader
        self.db_path = db_path
        self.threads = threads
        self.cmd = cmd
        # FUNC PARAMETERS END

        # check files - completed
        check_files_completed = [
            os.path.join(out_dir, f"{fileHeader}.kraken2.out"), 
            os.path.join(out_dir, f"{fileHeader}.kraken2.report")
        ]
        self.steps = [
            Step(
                step_name="Kraken2 read recruitment", 
                func=run_kraken2, 
                force=self.force,
                # FUNC PARAMETERS START
                input=self.input, 
                out_dir=self.out_dir, 
                tool_path=self.tool_path, 
                conda_dir=self.conda_dir,
                env_name=self.env_name,
                fileHeader=self.fileHeader, 
                db_path=self.db_path, 
                threads=self.threads, 
                cmd=self.cmd, 
                # FUNC PARAMETERS END

                check_files_dependencies=[],
                check_files_completed=check_files_completed
            )
        ]

        os.makedirs(out_dir, exist_ok=True)


class ReadExtractor(Workflow):
    def __init__(
        self, 
        kraken_out,
        kraken_report, 
        nohost_fq_gz, 

        out_dir, 
        fileHeader, 

        conda_dir, 
        env_name,
        split_script_path,

        reads_type, 
        threads=1,
        force=True
    ):
        
        self.taxid_dict = {
            "bacteria": 2,
            "viruses": 10239,
            "fungi": 4751
        }
        # Parameter validation
        assert isinstance(kraken_out, str) and os.path.exists(kraken_out), f"Kraken2 output file {kraken_out} not found."
        assert isinstance(kraken_report, str) and os.path.exists(kraken_report), f"Kraken2 report file {kraken_report} not found."
        assert isinstance(nohost_fq_gz, str) and os.path.exists(nohost_fq_gz), f"Input file {nohost_fq_gz} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(fileHeader, str), "fileHeader must be str"
        assert isinstance(reads_type, str) and reads_type.lower() in self.taxid_dict.keys(), f"Invalid reads_type: {self.reads_type}. Must be one of {list(self.taxid_dict.keys())}."
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"
        assert isinstance(split_script_path, str) and os.path.exists(split_script_path), f"Split script path {split_script_path} not found."

        super().__init__(
            wf_name="Read extraction", 
            force=force, 
            completed=False
        )
        # FUNC PARAMETERS START
        self.kraken_out = kraken_out
        self.kraken_report = kraken_report
        self.nohost_fq_gz = nohost_fq_gz
        self.out_dir = out_dir
        self.fileHeader = fileHeader
        self.reads_type = reads_type

        self.conda_dir = conda_dir
        self.env_name = env_name
        self.split_script_path = split_script_path
        # FUNC PARAMETERS END
        
        os.makedirs(out_dir, exist_ok=True)

        # check files - completed
        if self.reads_type.lower() == "bacteria":
            check_files_completed = [os.path.join(out_dir, f"{fileHeader}.bacteria.fq.gz")]
        elif self.reads_type.lower() == "viruses":
            check_files_completed = [os.path.join(out_dir, f"{fileHeader}.viruses.fq.gz")]
        elif self.reads_type.lower() == "fungi":    
            check_files_completed = [os.path.join(out_dir, f"{fileHeader}.fungi.fq.gz")]

        self.steps = [
            Step(
                step_name=f"Split {self.reads_type} reads from Kraken2 output", 
                func=split_reads, 
                force=self.force,
                # FUNC PARAMETERS START
                kraken_out=self.kraken_out,
                kraken_report=self.kraken_report,
                nohost_fq_gz=self.nohost_fq_gz,
                out_dir=self.out_dir,
                fileHeader=self.fileHeader,
                reads_type=self.reads_type,
                conda_dir=self.conda_dir,
                env_name=self.env_name,
                split_script_path=self.split_script_path,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=check_files_completed
            )
        ]

class TaxonomyBacteriaWorkflow(Workflow):
    """
    Taxonomic classification workflow using MetaPhlAn4, with stepwise timing and safety checks.
    """
    def __init__(
        self, 
        input, 
        k2_report_path, 
        mpa_report_path,
        bac_report_path,
        out_dir, 
        fileHeader, 

        tool_path, 
        metaphlan4_db_path, 
        conda_dir,
        env_name,

        threads=1, 
        RPM_threshold=1.0,
        cmd='',
        force=True
    ):
        # Parameter validation
        assert isinstance(input, str) and os.path.exists(input), f"Input file {input} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(tool_path, str) and os.path.exists(tool_path), f"MetaPhlAn4 path {tool_path} not found."
        assert isinstance(fileHeader, str), "fileHeader must be str"
        # assert isinstance(metaphlan4_db_path, str) and os.path.exists(metaphlan4_db_path), f"MetaPhlAn4 db path {metaphlan4_db_path} not found."
        assert isinstance(threads, int) and threads > 0, "threads must be positive int"
        assert isinstance(RPM_threshold, (int, float)) and RPM_threshold >= 0, "RPM_threshold must be non-negative number"
        assert isinstance(k2_report_path, str) and os.path.exists(k2_report_path), f"Kraken2 report path {k2_report_path} not found."
        assert isinstance(mpa_report_path, str), "mpa_report_path must be str"
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"
        assert isinstance(cmd, str), "cmd must be str"

        super().__init__(
            wf_name="Bacterial taxonomic classification by MetaPhlAn4", 
            force=force, 
            completed=False
        )
        # FUNC PARAMETERS START
        self.input = input
        self.k2_report_path = k2_report_path
        self.mpa_report_path = mpa_report_path
        self.bac_report_path = bac_report_path
        self.out_dir = out_dir
        self.fileHeader = fileHeader

        self.tool_path = tool_path
        self.metaphlan4_db_path = metaphlan4_db_path
        self.conda_dir = conda_dir
        self.env_name = env_name

        self.threads = threads
        self.RPM_threshold = RPM_threshold
        self.cmd = cmd
        self.force = force
        # FUNC PARAMETERS END

        os.makedirs(out_dir, exist_ok=True)

        self.steps = [
            Step(
                step_name="MetaPhlAn4 taxonomic classification", 
                func=run_metaphlan4, 
                force=self.force,
                # FUNC PARAMETERS START
                input=self.input,
                out_dir=self.out_dir,
                fileHeader=self.fileHeader,

                tool_path=self.tool_path,
                metaphlan4_db_path=self.metaphlan4_db_path,
                
                threads=self.threads,
                cmd=self.cmd,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.metaphlan4.bowtie2.out"),
                    os.path.join(out_dir, f"{fileHeader}.metaphlan4.report")
                ]
            ),
            Step(
                step_name="Generate bacteria report", 
                func=generate_mpa_report, 
                force=self.force,
                # FUNC PARAMETERS START
                # out_dir=self.out_dir,
                # fileHeader=self.fileHeader,
                k2_report_path=self.k2_report_path,
                bac_report_path=self.bac_report_path,
                mpa_report_path=self.mpa_report_path,
                RPM_threshold=self.RPM_threshold,
                # FUNC PARAMETERS END
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.bacteria.report")
                ]
            )
        ]


class TaxonomyBowtie2Workflow(Workflow):
    """
    Taxonomic mapping workflow using Bowtie2 and TaxonKit, with stepwise timing and safety checks.
    """
    def __init__(
            self,
            input, 
            k2_report_path, 
            lca_report_path,
            bt_report_path, 
            out_dir, 
            fileHeader, 

            bt_path, 
            taxonkit_path, 
            db_type, 
            acc2taxid_dir, 
            conda_dir, 
            env_name, 

            threads=1,
            RPM_threshold=1.0,
            cmd='',

            force=True
    ):
        # Parameter validation
        assert isinstance(input, str) and os.path.exists(input), f"Input file {input} not found."
        assert isinstance(k2_report_path, str) and os.path.exists(k2_report_path), f"Kraken2 report path {k2_report_path} not found."
        assert isinstance(lca_report_path, str), "lca_report_path must be str"
        assert isinstance(bt_report_path, str), "bt_report_path must be str"
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(bt_path, str) and os.path.exists(bt_path), f"Bowtie2 path {bt_path} not found."
        assert isinstance(taxonkit_path, str) and os.path.exists(taxonkit_path), f"TaxonKit path {taxonkit_path} not found."
        assert isinstance(fileHeader, str), "fileHeader must be str"
        assert db_type in ["viruses", "fungi"], f"db_type must be 'viruses' or 'fungi', got {db_type}"
        assert isinstance(acc2taxid_dir, str) and os.path.exists(acc2taxid_dir), f"acc2taxid_dir {acc2taxid_dir} not found."
        assert isinstance(threads, int) and threads > 0, "threads must be positive int"
        assert isinstance(RPM_threshold, (int, float)) and RPM_threshold >= 0, "RPM_threshold must be non-negative number"
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"
        assert isinstance(cmd, str), "cmd must be str"
        
        super().__init__(
            wf_name=f"{db_type.capitalize()} taxonomic classification by Bowtie2", 
            force=force, 
            completed=False
        )

        # FUNC PARAMETERS START
        self.input = input
        self.k2_report_path = k2_report_path
        self.lca_out_path = lca_report_path
        self.bt_report_path = bt_report_path
        self.out_dir = out_dir
        self.fileHeader = fileHeader

        self.bt_path = bt_path
        self.taxonkit_path = taxonkit_path
        self.acc2taxid_dir = acc2taxid_dir
        self.db_path_dict = {
            "viruses": envs.BT_VIRAL_DB_PATH,
            "fungi": envs.BT_FUNGAL_DB_PATH
        }
        self.acc2taxid_dict = {
            "viruses": os.path.join(acc2taxid_dir, "viruses.acc2taxid.txt"),
            "fungi": os.path.join(acc2taxid_dir, "EuPathDB46.seq2taxid.txt"),
        }
        if db_type not in self.db_path_dict:
            raise ValueError(f"Invalid db_type: {db_type}. Must be one of {list(self.db_path_dict.keys())}.")
        self.db_path = self.db_path_dict[db_type]
        self.acc2taxid_path = self.acc2taxid_dict[db_type]
        self.conda_dir = conda_dir
        self.env_name = env_name
        
        self.db_type = db_type
        self.threads = threads
        self.RPM_threshold = RPM_threshold
        self.cmd = cmd
        # FUNC PARAMETERS END

        self.steps = [
            Step(
                step_name="Bowtie2 mapping and extract uniq/multi mapping results", 
                func=run_mapping_and_extract, 
                force=self.force,
                # FUNC PARAMETERS START
                input=self.input, 
                out_dir=self.out_dir, 
                fileHeader=self.fileHeader, 

                bt_path=self.bt_path, 
                db_type=self.db_type, 
                db_path=self.db_path, 
                conda_dir=self.conda_dir,
                env_name=self.env_name,
                
                threads=self.threads,
                cmd=self.cmd,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.tmp"),
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.tmp")
                ]
            ),
            Step(
                step_name="Append taxid to uniq/multi mapping results", 
                func=append_taxid, 
                force=self.force,
                # FUNC PARAMETERS START
                out_dir=self.out_dir, 
                db_type=self.db_type, 
                acc2taxid_path=self.acc2taxid_path, 
                fileHeader=self.fileHeader,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.reads2acc2taxid.tmp"),
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.reads2acc2taxid.tmp")
                ]
            ),
            Step(
                step_name="Convert taxid to lineage path using TaxonKit", 
                func=convert_taxid_to_path, 
                force=self.force,
                # FUNC PARAMETERS START
                out_dir=self.out_dir, 
                db_type=self.db_type, 
                taxonkit_path=self.taxonkit_path, 
                fileHeader=self.fileHeader,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.lineage.tmp"),
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.lineage.tmp")
                ]
            ),
            Step(
                step_name="Calculate LCA for multi-mapping results", 
                func=run_lca, 
                force=self.force,
                # FUNC PARAMETERS START
                out_dir=self.out_dir, 
                db_type=self.db_type, 
                fileHeader=self.fileHeader,
                # FUNC PARAMETERS END
                check_files_dependencies=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.lineage.tmp")
                ],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.LCA.tmp")
                ]
            ),
            Step(
                step_name="Merge uniq/multi results and cleanup temporary files", 
                func=merge_and_cleanup, 
                force=self.force,
                # FUNC PARAMETERS START
                out_dir=self.out_dir, 
                db_type=self.db_type, 
                fileHeader=self.fileHeader,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.LCA.out")
                ]
            ),
            Step(
                step_name=f"Generate {self.db_type} taxonomic report", 
                func=generate_bt2_report, 
                force=self.force,
                # FUNC PARAMETERS START
                # out_dir=self.out_dir, 
                # fileHeader=self.fileHeader, 
                # db_type=self.db_type, 
                lca_out_path=self.lca_out_path,
                k2_report_path=self.k2_report_path, 
                RPM_threshold=self.RPM_threshold,
                bt_report_path=self.bt_report_path,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(out_dir, f"{fileHeader}.{db_type}.report")
                ]
            )
        ]

        os.makedirs(out_dir, exist_ok=True)


class ARGIdentification(Workflow):
    def __init__(
        self, 
        input_fq, 
        out_dir, 
        fileHeader, 

        tool_path, 
        db_path, 
        conda_dir,
        env_name,

        threads=1,
        cmd='',
        force=True
    ):
        # Parameter validation
        assert isinstance(input_fq, str) and os.path.exists(input_fq), f"Input file {input_fq} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(tool_path, str) and os.path.exists(tool_path), f"Tool path {tool_path} not found."
        assert isinstance(db_path, str) and os.path.exists(db_path), f"Database path {db_path} not found."
        assert isinstance(fileHeader, str), "fileHeader must be str"
        assert isinstance(threads, int) and threads > 0, "threads must be positive int"

        super().__init__(
            wf_name="Antibiotic resistance gene identification by RGI", 
            force=force, 
            completed=False
        )
        # FUNC PARAMETERS START
        self.input_fq = input_fq
        self.out_dir = out_dir
        self.fileHeader = fileHeader

        self.tool_path = tool_path
        self.db_path = db_path
        self.conda_dir = conda_dir
        self.env_name = env_name
        
        self.threads = threads
        # FUNC PARAMETERS END

        self.steps = [
            Step(
                step_name="RGI antibiotic resistance gene identification", 
                func=run_rgi, 
                force=self.force,
                # FUNC PARAMETERS START
                input_fq=self.input_fq,
                out_dir=self.out_dir,
                fileHeader=self.fileHeader,

                tool_path=self.tool_path,
                db_path=self.db_path,
                conda_dir=self.conda_dir,
                env_name=self.env_name,

                threads=self.threads,
                cmd=cmd,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(self.out_dir, f"{self.fileHeader}.rgi.out.overall_mapping_stats.txt")
                ]
            )
        ]

        os.makedirs(self.out_dir, exist_ok=True)

class VirusStrainProfiling(Workflow):
    def __init__(
        self, 
        input_fq_gz, 
        out_dir, 
        fileHeader, 

        virus_report_path, 
        tool_path, 
        db_path, 
        virstrain_list_path,
        conda_dir,
        env_name,

        cmd ='',
        threads=1,

        force=True
    ):
        # Parameter validation
        assert isinstance(input_fq_gz, str) and os.path.exists(input_fq_gz), f"Input file {input_fq_gz} not found."
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(fileHeader, str), "fileHeader must be str"
        assert isinstance(virus_report_path, str) and os.path.exists(virus_report_path), f"Report file {virus_report_path} not found."
        assert isinstance(tool_path, str) and os.path.exists(tool_path), f"Tool path {tool_path} not found."
        assert isinstance(db_path, str) and os.path.exists(db_path), f"Database path {db_path} not found."
        assert isinstance(virstrain_list_path, str) and os.path.exists(virstrain_list_path), f"VirStrain list file {virstrain_list_path} not found."
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(conda_dir, str) and os.path.exists(conda_dir), f"Conda path {conda_dir} not found."
        assert isinstance(env_name, str) and len(env_name)>0, "env_name must be str"
        assert isinstance(cmd, str), "cmd must be str"

        super().__init__(
            wf_name="Virus strain profiling by VirStrain", 
            force=force, 
            completed=False
        )

        # FUNC PARAMETERS START
        self.input_fq_gz = input_fq_gz
        self.out_dir = out_dir
        self.fileHeader = fileHeader

        self.virus_report_path = virus_report_path
        self.tool_path = tool_path
        self.db_path = db_path
        self.virstrain_list_path = virstrain_list_path
        self.conda_dir = conda_dir
        self.env_name = env_name
        self.threads = threads

        self.cmd = cmd
        # FUNC PARAMETERS END

        self.steps = [
            Step(
                step_name="Virus strain profiling by VirStrain", 
                func=run_virstrain, 
                force=self.force,
                # FUNC PARAMETERS START
                input_fq_gz=self.input_fq_gz,
                out_dir=self.out_dir,
                fileHeader=self.fileHeader,

                virus_report_path=self.virus_report_path,
                tool_path=self.tool_path,
                db_path=self.db_path,
                virstrain_list_path=self.virstrain_list_path,
                conda_dir=self.conda_dir,
                env_name=self.env_name,

                threads=self.threads,
                cmd=self.cmd,
                # FUNC PARAMETERS END
                check_files_dependencies=[],
                check_files_completed=[
                    os.path.join(self.out_dir, f"{self.fileHeader}.virstrain.report")
                ]
            )
        ]

        os.makedirs(out_dir, exist_ok=True)

class GenerateReport(Workflow):
    def __init__(
        self, 
        out_dir, 
        # FUNC PARAMETERS START
        k2_report_path, 
        mpa_report_path,
        vir_lca_out_path,
        fun_lca_out_path,
        bac_report_path, 
        vir_report_path, 
        fun_report_path,
        overall_report_path,
        bac_rpm=10.0,
        vir_rpm=1.0,
        fun_rpm=10.0,
        # FUNC PARAMETERS END
        bac=False,
        vir=False,
        fun=False,
        overall=True,
        all=False,

        force=False,
        threads=1,
    ):
        assert isinstance(out_dir, str), "out_dir must be str"
        assert isinstance(k2_report_path, str) and os.path.exists(k2_report_path), f"Kraken2 report {k2_report_path} not found."
        assert isinstance(mpa_report_path, str) and os.path.exists(mpa_report_path), f"MetaPhlAn4 report {mpa_report_path} not found."
        assert isinstance(bac_report_path, str), "bac_report_path must be str"
        assert isinstance(vir_report_path, str), "vir_report_path must be str"
        assert isinstance(fun_report_path, str), "fun_report_path must be str"
        assert isinstance(overall_report_path, str), "overall_report_path must be str"
        assert isinstance(bac_rpm, (int, float)) and bac_rpm >= 0, "bac_rpm must be non-negative number"
        assert isinstance(vir_rpm, (int, float)) and vir_rpm >= 0, "vir_rpm must be non-negative number"
        assert isinstance(fun_rpm, (int, float)) and fun_rpm >= 0, "fun_rpm must be non-negative number"
        assert isinstance(bac, bool), "bac must be bool"
        assert isinstance(vir, bool), "vir must be bool"
        assert isinstance(fun, bool), "fun must be bool"
        assert isinstance(overall, bool), "overall must be bool"
        assert isinstance(all, bool), "all must be bool"
        assert isinstance(force, bool), "force must be bool"
        assert isinstance(threads, int) and threads > 0, "threads must be positive int"

        super().__init__(
            wf_name="Generate Taxonomic Report", 
            force=force, 
            completed=False
        )

        self.out_dir = out_dir
        self.k2_report_path = k2_report_path
        self.mpa_report_path = mpa_report_path
        self.vir_lca_out_path = vir_lca_out_path
        self.fun_lca_out_path = fun_lca_out_path
        self.bac_report_path = bac_report_path
        self.vir_report_path = vir_report_path
        self.fun_report_path = fun_report_path
        self.overall_report_path = overall_report_path
        self.bac_rpm = bac_rpm
        self.vir_rpm = vir_rpm
        self.fun_rpm = fun_rpm
        self.bac = bac
        self.vir = vir
        self.fun = fun
        self.overall = overall
        self.all = all

        os.makedirs(out_dir, exist_ok=True)

        step_report_bac = Step(
            step_name="Generate bacteria report", 
            func=generate_mpa_report, 
            force=self.force,
            # FUNC PARAMETERS START
            k2_report_path=self.k2_report_path,
            bac_report_path=self.bac_report_path,
            mpa_report_path=self.mpa_report_path,
            RPM_threshold=self.bac_rpm,
            # FUNC PARAMETERS END
            check_files_completed=[
                self.bac_report_path
            ]
        )
        step_report_vir = Step(
            step_name=f"Generate viruses taxonomic report", 
            func=generate_bt2_report, 
            force=self.force,
            # FUNC PARAMETERS START
            lca_out_path=self.vir_lca_out_path,
            k2_report_path=self.k2_report_path, 
            RPM_threshold=self.vir_rpm,
            bt_report_path=self.vir_report_path,
            # FUNC PARAMETERS END
            check_files_dependencies=[],
            check_files_completed=[
                self.vir_report_path
            ]
        )
        step_report_fun = Step(
            step_name=f"Generate fungi taxonomic report", 
            func=generate_bt2_report, 
            force=self.force,
            # FUNC PARAMETERS START
            lca_out_path=self.fun_lca_out_path,
            k2_report_path=self.k2_report_path, 
            RPM_threshold=self.fun_rpm,
            bt_report_path=self.fun_report_path,
            # FUNC PARAMETERS END
            check_files_dependencies=[],
            check_files_completed=[
                self.fun_report_path
            ]
        )
        step_report_overall = Step(
            step_name="Generate overall report",
            func=generate_overall_report,
            # FUNC PARAMETERS START
            bacteria_report_path=self.bac_report_path,
            viruses_report_path=self.vir_report_path,
            fungi_report_path=self.fun_report_path,
            out_path=self.overall_report_path,
            bac_rpm=self.bac_rpm,
            vir_rpm=self.vir_rpm,
            fun_rpm=self.fun_rpm,
            # FUNC PARAMETERS END
            check_files_dependencies=[
                self.bac_report_path,
                self.vir_report_path,
                self.fun_report_path,
            ],
            check_files_completed=[
                self.overall_report_path,
            ],
            force=self.force
        )
        self.steps = []

        if self.bac or self.all:
            self.steps.append(step_report_bac)
        if self.vir or self.all:
            self.steps.append(step_report_vir)
        if self.fun or self.all:
            self.steps.append(step_report_fun)
        if self.overall or self.all:
            self.steps.append(step_report_overall)
        