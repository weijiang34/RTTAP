import os
import sys
import subprocess as sp
import logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

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

def safe_run_cmd(cmd: str, log_file: str = None):
    """
    Safely execute shell commands with automatic logging and error handling.
    Args:
        cmd (str): The shell command to execute, or a path to a shell script.
    """
    try:
        logging.debug(f"\tRunning script/command: {cmd}")
        if log_file:
            with open(log_file, 'w') as f:
                sp.run(cmd, shell=True, check=True, stdout=f, stderr=f)
        else:
            sp.run(cmd, shell=True, check=True)
    except sp.CalledProcessError as e:
        logging.error(f"\tScript/command failed: {cmd}\n{e}")
        raise

def safe_remove(path: str):
    """
    Safely remove files to avoid accidental deletion.
    Args:
        path (str): The file path to remove.
    """
    if os.path.exists(path):
        os.remove(path)
        # logging.info(f"Removed file: {path}")
    else:
        logging.warning(f"\tFile {path} not found for removal.")

# def check_file_exists(file_path, description):
#     """Check if file exists. If not, report error and exit."""
#     if not os.path.exists(file_path):
#         logging.error(f"{description} file not found: {file_path}")
#         sys.exit(1)
#     return True

# def check_step_completed(step_file, description):
#     """Check if steps are finished."""
#     if os.path.exists(step_file):
#         logging.info(f"{description} already completed, skipping...")
#         return True
#     return False