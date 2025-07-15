import os
import yaml
import argparse 
import datetime
import time
import sys
import logging
from subprocess import Popen, PIPE


def setup_logger(log_file):
    # Create a logger
    logger = logging.getLogger('main')
    logger.setLevel(logging.DEBUG)

    # Create a file handler
    filemode = 'a+'
    file_handler = logging.FileHandler(log_file, mode=filemode)
    file_handler.setLevel(logging.DEBUG)

    # Create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)

    # Create a formatter and add it to the file handler
    formatter = logging.Formatter('%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s','%Y-%m-%d:%H:%M:%S')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add the file handler to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def runShellCmd(cmdList):
    process = Popen(cmdList, stdout=PIPE, stderr=PIPE)
    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = process.poll()


def run_command(command, logger):
    process = Popen([*command, 'color=always'], stdout=PIPE, stderr=PIPE, text=True, bufsize=1)

    # Log and display output in real-time
    with process.stdout, process.stderr:
        for line in iter(process.stdout.readline, ''):
            sys.stdout.write(line)
            logger.debug(line.strip())

        for line in iter(process.stderr.readline, ''):
            sys.stderr.write(line)
            logger.error(line.strip())
    
    process.wait()

#from utilities import setup_logger, runShellCmd, run_command
datetime_tag = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
ListOfMainFuncs = ['CalibrateEvents',
                   'SelectEvents',
                   'ReduceEvents',
                   'PlotCutFlow',
                   'ProduceColumns',
                   'CreateHistograms','MergeHistograms',
                   'PlotVariables1D', 'PlotVariables2D', 'PlotVariablesPerProcess2D',
                   'PlotShiftedVariables1D', 'PlotShiftedVariables2D',
                   'UniteColumns']


if not os.path.exists("cmdlogs"): os.mkdir("cmdlogs")
else: print("logs exists")

parser = argparse.ArgumentParser(description='Pre-Processing')
parser.add_argument('-i',
                    '--inconfig',
                    type=str,
                    required=True,
                    help="cf.SelectEvents arguments")
parser.add_argument('-f',
                    '--func',
                    type=str,
                    required=True,
                    help="e.g. SelectEvents")
parser.add_argument('-l',
                    '--limited',
                    action='store_true',
                    required=False,
                    default=False,
                    help="full or limited")

pargs = parser.parse_args()

yml_config = None
yml_config_file = pargs.inconfig

with open(yml_config_file,'r') as conf:
    yml_config = yaml.safe_load(conf)

main_func = pargs.func
assert main_func in ListOfMainFuncs, f"Error: {main_func} is wrong"

islimited = pargs.limited
yml_config = yml_config.get("FULL") if not islimited else yml_config.get("LIMITED")

#from IPython import embed; embed()

era           = yml_config.get("era")
run           = yml_config.get("run")
postfix       = yml_config.get("postfix")
limit         = "limited" if islimited else "full"
nworkers      = yml_config.get("workers")
tasks_per_job = yml_config.get("tasks_per_job")
wrapper       = yml_config.get("wrapper")

main_args = yml_config.get("args")
version_tag = main_args.get("version")
if version_tag == "dummy":
    version_tag = f"Run{run}_{era}{postfix}_{limit}_{datetime_tag}"

logfile = os.path.join(os.getcwd(),"cmdlogs",f"cmd_{main_func}_{version_tag}.log")
logger  = setup_logger(logfile)

logger.info(f"Yaml     : {yml_config_file}")
logger.info(f"Execute  : cf.{main_func}")
logger.info(f"Era      : {era}")
logger.info(f"Run      : {run}")
logger.info(f"Postfix  : {postfix}")
logger.info(f"Limited? : {islimited}")
logger.info(f"nWorkers : {nworkers}")
logger.info(f"nTasks_per_job : {tasks_per_job}")
logger.info(f"Wrapper  : {wrapper}")

config    = main_args.get("config")
config    = config if not islimited else f"{config}_limited" 
logger.info(f"Config   : {config}")

dataset_list = main_args.get("datasets")
if len(dataset_list) > 1 and wrapper == False:
    logger.info(f"\tmore than one dataset, making wrapper True")
    wrapper = True
datasets  = ",".join(dataset_list)
logger.info(f"datasets : {datasets}")

processes = ",".join(main_args.get("processes"))
logger.info(f"processes: {processes}")

shiftList = main_args.get("shifts")
shifts = ",".join(shiftList) if shiftList else ""
logger.info(f"shifts: {shifts}")

workflow  = main_args.get("workflow")
logger.info(f"worklflow: {workflow}")

branch    = main_args.get("branch")
logger.info(f"branch   : {branch}")

version = version_tag
logger.info(f"version  : {version}")

extras = []
if "extras" in list(main_args.keys()):
    extras = [f"--{extra}" for extra in main_args.get("extras")]

categories = ",".join(main_args.get("categories")) if "categories" in list(main_args.keys()) else None
variables  = ",".join(main_args.get("variables")) if "variables" in list(main_args.keys()) else None
variables1D = ",".join([var for var in main_args.get("variables") if len(var.split('-')) < 2])
variables2D = ",".join([var for var in main_args.get("variables") if len(var.split('-')) == 2])

shifts_val = "nominal,{"+f"{shifts}"+"}_{up,down}" if shiftList else "nominal"

skip_wrapping = main_func.startswith("PlotVariables")
if wrapper:    
    if not skip_wrapping:
        cmd_list = [
            "law", "run", f"cf.{main_func}Wrapper",
            "--config", config,
            "--datasets", datasets,
            f"--shifts", shifts_val,
            f"--cf.{main_func}-workflow", workflow,
            #f"--cf.{main_func}-branch", branch,
            f"--cf.{main_func}-tasks-per-job", tasks_per_job,
            f"--cf.{main_func}-pilot",
            "--version", version,
            "--workers", nworkers,
            #f"--cf.{main_func}-log-file", jobfile,
        ]
        if main_func.startswith("CreateHistograms") or main_func.startswith("MergeHistograms"):
            cmd_list.append(f"--cf.{main_func}-variables")
            cmd_list.append(variables)
    else:
        _variables = variables1D if main_func == "PlotVariables1D" else variables2D
        logger.info(f"variables : {_variables}")
        cmd_list = [
            "law", "run", f"cf.{main_func}",
            "--config", config,
            "--datasets", datasets,
            "--processes", processes,
            f"--workflow", workflow,
            #f"--branch", branch,
            "--version", version,
            "--categories", categories,
            "--variables", _variables,
            "--workers", nworkers,
            "--pilot",
            #f"--log-file", jobfile,
        ] + extras

else:
    cmd_list = [
        "law", "run", f"cf.{main_func}",
        "--config", config,
        "--dataset", datasets,
        "--workflow", workflow,
        "--workers", nworkers,
        #"--branch", branch,
        "--version", version,
        "--pilot",
        #"&>", jobfile, '&'
    ]
    if branch > 0:
        cmd_list.append("--branch")
        cmd_list.append(branch)
    if main_func in ["CreateHistograms", "MergeHistograms", "PlotVariables1D", "PlotVariables2D"]:
        cmd_list.append("--variables")
        if main_func in ["PlotVariables1D", "PlotVariables2D"]:
            cmd_list.append(variables1D) if "1D" in main_func else cmd_list.append(variables2D)
            cmd_list.append("--categories")
            cmd_list.append(categories)
        else:
            cmd_list.append(variables)

    
        
cmd_list = [str(item) for item in cmd_list]
cmd = " ".join(cmd_list)

logger.info(f"CMD      : \n\n{cmd}\n")
logger.info(f"Copy this line above, paste it and hit enter : BABUSHCHA !")

#logger.info(f"writing CMD in {cmdfile}")
#cmdf = open(cmdfile, 'w')
#cmdf.write("#!/bin/sh\n\n")
#cmdf.write(cmd + '\n')
#cmdf.close()
