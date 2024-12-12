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
                   'CreateHistograms',
                   'PlotVariables1D', 'PlotVariables2D',
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



pargs = parser.parse_args()

yml_config = None
yml_config_file = pargs.inconfig

with open(yml_config_file,'r') as conf:
    yml_config = yaml.safe_load(conf)

#main_func = yml_config.get("main")
main_func = pargs.func
assert main_func in ListOfMainFuncs, f"Error: {main_func} is wrong"

era           = yml_config.get("era")
run           = yml_config.get("run")
postfix       = yml_config.get("postfix")
islimited     = yml_config.get("limited")
limit         = "limited" if islimited else "full"
nworkers      = yml_config.get("workers")
tasks_per_job = yml_config.get("tasks_per_job")

#cmdfile = os.path.join(os.getcwd(),"cmdlogs",f"cmdlog_selectEvents_{datetime_tag}.sh")
jobfile = os.path.join(os.getcwd(),"cmdlogs",f"joblog_run{run}_{era}{postfix}_{limit}_{main_func}_{datetime_tag}.log")
logfile = os.path.join(os.getcwd(),"cmdlogs",f"cmdlog_run{run}_{era}{postfix}_{limit}_{main_func}_{datetime_tag}.log")
logger  = setup_logger(logfile)


logger.info(f"Yaml     : {yml_config_file}")
logger.info(f"Execute  : cf.{main_func}")
logger.info(f"Era      : {era}")
logger.info(f"Run      : {run}")
logger.info(f"Postfix  : {postfix}")
logger.info(f"Limited? : {islimited}")
logger.info(f"nWorkers : {nworkers}")
logger.info(f"nTasks_per_job : {tasks_per_job}")

wrapper   = yml_config.get("wrapper")
logger.info(f"Wrapper  : {wrapper}")

main_args = yml_config.get("args")
config    = main_args.get("config")
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

version   = main_args.get('version')
if version.startswith("dummy"):
    version = f"Run{run}_{era}{postfix}_{limit}_{datetime_tag}"
logger.info(f"version  : {version}")

extras = []
if "extras" in list(main_args.keys()):
    extras = [f"--{extra}" for extra in main_args.get("extras")]

categories = ",".join(main_args.get("categories")) if "categories" in list(main_args.keys()) else None
variables  = ",".join(main_args.get("variables")) if "variables" in list(main_args.keys()) else None

cmd_list = [
    "law", "run", f"cf.{main_func}",
    "--config", config,
    "--dataset", datasets,
    "--workflow", workflow,
    #"--branch", branch,
    "--version", version,
    #"&>", jobfile, '&'
]

#skip_wrapping = ((main_func == "PlotVariables1D") or (main_func == "PlotVariables2D")
skip_wrapping = main_func.startswith("PlotVariables")
shifts_val = "nominal,{"+f"{shifts}"+"}_{up,down}" if shiftList else "nominal"
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
        if main_func.startswith("CreateHistograms"):
            cmd_list.append(f"--cf.{main_func}-variables")
            cmd_list.append(variables)
    else:
        cmd_list = [
            "law", "run", f"cf.{main_func}",
            "--config", config,
            "--datasets", datasets,
            "--processes", processes,
            f"--workflow", workflow,
            f"--branch", branch,
            "--version", version,
            "--categories", categories,
            "--variables", variables,
            "--workers", nworkers,
            "--pilot",
            f"--log-file", jobfile,
        ] + extras
    
cmd_list = [str(item) for item in cmd_list]
cmd = " ".join(cmd_list)

logger.info(f"CMD      : \n\n{cmd}\n")
logger.info(f"Copy this line above, paste it and hit enter : BABUSHCHA !")

#logger.info(f"writing CMD in {cmdfile}")
#cmdf = open(cmdfile, 'w')
#cmdf.write("#!/bin/sh\n\n")
#cmdf.write(cmd + '\n')
#cmdf.close()
