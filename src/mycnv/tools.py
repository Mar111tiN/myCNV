import os
import pandas as pd
from subprocess import PIPE, run
from io import StringIO
from datetime import datetime as dt
from yaml import CLoader as Loader, load

ansii_colors = {
    "magenta": "[1;35;2m",
    "green": "[1;32;2m",
    "red": "[1;31;1m",
    "cyan": "[1;36;1m",
    "gray": "[1;90;1m",
    "black": "[0m",
    "blue": "[1;94;1m",
}

colors = {
    "process": ansii_colors["green"],
    "time": ansii_colors["magenta"],
    "normal": ansii_colors["gray"],
    "warning": ansii_colors["red"],
    "success": ansii_colors["cyan"],
}


def show_output(text, color="normal", multi=False, time=False, **kwargs):
    """
    get colored output to the terminal
    """
    time = (
        f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ""
    )
    proc = f"\033{colors['process']}Process {os.getpid()}\033[0m : " if multi else ""
    text = f"\033{colors[color]}{text}\033[0m"
    print(time + proc + text, **kwargs)


def show_command(command, list=False, multi=True, **kwargs):
    """
    prints the command line if debugging is active
    """

    proc = f"\033[92mProcess {os.getpid()}\033[0m : " if multi else ""
    if list:
        command = f"\033[1m$ {' '.join(command)}\033[0m"
    else:
        command = f"\033[1m$ {command}\033[0m"
    print(proc + command, **kwargs)
    return


def run_cmd(cmd, show=False, **kwargs):
    if show:
        show_command(cmd, **kwargs)
    exit = run(cmd, shell=True, check=True)
    return exit == 0


def cmd2df(cmd, show=False, **kwargs):
    """
    wrapper for running shell commands directly into a dataframe
    optional output with show argument that passes kwargs to show_command
    """

    if show:
        show_command(cmd, **kwargs)
    cmd_df = pd.read_csv(
        StringIO(run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode("utf-8")),
        sep="\t",
    )
    return cmd_df


def get_CNVconfig(config_file="", local_config={}):
    """
    load the global CNVconfig updated with local config
    """

    with open(config_file, "r") as stream:
        added_config = load(stream, Loader=Loader)["CNV"]
        added_config.update(local_config)
    return added_config
