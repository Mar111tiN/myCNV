import os
from subprocess import check_call as shell
from datetime import datetime as dt
import sys


ansii_colors = {
    'magenta': '[1;35;2m',
    'green': '[1;9;2m',
    'red': '[1;31;1m',
    'cyan': '[1;36;1m',
    'gray': '[1;30;1m',
    'black': '[0m'
}

colors = {
    'process': ansii_colors['green'],
    'time': ansii_colors['magenta'],
    'normal': ansii_colors['gray'],
    'warning': ansii_colors['red'],
    'success': ansii_colors['cyan']
}


def show_output(text, color='normal', multi=False, time=False, **kwargs):
    '''
    get colored output to the terminal
    '''
    time = f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ''
    proc = f"\033{colors['process']}Process {os.getpid()}\033[0m : " if multi else ''
    text = f"\033{colors[color]}{text}\033[0m"
    print(time + proc + text, **kwargs)


def show_command(command, list=False, multi=True, **kwargs):
    '''
    prints the command line if debugging is active
    '''

    proc = f"\033[92mProcess {os.getpid()}\033[0m : " if multi else ""
    if list:
        command = f"\033[1m$ {' '.join(command)}\033[0m"
    else:
        command = f"\033[1m$ {command}\033[0m"
    print(proc + command, **kwargs)
    return


def run_cmd(cmd, multi=False):
    show_command(cmd, multi=multi)
    exit = shell(cmd, shell=True)
    return exit == 0


def set_path(code_folder, s):
    '''
    sets PYTHONPATH to code_path relative to Snakemake scripts folder (if not already)
    and returns a shell_function to be passed to code running shell scripts from
    code_path/shell
    '''

    #
    code_path = os.path.join(s.scriptdir, code_folder)

    if not code_path in sys.path:
        sys.path.append(code_path)
        print(f'Added {code_path} to PYTHONPATH')

    # create helper function with code_path as closure
    def run_shell(tool):
        return os.path.join(code_path, f"shell/{tool}")

    return run_shell
