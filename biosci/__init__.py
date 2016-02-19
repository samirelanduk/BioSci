import os
import subprocess

PYMOL_LOCATION = None

def run_pymol(path="", command_line=False, orders=[]):
    #Try to find pymol
    pymol = ""
    if PYMOL_LOCATION:
        pymol = PYMOL_LOCATION
    else:
        defpath = os.defpath
        if defpath[0] == ":":
            defpath = defpath[1:]
        for directory in defpath.split(":"):
            if "pymol" in os.listdir(directory):
                pymol = "pymol"

    #If pymol was found, run it
    if pymol:
        command = pymol
        if path:
            command += ' "%s"' % path
        if command_line or orders:
            command += " -"
        if command_line:
            command += "cq"
        if orders:
            command += 'd "%s"' %  ";".join(orders)
        print(command)
        subprocess.call(command, shell=True)
    else:
        print("Don't know where Pymol is - please set PYMOL_LOCATION first.")
