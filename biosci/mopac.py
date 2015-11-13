from __future__ import print_function
import subprocess
import datetime
import os

MOPAC_LOCATION = None

#A quick and dirty function for writing to the log
def write_to_log(log, line, send_back=True):
    line = datetime.datetime.strftime(datetime.datetime.now(), "%d %b %Y, %H:%M:%S") + ": " + line
    if log is not None:
        f = open(log, "a")
        f.write(line + "\n")
        f.close()
    if send_back or log is None:
        return line


#Assigns keywords to a .mop file (whatever the top line is to begin with will be removed)
def give_keywords(mop_path, keywords):
    f = open(mop_path, "r")
    lines = f.readlines()
    f.close()
    f = open(mop_path, "w")
    f.writelines([keywords + "\n"] + lines[1:])
    f.close()


#gets the keywords of a .mop file
def get_keywords(mop_path):
    f = open(mop_path, "r")
    lines = f.readlines()
    f.close()
    return lines[0].strip()


#Get path of file name (the location of the file without the file name itself)
def get_path(file_name, symbol_if_this=""):
    #Is there a path?
    if os.sep in file_name:
        #yes
        return os.sep.join(file_name.split(os.sep)[:-1]) + os.sep
    else:
        return symbol_if_this


#Get the file_name from a path
def get_file_name(path):
    if os.sep in path:
        return path.split(os.sep)[-1]
    else:
        return path


#Get all of the file path before the file extension
def get_pre_dot(file_path):
    file_name = get_file_name(file_path)
    if "." in file_name:
        file_name = ".".join(file_name.split(".")[:-1])
    return get_path(file_path) + file_name


#Check if MOPAC_LOCATION is set
def check_mopac_location(func):
    def new_func(*args, **kwargs):
        if MOPAC_LOCATION is not None:
            func(*args, **kwargs)
        else:
            print("You must set the MOPAC_LOCATION constant first.")

    return new_func


#Run MOPAC on a .mop file
@check_mopac_location
def run_mopac(mop_path, log=None, delete=False):
    #What files are present right now?
    start_files = os.listdir(get_path(mop_path, symbol_if_this="./"))

    #Run MOPAC
    print(write_to_log(log,
     "Running MOPAC (with %s as keywords) on %s..." % (get_keywords(mop_path), mop_path)))
    os.environ["LD_LIBRARY_PATH"] = MOPAC_LOCATION
    os.environ["MOPAC_LICENSE"] = MOPAC_LOCATION
    subprocess.call("%s/MOPAC2012.exe %s" % (MOPAC_LOCATION, mop_path), shell=True)
    print(write_to_log(log, "MOPAC run of %s complete" % mop_path))

    #Delete files
    if delete:
        for f in os.listdir(get_path(mop_path, symbol_if_this="./")):
            if f not in start_files and ".pdb" not in f and ".mop" not in f:
                print(write_to_log(log, "Removing %s." % (get_path(mop_path) + f)))
                os.remove(get_path(mop_path) + f)


#Run MOPAC on a PDB file (REQUIRES OPENBABEL CURRENTLY)
@check_mopac_location
def run_mopac_pdb(pdb_path, suffix="_", keywords="PDBOUT", log=None, delete=False):

    #Convert to MOP
    mop_name = get_pre_dot(pdb_path) + suffix + ".mop"
    print(write_to_log(log, "Converting %s to %s" % (pdb_path, mop_name)))
    subprocess.call("babel -ipdb %s -omop %s" % (pdb_path, mop_name), shell=True)

    #Give MOPAC keywords
    give_keywords(mop_name, keywords)

    #Run MOPAC
    run_mopac(mop_name, log=log, delete=delete)
