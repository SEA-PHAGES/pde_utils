import shlex
from subprocess import Popen, PIPE, DEVNULL

from pathlib import Path
from pdm_utils.classes.alchemyhandler import AlchemyHandler

#GLOBAL VARIABLES
#-----------------------------------------------------------------------------

#SEARCH TOOLS
#-----------------------------------------------------------------------------
def hhsearch(infile_path, db_path, hhr_path, blasttab_path=None, add_cons=False,
                                             M=50, e=1, qid=0, cov=0,
                                             verbose=False):
    command = (f"hhsearch -i {str(alignment_path)} -d {str(db_path)} " 
               f"-o {str(hhr_path)} -M {M} "
               f"-e {e} -qid {qid} -cov {cov}")

    if blasttab_path:
        command = " ".join([command, f"-blasttab {str(blasttab_path)}"])
        
    if add_cons:
        command = " ".join([command, "-add_cons"])

    if verbose:
        command = " ".join([command, "-v 2"])
    else:
        command = " ".join([command, "-v 0"])


    split_command = shlex.split(command)

    with Popen(args=split_command, stdout=PIPE) as process:
        out, errors = process.communicate()

def hhblitz():
    pass

def not_main():
    pass

if __name__ == "__main__":
    pass    
