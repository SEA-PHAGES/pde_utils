import shlex
from subprocess import Popen, PIPE


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------


# SEARCH TOOLS
# -----------------------------------------------------------------------------
def hhsearch(infile_path, db_path, hhr_path, blasttab_path=None,
             add_cons=False, M=50, e=1, qid=0, cov=0, verbose=False):
    verbose_num = 0

    command = (f"hhsearch -i {str(infile_path)} -d {str(db_path)} "
               f"-o {str(hhr_path)} -M {M} "
               f"-e {e} -qid {qid} -cov {cov}")

    if blasttab_path:
        command = " ".join([command, f"-blasttab {str(blasttab_path)}"])

    if add_cons:
        command = " ".join([command, "-add_cons"])

    command = " ".join([command, "-v", str(verbose_num)])

    split_command = shlex.split(command)

    with Popen(args=split_command, stdout=PIPE) as process:
        out, errors = process.communicate()


def hhblits(infile_path, db_path, hhr_path, blasttab_path=None,
            add_cons=False, e=1, qid=0, cov=0, n=2, verbose=False):
    verbose_num = 0

    command = (f"hhblits -i {str(infile_path)} -d {str(db_path)} "
               f"-o {str(hhr_path)} "
               f"-e {e} -qid {qid} -cov {cov} -n {n}")

    if blasttab_path:
        command = " ".join([command, f"-blasttab {str(blasttab_path)}"])

    if add_cons:
        command = " ".join([command, "-add_cons"])

    command = " ".join([command, "-v", str(verbose_num)])

    split_command = shlex.split(command)

    with Popen(args=split_command, stdout=PIPE) as process:
        out, errors = process.communicate()
    pass


if __name__ == "__main__":
    pass
