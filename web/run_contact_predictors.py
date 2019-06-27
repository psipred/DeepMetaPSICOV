import sys
from shutil import copyfile
import os
import subprocess
# sys.argv[1] the final .aln file
# sys.argv[2] psicov path
# sys.argv[3] freecontac path
# sys.argv[4] ccmpred path
# sys.argv[5] ccmpred output file
# sys.argv[6] tmpdir


def file_len(fname):
    with open(fname) as f:
        i = None
        for i, l in enumerate(f):
            pass
        if i:
            return i + 1
        else:
            return 0


def run_exe(args, name, outfile):
    """
        Function takes a list of command line args and executes the subprocess.
        Sensibly tries to catch some errors too!
    """
    code = 0
    print("Running "+name)
    try:
        print(' '.join(args))
        code = subprocess.Popen(' '.join(args), shell=True, stdout=subprocess.PIPE)
        out = code.communicate()[0]
        if outfile:
            output = open(outfile, "w")
            output.write(out.decode())
            output.close()
    except Exception as e:
        print(str(e))
        sys.exit(1)
    # if code.returncode != 0:
    #     print(name+" Non Zero Exit status: "+str(code))
    #     sys.exit(code)


aln_length = file_len(sys.argv[1])
file_id = sys.argv[1][:-6]
if aln_length >= 10:
    print("Gonna run psicov and shizzle")
    processPSICOV_args = ["timeout",
                          #"86400",
                          "7200",
                          sys.argv[2],
                          "-o",
                          "-d",
                          "0.03",
                          sys.argv[1],
                          ]
    run_exe(processPSICOV_args, "psicov", sys.argv[6]+"/"+sys.argv[1][:-4]+".psicov")

    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.dirname(sys.argv[3])))
    input_path = sys.argv[6]+"/"+sys.argv[1]
    processFreecontact_args = ["timeout",
                               #"86400",
                               "7200",
                               sys.argv[3],
                               "<",
                               input_path,
                               ]
    run_exe(processFreecontact_args, "freecontact", sys.argv[6]+"/"+sys.argv[1][:-4]+".evfold")
    os.chdir(cwd)

    processCcmpred_args = ["timeout",
                           #"86400",
                           "7200",
                           sys.argv[4],
                           "-t",
                           "24",
                           sys.argv[1],
                           sys.argv[5],
                           ]
    run_exe(processCcmpred_args, "ccmpred", None)

if not os.path.isfile(sys.argv[6]+"/"+sys.argv[1][:-4]+".psicov"):
    open(sys.argv[6]+"/"+sys.argv[1][:-4]+".psicov", 'a').close()
if not os.path.isfile(sys.argv[6]+"/"+sys.argv[1][:-4]+".evfold"):
    open(sys.argv[6]+"/"+sys.argv[1][:-4]+".evfold", 'a').close()
if not os.path.isfile(sys.argv[6]+"/"+sys.argv[1][:-4]+".ccmpred"):
    open(sys.argv[6]+"/"+sys.argv[1][:-4]+".ccmpred", 'a').close()
