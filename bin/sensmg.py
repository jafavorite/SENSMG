#!/usr/bin/python
# script to control first-order sensitivity calculation for partisn.

import os, sys, shutil, re, glob, string

def write_control(ictrl):
    ctrl = open("control", "w")
    ctrl.write("%d %d\n" % (ictrl, ITER))
    ctrl.write("%d %d %d %d %d %d %d %d %d %d\n" % (ISN, ISCT, NGROUP, ICHINORM, ISRCACC_NO, FISSDATA, IANGFLUX, AFLXFRM, IVER, IDBG))
    ctrl.write("%d %d %d %d %d\n" % (NM, NEL, NR, NZ, NRRR))
    ctrl.write(NCBC+"\n")
    ctrl.write("%d %d %d\n" % (IMISC, IALPHAN, NAG))
    ctrl.write("%d %d\n" % (NOFXUP, ISECORDER))
    ctrl.write(str(EPSI)+" "+str(EPSIG)+"\n")
    ctrl.write(IFILE+"\n")
    ctrl.write(PART_SHORT+"\n")
    ctrl.write("%d %d\n" % (RPLANE, ZPLANE))
    ctrl.close()

def exit_sens(ecode):
    if ecode == 3:
        print "error in input file."
        sys.stdout.flush()
# reset modules if they were set. (this is not necessary if only this process
# has the modified modules; i'm not sure.)
    if USE_EXISTING == "no" and MY_MODULES == "no":
        if LOADEDMODULES_org != None:
            pm=module("purge")
            for m in LOADEDMODULES_org.split(":"):
                pm=module("load", m)
# close log file if needed
    if log.closed == False:
        log.close()
# exit
    sys.exit(ecode)

def module(command, *arguments):
# load and manipulate modules using python.
    if MACH == "luna":
        commands = command+" "+string.join(arguments)
        python.module(commands)
    elif MACH == "sn" or MACH == "frost":
# this is how to load and manipulate modules using python on snow.
# adapted from /usr/share/lmod/lmod/init/env_modules_python.py
        commands = os.popen('/usr/share/lmod/lmod/libexec/lmod python %s %s'\
                    % (command, string.join(arguments))).read()
        exec commands

# set location of codes and sources4c/misc data.
sensmg_exe = "/usr/projects/data/nuclear/working/sensitivities/bin/sensmg"
sources_exe = "/usr/projects/data/nuclear/working/sensitivities/sources4c/bin/sources4c.jaf"
sources_dir = "/usr/projects/data/nuclear/working/sensitivities/sources4c/data"
misc_exe = "/usr/projects/data/nuclear/working/sensitivities/isc-1.3.0/bin/misc"
os.environ["ISCDATA"] = "/usr/projects/data/nuclear/working/sensitivities/isc-1.3.0/data"
os.environ["SENS_DATA"] = "/usr/projects/data/nuclear/working/sensitivities/data"

# setenv NDI_GENDIR_PATH /usr/projects/data/nuclear/ndi/2.0.20/share/gendir.all
# setenv NDI_GENDIR_PATH /usr/projects/data/nuclear/ndi/2.1.3/share/gendir.all
# setenv NDI_GENDIR_PATH /usr/projects/data/nuclear/ndi/devel/share/gendir.all
# these are special libraries for testing
# setenv NDI_GENDIR_PATH /users/fave/mt71x_h_sab/gendir
# setenv NDI_GENDIR_PATH /users/fave/mt71x_c_sab/gendir

# LIBNAME = "mt71x"
# LIBNAME = "mendf71x"
# LIBNAME = "mendf70x"
# LIBNAME = "mtmg01"
# LIBNAME = "kynea3"
# LIBNAME = "mendf6"
# LIBNAME = "mendf5"
# testing only
# LIBNAME = "vitb6"
# LIBNAME = "scale"
# LIBNAME = "acti"

#  edit           partisn_7_72                 partisn_8_x
# position  description          name    description          name
#           ("MENDF Library Edit...")    ("NDI Library Edit...")
#   1       fission spectrum     chi     fission spectrum     chi
#   2       eff. nusigf          nusigf  eff. nusigf          nusigf
#   3       total                total   total                total
#   4       absorption           abs     absorption           abs
#   5       chi                  ?       chi                  ?
#   6       elastic              mend1   elastic              (n,n)
#   7       inelastic            mend2   inelastic            (n,n')
#   8       (n,2n)               mend3   (n,2n)               (n,2n)
#   9       (n,3n)               mend4   (n,3n)               (n,3n)
#  10       (n,gamma)            mend5   (n,gamma)            (n,g)
#  11       (n,alpha)            mend6   (n,proton)           (n,p)
#  12       (n,proton)           mend7   (n,alpha)            (n,a)
#  13       (n,deuteron)         mend8   1st-chance fission   (n,f)
#  14       (n,triton)           mend9   2nd-chance fission   (n,n')f
#  15       (n,3He)              mend10  3rd-chance fission   (n,2n)f
#  16       sum of all fission   n-fiss  sum of all fission   (n,F)
#  17       prompt fission spec. mend12  prompt fission spec. chi_pr
#  18       total fission spec.  mend13  total fission spec.  chi_tot
#  19             --               --    (n,deuteron)         (n,d)
#  20             --               --    (n,triton)           (n,t)
# =====================================================================

# iangflux = 0/1 use only moments/use angular fluxes for sigt term
# for one-d spheres, this should be 1. it is automatically set to 0
# for cylinders.
IANGFLUX = 1

# normally this is 0.
# set to 1 to write verification input files for central diffs.
# set to 2 to write verification input files for a pert study.
IVER = 0

# NOFXUP = 0/1 negative flux fix-up/no fix-up
NOFXUP = 1

# normally this is 0.
# set to 1 to write some debug lines; not fully implemented.
IDBG = 0

# defaults for command-line variables.
IFILE = "sensmg_inp"
NP = 1
NGROUP = 30
ISN = 32
ISCT = 3
# EPSI is the convergence requirement for all partisn calculations.
# EPSIG is the convergence requirement for the successive approximations
# of the generalized adjoint functions.
EPSI = 1.e-6
EPSIG = 1.e-5
CHINORM = "full"
# FISSDATA = 0/1/2 fission transfer matrix/chi matrix/chi vector
# should be 0 for ndi. use 2 for debugging with central differences.
FISSDATA = 0
# AFLXFRM = 0/1 no/yes Use the angular flux formulation (ievt=2 & dsasrch=2 only).
# may help with subcritical alphas.
AFLXFRM = 0
# SRCACC_NO = none/for/adj/for+adj turns off source acceleration for 
# neither/forward/adjoint/both. source acceleration is always off for generalized
# adjoints. when source acceleration is on, it is the PARTISN default.
SRCACC_NO = "none"
USE_EXISTING = "no"
MY_MODULES = "no"

# ALPHA_N = yes/no use/ignore (alpha,n) neutron sources
ALPHA_N = "yes"
# USE_MISC = yes/no use MISC/use SOURCES to compute spontaneous fission sources
USE_MISC = "yes"
# NAG = number of alpha-particle energy groups for sources4c
NAG = 100

# RPLANE and ZPLANE are r and z planes (1 to nr, 0 to nz) to use angular fluxes.
# for older partisn versions (before 8_27).
RPLANE = -1
ZPLANE = -1

# TIMEDEP = 0/1 = don't/do time-dependent subcritical multiplication calculation
# must be 0 for now
TIMEDEP = 0

# SECORDER = yes/no do/don't compute 2nd-order sensitivities
SECORDER = "no"

# set here in case of input errors
LOADEDMODULES_org = None

# input parser. there must be an even number of arguments, but not more than 36.
IERROR = 0
rem = len(sys.argv) % 2 # remainder operator
if rem == 0 or len(sys.argv) > 38:
    print "error on command line. odd number of entries or too many."
    IERROR = 1
i = 2
while  i < len(sys.argv):
    if sys.argv[i-1] == "-i":
        IFILE = sys.argv[i]
    elif sys.argv[i-1] == "-np":
        NP = int(sys.argv[i])
    elif sys.argv[i-1] == "-ngroup":
        NGROUP = int(sys.argv[i])
    elif sys.argv[i-1] == "-isn":
        ISN = int(sys.argv[i])
    elif sys.argv[i-1] == "-isct":
        ISCT = int(sys.argv[i])
    elif sys.argv[i-1] == "-epsi":
        EPSI = sys.argv[i]
    elif sys.argv[i-1] == "-epsig":
        EPSIG = sys.argv[i]
    elif sys.argv[i-1] == "-chinorm":
        CHINORM = sys.argv[i]
    elif sys.argv[i-1] == "-fissdata":
        FISSDATA = int(sys.argv[i])
    elif sys.argv[i-1] == "-aflxfrm":
        AFLXFRM = int(sys.argv[i])
    elif sys.argv[i-1] == "-alpha_n":
        ALPHA_N = sys.argv[i]
    elif sys.argv[i-1] == "-misc":
        USE_MISC = sys.argv[i]
    elif sys.argv[i-1] == "-nag":
        NAG = int(sys.argv[i])
    elif sys.argv[i-1] == "-srcacc_no":
        SRCACC_NO = sys.argv[i]
    elif sys.argv[i-1] == "-use_existing":
        USE_EXISTING = sys.argv[i]
    elif sys.argv[i-1] == "-my_modules":
        MY_MODULES = sys.argv[i]
    elif sys.argv[i-1] == "-rplane":
        RPLANE = int(sys.argv[i])
    elif sys.argv[i-1] == "-zplane":
        ZPLANE = int(sys.argv[i])
    elif sys.argv[i-1] == "-2nd_order":
        SECORDER = sys.argv[i]
    else:
        print "error on command line. illegal entry:", sys.argv[i-1]
        IERROR = 1
    i = i+2

# echo paramaters
log = open("sensmg.log", "w")
log.write("PARAMETERS:\n")
comm_line = " ".join(sys.argv)
log.write("  COMMAND LINE="+comm_line+"\n")
log.write("  CONTROL SCRIPT="+sys.argv[0]+"\n") 
log.write("  SENSMG_CODE="+sensmg_exe+"\n")
log.write("  INPUT="+IFILE+"\n")
log.write("  NGROUP="+str(NGROUP)+"\n")
log.write("  ISN="+str(ISN)+"\n")
log.write("  ISCT="+str(ISCT)+"\n")
log.write("  EPSI="+str(EPSI)+"\n")
log.write("  EPSIG="+str(EPSIG)+"\n")
log.write("  NOFXUP="+str(NOFXUP)+"\n")
log.write("  SECORDER="+str(SECORDER)+"\n")
log.write("  IVER="+str(IVER)+"\n")
log.write("  CHINORM="+CHINORM+"\n")
log.write("  FISSDATA="+str(FISSDATA)+"\n")
log.write("  AFLXFRM="+str(AFLXFRM)+"\n")
log.write("  SRCACC_NO="+SRCACC_NO+"\n")
log.write("  IANGFLUX="+str(IANGFLUX)+"\n")
log.write("  USE_EXISTING="+USE_EXISTING+"\n")
log.write("  MY_MODULES="+MY_MODULES+"\n")
log.write("  IDBG="+str(IDBG)+"\n")
log.write("  ALPHA_N="+str(ALPHA_N)+"\n")
log.write("  USE_MISC="+str(USE_MISC)+"\n")
log.write("  RPLANE="+str(RPLANE)+"\n")
log.write("  ZPLANE="+str(ZPLANE)+"\n")
log.write("  TIMEDEP="+str(TIMEDEP)+"\n")
if USE_MISC == "yes":
    log.write("  MISC="+misc_exe+"\n")
    log.write("  ISCDATA="+os.environ.get("ISCDATA")+"\n")
if USE_MISC == "no" or ALPHA_N == "yes":
    log.write("  SOURCES="+sources_exe+"\n")
    log.write("  SOURCES_DIR="+sources_dir+"\n")
    log.write("  NAG="+str(NAG)+"\n")
sys.stdout.flush()
log.flush()

# ensure OMP_NUM_THREADS is 1 due to bug in SENSMG.
if os.environ.get("OMP_NUM_THREADS") != "1":
    os.environ["OMP_NUM_THREADS"] = "1"
log.write("  OMP_NUM_THREADS="+os.environ.get("OMP_NUM_THREADS")+"\n")
sys.stdout.flush()
log.flush()

# HOSTNAME is needed for partisn executable and modules.
HN = os.environ.get("HOSTNAME")[0:2]
# print HN
if HN == "sn" or HN == "tt" or HN == "tr":
    MACH = HN
elif HN == "fi" or HN == "ic":
    MACH = "frost"
elif HN == "lu":
    MACH = "luna"
else:
    MACH = "None"
    print "warning. unrecognized HOSTNAME. no support for default partisn or for modules."
    log.write("warning. unrecognized HOSTNAME. no support for default partisn or for modules.\n")
sys.stdout.flush()
log.flush()

if MY_MODULES != "no" and MY_MODULES != "yes":
    print "error on command line. -my_modules=", MY_MODULES
    MY_MODULES == "no"
    IERROR = 1

# set partisn executable.
IERRORP = 0
SENS_PARTISN = os.environ.get("SENS_PARTISN")
if SENS_PARTISN == None:
    print "warning. SENS_PARTISN not set in environment; using default."
    log.write("warning. SENS_PARTISN not set in environment; using default.\n")
    if MACH == "sn":
        SENS_PARTISN = "/usr/projects/lindet/rel8_27/8_27_15/snow-intel-17.0.4-openmpi-2.1.2/partisn"
    elif MACH == "frost":
        SENS_PARTISN = "/usr/projects/lindet/rel8_27/8_27_15/frost-intel-17.0.4-openmpi-2.1.2/partisn"
    elif MACH == "luna":
        SENS_PARTISN = "/usr/projects/lindet/rel8_27/8_27_15/luna-intel-17.0.4-openmpi-2.1.2/partisn"
    else:
        SENS_PARTISN = "None"
        print "error. no SENS_PARTISN for this machine."
        log.write( "error. no SENS_PARTISN for this machine.\n")
        IERRORP = 1

if "partisn_5_97" in SENS_PARTISN:
    PARTISN_EXE = SENS_PARTISN
else:
    PARTISN_EXE = "mpirun -np "+str(NP)+" "+SENS_PARTISN
log.write("  NP="+str(NP)+"\n")
log.write("  PARTISN_EXE="+PARTISN_EXE+"\n")
sys.stdout.flush()
log.flush()
# TODO list of supported partisn appears in two places; consolidate.
# use a dictionary?
if "partisn_7_72" in PARTISN_EXE:
    PART_SHORT = "7_72"
elif "partisn_8_03" in PARTISN_EXE:
    PART_SHORT = "8_03"
elif "partisn_8_20a" in PARTISN_EXE:
    PART_SHORT = "8_20a"
elif "8_27_15" in PARTISN_EXE:
    PART_SHORT = "8_27_15"
elif "partisn_5_97" in PARTISN_EXE:
    PART_SHORT = "5_97"
else:
    PART_SHORT = "None"
    print "error. the only versions of partisn supported are 5_97, 7_72, 8_03, 8_20a, and 8_27_15."
    log.write("error. the only versions of partisn supported are 5_97, 7_72, 8_03, 8_20a, and 8_27_15.\n")
    IERRORP = 1

# check partisn executable
if os.path.exists(SENS_PARTISN) == False:
    print "error. "+SENS_PARTISN+" does not exist."
    log.write("error. "+SENS_PARTISN+" does not exist.\n")
    IERRORP = 1
else:
    if os.access(SENS_PARTISN, os.X_OK) == False: 
        print "error. "+SENS_PARTISN+" cannot be executed."
        log.write("error. "+SENS_PARTISN+" cannot be executed.\n")
        IERRORP = 1
sys.stdout.flush()
log.flush()

# load modules for partisn. also check for a supported version of partisn.
if IERROR == 0 and USE_EXISTING == "no":
# TODO catch error on module load, e.g. if the module does not exist.
# not necessary if this will be more general?
    if MY_MODULES == "no":
        LOADEDMODULES_org = os.environ.get("LOADEDMODULES")
# import python needed on moonlight and luna.
        if MACH == "luna":
            sys.path.append("/usr/share/Modules/init")
            import python
        pm=module("purge")
        if "partisn_7_72" in PARTISN_EXE:
            pm=module("load", "intel/13.1.0 openmpi/1.6.3")
        elif "partisn_8_03" in PARTISN_EXE:
            pm=module("load", "intel/15.0.5 openmpi/1.6.5 quo/1.2.4")
        elif "partisn_8_20a" in PARTISN_EXE:
            pm=module("load", "friendly-testing")
            pm=module("load", "user_contrib")
            pm=module("load", "python/2.7-anaconda-2.1.0")
            pm=module("load", "intel/16.0.3 openmpi/1.10.5 quo/1.2.9")
        elif "8_27_15" in PARTISN_EXE:
            pm=module("load", "friendly-testing")
            pm=module("load", "user_contrib")
            pm=module("load", "python/2.7-anaconda-4.1.1")
            pm=module("load", "intel/17.0.4 openmpi/2.1.2 quo/1.3")
        elif "partisn_5_97" in PARTISN_EXE:
            pm=module("load", "intel/17.0.4")
sys.stdout.flush()
log.flush()

# check sensmg executable
IERRORS = 0
if os.path.exists(sensmg_exe) == False:
    print "error. executable "+sensmg_exe+" does not exist."
    log.write("error. executable "+sensmg_exe+" does not exist.\n")
    IERRORS = 1
else:
    if os.access(sensmg_exe, os.X_OK) == False: 
        print "error. "+sensmg_exe+" cannot be executed."
        log.write("error. "+sensmg_exe+" cannot be executed.\n")
        IERRORS = 1
if IERRORS == 1:
    print "  check file system or reset variable sensmg_exe in "+sys.argv[0]+"."
    log.write("  check file system or reset variable sensmg_exe in "+sys.argv[0]+".\n")
    IERRORP = 1
sys.stdout.flush()
log.flush()

if FISSDATA  < 0 or FISSDATA > 2:
    print "error on command line. -fissdata=", FISSDATA
    IERROR = 1

if AFLXFRM  < 0 or AFLXFRM > 1:
    print "error on command line. -aflxfrm=", AFLXFRM
    IERROR = 1

if ALPHA_N == "yes":
    IALPHAN = 1
elif ALPHA_N == "no":
    IALPHAN = 0
else:
    print "error on command line. -alpha_n=", ALPHA_N
    IERROR = 1

if USE_MISC == "yes":
    IMISC = 1
elif USE_MISC == "no":
    IMISC = 0
else:
    print "error on command line. -misc=", USE_MISC
    IERROR = 1

if CHINORM == "none":
    ICHINORM = 0
elif CHINORM == "full":
    ICHINORM = 1
elif CHINORM == "partial":
    ICHINORM = 2
else:
    print "error on command line. -chinorm=", CHINORM
    IERROR = 1

if SRCACC_NO == "none":
    ISRCACC_NO = 0
elif SRCACC_NO == "for":
    ISRCACC_NO = 1
elif SRCACC_NO == "adj":
    ISRCACC_NO = 2
elif SRCACC_NO == "for+adj" or SRCACC_NO == "adj+for":
    ISRCACC_NO = 3
else:
    print "error on command line. -srcacc_no=", SRCACC_NO
    IERROR = 1

if USE_EXISTING != "no" and USE_EXISTING != "yes":
    print "error on command line. -use_existing=", USE_EXISTING
    IERROR = 1

if NAG  < 1:
    print "error on command line. -nag=", NAG
    IERROR = 1

if SECORDER == "yes":
    ISECORDER = 1
elif SECORDER == "no":
    ISECORDER = 0
else:
    print "error on command line. -2nd_order=", SECORDER
    IERROR = 1

# exit if there is an input error.
if IERROR != 0:
    print
    print "usage: sensmg.py [-i <inputfile>] [-np <# procs>] [-ngroup <ngroup>] [-isn <isn>] [-isct <isct>]"
    print "         [-epsi <epsi>] [-epsig <epsig>] [-srcacc_no <none, for, adj, or for+adj>]"
    print "         [-fissdata <0, 1, or 2>] [-chinorm <none, full, or partial>]"
    print "         [-aflxfrm <0 or 1>]"
    print "         [-misc <yes or no>] [-alpha_n <yes or no>] [-nag <# alpha groups>]"
    print "         [-use_existing <no or yes>] [-my_modules <no or yes>]"
    print "         [-rplane <rplane>] [-zplane <zplane>]"
    print "  default inputfile is sensmg_inp"
    print "  default np is 1"
    print "  default ngroup is 30"
    print "  default isn is 32"
    print "  default isct is 3"
    print "  default epsi is 1.e-6"
    print "  default epsig is 1.e-5"
    print "  default srcacc_no is none."
    print "  default fissdata is 0."
    print "  default aflxfrm is 0."
    print "  default chinorm is full."
    print "  default misc is yes."
    print "  default alpha_n is yes."
    print "  default nag is 100."
    print "  default USE_EXISTING is no; use yes to use existing output files."
    print "  default MY_MODULES is no; use yes to skip automatic module loading (use your own);"
    print "          also use yes if your system does not use modules."
    print "  default rplane is -1 (i.e. to not use this feature)."
    print "  default zplane is -1 (i.e. to not use this feature)."
    log.write( "usage: sensmg.py [-i <inputfile>] [-np <# procs>] [-ngroup <ngroup>] [-isn <isn>] [-isct <isct>]\n")
    log.write( "         [-epsi <epsi>] [-epsig <epsig>] [-srcacc_no <none, for, adj, or for+adj>]\n")
    log.write( "         [-fissdata <0, 1, or 2>] [-chinorm <none, full, or partial>]\n")
    log.write( "         [-aflxfrm <0 or 1>]\n")
    log.write( "         [-misc <yes or no>] [-alpha_n <yes or no>] [-nag <# alpha groups>]\n")
    log.write( "         [-use_existing <no or yes>] [-my_modules <no or yes>]\n")
    log.write( "         [-rplane <rplane>] [-zplane <zplane>]\n")
    log.write( "  default inputfile is sensmg_inp\n")
    log.write( "  default np is 1\n")
    log.write( "  default ngroup is 30\n")
    log.write( "  default isn is 32\n")
    log.write( "  default isct is 3\n")
    log.write( "  default epsi is 1.e-6\n")
    log.write( "  default epsig is 1.e-5\n")
    log.write( "  default srcacc_no is none.\n")
    log.write( "  default fissdata is 0.\n")
    log.write( "  default aflxfrm is 0.\n")
    log.write( "  default chinorm is full.\n")
    log.write( "  default misc is yes.\n")
    log.write( "  default alpha_n is yes.\n")
    log.write( "  default nag is 100.\n")
    log.write( "  default USE_EXISTING is no; use yes to use existing output files.\n")
    log.write( "  default MY_MODULES is no; use yes to skip automatic module loading (use your own);\n")
    log.write( "          also use yes if your system does not use modules.\n")
    log.write( "  default rplane is -1 (i.e. to not use this feature.)\n")
    log.write( "  default zplane is -1 (i.e. to not use this feature.)\n")
    es=exit_sens(2)
if IERRORP != 0:
    es=exit_sens(1)

# test input file and get NM, NR, NZ, NRRR, LIBNAME.
if os.path.exists(IFILE) == False:
    print "\nerror. input file not found: "+IFILE
    log.write("\nerror. input file not found: "+IFILE+"\n")
    es=exit_sens(1)
inpf = open(IFILE, "r")
line = inpf.readline() # title
line = inpf.readline() # keywords
if "cyl" in line:
    NDIMS = 2
else:
    NDIMS = 1
IFS = 0
IFEYNY = 0
if "lkg" in line or "feyny" in line or "sm2" in line:
    IFS = 1
    if "feyny" in line or "sm2" in line:
        IFEYNY = 1
line = inpf.readline() # LIBNAME
LIBNAME = line.split()[0]
log.write("  LIBNAME="+LIBNAME+"\n")
line = inpf.readline() # NM
try:
    NM = int(line.split()[0])
except:
    es=exit_sens(3)
n = 0
NCB = []
while n < NM:
    line = inpf.readline() # materials
    line = line.split("/")[0]
    NCB.append( ( len(line.split()) -1)/2 )
    n = n + 1
NEL = sum(NCB)
NCBC = " ".join(map(str, NCB)) # string with the list of integers
line = inpf.readline() # densities
line = inpf.readline() # NR, NZ
try:
    NR = int(line.split()[0])
except:
    es=exit_sens(3)
if NDIMS == 2:
    try:
        NZ = int(line.split()[1])
    except:
        es=exit_sens(3)
else:
    NZ = 1
line = inpf.readline() # radii
if NDIMS == 2:
    line = inpf.readline() # heights
n = 0
while n < NZ:
    line = inpf.readline() # material assignments
    n = n + 1
line = inpf.readline() # coarse mesh for reaction rate
line = inpf.readline() # NRRR
try:
    NRRR = int(line.split()[0])
except:
    es=exit_sens(3)
inpf.close()

# NM is limited to be less than 1000 in this script.
# NRRR is limited to be less than 100 in this script.
IERROR = 0
if NM > 999:
    print "error. number of materials must be less than 1000."
    print "  NM="+str(NM)
    log.write("error. number of materials must be less than 1000.\n")
    log.write("  NM="+str(NM)+"\n")
    IERROR = 1
if NRRR > 99:
    print "error. number of reaction-rate ratios must be less than 100."
    print "  NRRR="+str(NRRR)
    log.write("error. number of reaction-rate ratios must be less than 100.\n")
    log.write("  NRRR="+str(NRRR)+"\n")
    IERROR = 1
if IERROR != 0:
    es=exit_sens(4)

# ensure consistency among LIBNAME and NGROUP.
# don't check the regular NDI libraries as there are many collapse options. see
# https://xweb.lanl.gov/projects/data/nuclear/ndi_data/transport/group_structure.html
IERROR = 0
NG0 = NGROUP
if LIBNAME == "mendf6":
    if NGROUP != 30:
        NGROUP = 30
        IERROR = 1
elif LIBNAME == "kynea3":
    if NGROUP != 79:
        NGROUP = 79
        IERROR = 1
elif LIBNAME == "acti":
    if NGROUP != 130:
        NGROUP = 130
        IERROR = 1
elif LIBNAME == "vitb6":
    if NGROUP != 241:
        NGROUP = 241
        IERROR = 1
elif LIBNAME == "scale":
    if NGROUP != 44:
        NGROUP = 44
        IERROR = 1
if IERROR != 0:
    print "warning. NGROUP and LIBNAME were not consistent; NGROUP was modified."
    print "  old NGROUP="+str(NG0)+"  new NGROUP="+str(NGROUP)+"  LIBNAME="+LIBNAME
    log.write("warning. NGROUP and LIBNAME were not consistent; NGROUP was modified.\n")
    log.write("  old NGROUP="+str(NG0)+"  new NGROUP="+str(NGROUP)+"  LIBNAME="+LIBNAME+"\n")

NDI_GENDIR_PATH = os.environ.get("NDI_GENDIR_PATH")
if NDI_GENDIR_PATH == None:
    print "warning. NDI_GENDIR_PATH not set; using default."
    log.write("warning. NDI_GENDIR_PATH not set; using default.\n")
    NDI_GENDIR_PATH = "/usr/projects/data/nuclear/ndi/2.1.3/share/gendir.all"
    os.environ["NDI_GENDIR_PATH"] = NDI_GENDIR_PATH
log.write("  NDI_GENDIR_PATH="+NDI_GENDIR_PATH+"\n")
sys.stdout.flush()
log.flush()

# if USE_EXISTING, check for for/for_out
if USE_EXISTING == "yes":
    if os.path.exists("for/for_out") == False:
        print "\nerror. use_existing=yes but for/for_out not found."
        log.write("\nerror. use_existing=yes but for/for_out not found.\n")
        es=exit_sens(1)
# if not, remove old directories and recreate.
else:
    for d in [ "for", "adj", "xs1" ]:
        if os.path.exists(d) == True:
            shutil.rmtree(d)
        os.mkdir(d)
    for d in [ "sources", "misc" ]:
        if os.path.exists(d) == True:
            shutil.rmtree(d)
        if IFS == 1:
            os.mkdir(d)
    for d in [ "smf", "sma" ]:
        if os.path.exists(d) == True:
            shutil.rmtree(d)
        if IFEYNY == 1:
            os.mkdir(d)

# add sources4c data.
    if IFS == 1:
        if IMISC == 0 or IALPHAN == 1:
            os.chdir("sources")
            for t in [ "tape2", "tape3", "tape4", "tape5" ]:
                if os.path.lexists(t) == False:
                    os.symlink(sources_dir+"/"+t+".txt", t)
            os.chdir("..")

# if using kynea3, vitb6, or scale, link to data.
    SENS_DATA = os.environ.get("SENS_DATA")
    if LIBNAME == "kynea3":
        for d in [ "for", "xs1" ]:
            lnk = d+"/kynea3"
            if os.path.exists(lnk) == False:
                os.symlink(SENS_DATA+"/sandia.neutronxs.kynea3.bxs", lnk)
    elif LIBNAME == "vitb6":
        for d in [ "for", "xs1" ]:
            lnk = d+"/vitb6"
            if os.path.exists(lnk) == False:
                os.symlink(SENS_DATA+"/vitb6.xslib", lnk)
    elif LIBNAME == "scale":
        for d in [ "for", "xs1" ]:
            lnk = d+"/scale"
            if os.path.exists(lnk) == False:
                os.symlink(SENS_DATA+"/bxslib.scale.44", lnk)
    elif LIBNAME == "special":
        for l in [ "macrxs", "snxedt", "ndxsrf", "znatdn" ]:
            f = "for_"+l
            lnk = "for/"+l
            if os.path.exists(lnk) == False:
                os.symlink(SENS_DATA+"/"+f, lnk)
            f = "xs1_"+l
            lnk = "xs1/"+l
            if os.path.exists(lnk) == False:
                os.symlink(SENS_DATA+"/"+f, lnk)

# remove old files.
    for to_rm in [ "control", "stopconverged" ]:
        if os.path.lexists(to_rm):
            os.remove(to_rm)
# a01-a99 are 1st-LASS adjoint directories for reaction-rate ratios
    to_rm = glob.glob("a[0-9][0-9]")
    for f in to_rm:
        shutil.rmtree(f)
    for i in range (1,NRRR+1):
        d = "a"+"%02d" % (i)
        os.mkdir(d)
sys.stdout.flush()
log.flush()

# remove old files whether or not USE_EXISTING.
for to_rm in [ "stoponerror", "sens_l_x", "sens_k_x", "sens_a_x", "sens_rr_x",
               "sens_l_r", "sens_k_r", "sens_a_r", "sens_rr_r",
               "senslx",  "senssm" ]:
    if os.path.lexists(to_rm):
        os.remove(to_rm)

IERROR = 0
if LIBNAME == "mt71x" and ( "mt71x" not in NDI_GENDIR_PATH and "devel" not in NDI_GENDIR_PATH ):
    version=[2,1,2] # need NDI version greater than or equal to 2.1.2
    for m in NDI_GENDIR_PATH.split("/"):
        if '.' in m:
            for index, n in enumerate(m.split('.')):
                if int(n) > version[index]:
                    break
                elif int(n) == version[index]:
                    continue
                else:
                    IERROR = 1
            break
elif ( LIBNAME != "mt71x" and LIBNAME != "kynea3" and LIBNAME != "vitb6" and LIBNAME != "scale" ) and "mt71x" in NDI_GENDIR_PATH:
    IERROR = 1
elif ( LIBNAME == "mendf71x" ):
    version=[2,0,20] # need NDI version greater than or equal to 2.0.20
    for m in NDI_GENDIR_PATH.split("/"):
        if '.' in m:
            for index, n in enumerate(m.split('.')):
                if int(n) > version[index]:
                    break
                elif int(n) == version[index]:
                    continue
                else:
                    IERROR = 1
            break
if IERROR != 0:
    print "error. wrong NDI_GENDIR_PATH for LIBNAME="+LIBNAME
    print "  NDI_GENDIR_PATH="+NDI_GENDIR_PATH
    log.write("error. wrong NDI_GENDIR_PATH for LIBNAME="+LIBNAME+"\n")
    log.write("  NDI_GENDIR_PATH="+NDI_GENDIR_PATH+"\n")
    es=exit_sens(1)

ITER = -1

# start calculations.
# for fixed-source, read input file, write misc and/or sources4c input files,
# and run misc and/or sources4c.
if IFS == 1:
    wc=write_control(1)
    log.close()
    os.system(sensmg_exe)
    if os.path.exists("stoponerror") == True:
        es=exit_sens(1)
    log = open("sensmg.log", "a")
    if IMISC == 1 and USE_EXISTING == "no":
        print "running misc for each material...."
        log.write("running misc for each material....\n")
        sys.stdout.flush()
        log.flush()
# MY_MODULES is for the partisn modules. this needs to be done for misc.
# this logic should work when not on LANL ICN.
        LOADEDMODULES_pre_misc = os.environ.get("LOADEDMODULES")
        if LOADEDMODULES_pre_misc != None:
            pm=module("purge")
            pm=module("load", "gcc/5.3.0")
        os.chdir("misc")
        mat_files = glob.glob("m[0-9][0-9][0-9]_misc.inp")
        mat_files.sort()
        for n2 in  mat_files:
            mxx = n2[0:4] # m001, m002, m003, etc.
            print "material "+mxx
            inp = mxx+"_misc.inp"
            out = mxx+"_misc.out"
            src = mxx+"_misc.src"
# run misc
            os.system(misc_exe+" "+inp)
            if os.path.exists(out):
                if not os.path.exists(src):
# or use error code 2, which means no sources; how to get exit codes in python?
                    srcf = open(src, "w")
                    srcf.write("no sources\n")
                    srcf.close()
        if LOADEDMODULES_pre_misc != None:
            pm=module("purge")
            for m in LOADEDMODULES_pre_misc.split(":"):
                pm=module("load", m)
#           pm=module("list")
        os.chdir("..")
    if ( IMISC == 0 or IALPHAN == 1 ) and USE_EXISTING == "no":
        print "running sources4c for each material...."
        log.write("running sources4c for each material....\n")
        sys.stdout.flush()
        log.flush()
        os.chdir("sources")
        for to_rm in [ "outp", "outp2", "tape6", "pdata", "tape7", "tape8", "tape9" ]:
            if os.path.lexists(to_rm):
                os.remove(to_rm)
        mat_files = glob.glob("m[0-9][0-9][0-9]_tape1")
        mat_files.sort()
        for n2 in  mat_files:
            mxx = n2[0:4] # m001, m002, m003, etc.
            print "material "+mxx
            if os.path.lexists("tape1"):
                os.remove("tape1")
            os.symlink(n2, "tape1")
# run sources4c
            os.system(sources_exe)
            if os.path.exists("outp2"):
                out=mxx+"_outp2"
                shutil.move("outp2", out)
                out=mxx+"_tape6"
                shutil.move("tape6", out)
                out=mxx+"_tape7"
                shutil.move("tape7", out)
                out=mxx+"_tape8"
                shutil.move("tape8", out)
                out=mxx+"_outp"
                shutil.move("outp", out)
# pdata file is only in Favorite's version
            if os.path.exists("pdata"):
                out=mxx+"_pdata"
                shutil.move("pdata", out)
        os.chdir("..")
    sys.stdout.flush()
    log.flush()

# read input file, sources4c or misc output file (for lkg, feyny, or sm2);
# write forward and cross-section files for partisn.
wc=write_control(2)
log.close()
os.system(sensmg_exe)
if os.path.exists("stoponerror") == True:
    es=exit_sens(1)
log = open("sensmg.log", "a")
# run partisn forward and cross-section files.
for f in [ "for", "xs1" ]:
    if f == "xs1" and LIBNAME == "special":
        print "don't run partisn for "+f+"_inp...."
        log.write("don't run partisn for "+f+"_inp....\n")
        break
    print "running partisn for "+f+"_inp...."
    log.write("running partisn for "+f+"_inp....\n")
    sys.stdout.flush()
    log.flush()
    os.chdir(f)
    inp=f+"_inp"
    out=f+"_out"
    if USE_EXISTING == "no":
        os.system(PARTISN_EXE+" "+inp+" "+out)
    outf = open(out, "r")
    IERRSPHERE = 0
    for line in outf:
        if re.search("tinp213d", line):
            IERRSPHERE = 1
        if re.search("fatal input error", line):
            print "\nerror. fatal input errors for "+inp+"."
            log.write("\nerror. fatal input errors for "+inp+".\n")
            if IERRSPHERE == 1:
                print "for spheres, use \"-np 1\"."
                log.write( "for spheres, use \"-np 1\".\n")
            outf.close()
            es=exit_sens(1)
        elif re.search("not converged", line):
            print "error. partisn not converged for "+inp+"."
            log.write("error. partisn not converged for "+inp+".\n")
        elif re.search("alpha too negative", line):
            print "error. alpha for "+inp+" is too negative for partisn."
            log.write("error. alpha for "+inp+" is too negative for partisn.\n")
            outf.close()
            es=exit_sens(1)
        elif re.search("processor mesh .gt. problem mesh", line):
            print "error. bad number of processors for "+inp+"."
            log.write("error. bad number of processors for "+inp+".\n")
            outf.close()
            es=exit_sens(1)
        elif re.search(" nocore", line):
            print "error. "+inp+" is too large for partisn."
            log.write("error. "+inp+" is too large for partisn.\n")
            outf.close()
            es=exit_sens(1)
    outf.close()
# clean up unused partisn files.
    for p in [ "altinp", "asgmat", "editit", "fissrc", "geodst", "redoin",
               "rtflux", "sncons", "solinp", "summry" ]:
        if os.path.exists(p) == True:
            os.remove(p)
    os.chdir("..")
    sys.stdout.flush()
    log.flush()

# read partisn forward and cross-section files.
ITER = ITER + 1
wc=write_control(3)
log.close()
os.system(sensmg_exe)
if os.path.exists("stoponerror") == True:
    es=exit_sens(1)
log = open("sensmg.log", "a")

# run partisn regular adjoint file.
for f in [ "adj" ]:
    print "running partisn for "+f+"_inp...."
    log.write("running partisn for "+f+"_inp....\n")
    sys.stdout.flush()
    log.flush()
    os.chdir(f)
    inp=f+"_inp"
    out=f+"_out"
    for l in [ "macrxs", "snxedt", "ndxsrf", "znatdn" ]:
        if os.path.exists(l) == False:
            os.symlink("../for/"+l, l)
    if USE_EXISTING == "no":
        os.system(PARTISN_EXE+" "+inp+" "+out)
    outf = open(out, "r")
    for line in outf:
        if re.search("not converged", line):
            print "error. partisn not converged for "+inp+"."
            log.write("error. partisn not converged for "+inp+".\n")
        elif re.search("alpha too negative", line):
            print "error. alpha for "+inp+" is too negative for partisn."
            log.write("error. alpha for "+inp+" is too negative for partisn.\n")
            outf.close()
            es=exit_sens(1)
        elif re.search(" nocore", line):
            print "error. "+inp+" is too large for partisn."
            log.write("error. "+inp+" is too large for partisn.\n")
            outf.close()
            es=exit_sens(1)
    outf.close()
# clean up unused partisn files.
    for p in [ "adjmac", "altinp", "asgmat", "atflux", "editit",
               "fissrc", "geodst", "redoin", "sncons", "solinp", "summry" ]:
        if os.path.exists(p) == True:
            os.remove(p)
    os.chdir("..")
    sys.stdout.flush()
    log.flush()

# if there are no reaction rates or if USE_EXISTING is yes.
# after this, will skip the generalized adjoint section.
if os.path.exists("stopconverged") == True or USE_EXISTING == "yes":
    wc=write_control(4)
    log.close()
    os.system(sensmg_exe)
    if os.path.exists("stoponerror") == True:
        es=exit_sens(1)
    log = open("sensmg.log", "a")

# run partisn generalized adjoint files.
NADJ = 0
while os.path.exists("stopconverged") == False:
    if ITER < 100:
        i0 = "%02d" % (ITER)
    else:
        print "iterations reached 100."
        log.write("iterations reached 100.\n")
        break
    adj_files = glob.glob("a[0-9][0-9]")
    adj_files.sort()
    if NADJ == 0:
        NADJ = len(adj_files)
        IADJCONV = range(NADJ)
        IADJCONV = [0] * NADJ
    NA = -1
    for axx in  adj_files:
        NA = NA + 1
        os.chdir(axx)
        for l in [ "macrxs", "snxedt", "ndxsrf", "znatdn" ]:
            if os.path.exists(l) == False:
                os.symlink("../for/"+l, l)
        src = axx+"_"+i0+"_fixsrc"
        inp = axx+"_"+i0+"_inp"
        out = axx+"_"+i0+"_out"
        if IADJCONV[NA] == 0:
            shutil.move(axx+"_inp", inp)
            inpf = open(inp, "r")
            for line in inpf:
                if re.search("generalized adjoint converged", line):
                    IADJCONV[NA] = 1
            inpf.close()
        if IADJCONV[NA] == 1:
            print "generalized adjoint converged for "+inp+"."
            log.write("generalized adjoint converged for "+inp+".\n")
            os.chdir("..")
        else:
            shutil.move(axx+"_fixsrc", src)
            os.symlink(src, "fixsrc")
            print "running partisn for "+inp+"...."
            log.write("running partisn for "+inp+"....\n")
            sys.stdout.flush()
            log.flush()
            if USE_EXISTING == "no":
                os.system(PARTISN_EXE+" "+inp+" "+out)
            outf = open(out, "r")
            for line in outf:
                if re.search("not converged", line):
                    print "error. partisn not converged for "+inp+"."
                    log.write("error. partisn not converged for "+inp+".\n")
                elif re.search(" nocore", line):
                    print "error. "+inp+" is too large for partisn."
                    log.write("error. "+inp+" is too large for partisn.\n")
                    outf.close()
                    es=exit_sens(1)
            outf.close()
# clean up unused partisn files.
            for p in [ "adjmac", "altinp", "asgmat", "atflux", "editit", "fissrc", "fixsrc",
                       "geodst", "redoin", "rtflux", "sncons", "solinp", "summry" ]:
                if os.path.exists(p) == True:
                    os.remove(p)
            os.chdir("..")
        sys.stdout.flush()
        log.flush()

    ITER = ITER + 1
    wc=write_control(3)
    log.close()
    os.system(sensmg_exe)
    if os.path.exists("stoponerror") == True:
        es=exit_sens(1)
    log = open("sensmg.log", "a")

# feyny or sm2 sensitivities. run partisn.
if IFEYNY == 1:
    for d in [ "smf", "sma" ]:
        os.chdir(d)
        for l in [ "macrxs", "snxedt", "ndxsrf", "znatdn" ]:
            if os.path.exists(l) == False:
                os.symlink("../for/"+l, l)
        o2_files = glob.glob("[0-9][0-9]_inp")
        o2_files.sort()
        for n2 in o2_files:
            inp = n2
            out = n2[0:2]+"_out"
            src = n2[0:2]+"_fixsrc"
            if os.path.exists(src):
                os.symlink(src, "fixsrc")
            print "running partisn for "+d+"/"+inp+"...."
            log.write("running partisn for "+d+"/"+inp+"....\n")
            sys.stdout.flush()
            log.flush()
            if USE_EXISTING == "no":
                os.system(PARTISN_EXE+" "+inp+" "+out)
            for l in [ "rmflux", "amflux", "raflxm", "aaflxm", "asfluxx", "asfluxy" ]:
                if os.path.exists(l) == True:
                    shutil.move(l, n2[0:2]+"_"+l)
            if os.path.exists("fixsrc") == True:
                os.remove("fixsrc")
            outf = open(out, "r")
            for line in outf:
                if re.search("not converged", line):
                    print "error. partisn not converged for "+inp+"."
                    log.write("error. partisn not converged for "+inp+".\n")
                    log.flush()
                elif re.search(" nocore", line):
                    print "error. "+inp+" is too large for partisn."
                    log.write("error. "+inp+" is too large for partisn.\n")
                    outf.close()
                    es=exit_sens(1)
            outf.close()
# clean up unused partisn files.
        for p in [ "aaflxm", "adjmac", "altinp", "amflux", "asgmat", "atflux", "editit", "fissrc",
                   "fixsrc", "geodst", "raflxm", "redoin", "rmflux", "rtflux", "sncons", "solinp",
                   "summry", "asfluxx", "asfluxy" ]:
            if os.path.exists(p) == True:
                os.remove(p)
        os.chdir("..")
        sys.stdout.flush()
        log.flush()

    wc=write_control(5)
    log.close()
    os.system(sensmg_exe)
    if os.path.exists("stoponerror") == True:
        es=exit_sens(1)
    log = open("sensmg.log", "a")

print "end of sensmg script"
log.write("end of sensmg script\n")
es=exit_sens(0)

