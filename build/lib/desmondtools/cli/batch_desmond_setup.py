import os
import sys
import argparse
from desmondtools import Multisim


def batch_setup():
    SCHRODINGER = os.environ["SCHRODINGER"]
    multisim = "{}/utilities/multisim".format(SCHRODINGER)
    setup_msj = Multisim(template="desmond-setup.msj")
    
    parser = argparse.ArgumentParser(description="batch gdesmond md system setup",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--conc', dest="conc", type=float, default=0.15, help="salt concentration in M")
    parser.add_argument('-d','--dist', dest="dist", type=float, default=10.0, help="buffer distance in A")
    parser.add_argument('-l','--lipid',  dest="lipid", default="", help="lipid bilayer")
    parser.add_argument('--cpp', dest="cpp", default=False, action="store_true", help="CPP simulation setup")
    parser.add_argument('-s','--solvent',  dest="solvent", default="TIP3P", help="solvent model")
    parser.add_argument('-i','--counterion',  dest="counterion", default="Na", help="neutralizing ion")
    parser.add_argument('-n','--negative', dest="neg", default="Cl", help="negative salt ion")
    parser.add_argument('-p','--positive', dest="pos", default="Na", help="positive salt ion")
    parser.add_argument('-f','--forcefield', dest="forcefield", default="S-OPLS", help="forcefield")
    parser.add_argument('-j','--jobfile', dest="job_file", default="desmond_setup_1.sh", help="job filename")
    parser.add_argument('-m','--msjfile', dest="msj_file", default="desmond_setup_1.msj", help="msj filename")
    parser.add_argument('-a','--appendix', dest="appendix", default="", help="job name appendix")
    parser.add_argument('mae', nargs="+", help="desmond mae file")
    args = parser.parse_args()

    if args.appendix:
        msj_file = args.msj_file[:-4] + "_" + args.appendix + ".msj"
    else:
        msj_file = args.msj_file

    if args.appendix:
        job_file = args.job_file[:-3] + "_" + args.appendix + ".sh"
    else:
        job_file = args.job_file

    # job file (.sh) and msj file (.msj) should match
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"

    while os.path.exists(msj_file):
        splited = msj_file.replace(".msj","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        msj_file = "_".join(splited) + ".msj"

    with open(msj_file, "w") as msj:
        setup_msj.dot.build_geometry.solvent = str(args.solvent)
        setup_msj.dot.build_geometry.add_counterion.ion = str(args.counterion)
        setup_msj.dot.build_geometry.salt.concentration = str(args.conc)
        setup_msj.dot.build_geometry.box.size = f"[ {args.dist} {args.dist} {args.dist} ]"
        setup_msj.dot.build_geometry.salt.negative_ion = str(args.neg)
        setup_msj.dot.build_geometry.salt.positive_ion = str(args.pos)
        if args.cpp:
            setup_msj.dot.build_geometry.rezero_system = "False"
            args.lipid = "POPC"
        if args.lipid:
            setup_msj.dot.build_geometry.membrane_box.lipid = str(args.lipid)
        else:
            # remove membrane_box block
            setup_msj.dot.build_geometry.pop("membrane_box")
        setup_msj.dot.build_geometry.override_forcefield = str(args.forcefield)
        setup_msj.dot.assign_forcefield.forcefield = str(args.forcefield)
        setup_msj.dot.assign_forcefield.water = str(args.solvent)
        setup_msj.write(msj)


    with open("README","a") as readme, open(job_file,"w") as job:
        cmd_echo = ""
        for argv in sys.argv:
            if cmd_echo:
                cmd_echo += " "
            if " " in argv:
                cmd_echo += f'"{argv}"'
            else:
                cmd_echo += f'{argv}'
        readme.write(f"{cmd_echo}\n\n")
        readme.write(f"Force Field   = {args.forcefield}\n")
        readme.write(f"Solvent       = {args.solvent}\n")
        readme.write(f"Counter Ion   = {args.counterion}\n")
        readme.write(f"Positive Ion  = {args.pos}\n")
        readme.write(f"Negative Ion  = {args.neg}\n")
        readme.write(f"Concentration = {args.conc} (M)\n")
        readme.write(f"Size          = {args.dist} (A)\n")
        readme.write(f"Lipid         = {args.lipid}\n")
        readme.write(f"msjfile       = {msj_file}\n")
        readme.write(f"Jobfile       = {job_file}\n")
        readme.write(f"Input structure(s):\n")

        for i, infile in enumerate(args.mae):
            prefix = os.path.basename(infile).split(".")[0]
            job_name = f"desmond_setup-{prefix}"
            if args.appendix:
                job_name += f"-{args.appendix}"
            cms_file = f"{job_name}-out.cms"
            job.write(f"if [ ! -f {cms_file} ]\n")
            job.write(f"then\n")
            job.write(f"{multisim} \\\n")
            job.write(f"  -JOBNAME {job_name} \\\n")
            job.write(f"  -m {msj_file} {os.path.abspath(infile)} \\\n") 
            job.write(f"  -o {cms_file} \\\n")
            job.write(f"  -HOST localhost:20 -maxjob 20 -WAIT\n")
            job.write(f"fi\n")
            readme.write(f"[{i+1:02d}] {infile}\n")
        readme.write("\n\n")

    os.chmod(job_file, 0o777)