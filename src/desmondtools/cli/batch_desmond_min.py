import os
import sys
import argparse
import random
import re

from desmondtools import Multisim


def batch_min():
    SCHRODINGER = os.environ["SCHRODINGER"]
    multisim = "{}/utilities/multisim".format(SCHRODINGER)
    
    min_msj = Multisim(template="desmond-min.msj")
    min_cfg = Multisim(template="desmond-min.cfg")

    opt = '-HOST localhost -maxjob 1 -cpu 1 -mode umbrella '
    opt += '-lic "DESMOND_GPGPU:16" '
    opt += '-description "minimization"'

    parser = argparse.ArgumentParser(description="batch gdesmond minimization jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', dest="gpu_device", type=int, default=0, 
                        help="gpu device id")
    parser.add_argument('-t', dest="simulation_time", type=float, default=100.0, 
                        help="simulation time in ps")
    parser.add_argument('-p', dest="prefix", type=str, default="r", 
                        help="directory prefix")
    parser.add_argument('-s', dest="start", type=int, default=1, 
                        help="directory start")
    parser.add_argument('-r', dest="repeat", type=int, default=1, 
                        help="number of repeats")
    parser.add_argument('-j', dest="job_file", type=str, default="desmond_min_job_1.sh", 
                        help="job filename")
    parser.add_argument('cms', nargs="+", help="desmond cms file")
    args = parser.parse_args()

    try:
        cms_files = [os.path.abspath(f) for f in args.cms]
        assert(len(cms_files) > 0)
    except:
        print(".cms file(s) not found")
        sys.exit(0)


    job_file = args.job_file
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"


    with open("README","a") as readme, open(job_file,"w") as job:
        print("\n" + job_file + "\n")
        outdir_nums = list(range(args.start, args.start+args.repeat))
        outdirs = [f"{args.prefix}{num:02d}" for num in outdir_nums]
        cmd_echo = ""
        for argv in sys.argv:
            if cmd_echo:
                cmd_echo += " "
            if " " in argv:
                cmd_echo += f'"{argv}"'
            else:
                cmd_echo += f'{argv}'
        readme.write(f"{cmd_echo}\n\n")
        readme.write(f"GPU device              = {args.gpu_device}\n")
        readme.write(f"Simulation Time (ns)    = {args.simulation_time}\n")
        readme.write(f"Repeat                  = {args.repeat}\n")
        readme.write( "Directory               = %s\n" % " ".join(outdirs))
        readme.write(f"Jobfile                 = {job_file}\n\n")

        job.write(f'export CUDA_VISIBLE_DEVICES="{args.gpu_device}"\n\n')

        for i, infile in enumerate(cms_files):
            info = f"[{i+1}] {infile}"
            print(info)
            readme.write(info+"\n")
        print()
        readme.write("\n")

        for (num, outdir) in zip(outdir_nums, outdirs):
            outdir_abspath = os.path.abspath(outdir)
            job.write(f"cd {outdir_abspath}/\n\n")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            cfg_file = f"{outdir}/desmond_min_job_{num:02d}.cfg"
            msj_file = f"{outdir}/desmond_min_job_{num:02d}.msj"
            cfg_file_basename= os.path.basename(cfg_file)

            with open(cfg_file,"w") as cfg, open(msj_file,"w") as msj:
                min_cfg.dot.randomize_velocity.seed = str(random.randint(1000,9999))
                min_cfg.dot.time = str(args.simulation_time) # (ps)
                min_cfg.write(cfg)

                min_msj.dot.simulate.cfg_file = cfg_file_basename
                min_msj.write(msj)
            
            # append to job file
            for infile in cms_files:
                prefix_ = re.sub(r'desmond_setup[-_]','', os.path.basename(infile))
                prefix  = prefix_.replace("-out.cms","")
                job_name = f"desmond_min_job_{num:02d}_{prefix}"
                job.write('{} -JOBNAME {} -m {} -c {} {} {} -o {} -WAIT\n\n'.format(
                    multisim, 
                    job_name,
                    os.path.basename(msj_file),
                    os.path.basename(cfg_file),
                    opt,
                    infile,
                    f"{job_name}-out.cms",
                ))

    os.chmod(job_file, 0o777)