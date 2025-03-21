import os
import sys
import argparse
import random
import re

from desmondtools import Multisim


def batch_md():
    SCHRODINGER = os.environ["SCHRODINGER"]
    multisim = "{}/utilities/multisim".format(SCHRODINGER)
    
    md_msj = Multisim(template="desmond-md.msj")
    md_cfg = Multisim(template="desmond-md.cfg")

    parser = argparse.ArgumentParser(description="batch gdesmond md jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', dest="gpu_device", type=int, default=0, 
                        metavar="gpu_device", help="gpu device id")
    parser.add_argument('-T', dest="temp", nargs="+", type=float, default=[300.0,], 
                        metavar="temperature", help="temperature in K")
    parser.add_argument('-R', dest="ramp", nargs="+", type=float, default=[10.0,],
                        metavar="tempramp", help="heat and cool ramps in ps/K")
    parser.add_argument('-t', dest="simulation_time", nargs="+", type=float, default=[100.0,],
                        metavar="simulation_time", help="simulation time in ns")
    parser.add_argument('-i', dest="interval", type=float, default=100.0,
                        metavar="interval", help="frame interval in ps")
    parser.add_argument('--posres', dest="posres", default="res. UNK",
                        metavar="posres", help="ASL for positional restraints during prod. sim.")
    parser.add_argument('--posres-force', dest="posres_force", type=float, default=0.0,
                        metavar="posres_force", help="positional restraints force constant (kcal/mol/A**2)")
    parser.add_argument('-p', dest="prefix", default="r", help="directory prefix")
    parser.add_argument('-s', dest="start", type=int, default=1, help="directory start")
    parser.add_argument('-r', dest="repeat", type=int, default=1, help="number of repeats")
    parser.add_argument('-j', dest="job_file", default="desmond_md_job_1.sh", help="job filename")
    parser.add_argument('cms', nargs="+", help="desmond cms file")
    args = parser.parse_args()

    try:
        cms_files = [os.path.abspath(f) for f in args.cms]
        assert(len(cms_files) > 0)
    except:
        print(".cms file(s) not found")
        sys.exit(0)

    opt  = '-HOST localhost -maxjob 1 -cpu 1 -mode umbrella '
    opt += '-set stage[1].set_family.md.jlaunch_opt=["-gpu"] -lic "DESMOND_GPGPU:16"'

    job_file = args.job_file
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"

    if len(args.temp) > 1:
        try:
            assert len(args.temp) == len(args.simulation_time)
            assert len(args.temp) == (len(args.ramp) + 1)
        except:
            print("For a variable temperature simulaton, the number of temperatures and simulation times ")
            print("should match and temperature ramp(s) (ps/K) should be given between temperatures.")
            print("Please check -T, -t and -R options")
            sys.exit(0)
    else:
        args.ramp = [] # not used

    # Note: if not specified, 
    # all times are in the unit of ps and 
    # energy is in the unit of kcal/mol.

    with open("README","a") as readme, open(job_file,"w") as job:
        print(f"Job file = {job_file}")

        dirs = [ f"{args.prefix}{n:02d}" for n in range(args.start, args.start+args.repeat) ]

        t_schedule = [ [ args.temp[0], 0.0 ], ]

        if len(args.temp) > 1 :

            elapsed = 0
            prev_temp = args.temp[0]
            idx = 0

            for (temp,simulation_time) in zip(args.temp, args.simulation_time):


                deltaT = abs(temp - prev_temp)

                if deltaT > 0.001: # ramp
                    elapsed += args.ramp[idx] * deltaT # ramp: ps/K
                    t_schedule.append([temp, elapsed])
                    idx += 1

                elapsed += simulation_time * 1000. # ns -> ps
                t_schedule.append([temp, elapsed]) # ns -> ps
                prev_temp = temp
            total_simulation_time = elapsed # ps
        else:
            total_simulation_time = sum(args.simulation_time) * 1000.0 # ps

        if not args.interval:
            # args.simulation_time in ns and args.interval in ps
            # default: make 1000 frames
            args.interval = total_simulation_time / 1000.
        
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
        readme.write(f"Temperature (K)         = {args.temp}\n")
        readme.write(f"Temperature Ramp (ps/K) = {args.ramp}\n")
        readme.write(f"Simulation Time (ns)    = {args.simulation_time}\n")
        readme.write(f"Temperature schedule    = {str(t_schedule)}\n")
        readme.write(f"Total Sim. Time (ns)    = {total_simulation_time/1000.0}\n")
        readme.write(f"Trajectory interval(ps) = {args.interval}\n")
        readme.write(f"Repeat                  = {args.repeat}\n")
        readme.write( "Directory               = %s\n" % " ".join(dirs))
        readme.write(f"Jobfile                 = {job_file}\n\n")
        
        for i, infile in enumerate(cms_files):
            info = f"[{i+1:02d}] {infile}"
            print(info)
            readme.write(info+"\n")
        readme.write("\n\n")

        job.write(f'export CUDA_VISIBLE_DEVICES="{args.gpu_device}"\n\n')

        for n in range(args.start, args.start+args.repeat): 
            outdir = f"{args.prefix}{n:02d}"
            outdir_abspath = os.path.abspath(outdir)
            job.write(f"cd {outdir_abspath}/\n\n")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            msj_file = f"{outdir}/desmond_md_job_{n:02d}.msj"
            cfg_file = f"{outdir}/desmond_md_job_{n:02d}.cfg"
            cfg_file_basename= os.path.basename(cfg_file)
            
            with open(msj_file,"w") as msj:
                # modify desmond msj template
                md_msj.dot.simulate[-1].cfg_file = cfg_file_basename
                
                # Setting up restraints using the restraints keyword:
                # https://www.schrodinger.com/kb/332119
                if args.posres_force > 0.0:
                    # print the restraints in the multisim log file
                    md_msj.dot.simulate[-1].print_restraint = 'true'

                    # add the new terms defined in "restraints.new" to existing restraints.
                    # The default is restraints.existing = ignore which will 
                    # delete existing terms before adding any new ones.
                    # md_msj.dot.simulate[-1].restraints.existing = 'retain'

                    md_msj.dot.simulate[-1].restraints.new = [
                        {
                            'name'              : 'posre_harm',
                            'atoms'             : [ f'"{args.posres}"' ],
                            'force_constants'   : [ args.posres_force, ] * 3,
                        }
                        ]
                    # force constants in the x, y, and z direction

                md_msj.write(msj)

            with open(cfg_file,"w") as cfg:
                # modify desmond cfg template
                md_cfg.dot.randomize_velocity.seed = random.randint(1000, 9999)
                md_cfg.dot.time = total_simulation_time
                md_cfg.dot.temperature = t_schedule
                md_cfg.dot.trajectory.interval = args.interval
                md_cfg.write(cfg)

            for infile in cms_files:
                prefix_ = re.sub(r'desmond_setup[-_]','', os.path.basename(infile))
                prefix  = prefix_.replace("-out.cms","")
                
                job_name = f"desmond_md_job_{n:02d}_{prefix}"
                job.write('if [ ! -f {}/{} ]\n'.format( outdir_abspath, f"{job_name}-out.cms",))
                job.write('then\n')
                job.write('{} -JOBNAME {} -m {} -c {} -description "{}" {} {} -o {} -WAIT\n'.format(
                    multisim, 
                    job_name, 
                    os.path.basename(msj_file),
                    os.path.basename(cfg_file),
                    "GPU desmond MD",
                    opt,
                    os.path.join("..",infile),
                    f"{job_name}-out.cms",
                ))
                job.write('fi\n\n')

    os.chmod(job_file, 0o777)