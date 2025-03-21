import os
import sys
import argparse
import random
import re
import copy

from desmondtools import Multisim



def index_or_asl(expr:str) -> str:
    if len(expr.split()) == 1:
        return expr
    else:
        return f'"{expr}"'
    

def batch_metad():
    """Metadynamics Simulations.
    """

    SCHRODINGER = os.environ["SCHRODINGER"]
    multisim = "{}/utilities/multisim".format(SCHRODINGER)

    metad_msj = Multisim(template="desmond-metad.msj")
    metad_cfg = Multisim(template="desmond-metad.cfg")

    opt = '-HOST localhost -maxjob 1 -cpu 1 -mode umbrella -lic "DESMOND_GPGPU:16" '
    opt += '-description "metadynamics"'

    parser = argparse.ArgumentParser(description="batch gdesmond metadynamics jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--meta-height', dest="meta_height", type=float, default=0.03, help="Height of the repulsive Gaussian potential in kcal/mol (0.03)")
    parser.add_argument('--meta-interval', dest="meta_interval", type=float, default=0.09, help="Interval in ps at which Gaussian potentials are added (0.09)")
    parser.add_argument('--meta-first', dest="meta_first", type=float, default=0.0, help="Time in ps at which the Gaussian potentials are first added (0)")
    parser.add_argument('--meta-last', dest="meta_last", type=float, default=-1.0, help="Time in ps at which the Gaussian potentials are last added (simulation time)")
    parser.add_argument('--meta-kTemp', dest="meta_kTemp", type=float, default=2.4, help="Perform the well-tempered metadynamics (2.4)")
    parser.add_argument('--cv-dist', dest='cv_dist', nargs=4, action='append', default=[], help="atom1 atom2 width(0.05) wall")
    parser.add_argument('--cv-angle', dest='cv_angle', nargs=4, action='append', default=[], help="atom1 atom2 atom3 width(2.5)")
    parser.add_argument('--cv-dihedral', dest='cv_dihedral', nargs=5, action='append', default=[], help="atom1 atom2 atom3 atom4 width(5.0)")
    parser.add_argument('--cv-rgyr', dest='cv_rgyr', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-rgyr-mass', dest='cv_rgyr_mass', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-rmsd', dest='cv_rmsd', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-rmsd-symm', dest='cv_rmsd_symm', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-zdist', dest='cv_zdist', nargs=3, action='append', default=[], help="atom width wall(0.05)")
    parser.add_argument('--cv-zdist0', dest='cv_zdist0', nargs=3, action='append', default=[], help="atom width wall(0.1)")
    parser.add_argument('-g', dest="gpu_device", type=int, default=0, help="gpu device id")
    parser.add_argument('-T', dest="temperature", type=float, default=300.0, help="temperature in Kelvin")
    parser.add_argument('-t', dest="simulation_time", type=float, default=40.0, help="simulation time in ns")
    parser.add_argument('-i', dest="interval", type=float, default=40.0, help="frame interval in ps")
    parser.add_argument('-p', dest="prefix", type=str, default="m", help="directory prefix")
    parser.add_argument('-s', dest="start", type=int, default=1, help="directory start")
    parser.add_argument('-r', dest="repeat", type=int, default=1, help="number of repeats")
    parser.add_argument('-j', dest="job_file", type=str, default="desmond_metadynamics_job_1.sh", help="job filename")
    parser.add_argument('cms', nargs="+", help="desmond cms file")
    args = parser.parse_args()

    try:
        cms_files = [os.path.abspath(f) for f in args.cms]
        assert(len(cms_files) > 0)
    except:
        print(".cms file(s) not found")
        sys.exit(0)

    if not args.interval:
        args.interval = args.simulation_time # it will save 1000 frames

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
        readme.write(f"Temperature (K)         = {args.temperature}\n")
        readme.write(f"Simulation Time (ns)    = {args.simulation_time}\n")
        readme.write(f"Trajectory interval(ps) = {args.interval}\n")
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
            cfg_file = f"{outdir}/desmond_metadynamics_job_{num:02d}.cfg"
            msj_file = f"{outdir}/desmond_metadynamics_job_{num:02d}.msj"
            cfg_file_basename= os.path.basename(cfg_file)

            with open(cfg_file,"w") as cfg, open(msj_file,"w") as msj:
                metad_cfg.dot.randomize_velocity.seed = str(random.randint(1000,9999))
                metad_cfg.dot.time = str(args.simulation_time*1000.0)
                metad_cfg.dot.trajectory.interval = str(args.interval)
                metad_cfg.dot.temperature = f"[ [{args.temperature} 0] ]"
                metad_cfg.write(cfg)

                metad_msj.dot.simulate[-1].cfg_file = cfg_file_basename
                
                # metadynamics
                cv_list = []


                for atom1, atom2, width, wall in args.cv_dist:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    atom2 = index_or_asl(atom2)
                    cv.atom = f'[{atom1} {atom2}]'
                    cv.type = 'dist'
                    cv.width = str(width)
                    if not wall.startswith('-'):
                        cv.wall = str(wall)
                    cv_list.append(cv)
                
                for atom1, atom2, atom3, width in args.cv_angle:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    atom2 = index_or_asl(atom2)
                    atom3 = index_or_asl(atom3)
                    cv.atom = f'[{atom1} {atom2} {atom3}]'
                    cv.type = 'angle'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, atom2, atom3, atom4, width in args.cv_dihedral:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    atom2 = index_or_asl(atom2)
                    atom3 = index_or_asl(atom3)
                    atom4 = index_or_asl(atom4)
                    cv.atom = f'[{atom1} {atom2} {atom3} {atom4}]'
                    cv.type = 'dihedral'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rmsd:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rmsd'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rmsd_symm:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rmsd_symm'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rgyr:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rgyr'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rgyr_mass:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rgyr_mass'
                    cv.width = str(width)
                    cv_list.append(cv)

                for atom1, width, wall in args.cv_zdist:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'zdist'
                    cv.width = str(width)
                    if not wall.startswith('-'):
                        cv.wall = str(wall)
                    cv_list.append(cv)
                
                for atom1, width, wall in args.cv_zdist0:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'zdist0'
                    cv.width = str(width)
                    if not wall.startswith('-'):
                        cv.wall = str(wall)
                    cv_list.append(cv)

                metad_msj.dot.simulate[-1].meta.cv = cv_list
                metad_msj.dot.simulate[-1].meta.height = str(args.meta_height)
                metad_msj.dot.simulate[-1].meta.interval = str(args.meta_interval)
                metad_msj.dot.simulate[-1].meta.first = str(args.meta_first)

                if args.meta_last > 0.0:
                    metad_msj.dot.simulate[-1].meta.last = str(args.meta_last)
                
                # well-tempered metadynamics
                if args.meta_kTemp > 0.0:
                    metad_msj.dot.simulate[-1].meta.kTemp = str(args.meta_kTemp)
                
                metad_msj.write(msj)
            
            # append to job file
            for infile in cms_files:
                prefix_ = re.sub(r'desmond_setup[-_]','', os.path.basename(infile))
                prefix  = prefix_.replace("-out.cms","")
                job_name = f"desmond_metadynamics_job_{num:02d}_{prefix}"
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