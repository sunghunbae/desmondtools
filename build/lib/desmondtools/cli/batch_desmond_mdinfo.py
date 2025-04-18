import os
import re
import sys
import glob
import argparse

from datetime import timedelta, datetime    


def batch_mdinfo():
    SCHRODINGER = os.environ["SCHRODINGER"]
    USER=os.environ["USER"]
    schrodinger_hosts = "{}/schrodinger.hosts".format(SCHRODINGER)

    simulation_time = re.compile(r'last_time = "(?P<t>[.0-9]+)"')
    #    last_time = "100000.0"
    progress = re.compile(r'Chemical time:\s+(?P<finished>[.0-9]+) ps, Step:\s+[0-9]+, ns/day:\s+(?P<rate>[.0-9]+)')
    #Chemical time:         33507.6000 ps, Step: 5584600, ns/day:      296.231

    parser = argparse.ArgumentParser(description="desmond MD info.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tmpdir', dest='tmpdir', default=None, help='desmond temporary directory')
    parser.add_argument('-l', '--logfile', dest='logfile', default=None)
    args = parser.parse_args()

    if args.tmpdir is None:
        with open(schrodinger_hosts, "r") as f:
            for line in f:
                if line.strip().startswith("#"):
                    continue
                c = line.strip().split(":")
                if len(c) != 2: 
                    continue
                k = c[0].strip()
                v = c[1].strip()
                if k == 'tmpdir':
                    args.tmpdir = os.path.join(v, USER)
            if not os.path.exists(args.tmpdir):
                sys.stderr.write("Use -t or --tmpdir to give temporary directory\n")
                sys.exit(1)


    if args.logfile is None:
        args.logfile = args.tmpdir+"/desmond_*/*.log"
        logfile = None
        for f in glob.glob(args.logfile):
            if not "multisim" in f:
                logfile = f
                break
        if not logfile:
            sys.stderr.write(f"desmond log not found in {args.logfile}\n")
            sys.stderr.write("Use -l or --logfile to give log file path\n")
            sys.stderr.write("Use -t or --tmpdir to give temporary directory\n")
            sys.exit(2)


    total = None
    rate = []
    with open(logfile,"r") as f:
        for line in f:
            m = simulation_time.search(line)
            n = progress.search(line)
            if not total and m and m.group("t"):
                total = float(m.group("t"))*0.001
            if n and n.group("finished") and n.group("rate"):
                rate.append(float(n.group("rate")))
                finished = float(n.group("finished"))*0.001

        n = len(rate)
        print("-"*80)
        print(f"| SCHRODINGER= {SCHRODINGER}")
        print(f"| tmpdir: {args.tmpdir}")
        print(f"|")
        print(f"| Desmond MD Timing Info")
        print(f"| ----------------------")
        print(f"|")
        print(f"| {os.path.dirname(logfile)}/")
        print(f"| {os.path.basename(logfile)}")
        print(f"|")

        if n == 0:
            if total:
                print(f"| Total     {total:9.2f} ns")
            print(f"| Timing data not available yet")
            print("-"*80)
            sys.exit(0)

        remaining = total-finished
        avg_rate = sum(rate[-n:])/n
        eta = 24.0*remaining/avg_rate # hours
        eta_time = datetime.now() + timedelta(hours=eta)
        print(f"| Total     {total:9.2f} ns")
        print(f"| Completed {finished:9.2f} ns ({100.0*finished/total:5.1f}%)")
        print(f"| Remaining {remaining:9.2f} ns")
        print(f"| ")
        print(f"| Average timings")
        print(f"|    ns/day = {avg_rate:9.2f}")
        print(f"| ")
        print(f"| Estimated time")
        print(f"|    remaining = {eta:.2f} hours")
        print(f"|    completed = {eta_time.ctime()}")
        print("-"*80)