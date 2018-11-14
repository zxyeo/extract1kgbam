#!/usr/bin/env python

import sys
import os
import errno
import argparse
import time
import subprocess
import tenacity
import multiprocessing
from functools import partial
import logging
import re

'''
input path of bams
input target files containing chromosome:start-stop
input working directory

./extract1kgBAM.py 
--bamlist bam_batch.list 
--target primary_targets.list 
--workdir ./outputs
'''


# parse arguments: input 1000genome bams, region of interest
def parse_input():
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(description='Download 1000genomes BAM by regions (via FTP or AWS S3)')
    parser.add_argument('--bam', help='path of bam', type=str)
    parser.add_argument('--bamlist', help='path of bam list', type=str)
    parser.add_argument('--target', help='path to target file', type=str, required=True)
    parser.add_argument('--workdir', help='working directory', type=str, required=True)
    parser.add_argument('--force', help='force overwrite even if file exists', type=str2bool, nargs='?',
                        const=True, default=False)
    parser.add_argument('--version', action='version', version='%(prog)s (version 0.1)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    return args


def check_mkdir(outdir):
    try:
        os.makedirs(outdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def extract1kgBAM(target, wkdir, force, bam):
    def printlog(log_info, logger):
        logger.info(log_info)

    @tenacity.retry
    def extractbam(bam, pre, outdir, region, force):
        def create_bamlist(outbam):
            dn = os.path.dirname(outbam)
            path_bam = os.path.join(dn, "bam.list")

            if os.path.exists(path_bam):
                lines = [i.rstrip() for i in open(path_bam).readlines()]
                if outbam not in lines:
                    with open(path_bam, "a") as lines:
                        lines.write("{}\n".format(outbam))
                else:
                    print("Already in bam.list: {}".format(outbam))
            else:
                with open(path_bam, "w") as lines:
                    lines.write("{}\n".format(outbam))

        # start time
        tx = time.time()
        # create output directory if required
        outdir_sub = os.path.join(outdir, "subbam")
        check_mkdir(outdir_sub)
        outbam = "{}/{}.{}.bam".format(outdir_sub, region, pre)
        skip = ""
        cmd = "samtools view {0} {1} -o {2} -O BAM".format(bam, region, outbam)
        if force:
            subprocess.check_call(cmd, shell=True, executable='/bin/bash')
            create_bamlist(outbam)
        else:
            # run if not exist
            if not os.path.isfile(outbam):
                subprocess.check_call(cmd, shell=True, executable='/bin/bash')
                create_bamlist(outbam)
            else:
                skip = " SKIP: same file found"

        # end time: per bam
        txrun = time.time()
        runtime = str(round(txrun - tx, 6)) + "(s)" + skip
        log = '{}\t{}'.format(region, runtime)
        printlog(log, perbam_logger)

    def extractbams(bam, pre, outdir, target, force):
        with open(target, 'r') as f:
            regions = [line.rstrip('\n') for line in f]

        log = "Extract:\t{}".format(pre)
        printlog(log, perbam_logger)
        [extractbam(bam, pre, outdir, region, force) for region in regions]

    def mergebam(outdir, pre, force):
        # start time
        tx = time.time()

        skip = ""
        mergedbam = os.path.join(outdir, "{}.targets.bam".format(pre))
        cmd1 = "samtools merge -c -p -b {}/subbam/bam.list {}".format(outdir, mergedbam)
        cmd2 = "samtools index {}".format(mergedbam)

        if force:
            subprocess.check_call(cmd1, shell=True, executable='/bin/bash')
            subprocess.check_call(cmd2, shell=True, executable='/bin/bash')
        else:
            if not os.path.isfile(mergedbam):
                subprocess.check_call(cmd1, shell=True, executable='/bin/bash')
                subprocess.check_call(cmd2, shell=True, executable='/bin/bash')
            else:
                skip = " SKIP: same file found"

        checkrun = "merged bam"
        txrun = time.time()
        runtime = str(round(txrun - tx, 6)) + "(s)" + skip
        log = '{}\t{}'.format(checkrun, runtime)
        printlog(log, perbam_logger)

    # create output directory if required
    check_mkdir(wkdir)

    # get path
    bn = os.path.basename(bam)
    # remove extension
    pre = bn.rstrip('.bam')
    outdir = os.path.join(wkdir, pre)

    check_mkdir(outdir)
    # logging
    perbam_logger = setup_logger('perbam_logger', '{}/{}.log'.format(outdir, pre))

    extractbams(bam, pre, outdir, target, force)
    mergebam(outdir, pre)


def extract1kgBAM_batch(target, wkdir, force, bamlist):
    def setup_multiprocess():
        try:
            cpus = multiprocessing.cpu_count()
        except NotImplementedError:
            cpus = 2  # arbitrary default

        pool = multiprocessing.Pool(processes=cpus)
        return pool

    # generate list of input bams
    with open(bamlist, 'r') as f:
        bam_list = [line.rstrip('\n') for line in f]

    pool = setup_multiprocess()
    func = partial(extract1kgBAM, target, wkdir, force)
    pool.map(func, bam_list)
    pool.close()
    pool.join()


def setup_logger(name, log_file, level=logging.INFO):
    formatter = logging.Formatter('%(asctime)s %(message)s')
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def main():
    def runprogram():
        if (bam is None) and (bamlist is not None):
            extract1kgBAM_batch(target, wkdir, force, bamlist)
        if (bam is not None) and (bamlist is None):
            extract1kgBAM(target, wkdir, force, bam)

    def printlog(log_info, logger):
        print(log_info)
        logger.info(log_info)

    args = parse_input()
    bamlist = args.bamlist
    bam = args.bam
    target = args.target
    workdir = args.workdir
    force = args.force
    if (bam is not None) or (bamlist is not None):
        try:
            exec_dir = os.getcwd()
            wkdir = os.path.abspath(workdir)
            check_mkdir(wkdir)

            main_logger = setup_logger('main_logger', '{}/extract1kgBAM.log'.format(wkdir))
            printlog('Preparing processes ......\n Current directory: {}'.format(exec_dir), main_logger)
            # cleaning
            purge(exec_dir, ".bam.bai")

            # start time
            tx = time.time()

            printlog('Started ......', main_logger)
            runprogram()

            # cleaning
            purge(exec_dir, ".bam.bai")

            # end time
            txrun = time.time()
            runtime = str(round(txrun - tx, 6)) + "(s)"
            printlog('Completed\n Output directory {}\nTotal run time\t{}'.format(wkdir, runtime), main_logger)
        except IOError:
            print("ERROR! Check inputs!")
    else:
        print("Please provide input bam!")


if __name__ == "__main__":
    main()


