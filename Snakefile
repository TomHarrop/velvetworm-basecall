#!/usr/bin/env python3

import os
import pathlib2
import tarfile
import warnings
import multiprocessing


#############
# FUNCTIONS #
#############

def aggregate_untar(wildcards):
    checkpoint_output = checkpoints.untar.get(**wildcards).output['tardir']
    return expand(
        'output/01_untar/{fc}/{tar_file}/{f5}.fast5',
        fc=wildcards.fc,
        tar_file=wildcards.tar_file,
        f5=glob_wildcards(os.path.join(checkpoint_output, '{f5}.fast5')).f5)


def find_all_by_ext(read_dir, extension):
    my_read_path = pathlib2.Path(read_dir)
    my_files = []
    for dirpath, dirnames, filenames in os.walk(my_read_path):
        for fn in filenames:
            if fn.endswith(extension):
                my_files.append(
                    str(pathlib2.Path(dirpath, fn).resolve()))
    return my_files


def get_basecall_input(wildcards):
    if wildcards.fc in fcs_with_tar:
        all_tarnames = [pathlib2.Path(x).stem for x in fcs_with_tar[wildcards.fc]]
        return expand('output/01_untar/flags/{{fc}}_{tar_file}.untarred',
                      tar_file=all_tarnames)
    if wildcards.fc in fcs_with_fast5:
        return str(pathlib2.Path(read_dir, '{fc}'))


def get_basecall_wd(fc):
    if fc in fcs_with_tar:
        return 'output/01_untar/{}'.format(fc)
    if fc in fcs_with_fast5:
        return str(pathlib2.Path(read_dir, fc))


def get_untar_input(wildcards):
    my_tarfile = [x for x in fcs_with_tar[wildcards.fc]
                  if pathlib2.Path(x).stem == wildcards.tar_file]
    return(my_tarfile)


def list_fast5_in_tar(tar_file):
    with tarfile.open(tar_file) as my_tar:
        warnings.warn(
            'Scanning {} for fast5, will take a while'.format(tar_file))
        my_members = my_tar.getmembers()
    my_tarpaths = [pathlib2.Path(x.path).name for x in my_members 
                   if x.path.endswith('.fast5')]
    return(my_tarpaths)


def match_filename_to_flowcell(filename, fc_list):
    my_parts = pathlib2.Path(filename).parts
    my_fc_match = [x for x in my_parts if x in fc_list][0]
    if len(my_fc_match) == 0:
        raise ValueError('No match for {}'.format(filename))
    return(my_fc_match)


###########
# GLOBALS #
###########

read_dir = 'data/raw'
fc_list = [
    '20190211-NPL0612-P2-A5-D5',
    '20190215-NPL0612-P2-E1-H1',
    '20190216-NPL0612-P1-A9-D9',
    '20190216-NPL0612-P1-E9-H9',
    '20190217-NPL0612-P1-A11-D11']

pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
guppy_container = 'shub://TomHarrop/singularity-containers:guppy_2.3.7'
minionqc_container = 'shub://TomHarrop/singularity-containers:minionqc_1.4.1'

########
# MAIN #
########

all_tarfiles = find_all_by_ext(read_dir, '.tar')
all_fast5 = find_all_by_ext(read_dir, '.fast5')

# find the flowcells that have fast5 files
flowcell_to_fast5 = {key: [] for key in fc_list}
for fast5 in all_fast5:
    my_fc = match_filename_to_flowcell(fast5, fc_list)
    flowcell_to_fast5[my_fc].append(fast5)

fcs_with_fast5 = {k: v for k, v in flowcell_to_fast5.items() if len(v) > 0}

# find the flowcells that have tar files
flowcell_to_tar = {key: [] for key in fc_list}
for tar in all_tarfiles:
    my_fc = match_filename_to_flowcell(tar, fc_list)
    flowcell_to_tar[my_fc].append(tar)

fcs_with_tar = {k: v for k, v in flowcell_to_tar.items() if len(v) > 0}

#########
# RULES #
#########

rule target:
    input:
        expand('output/02_basecalled/{fc}/sequencing_summary.txt',
               fc=list(fcs_with_tar.keys()) + list(fcs_with_fast5.keys())),
        'output/03_minionqc/combinedQC/summary.yaml'

# run minion qc on output
rule minionqc:
    input:
        expand('output/02_basecalled/{fc}/sequencing_summary.txt',
               fc=fc_list)
    output:
        'output/03_minionqc/combinedQC/summary.yaml'
    params:
        search_dir = 'output/02_basecalled',
        outdir = 'output/03_minionqc'
    threads:
        min(len(fc_list), multiprocessing.cpu_count())
    priority:
        1
    singularity:
        minionqc_container
    log:
        'output/logs/minionqc.log'
    shell:
        'MinIONQC.R '
        '--processors={threads} '
        '--input={params.search_dir} '
        '--outputdirectory={params.outdir} '
        '&> {log}'


# test aggregate fast5 files
rule basecall:
    input:
        get_basecall_input
    output:
        'output/02_basecalled/{fc}/sequencing_summary.txt'
    params:
        wd = lambda wildcards: get_basecall_wd(wildcards.fc),
        outdir = 'output/02_basecalled/{fc}'
    log:
        'output/logs/{fc}_guppy.log'
    resources:
        gpu = 1
    priority:
        1
    singularity:
        guppy_container
    shell:
        'guppy_basecaller '
        '--flowcell FLO-PRO002 '
        '--kit SQK-LSK109 '
        '--input {params.wd} '
        '--save_path {params.outdir} '
        '--recursive '
        '--device auto '
        '--num_callers 16 '
        '--chunks_per_runner 96 '
        '&> {log}'


# test aggregate untar
rule aggregate:
    input:
        aggregate_untar
    output:
        'output/01_untar/flags/{fc}_{tar_file}.untarred'
    wildcard_constraints:
        tar_file = '\d+'
    shell:
        'touch {output}'

# generic untar rule
checkpoint untar:
    input:
        get_untar_input
    output:
        tardir = directory('output/01_untar/{fc}/{tar_file}')
    params:
        tardir = 'output/01_untar/{fc}'
    wildcard_constraints:
        tar_file = '\d+'
    log:
        'output/logs/{fc}_{tar_file}.untar.log'
    singularity:
        pigz_container
    shell:
        'tar -x -f "{input}" -C {params.tardir} &> {log}'
