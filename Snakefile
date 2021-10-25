#!/usr/bin/env python3

import fileinput
import os
import pathlib
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
    my_read_path = pathlib.Path(read_dir)
    my_files = []
    for dirpath, dirnames, filenames in os.walk(my_read_path):
        for fn in filenames:
            if fn.endswith(extension):
                my_files.append(
                    str(pathlib.Path(dirpath, fn).resolve()))
    return my_files


def get_basecall_input(wildcards):
    if wildcards.fc in fcs_with_tar:
        all_tarnames = [pathlib.Path(x).stem for x in fcs_with_tar[wildcards.fc]]
        return expand('output/01_untar/flags/{{fc}}_{tar_file}.untarred',
                      tar_file=all_tarnames)
    if wildcards.fc in fcs_with_fast5:
        return str(pathlib.Path(read_dir, '{fc}'))


def get_basecall_wd(fc):
    if fc in fcs_with_tar:
        return 'output/01_untar/{}'.format(fc)
    if fc in fcs_with_fast5:
        return str(pathlib.Path(read_dir, fc))

def get_fastq_files(wildcards):
    glob_path = 'output/02_basecalled/{}/'.format(wildcards.fc)
    fq_list = glob_wildcards(os.path.join(glob_path, '{fq}.fastq')).fq
    return(expand(os.path.join(glob_path, '{fq}.fastq'),
                  fq=fq_list))

def get_untar_input(wildcards):
    my_tarfile = [x for x in fcs_with_tar[wildcards.fc]
                  if pathlib.Path(x).stem == wildcards.tar_file]
    return(my_tarfile)


def list_fast5_in_tar(tar_file):
    with tarfile.open(tar_file) as my_tar:
        warnings.warn(
            'Scanning {} for fast5, will take a while'.format(tar_file))
        my_members = my_tar.getmembers()
    my_tarpaths = [pathlib.Path(x.path).name for x in my_members 
                   if x.path.endswith('.fast5')]
    return(my_tarpaths)


def match_filename_to_flowcell(filename, fc_list):
    my_parts = pathlib.Path(filename).parts
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
guppy_container = 'docker://ghcr.io/tomharrop/container-guppy:5.0.16'
minionqc_container = 'shub://TomHarrop/ont-containers:minionqc_1.4.1'
bbduk_container = 'docker://ghcr.io/deardenlab/container-bbmap:bbmap_38.90'

# get the host from the cli
# snakemake config host=host.address:port
# host = config['host']

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
               fc=fc_list),
        'output/03_minionqc/combinedQC/summary.yaml',
        'output/04_filtered/all_passed_reads.fastq'


# combine and filter sequencing_summary files
rule join_passed_reads:
    input:
        expand('output/04_filtered/{fc}/kept_reads.fastq',
               fc=fc_list)
    output:
        'output/04_filtered/all_passed_reads.fastq'
    shell:
        'cat {input} > {output}'

rule filter_by_name:
    input:
        fq = 'output/04_filtered/{fc}/all_reads.fastq',
        names = 'output/04_filtered/{fc}/reads_to_keep.txt'
    output:
        temp('output/04_filtered/{fc}/kept_reads.fastq')
    log:
        'output/logs/filter_by_name_{fc}.log'
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input.fq} '
        'names={input.names} '
        'include=t '
        'out={output} '
        '2> {log} '

rule merge_by_flowcell:
    input:
        'output/02_basecalled/{fc}/sequencing_summary.txt'
    output:
        fq = temp('output/04_filtered/{fc}/all_reads.fastq')
    params:
        fq_list = lambda wildcards: get_fastq_files(wildcards)
    run:
        with open(output.fq, 'wt') as f:
            for line in fileinput.input(params.fq_list):
                f.write(line)

rule choose_reads_to_keep:
    input:
        'output/02_basecalled/{fc}/sequencing_summary.txt'
    output:
        'output/04_filtered/{fc}/reads_to_keep.txt'
    log:
        'output/logs/choose_reads_to_keep_{fc}.log'
    singularity:
        minionqc_container
    script:
        'src/choose_reads_to_keep.R'

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
        outdir = 'output/02_basecalled/{fc}',
        ont_config = '/opt/ont/guppy/data/dna_r9.4.1_450bps_sup_prom.cfg'
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
        '--device auto '
        '--input_path {params.wd} '
        '--save_path {params.outdir} '
        '--recursive '
        '--config {params.ont_config} '
        '&> {log}'


# test aggregate untar
rule aggregate:
    input:
        'output/01_untar/{fc}/{tar_file}'
    output:
        'output/01_untar/flags/{fc}_{tar_file}.untarred'
    wildcard_constraints:
        tar_file = '\d+'
    shell:
        'touch {output}'

# generic untar rule
rule untar:
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
