configfile: 'config.yaml'

from pathlib import Path
import datetime as dt
import os

DATAPATH = 'Data/CPTAC'
INTPATH = 'Intermediate_Data/CPTAC'
DATE = dt.datetime.today().strftime("%Y_%m_%d")
DAY = dt.datetime.today().strftime("%d")

CHROMOSOMES = ['chr'+str(n) for n in range(1,23)]

PARTONE = os.listdir(DATAPATH)
PARTTWO = {(a, [b for b in os.listdir(DATAPATH+'/'+a) if b.endswith('bam')][0][:-4]) for a in PARTONE}
#Note: change back index to -18 to take off full _wgs_ tag
SAMPLEPATHS = ['/'.join(a) for a in PARTTWO]

#Test rule (shows that snakemake is working)
rule checkVars:
  params:
    bams=expand('{intpath}/{samplepaths}.{chromosomes}.bam', intpath = INTPATH, samplepaths=SAMPLEPATHS[0], chromosomes=CHROMOSOMES)
  shell:
    '''
    echo {params.bams}
    '''

rule test:
  input:
    bams=expand('{intpath}/{samplepaths}.{chromosomes}.bam', intpath = INTPATH, samplepaths=SAMPLEPATHS, chromosomes=CHROMOSOMES)

rule bamsplit:
  input:
    bam=expand('{datapath}/{{code1}}/{{code2}}_wgs_gdc_realn.bam', datapath=DATAPATH),
    bai=expand('{datapath}/{{code1}}/{{code2}}_wgs_gdc_realn.bai', datapath=DATAPATH)
  output:
    bam=expand('{intpath}/{{code1}}/{{code2}}_wgs_gdc_realn.{{chrN}}.bam', intpath=INTPATH)
  params:
    chrP='{chrN}'
  shell:
    '''
    module load samtools
    samtools view {input.bam} {params.chrP} > {output.bam}
    '''

