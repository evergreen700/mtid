configfile: 'config.yaml'

from pathlib import Path
import datetime as dt
import os
import json

MANIFEST = 'Intermediate_Data/CPTAC/manifest.json'
DATAPATH = 'Data/CPTAC'
INTPATH = 'Intermediate_Data/CPTAC'
DATE = dt.datetime.today().strftime("%Y_%m_%d")
DAY = dt.datetime.today().strftime("%d")
CASES = {}
with open(MANIFEST, 'r') as manifest:
  for case in manifest.readlines():
    i = json.loads(case)
    CASES[i['case_id']] = i

CHROMOSOMES = ['chr'+str(n) for n in range(1,23)]

PINSUFFIX = ['BP', 'CloseEndMapped', 'D', 'INT_final', 'INV', 'LI', 'RP', 'SI', 'TD']
GDCH38 = 'Data/GDC_h38/GRCh38.d1.vd1.fa'
GDCH38INFO = 'GRCh38.d1.vd1'

PARTONE = os.listdir(DATAPATH)
FCASES = [CASES[a] for a in CASES.keys() if CASES[a]['Primary Tumor'][0].split('/')[0] in PARTONE]
CASESIMP = {a['case_id']:{'Tumor':a['Primary Tumor'][0].split('/')[0], 'Normal':a['Blood Derived Normal'][0].split('/')[0] if a['Blood Derived Normal'] is not None else a['Solid Tissue Normal'].split('/')[0]} for a in FCASES}
CASENAMES = [a['case_id'] for a in FCASES]
TUMORS = [a['Primary Tumor'][0] for a in FCASES]
NORMALS = [a['Blood Derived Normal'][0] if a['Blood Derived Normal'] is not None else a['Solid Tissue Normal'][0] for a in FCASES]
PARTTWO = {(a, [b for b in os.listdir(DATAPATH+'/'+a) if b.endswith('bam')][0][:-4]) for a in PARTONE}
#Note: change back index to -18 to take off full _wgs_ tag
SAMPLEPATHS = ['/'.join(a) for a in PARTTWO]


#Test rule (shows that snakemake is working)
rule checkVars:
  input:
    bams=[a['case_id'] for a in FCASES]
  shell:
    '''
    echo {params.bams}
    '''

rule week1Benchmark:
  input:
    pindel=expand([INTPATH+'/'+a['case_id']+'/'+CASESIMP[a['case_id']]['Tumor']+'_AGAINST_'+CASESIMP[a['case_id']]['Normal']+'_somatic_pindel.{chrN}.vcf.gz' for a in FCASES], chrN='chr22'),
    abra2=expand([INTPATH+'/'+a['case_id']+'/'+CASESIMP[a['case_id']]['Tumor']+'_AGAINST_'+CASESIMP[a['case_id']]['Normal']+'_somatic_abra.{chrN}.vcf.gz' for a in FCASES], chrN='chr22'),
    rufus=expand([INTPATH+'/'+a['case_id']+'/'+CASESIMP[a['case_id']]['Tumor']+'_AGAINST_'+CASESIMP[a['case_id']]['Normal']+'_somatic_rufus.{chrN}.vcf.gz' for a in FCASES], chrN='chr22'),
    platypus=expand([INTPATH+'/'+a['case_id']+'/'+CASESIMP[a['case_id']]['Tumor']+'_AGAINST_'+CASESIMP[a['case_id']]['Normal']+'_somatic_platypus.{chrN}.vcf.gz' for a in FCASES], chrN='chr22')

rule week2Benchmark:
  input:
    graph=expand([INTPATH+'/'+a['case_id']+'/'+CASESIMP[a['case_id']]['Tumor']+'_AGAINST_'+CASESIMP[a['case_id']]['Normal']+'_somatic_intersections.{chrN}.tsv' for a in FCASES], chrN=CHROMOSOMES),

rule upSetPlot:
  input:
    tsv=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_intersections.{chrN}.tsv'
  output:
    png=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_UPSETPLOT.{chrN}.png'
  shell:
    '''
    module load python/3.11
    python upset_plot.py {input.tsv} {output.png}
    '''

rule intersect:
   input:
     vcfs=expand(INTPATH+'/{{case_id}}/{{id1}}_AGAINST_{{id2}}_somatic_{tool}.{{chrN}}.vcf.gz', tool=['pindel','abra','rufus','platypus'])
   output:
     tsv=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_intersections.{chrN}.tsv'
   params:
     inpath=INTPATH+'/{case_id}/',
     chrP='{chrN}'
   shell:
     '''
     ./compare.sh {params.inpath} {params.chrP} > {output.tsv}
     '''

rule somaticFromGermlinePindel:
  input:
    Tvcf=INTPATH+'/{case_id}/{id1}.{chrN}.pindel.vcf',
    Nvcf=INTPATH+'/{case_id}/{id2}.{chrN}.pindel.vcf'
  output:
    Svcf=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_pindel.{chrN}.vcf'
  shell:
    '''
    module load bedtools
    bedtools intersect -a {input.Tvcf} -b {input.Nvcf} -v -header > {output.Svcf}
    '''

rule somaticFromGermlinePlatypus:
  input:
    Tvcf=INTPATH+'/{case_id}/{id1}.{chrN}.platypus.vcf',
    Nvcf=INTPATH+'/{case_id}/{id2}.{chrN}.platypus.vcf'
  output:
    Svcf=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_platypus.{chrN}.vcf'
  shell:
    '''
    module load bedtools
    bedtools intersect -a {input.Tvcf} -b {input.Nvcf} -v -header > {output.Svcf}
    '''

rule test:
  input:
    bams=expand('{intpath}/{samplepaths}_abra.{chromosomes}.bai', intpath = INTPATH, samplepaths=SAMPLEPATHS, chromosomes=CHROMOSOMES)

rule test2:
  input:
    pindels=expand('{intpath}/{samplepaths}.{chrN}.pindel.vcf', intpath=INTPATH, samplepaths=SAMPLEPATHS, chrN="chr22", suf=PINSUFFIX)

rule test3:
  input:
    pindels=expand('{intpath}/{samplepaths}.{chrN}_{suf}', intpath=INTPATH, samplepaths=SAMPLEPATHS, chrN=CHROMOSOMES, suf='COMBINED')

#pindel rules here

rule pindelVCF:
  input:
    pindel='{path}_COMBINED',
    pin2vcf='Tools/Pindel/pindel2vcf'
  output:
    vcf='{path}.pindel.vcf'
  shell:
    '''
    ./{input.pin2vcf} -p {input.pindel} -r {GDCH38} -R {GDCH38INFO} -d 2023 -v {output.vcf}
    '''

rule pindelJoin:
  input:
    pindels=expand('{{path}}_{suf}', suf=PINSUFFIX)
  output:
    pindel='{path}_COMBINED'
  shell:
    '''
    cat {input.pindels} > {output.pindel}
    '''

rule pindel:
  input:
    bam=expand('{intpath}/{{case}}/{{id}}.{{chrN}}.bam', intpath=INTPATH),
    bai=expand('{intpath}/{{case}}/{{id}}.{{chrN}}.bai', intpath=INTPATH),
    pinPath='Tools/Pindel/pindel',
    ref=expand('{refpath}', refpath=GDCH38)
  output:
    pindels=expand('{intpath}/{{case}}/{{id}}.{{chrN}}_{suf}', intpath=INTPATH, suf=PINSUFFIX)
  params:
    bconfig=expand('{intpath}/{{case}}/{{id}}.{{chrN}}.config.txt', intpath=INTPATH),
    prefix=expand('{intpath}/{{case}}/{{id}}.{{chrN}}', intpath=INTPATH),
    chrP='{chrN}'
  shell:
    '''
    echo "{input.bam}  250  pindel" > {params.bconfig}
    ./{input.pinPath} -f {input.ref} -i {params.bconfig} -c {params.chrP} -A 30 -M 8 -o {params.prefix} 
    '''

#Insert size of 250 is entirely arbitrary and I don't know how to get a better one. At least this one works without raising errors.

#ABRA2 rules here

rule CadabraVCF:
  input:
    Tbam=INTPATH+'/{case}/{Tid}.WITH_NORMAL.{Nid}_abra.{chrN}.bam',
    Tbai=INTPATH+'/{case}/{Tid}.WITH_NORMAL.{Nid}_abra.{chrN}.bai',
    Nbam=INTPATH+'/{case}/{Nid}.WITH_TUMOR.{Tid}_abra.{chrN}.bam',
    Nbai=INTPATH+'/{case}/{Nid}.WITH_TUMOR.{Tid}_abra.{chrN}.bai',
    ref=GDCH38,
    abra='Tools/Abra/abra2-2.23.jar'
  output:
    vcf=INTPATH+'/{case}/{Tid}_AGAINST_{Nid}_somatic_abra.{chrN}.vcf'
  shell:
    '''
    java -Xmx4G -cp {input.abra} abra.cadabra.Cadabra --ref {input.ref} --normal {input.Nbam} --tumor {input.Tbam} > {output.vcf}
    '''

rule ABRAlignSomatic:
  input:
    Tbam=expand('{intpath}/{{case}}/{{Tid}}.{{chrN}}.bam', intpath=INTPATH),
    Tbai=expand('{intpath}/{{case}}/{{Tid}}.{{chrN}}.bai', intpath=INTPATH),
    Nbam=expand('{intpath}/{{case}}/{{Nid}}.{{chrN}}.bam', intpath=INTPATH),
    Nbai=expand('{intpath}/{{case}}/{{Nid}}.{{chrN}}.bai', intpath=INTPATH),
    abra='Tools/Abra/abra2-2.23.jar'
  output:
    Tbam=INTPATH+'/{case}/{Tid}.WITH_NORMAL.{Nid}_abra.{chrN}.bam',
    Nbam=INTPATH+'/{case}/{Nid}.WITH_TUMOR.{Tid}_abra.{chrN}.bam'
  shell:
    '''
    java -Xmx32G -jar {input.abra} --in {input.Nbam},{input.Tbam} --out {output.Nbam},{output.Tbam} --ref {GDCH38} --threads 8
    '''

rule RUFUS:
  input:
    Tbam=INTPATH+'/{case}/{Tid}.{chrN}.bam',
    Tbai=INTPATH+'/{case}/{Tid}.{chrN}.bai',
    Nbam=INTPATH+'/{case}/{Nid}.{chrN}.bam',
    Nbai=INTPATH+'/{case}/{Nid}.{chrN}.bai',
    ref=GDCH38,
    rufus='Tools/Rufus/RUFUS/runRufus.sh'
  output:
    vcf=INTPATH+'/{case}/{Tid}_AGAINST_{Nid}_somatic_rufus.{chrN}.vcf.gz'
  params:
    ksize=25,
    rawOutput='{Tid}.{chrN}.bam.generator.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz'
  shell:
    '''
    module load samtools
    module load bedtools
    ./{input.rufus} --subject {input.Tbam} --controls {input.Nbam} --kmersize {params.ksize} --threads 41 --ref {input.ref}
    mv {params.rawOutput} {output.vcf}
    '''

rule platypus:
  input:
    bam=INTPATH+'/{case}/{id}.{chrN}.bam',
    bai=INTPATH+'/{case}/{id}.{chrN}.bai',
    ref=GDCH38,
    platypus='Tools/Platypus_0.8.1/Platypus.py'
  output:
    vcf=INTPATH+'/{case}/{id}.{chrN}.platypus.vcf'
  shell:
    '''
    module load python/2.7
    python {input.platypus} callVariants --refFile {input.ref} --bamFiles {input.bam} --minReads 8 --output {output.vcf}
    '''  
#minReads 8 is too friendly

rule index:
  input:
    bam='{path}.bam'
  output:
    bai='{path}.bai'
  shell:
    '''
    module load samtools
    samtools index {input.bam} -b -o {output.bai}
    '''

rule bamsplit:
  input:
    bamdir=expand('{datapath}/{{id}}/', datapath=DATAPATH),
  output:
    bam=expand('{intpath}/{{case}}/{{id}}.{{chrN}}.bam', intpath=INTPATH)
  params:
    chrP='{chrN}'
  shell:
    '''
    module load samtools
    samtools view -b {input.bamdir}/*.bam {params.chrP} > {output.bam}
    '''

rule vcfGzip:
  input: '{path}_somatic_{tool}.{chrN}.vcf'
  output:'{path}_somatic_{tool}.{chrN}.vcf.gz'
  wildcard_constraints:
    tool="pindel|abra|platypus"
  shell:
    '''
    gzip -c {input} > {output}
    '''
