configfile: 'config.yaml'

from pathlib import Path
import datetime as dt
import os
import json

MANIFEST = 'Intermediate_Data/CPTAC/newManifest.json'
DATAPATH = 'Data/CPTAC_DOWNLOADS'
INTPATH = 'Intermediate_Data/CPTAC'
DATE = dt.datetime.today().strftime("%Y_%m_%d")
DAY = dt.datetime.today().strftime("%d")
CASES = []
with open(MANIFEST, 'r') as manifest:
  for case in manifest.readlines():
    i = json.loads(case)
    CASES.append(i)

CHROMOSOMES = ['chr'+str(n) for n in range(1,23)]

PINSUFFIX = ['BP', 'CloseEndMapped', 'D', 'INT_final', 'INV', 'LI', 'RP', 'SI', 'TD']
GDCH38 = 'Data/GDC_h38/GRCh38.d1.vd1.fa'
GDCH38INFO = 'GRCh38.d1.vd1'

#INCLUDEDCASES = os.listdir(DATAPATH)
FCASES = [CASES[0]]#[row for row in CASES if row["case_id"] in INCLUDEDCASES]

CASEBAMS = {a['case_id']:{'Tumor':a['Aligned Reads- Primary Tumor'][0]["ID"], 'Normal':a['Aligned Reads- Blood Derived Normal'][0]["ID"] if a['Aligned Reads- Blood Derived Normal'] is not None else a['Aligned Reads- Solid Tissue Normal'][0]["ID"]} for a in FCASES}
#CASEVCFS = {a['case_id']: [].extend(a["Raw Simple Somatic Mutation"]).extend(a["Structural Rearrangement"]) for a in FCASES}

##For downloading without snakemake
#with open("symDirs.sh", "w") as outSymManifest:
#  with open("vcfLog.txt", "w") as outManifest:
#    for c in CASES:
#      outManifest.write("mkdir -p "+DATAPATH+"/"+c["case_id"]+"\n")
#      outSymManifest.write("ln -sfnr "+DATAPATH+"/"+c["case_id"]+ " Data/CASE_SYMLINKS/" + c["case_submitter_id"] + "\n")
#      if c["Raw Simple Somatic Mutation"]!=None:
#        for v in c["Raw Simple Somatic Mutation"]:
#          outManifest.write("./Tools/gdc-client download "+v["ID"]+" -t Data/New/gdc-user-token.2023-09-26T20_07_40.822Z.txt -d "+DATAPATH+"/"+c["case_id"]+"\n")
#      if c["Structural Rearrangement"]!=None:
#        for v in c["Structural Rearrangement"]:
#          outManifest.write("./Tools/gdc-client download "+v["ID"]+" -t Data/New/gdc-user-token.2023-09-26T20_07_40.822Z.txt -d "+DATAPATH+"/"+c["case_id"]+"\n")
#  


CASENAMES = [a['case_id'] for a in FCASES]
#TUMORS = [a['Primary Tumor'][0] for a in FCASES]
#NORMALS = [a['Blood Derived Normal'][0] if a['Blood Derived Normal'] is not None else a['Solid Tissue Normal'][0] for a in FCASES]
#PARTTWO = {(a, [b for b in os.listdir(DATAPATH+'/'+a) if b.endswith('bam')][0][:-4]) for a in PARTONE}
#Note: change back index to -18 to take off full _wgs_ tag
#SAMPLEPATHS = ['/'.join(a) for a in PARTTWO]
#CASEVCFPATHS = []
#for a in CASENAMES:
#  for b in CASEVCFS[a]:
#    CASEVCFPATHS.append(

localrules: downloadNewFiles

wildcard_constraints:
  file_id="[\w-]+",
  Tfile_id="[\w-]+",
  Nfile_id="[\w-]+",
  chrN="chr[\d]+",
  gTool="pindel|platypus"


#Test rule (shows that snakemake is working)
rule checkVars:
  input:
    bams=[a['case_id'] for a in FCASES]
  shell:
    '''
    echo {params.bams}
    '''

rule downloadVCFs:
  input:
    somaticMutation=[[DATAPATH+'/'+a['case_id'] +'/'+b["ID"] for b in a["Raw Simple Somatic Mutation"]] for a in FCASES if a["Raw Simple Somatic Mutation"]!=None],
    structuralVariants=[[DATAPATH+'/'+a['case_id'] +'/'+b["ID"] for b in a["Structural Rearrangement"]] for a in FCASES if a["Structural Rearrangement"]!=None]



rule week1Benchmark:
  input:
    pindel=expand([INTPATH+'/'+a['case_id']+'/'+CASEBAMS[a['case_id']]['Tumor']+'_AGAINST_'+CASEBAMS[a['case_id']]['Normal']+'_somatic_pindel.{chrN}.vcf.gz' for a in FCASES], chrN='chr22'),
    abra2=expand([INTPATH+'/'+a['case_id']+'/'+CASEBAMS[a['case_id']]['Tumor']+'_AGAINST_'+CASEBAMS[a['case_id']]['Normal']+'_somatic_abra.{chrN}.vcf.gz' for a in FCASES], chrN='chr22'),
    rufus=expand([INTPATH+'/'+a['case_id']+'/'+CASEBAMS[a['case_id']]['Tumor']+'_AGAINST_'+CASEBAMS[a['case_id']]['Normal']+'_somatic_rufus.{chrN}.vcf.gz' for a in FCASES], chrN='chr22'),
    platypus=expand([INTPATH+'/'+a['case_id']+'/'+CASEBAMS[a['case_id']]['Tumor']+'_AGAINST_'+CASEBAMS[a['case_id']]['Normal']+'_somatic_platypus.{chrN}.vcf.gz' for a in FCASES], chrN='chr22')

rule week2Benchmark:
  input:
    graph=expand([INTPATH+'/'+a['case_id']+'/'+CASEBAMS[a['case_id']]['Tumor']+'_AGAINST_'+CASEBAMS[a['case_id']]['Normal']+'_somatic_intersections.{chrN}.tsv' for a in FCASES], chrN="chr1"),

rule week3Benchmark:
  input:
    table=[INTPATH+'/'+a['case_id']+'/'+CASEBAMS[a['case_id']]['Tumor']+'_AGAINST_'+CASEBAMS[a['case_id']]['Normal']+'_somatic_intersections.tsv' for a in FCASES]

rule downloadNewFiles:
  priority: 0
  input:
    gdcClient="Tools/gdc-client",
    token="Data/New/gdc-user-token.2023-09-26T20_07_40.822Z.txt"
  output:
    cdir=directory(DATAPATH+"/{case_id}/{file_id}"),
  resources:
    gdc_download_limit=1
  params:
    casePath=DATAPATH+"/{case_id}",
    fileID="{file_id}"
  shell:
    '''
    mkdir -p {params.casePath}
    ./{input.gdcClient} download {params.fileID} -d {params.casePath} -t {input.token}
    '''

rule upSetPlot:
  input:
    tsv=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_intersections.tsv'
  output:
    png=INTPATH+'/{case_id}/{id1}_AGAINST_{id2}_somatic_UPSETPLOT.png'
  shell:
    '''
    module load python/3.11
    python upset_plot.py {input.tsv} {output.png}
    '''

rule intersect:
   priority: 2
   input:
     vcfs=expand(INTPATH+'/{{case_id}}/{{Tfile_id}}_AGAINST_{{Nfile_id}}_somatic_{tool}.{{chrN}}.vcf.gz', tool=['pindel','abra','rufus','platypus'])
   output:
     tsv=INTPATH+'/{case_id}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_intersections.{chrN}.tsv'
   params:
     inpath=INTPATH+'/{case_id}/',
     chrP='{chrN}'
   shell:
     '''
     ./compare.sh {params.inpath} {params.chrP} > {output.tsv}
     '''

rule intersectAll:
   priority: 2
   input:
     vcfs=expand(INTPATH+'/{{case_id}}/{{Tfile_id}}_AGAINST_{{Nfile_id}}_somatic_{tool}.vcf.gz', tool=['pindel','abra','rufus','platypus'])
   output:
     tsv=INTPATH+'/{case_id}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_intersections.tsv'
   params:
     inpath=INTPATH+'/{case_id}/',
     Tdatapath=DATAPATH+'/{case_id}/{Tfile_id}',
     Ndatapath=DATAPATH+'/{case_id}/{Nfile_id}'
   shell:
     '''
     ./compareAll.sh {params.inpath} > {output.tsv}
     rm {params.inpath}*.bam
     rm {params.inpath}*.bai
     rm -rf {params.Tdatapath}
     rm -rf {params.Ndatapath}
     '''

rule wholeGenomeVCF:
  priority: 2
  input:
    vcfs=expand(INTPATH+'/{{case_id}}/{{Tfile_id}}_AGAINST_{{Nfile_id}}_somatic_{{tool}}.{chrN}.vcf.gz', chrN=CHROMOSOMES)
  output:
    vcf=INTPATH+'/{case_id}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_{tool}.vcf.gz'
  params:
    first=INTPATH+'/{case_id}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_{tool}.chr1.vcf.gz'
  shell:
    '''
    cat {params.first} | grep "^#" > {output.vcf}
    cat {input.vcfs} | grep -v "^#" >> {output.vcf}
    '''

rule somaticFromGermline:
  priority: 2
  input:
    Tvcf=INTPATH+'/{case_id}/{Tfile_id}.{chrN}.{gTool}.vcf',
    Nvcf=INTPATH+'/{case_id}/{Nfile_id}.{chrN}.{gTool}.vcf'
  output:
    Svcf=INTPATH+'/{case_id}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_{gTool}.{chrN}.vcf'
  shell:
    '''
    module load bedtools
    bedtools intersect -a {input.Tvcf} -b {input.Nvcf} -v -header > {output.Svcf}
    rm {input.Tvcf} {input.Nvcf}
    '''

#rule somaticFromGermlinePlatypus:
#  priority: 2
#  input:
#    Tvcf=INTPATH+'/{case_id}/{Tfile_id}.{chrN}.platypus.vcf',
#    Nvcf=INTPATH+'/{case_id}/{Nfile_id}.{chrN}.platypus.vcf'
#  output:
#    Svcf=INTPATH+'/{case_id}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_platypus.{chrN}.vcf'
#  shell:
#    '''
#    module load bedtools
#    bedtools intersect -a {input.Tvcf} -b {input.Nvcf} -v -header > {output.Svcf}
#    '''

#rule test:
#  input:
#    bams=expand('{intpath}/{samplepaths}_abra.{chromosomes}.bai', intpath = INTPATH, samplepaths=SAMPLEPATHS, chromosomes=CHROMOSOMES)
#
#rule test2:
#  input:
#    pindels=expand('{intpath}/{samplepaths}.{chrN}.pindel.vcf', intpath=INTPATH, samplepaths=SAMPLEPATHS, chrN="chr22", suf=PINSUFFIX)
#
#rule test3:
#  input:
#    pindels=expand('{intpath}/{samplepaths}.{chrN}_{suf}', intpath=INTPATH, samplepaths=SAMPLEPATHS, chrN=CHROMOSOMES, suf='COMBINED')
#
#pindel rules here

rule pindelVCF:
  priority: 2
  input:
    pindel='{path}_COMBINED',
    pin2vcf='Tools/Pindel/pindel2vcf'
  output:
    vcf='{path}.pindel.vcf'
  shell:
    '''
    ./{input.pin2vcf} -p {input.pindel} -r {GDCH38} -R {GDCH38INFO} -d 2023 -v {output.vcf}
    rm {input.pindel}
    '''

rule pindelJoin:
  priority: 2
  input:
    pindels=expand('{{path}}_{suf}', suf=PINSUFFIX)
  output:
    pindel='{path}_COMBINED'
  shell:
    '''
    cat {input.pindels} > {output.pindel}
    rm {input.pindels}
    '''

rule pindel:
  priority: 2
  input:
    bam=expand('{intpath}/{{case}}/{{file_id}}.{{chrN}}.bam', intpath=INTPATH),
    bai=expand('{intpath}/{{case}}/{{file_id}}.{{chrN}}.bai', intpath=INTPATH),
    pinPath='Tools/Pindel/pindel',
    ref=expand('{refpath}', refpath=GDCH38)
  output:
    pindels=expand('{intpath}/{{case}}/{{file_id}}.{{chrN}}_{suf}', intpath=INTPATH, suf=PINSUFFIX)
  params:
    bconfig=expand('{intpath}/{{case}}/{{file_id}}.{{chrN}}.config.txt', intpath=INTPATH),
    prefix=expand('{intpath}/{{case}}/{{file_id}}.{{chrN}}', intpath=INTPATH),
    chrP='{chrN}'
  shell:
    '''
    echo "{input.bam}  250  pindel" > {params.bconfig}
    ./{input.pinPath} -f {input.ref} -i {params.bconfig} -c {params.chrP} -A 30 -M 8 -o {params.prefix} 
    '''

#Insert size of 250 is entirely arbitrary and I don't know how to get a better one. At least this one works without raising errors.

#ABRA2 rules here

rule CadabraVCF:
  priority: 3
  input:
    Tbam=INTPATH+'/{case}/{Tfile_id}.WITH_NORMAL.{Nfile_id}_abra.{chrN}.bam',
    Tbai=INTPATH+'/{case}/{Tfile_id}.WITH_NORMAL.{Nfile_id}_abra.{chrN}.bai',
    Nbam=INTPATH+'/{case}/{Nfile_id}.WITH_TUMOR.{Tfile_id}_abra.{chrN}.bam',
    Nbai=INTPATH+'/{case}/{Nfile_id}.WITH_TUMOR.{Tfile_id}_abra.{chrN}.bai',
    ref=GDCH38,
    abra='Tools/Abra/abra2-2.23.jar'
  output:
    vcf=INTPATH+'/{case}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_abra.{chrN}.vcf'
  shell:
    '''
    java -Xmx4G -cp {input.abra} abra.cadabra.Cadabra --ref {input.ref} --normal {input.Nbam} --tumor {input.Tbam} > {output.vcf}
    '''

rule ABRAlignSomatic:
  priority: 3
  input:
    Tbam=expand('{intpath}/{{case}}/{{Tfile_id}}.{{chrN}}.bam', intpath=INTPATH),
    Tbai=expand('{intpath}/{{case}}/{{Tfile_id}}.{{chrN}}.bai', intpath=INTPATH),
    Nbam=expand('{intpath}/{{case}}/{{Nfile_id}}.{{chrN}}.bam', intpath=INTPATH),
    Nbai=expand('{intpath}/{{case}}/{{Nfile_id}}.{{chrN}}.bai', intpath=INTPATH),
    abra='Tools/Abra/abra2-2.23.jar'
  output:
    Tbam=INTPATH+'/{case}/{Tfile_id}.WITH_NORMAL.{Nfile_id}_abra.{chrN}.bam',
    Nbam=INTPATH+'/{case}/{Nfile_id}.WITH_TUMOR.{Tfile_id}_abra.{chrN}.bam'
  shell:
    '''
    java -Xmx16G -jar {input.abra} --in {input.Nbam},{input.Tbam} --out {output.Nbam},{output.Tbam} --ref {GDCH38} --threads 8
    '''

rule RUFUS:
  priority: 3
  shadow: "shallow"
  input:
    Tbam=INTPATH+'/{case}/{Tfile_id}.{chrN}.bam',
    Tbai=INTPATH+'/{case}/{Tfile_id}.{chrN}.bai',
    Nbam=INTPATH+'/{case}/{Nfile_id}.{chrN}.bam',
    Nbai=INTPATH+'/{case}/{Nfile_id}.{chrN}.bai',
    ref=GDCH38,
    rufus='Tools/Rufus/RUFUS/runRufus.sh'
  output:
    vcf=INTPATH+'/{case}/{Tfile_id}_AGAINST_{Nfile_id}_somatic_rufus.{chrN}.vcf.gz'
  params:
    ksize=25,
    rawOutput='{Tfile_id}.{chrN}.bam.generator.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz',
  shell:
    '''
    module load samtools
    module load bedtools
    ./{input.rufus} --subject {input.Tbam} --controls {input.Nbam} --kmersize {params.ksize} --min 8 --threads 41 --ref {input.ref}
    mv {params.rawOutput} {output.vcf}
    '''

rule platypus:
  priority: 2
  input:
    bam=INTPATH+'/{case}/{file_id}.{chrN}.bam',
    bai=INTPATH+'/{case}/{file_id}.{chrN}.bai',
    ref=GDCH38,
    platypus='Tools/Platypus_0.8.1/Platypus.py'
  output:
    vcf=INTPATH+'/{case}/{file_id}.{chrN}.platypus.vcf'
  shell:
    '''
    module load python/2.7
    python {input.platypus} callVariants --refFile {input.ref} --bamFiles {input.bam} --minReads 8 --output {output.vcf}
    '''  
#minReads 8 is too friendly

rule index:
  priority: 2
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
  priority: 1
  input:
    bamdir=expand('{datapath}/{{case}}/{{file_id}}/', datapath=DATAPATH),
  output:
    bam=expand('{intpath}/{{case}}/{{file_id}}.{{chrN}}.bam', intpath=INTPATH)
  params:
    chrP='{chrN}'
  shell:
    '''
    module load samtools
    samtools view -b {input.bamdir}/*.bam {params.chrP} > {output.bam}
    '''

rule vcfGzip:
  priority: 2
  input: '{path}_somatic_{tool}.{chrN}.vcf'
  output:'{path}_somatic_{tool}.{chrN}.vcf.gz'
  wildcard_constraints:
    tool="pindel|abra|platypus"
  shell:
    '''
    gzip {input}
    '''
