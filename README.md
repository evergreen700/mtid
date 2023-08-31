# mtid
Multi-tool indel detection

Current snakemake rules:
  - week1Benchmark: tests current progress. Currently tries to generate Pindel and Abra somatic vcfs. Not all Abra rules are fully working but it's useful for making sure all rules eventually feed into each other. Example run: 
    ```snakemake -n week1Benchmark --cluster="sbatch --time 4:00:00 --mem 32G --output Log/slurm-%j.out" --jobs=10 --rerun-incomplete```

Current progress:
  - bamsplit: splits initial files by chromosome to make files more palatable.
  - index: indexes newly created bam files from samtools and abra2
  - pindel: runs pindel on one bam/bai combo. Generates a set of output files
  - pindeljoin: concatonates the pindel files into one big file
  - pindelvcf: takes the file from pindeljoin and converts it into a VCF
  - somaticFromGermlinePindel/Platypus: takes Tumor/Normal vcf pairs and returns a vcf of mutations not seen in normal tissue
  - ABRAlignSomatic: realigns bam files in tumor/normal pairs
  - CadabraVCF: somatic mutation caller that uses the flags created by ABRA
  - RUFUS: somatic mutation caller that uses k-mers
  - platypus: other mutation caller that Ohio State people recommended

MISC NOTES:
  - Find read length of file: ```samtools view file.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | less -S```
