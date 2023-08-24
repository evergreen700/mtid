# mtid
Multi-tool indel detection

Current snakemake rules:
  - test: tests current progress. Example run: 
    ```snakemake -n test --cluster="sbatch --time 4:00:00 --mem 32G --output Log/slurm-%j.out" --jobs=10 --rerun-incomplete```

Current progress:
  - bamsplit: splits initial files by chromosome to make files more palatable.
  - index: indexes newly created bam files from samtools and abra2
