

#before running this program, the following commands may need to be run on the files.
#tar -xf biospecimen.cases_selection.*.tar.gz
#tar -xf clinical.cases_selection.*.tar.gz 

import polars as pl
import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--in_directory", help="path to manifest directory", required=True)
parser.add_argument("-o", "--out_file", help="output location. Otherwise will write to stdout")

args = parser.parse_args()

inDirectory=args.in_directory
if inDirectory[-1] =="/":
  inDirectory = indirectory[:-1]

outPath=args.out_file

if inDirectory == None:
  print("Path to manifest directory required")
  parser.print_help()
  sys.exit(1)
try:
  aliquots = pl.read_csv(inDirectory+'/aliquot.tsv', separator='\t')
  samples = pl.read_csv(inDirectory+'/sample.tsv', separator='\t', infer_schema_length=10000)
  Ccases = pl.read_csv(inDirectory+'/clinical.tsv', separator='\t', infer_schema_length=10000)
  manifest = pl.read_csv(inDirectory+'/gdc_sample_sheet*.tsv', separator='\t')
except:
  print("This program requires the existance of the following files within the manifest directory (-o flag)")
  print(" - aliquot.tsv")
  print(" - sample.tsv")
  print(" - clinical.tsv")
  print(" - a gdc_sample_sheet*.tsv file. Only one manifest file allowed in the directory")
  sys.exit(1)

CcasesFinal = Ccases[["case_id", "case_submitter_id","primary_diagnosis", "tissue_or_organ_of_origin"]]

#aExtended = aliquots.join(samples, [a for a in aliquots.columns if a in samples.columns])

manEx = manifest.with_columns(pl.col('Case ID').str.split(by=", ").list.first().alias("case_submitter_id"))

aliEx = manEx.join(aliquots, [a for a in aliquots.columns if a in manEx.columns])
samEx = manEx.join(CcasesFinal, [a for a in manEx.columns if a in CcasesFinal.columns])
samReduced = samEx.with_columns(pl.when(pl.col("Data Type")=="Aligned Reads").then("Aligned Reads- "+pl.col("Sample Type")).otherwise(pl.col("Data Type")).alias("fulltype"))[["File ID", "File Name", "case_id","fulltype", "primary_diagnosis","tissue_or_organ_of_origin"]]


cases = samEx[['id','case_id', 'path', 'sample_type']].groupby(['case_id', 'sample_type']).all().pivot(index='case_id', columns='sample_type', values='path')

Tcases = cases.filter(pl.col('Primary Tumor').is_not_null())

casesFinal = Tcases.join(CcasesFinal, "case_id")

if outPath!=None:
  casesFinal.write_ndjson(outPath)
else:
  print(casesFinal.write_ndjson())

