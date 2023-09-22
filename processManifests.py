#tar -xf biospecimen.cases_selection.2023-08-25.tar.gz
#tar -xf clinical.cases_selection.2023-08-25.tar.gz 

import polars as pl

aliquots = pl.read_csv('Data/Manifests/aliquot.tsv', separator='\t')
samples = pl.read_csv('Data/Manifests/sample.tsv', separator='\t')
Ccases = pl.read_csv('Data/Manifests/clinical.tsv', separator='\t')
manifest = pl.read_csv('Data/Manifests/gdc_manifest.2023-08-25.txt', separator='\t')

CcasesFinal = Ccases[["case_id","primary_diagnosis", "tissue_or_organ_of_origin"]]

aExtended = aliquots.join(samples, [a for a in aliquots.columns if a in samples.columns])

manEx = manifest.with_columns((pl.col('id')+'/'+pl.col('filename')).alias('path'), pl.col('filename').str.extract(r"^(.*)_wgs", 1).alias('aliquot_id'))

aliEx = manEx.join(aliquots, [a for a in aliquots.columns if a in manEx.columns])
samEx = aliEx.join(samples, [a for a in aliEx.columns if a in samples.columns])
cases = samEx[['case_id', 'path', 'sample_type']].groupby(['case_id', 'sample_type']).all().pivot(index='case_id', columns='sample_type', values='path')

Tcases = cases.filter(pl.col('Primary Tumor').is_not_null())

casesFinal = Tcases.join(CcasesFinal, "case_id")

casesFinal.write_ndjson('Intermediate_Data/CPTAC/manifest.json')
