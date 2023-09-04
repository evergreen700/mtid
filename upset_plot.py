import polars as pl
import upsetplot as up
from matplotlib import pyplot

tools = ["PINDEL", "ABRA", "RUFUS", "PLATYPUS"]

rawTable = pl.read_csv("chr22.tsv", separator="\t", skip_rows=2)
midTable = rawTable.with_columns(pl.all().exclude(["UNIQUE_ENTRIES", "SHARED_ENTRIES"]) != 0)

for t in tools:
   rawTable=rawTable.with_columns(pl.when(pl.col(t)).then(pl.lit(t)).otherwise(pl.lit(None)).alias(t))

finalTable=midTable.group_by(pl.all().exclude(["UNIQUE_ENTRIES", "SHARED_ENTRIES"])).mean().with_columns(pl.col('SHARED_ENTRIES').cast(int))

categories=[]
counts=[]
for row in finalTable.iter_rows():
  labels=[]
  for i in range(0,len(tools)):
    if row[i]:
      labels.append(tools[i])
  categories.append(labels)
  counts.append(row[-1])

data = up.from_memberships(categories, data=counts)
up.plot(data)

pyplot.show()
