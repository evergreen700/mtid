import polars as pl
import upsetplot as up
from matplotlib import pyplot
import sys

inTSV = sys.argv[1]
outPNG = sys.argv[2]

tools = ["PINDEL", "ABRA", "RUFUS", "PLATYPUS"]

rawTable = pl.read_csv(inTSV, separator="\t", skip_rows=2)
midTable = rawTable.with_columns(pl.all().exclude(["UNIQUE_ENTRIES", "SHARED_ENTRIES", "EXCLUSIVE_ENTRIES"]) != 0)
finalTable=midTable.group_by(pl.all().exclude(["UNIQUE_ENTRIES", "SHARED_ENTRIES", "EXCLUSIVE_ENTRIES"])).mean().with_columns(pl.col('SHARED_ENTRIES').cast(int))

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
up.plot(data, show_counts=True, min_subset_size=1)

pyplot.savefig(outPNG)
