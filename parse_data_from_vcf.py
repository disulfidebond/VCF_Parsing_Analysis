import pandas as pd
import numpy as np

vcfDF = pd.read_csv('Cyno.box13.snps.hg38.geneListFiltered.vcf.tsv', sep='\t')
vcfDF.head()

# create list with column IDs
colIDS = ['FORMAT']
for i in range(40156, 40226):
    iString = str(i)
    colIDS.append(iString)
print(colIDS)
vcf_samplesDF = pd.DataFrame(vcfDF, columns=colIDS)
# preview data
vcf_samplesDF.head()

# filter out only the Genotype Tag:
rVals = []
for index,row in vcf_samplesDF.iterrows():
    rVal = []
for cIdx in colIDS:
    rValueRow = row[cIdx]
    rValueList = rValueRow.split(':')
    rValue = rValueList[0]
    rVal.append(rValue)
rVals.append(rVal)

parsed_samplesDF_vcf = pd.DataFrame(rVals, columns=colIDS)
# preview data
parsed_samplesDF_vcf.head()

# check the data:
parsed_samplesDF_vcf['40156'].value_counts()
# outputs 
0/0    505
0/1    168
1/1    152
1/0     62
./.      7
Name: 40156, dtype: int64

# check again with Bash from Jupyter Notebook, remove the '!' to run in a Bash session
!gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | tail -n 895 | cut -d$'\t' -f10 | grep -c '^0/0'
!gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | tail -n 895 | cut -d$'\t' -f10 | grep -c '^0/1'
!gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | tail -n 895 | cut -d$'\t' -f10 | grep -c '^1/1'
!gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | tail -n 895 | cut -d$'\t' -f10 | grep -c '^1/0'

# outputs
505
168
152
62
