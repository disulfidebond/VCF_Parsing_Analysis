# first, create a list of column identifiers. In the example VCF file, the sample column names are consecutive integers greater than 40227 and less than 40297
colIDS = ['FORMAT']
for i in range(40227, 40297):
    iString = str(i)
    colIDS.append(iString)
print(colIDS)
# create a dataframe from vcfDF dataframe, only use the column identifiers in the colIDS list
vcf_samplesDF = pd.DataFrame(vcfDF, columns=colIDS)
# preview the dataframe
vcf_samplesDF.head()

# then, pull the genotype data from the column
# cols_rawData will be a list of lists with the format [('1/1', N),('1/0', N),('0/1', N),('0/0', N),('./.', N)]
cols_rawData = [] 
for column in parsed_samplesDF_vcf:
    if column == 'FORMAT': # skip this column
        continue
    k = parsed_samplesDF_vcf[column].value_counts().keys().to_list()
    v = parsed_samplesDF_vcf[column].value_counts().to_list()
    zippedItems = list(zip(k,v))
    zippedItems.sort(key=lambda tup: tup[0], reverse=True) # keep a consistent order of tuples within each list
    cols_rawData.append(zippedItems)

# grab all phased '1/1' results
filtered_controlList_1_1 = []
for v in filtered_controlList:
    tmp = filter_from_listOfTuples(v, '1/1')
    filtered_controlList_1_1.append(tmp)
 
phased_genotype = [x[1] for x in filtered_controlList_1_1])
# count the number of phased genotypes, and get mean, etc.
reduce(lambda x, y: x + y, phased_genotype) / len(phased_genotype) 
