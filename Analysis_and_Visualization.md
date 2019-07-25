## Analysis and Visualization

Start by importing the VCF file into a pandas dataframe:


        import pandas as pd
        import numpy as np
        
        vcfDF = pd.read_csv('Cyno.box13.snps.hg38.geneListFiltered.vcf.tsv', sep='\t')
        vcfDF.head() # optional, preview data frame
        
        colIDS = ['FORMAT']
        for i in range(40156, 40226):
            iString = str(i)
            colIDS.append(iString)
        vcf_samplesDF = pd.DataFrame(vcfDF, columns=colIDS)
        
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
        
Then, start to get statistics on the data:

        # cols_rawData will be a list of tuples with the format [('1/1', N),('1/0', N),('0/1', N),('0/0', N),('./.', N)]
        cols_rawData = [] 
        for column in parsed_samplesDF_vcf:
            if column == 'FORMAT': # skip this column
                continue
            k = parsed_samplesDF_vcf[column].value_counts().keys().to_list()
            v = parsed_samplesDF_vcf[column].value_counts().to_list()
            zippedItems = list(zip(k,v))
            zippedItems.sort(key=lambda tup: tup[0], reverse=True)
            cols_rawData.append(zippedItems)
