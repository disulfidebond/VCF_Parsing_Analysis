# start with terms that are already present in the header line
baylor14_IDs = {'#CHROM': 'CHROM', 'POS':'POS', 'ID':'ID', 'REF':'REF', 'ALT':'ALT', 'QUAL':'QUAL', 'FILTER':'FILTER', 'INFO':'INFO', 'FORMAT':'FORMAT'}

# then create a list of identifiers. The format for the example file is
# GroupInt,VCF_column_name
# 0,22074
# 1,22075

# finally, import the data from the list, and modify the terms appropriately
with open('baylor_14_mappedIDs.csv', 'r') as fOpen:
    headerLine = True
    for i in fOpen:
        i = i.rstrip('\r\n')
        iSplit = i.split(',')
        if headerLine:
            headerLine = False
            continue
        else:
            valInt = int(iSplit[0])
            if valInt == 0:
                baylor14_IDs[str(iSplit[1])] = 'Control-' + iSplit[1]
            elif valInt == 1:
                baylor14_IDs[str(iSplit[1])] = 'Treated_LD10-' + iSplit[1] 
            elif valInt == 2:
                baylor14_IDs[str(iSplit[1])] = 'Treated_LD10_LD50-' + iSplit[1]
            elif valInt == 3:
                baylor14_IDs[str(iSplit[1])] = 'Treated_LD50_LD90-' + iSplit[1]
            elif valInt == 4:
                baylor14_IDs[str(iSplit[1])] = 'Treated_gt_LD90-' + iSplit[1]
            else:
                print('Warning, unknown value')
                print(i)
                break

# finally, rename the columns
vcfDF = vcfDF.rename(columns=baylor14_IDs)
