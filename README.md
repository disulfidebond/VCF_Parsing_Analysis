## Overview
Before beginning, it is **strongly** recommended to review the [VCF 4.0 specifications](https://samtools.github.io/hts-specs/VCFv4.2.pdf). It is also **strongly** recommended that the reader has a basic grasp of Biological concepts, and/or has [read through this guide](https://github.com/disulfidebond/Biology-for-ComputerScientists). The comments and code contained in this writeup will assume a basic understanding of the VCF format, and the data that is reported in VCF files.

Before continuing, I also should note that [GATK Tools](https://software.broadinstitute.org/gatk/download/) or [picard](https://broadinstitute.github.io/picard/), or the [IGV Viewer](https://software.broadinstitute.org/software/igv/download) from the Broad Institute may be completely sufficient for your needs. A description of the various commands for these software applications is well beyond the scope of this writeup, but along with Python and Bash commands, they are periodically mentioned to provide an in-depth description of what information can be contained within a VCF file. 

## Methods
The VCF format is similar to the Library of Congress: there's a wealth of information contained therein, but you need to understand how the [Classification System works](https://www.loc.gov/catdir/cpso/lcco/), or you're going to have an extremely difficult time finding anything. A Variant Caller Format (VCF) file contains numerous information lines (called 'meta-information' lines), which describe terms that will be used in the body of the VCF file. 

### Methods Introduction
At a very basic level, a VCF file has the following format:

**Description of Stuff** <- Meta-information lines

**Header line** <- Header Line

**Data that is formatted using the syntax described in the Description of Stuff** <- The actual data

The Description of Stuff section must follow the framework described in the VCF Standards, and should flow linearly down the document, which means a VCF file can contain information from multiple analyses:

* **Description of Stuff**
  * Description of what the Stuff means
  * Description of Stuff from First Analysis
  * Description of Stuff from Second Analysis
  * Description of Stuff from Third Analysis
* **Header line**
* **Data that is formatted using the syntax described in the Description of Stuff**


There is no limit to the number of these lines, and in the Cyno.box13.indels.hg38.vcf.gz gzipped VCF file, there are 506 lines of meta-information, and one header line. For example, in the vcf file Cyno.box13.snps.hg38.vcf.gz, this means that the first 507 lines of the document will contain information on how to read the rest of the document. VCF files are technically human readable, but realistically, you should write or [use a program](https://software.broadinstitute.org/software/igv/download) to read them for you.

If you want to remove the metadata and header lines, you could do:


        # bash
        DOCLEN=$(gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | wc -l )
        # find position where metadata and header end
        gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | grep -n 'chr1' 
        # TRIMSTART=((DOCLEN-HEADERSTOP))
        # in this case, TRIMSTART == line 895, where you will set the tail command:
        gzcat Cyno.box13.snps.hg38.geneListFiltered.vcf.gz | tail -n 895 | less
        
For filtering based on specific criteria, [you can also use GATK Tools](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/), [as well as picard](https://broadinstitute.github.io/picard/javadoc/picard/index.html?picard/vcf/filter/FilterVcf.html).

### Methods I
For something more concrete, consider the following segment from a VCF file.

**Metadata Lines**, More on this in a moment.

         ##fileformat=VCFv4.2
         ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
         ##FILTER=<ID=LowQual,Description="Low quality">
         ##FILTER=<ID=PASS,Description="All filters passed">
         ##FILTER=<ID=hardFilter.snp,Description="QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">
         ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
         ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
         ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
         ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
         ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
         ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
         ##GVCFBlock=minGQ=0(inclusive),maxGQ=1(exclusive)
         ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
         ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
         ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
         ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
         ##INFO=<ID=CCC,Number=1,Type=Integer,Description="Number of called chromosomes">
         ##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
         ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
         ##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
         ##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
         ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
         ##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of all GQ values">
         ##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description="Standard deviation of all GQ values">
         ##INFO=<ID=HWP,Number=1,Type=Float,Description="P value from test of Hardy Weinberg Equilibrium">       
         ##contig=<ID=chr11,length=135086622>
         ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">

**Header Line**
        
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  40156   40157   40158   40159

**Data**

        chr11   108227883       14:104278100:T:G        T       G       25611.30        PASS    AC=16;AF=0.114;AN=140;BaseQRankSum=1.14;ClippingRankSum=-0.328;DP=3880;FS=0;GQ_MEAN=415.56;GQ_STDDEV=590.36;InbreedingCoeff=-0.129;MLEAC=16;MLEAF=0.114;MQ=60;MQ0=0;MQRankSum=-0.259;NCC=0;QD=14.84;ReadPosRankSum=-0.704;SOR=0.72;CSQ=G|synonymous_variant|LOW||CCDS31669.1|Transcript|CCDS31669.1|protein_coding|2/62||||180|180|60|V|gtT/gtG|rs786201375||1||SNV|||YES|||CCDS31669.1|CCDS31669.1||||||||||||||||||||||||||||likely_benign||1|||||,G|synonymous_variant|LOW|ATM|ENSG00000149311|Transcript|ENST00000278616|protein_coding|3/63||||565|180|60|V|gtT/gtG|rs786201375||1||SNV|HGNC|HGNC:795|YES|5|P1|CCDS31669.1|ENSP00000278616|Q13315|A0A024R3C7|UPI000016B511|1|||PDB-ENSP_mappings:5np0.A&PDB-ENSP_mappings:5np0.B&PDB-ENSP_mappings:5np1.A&Pfam_domain:PF11640&hmmpanther:PTHR11139&hmmpanther:PTHR11139:SF96&SMART_domains:SM01342|||||||||||||||||||||likely_benign||1|||       GT:AD:DP:GQ:PL  0/0:37,0:37:99:0,102,1530       0/1:69,51:120:99:1619,0,2107    0/1:59,58:117:99:1773,0,1486    0/1:45,46:91:99:1404,0,1497
        
Even though this is only a snippet of a VCF file, at first glance, understanding what is happening can be quite arduous. Generally speaking, always start at the header line, which should start with '#CHROM'. This will tell you what is in the document, and how to decipher it. The example above has CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT fields.  

The numbers after 'FORMAT' are sample identifiers that are included in the VCF file, which can be a string, int, or float, provided each are unique. The 'CHROM' identifier denotes the chromosome number. For a quick primer on basic Biological concepts, [see this writeup](https://github.com/scienceystuff/Biology-for-ComputerScientists/blob/master/Biology-key-concepts.md#dna). The POS describes the position on the Chromosome. The 'FILTER' indicates whether the software that was used for the Variant Calling passed or failed the test for validity. 

The primary purpose of this data is to determine whether a given sample matches the string/character at that position in the reference file, or whether it has an alternative nucleotide at that position. 'REF' and 'ALT' indicate what the nucleotide at that position was recorded as in the reference (REF) file, and what an alternative (ALT) nucleotide that is contained in at least one sample could be at tha position. As an exercise, look up what the REF, ALT, INFO, and QUAL identifiers indicate in the VCF specifications file.

The 'INFO' contains a lot of information that may seem garbled at first glance, but it is actually very ordered. Each of the semicolon-delimited fields contains information on the variant call for the REF vs ALT nucleotide, and each is described in the INFO field. This can be followed by additional information, which in this case is delimited by '|'. If you look through the semicolondelimited INFO fields above (after the '##GVCFBlock line), you can trace what each means. 

* 'AC' -> Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed"
  * Number can be Any, the type is Integer, and the description is shown
* 'AF' -> Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed"
  * Number can be Any, the type is Float, and the description is shown
* 'AN' -> Number=1,Type=Integer,Description="Total number of alleles in called genotypes"
  * Only a single integer is allowed, the description is shown
* etc...

Looking at the data line above, it has an allele count of 16, an alleleic frequency of 0.114, and a total of 140 alleles.  As an exercise, review the remaining fields for the INFO.

The '|' delimited fields are also described. Look for the last metadata line, that starts with:

         ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: 
         
This line indicates a second analysis, Variant Effector Prediction (VEP), was completed on this dataset as well. The fields are described, and start with the Allele, then Consequence. An important and slightly confusing part here is this info line describes the predictions and effects of the **listed allele** at the beginning of this line ('T' changing to 'G'). It begins with the effects on the protein coding region of Gene CCDS31669.1, then continues with the affects on Gene ENSG00000149311. As a sidenote, there are different online data repositories for human/other genes, and the information contained in one may not be equivalent to another, thus requiring multiple descriptions. For a recap on how a single position can affect multiple genes and eukaryotic splicing, [see this guide](https://github.com/disulfidebond/Biology-for-ComputerScientists).

The last set of fields is also highly ordered. 

        GT:AD:DP:GQ:PL  0/0:37,0:37:99:0,102,1530       0/1:69,51:120:99:1619,0,2107    0/1:59,58:117:99:1773,0,1486    0/1:45,46:91:99:1404,0,1497
        
The first entry, 'GT:AD:DP:GQ:PL', indicates what data is contained and for each sample.  For the samples 40156, 40157,40158, and 40159, the following information is provided:

Genotype:AlleleicDepth:ReadDepth:GenotypeQuality:PhredLikelihood

See the INFO lines above for more inforamtion on each. On a final note, the Genotype will be a string with the format 'N/N', where N is an integer (usually 0 or 1), which describes the predicted genotype of that sample.

### Methods II
To view a VCF file, you could use less:

        gzcat some_compressed_vcf.vcf.gz | less

You could also import it into pandas:

        import pandas as pd
        import numpy as np
        
        vcfDF = pd.read_csv(some_uncompressed_vcf.vcf', sep='\t')
        vcfDF.head()
        
Note that if you import it into pandas, the metadata and header lines may cause the table to be malformed.  These lines can be removed with something similar to:

        # bash
        DOCLEN=$(gzcat some_compressed_vcf.vcf.gz | wc -l)
        # find position where metadata and header end
        gzcat some_compressed_vcf.vcf.gz | grep -n 'chr1' 
        # in this case, DOCLEN-HEADERSTOP == line 895, where you will set the tail command:
        gzcat some_compressed_vcf.vcf.gz | tail -n 895 | less
        
If you only wanted to view the INFO column, you could do something similar to:

        # bash
        gzcat some_compressed_vcf.vcf.gz | tail -n 895 | cut -d$'\t' -f8 | less
        
### Methods III
To parse the data, you could use python pandas, and verify it with Bash:

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

From here, continue to the [Analysis and Visualization]() writeup.
