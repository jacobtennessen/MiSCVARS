# MiSCVARS
Miscellaneous Scripts for Calling Variants by Aligning to Reference Sequence

MakeOnemapFormatFromVcfNoParentsGeneralFormat.pl - reads a vcf file and converts it to OneMap format. Designates all sites as heterozygous in one parent (the same parent); other site types may be included, but these will need to be changed by editing the OneMap output file (e.g. based on known parent genotypes)

MakeOnemapFormatFromVcfNoParentsTwoGenos.pl - reads a vcf file and converts it to OneMap format, only for sites inferred to be heterozygous in one parent (designed for a mapping situation in which offspring share maternal parent but multiple unknown fathers)

FilterSplitFastqs.pl - reads fastq files that consists of multiple samples identified in read title. Splits samples into individual files and renames them. Designed for output of dbcAmplicons from IBEST

MakeHTMLfromOneMapSex.pl- reads OneMap format file and coverts to color-coded HTML for easy viewing. Optionallly reads a file of sex phenotypes and adds in that information.
