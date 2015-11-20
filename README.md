# MiSCVARS
Miscellaneous Scripts for Calling Variants by Aligning to Reference Sequence

MakeOnemapFormatFromVcfNoParentsTwoGenos.pl - reads a vcf file and converts it to OneMap format, only for sites inferred to be heterozygous in one parent (designed for a mapping situation in which offspring share maternal parent but multiple unknown fathers)

FilterSplitFastqs.pl - reads fastq files that consists of multiple samples identified in read title. Splits samples into individual files and renames them. Designed for output of dbcAmplicons from IBEST
