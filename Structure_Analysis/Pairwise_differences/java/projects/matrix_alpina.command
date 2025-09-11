#!/bin/bash


###
#		Transform a VCF into a SNP matrix
#		Run with chromosome number from standard input
###
vcf=/full/path/to/vcf/file
matrix=/define/a/full/path/to/the/SNP/matrix/file/you/want/to/create

java -Xmx4G -classpath ~/java/lib/junit.jar:~/java/lib/jbzip2-0.9.1.jar:~/java/lib/args4j-2.0.12.jar:~/java/lib/commons-compress-1.0.jar:~/java/lib/gson-1.6-javadoc.jar:~/java/lib/gson-1.6-sources.jar:~/java/lib/gson-1.6.jar:bin c.e.data_processing.VcfCombined_to_snpMatrix_alpina ${vcf} ${1} ${matrix}
