#!/bin/bash

###
#		Call pairwise differences on many samples
###

matrix=/full/path/to/SNP/matrix
results=/full/path/to/results/folder
mkdir -p ${results}
mask=/full/path/to/mask/file
java -Xmx4G -classpath ~/java/lib/junit.jar:~/java/lib/jbzip2-0.9.1.jar:~/java/lib/args4j-2.0.12.jar:~/java/lib/commons-compress-1.0.jar:~/java/lib/gson-1.6-javadoc.jar:~/java/lib/gson-1.6-sources.jar:~/java/lib/gson-1.6.jar:bin c.e.data_processing.Pairwise_shore_clean ${matrix} ${1} ${results} ${mask}

