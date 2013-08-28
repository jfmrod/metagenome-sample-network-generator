metagenome-sample-network-generator
===================================

A tool to generate a network of metagenomic samples with link lengths adjusted according to the number of taxa shared.
This tool takes as input the sets of OTUs identified in one or more samples. From this data the tool computes
the fraction of shared OTUs between all samples and generates a (core) network in which the link lengths are adjusted to reflect
the fraction of shared OTUs, shorter link lengths representing more similar samples. Optionally, a set of reference samples
can be given and a network generated which includes the previous core network plus links and samples that are found to be
sufficiently similar to the core samples. In this case, no links are generated between reference samples.


REQUIREMENTS
=======================

This tool requires the neato program found in the graphviz package.
This package can be downloaded from: http://www.graphviz.org/



CORE SAMPLE FILE
=======================

The set of OTUs observed in each sample should be provided on a single file per sample, and as one sequence and OTU per line.
The format of a single line should be:

SEQID <TAB> OTUID  

This should be prepared previously from the sequence read data, the first column should represent the sequence read ID followed
by a tab character and the assigned OTU. The OTU id can be anything from a numerical id to the inferred taxonomy of the sequence.

For example, the contents of the example/sample1.otu file are:
seq1    OTU1  
seq2    OTU2  
seq3    OTU2  
seq4    OTU1  
seq5    OTU3  
...  


REFERENCE SAMPLES FILE
=======================

The set of reference samples should be provided in the following format. Where lines indicating the start of a new OTU begin with the '>'
character. And every line that does not begin with the '>' character indicates the sequence ID and the sample it belongs too.

>OTUID  
SEQID <TAB> SAMPLEID  
SEQID <TAB> SAMPLEID  
SEQID <TAB> SAMPLEID  
>OTUID  
SEQID <TAB> SAMPLEID  
SEQID <TAB> SAMPLEID  
SEQID <TAB> SAMPLEID  

For example, the contents of the example/refsamples.otu file are:
>OTU1  
seq1.1	osample1  
seq2.1	osample2  
seq3.1	osample3  
>OTU2  
seq1.2	osample1  
seq1.3	osample1  
seq2.2	osample2  
seq2.3	osample2  
...  

It is important that the OTU IDs used in the reference file are the same as the ones used in the classification of the (core) sample files.


ATTRIBUTES FILE FORMAT
=======================

The attributes file is used to customize the appearance of the nodes in the generated network. For example, this can be used
to change the color, shape or size of the nodes according to the sample types or any property of the data associated
with the samples. This applies equally to core as well as reference samples.

Any property of the nodes can be changed provided that it is a recognized node property of the neato program.
Please consult the graphviz documentation for the relevant information.

The file format of the attributes file is:

>SAMPLEID  
property=value  
property=value  
property=value  
>SAMPLEID  
property=value  
property=value  
property=value  

Where SAMPLEID is the filename of the sample provided to neato if it is a core sample, or, in the case of a reference sample, it is simply the
SAMPLEID provided in the file.

For example, the contents of the example/samples.attr file are:
>sample1.otu  
shape=box  
color="#33FF33"  
label="sample1"  
>sample2.otu  
shape=box  
fillcolor="#3333FF"  
label="sample2"  

These properties are set for each sample and passed along to neato when the script generates the .dot file.


USAGE
=======================

syntax:  
  make-sample-network.sh [-a samples.attrs] [-r refsamples.otu] <sample1.otu> [sample2.otu] [...]  

The core samples can be provided as filenames after the command "make-sample-network.sh".
Optionally, the reference samples file can be provided using the "-r" argument.
Another optional argument is the "-a" which is used to provide the customized sample attributes.


EXAMPLE
=======================

To produce the example network file one should run the following commands:
./make-sample-network.sh -a example/samples.attr -r example/refsamples.otu example/sample*.otu > example/network.dot  

This generates the .dot file that can then be used with neato to generate a pdf file, using the command:

neato -Tpdf -O example/network.dot  

The generated network can then be found in the file example/network.dot.pdf


