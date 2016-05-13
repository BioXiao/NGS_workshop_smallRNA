# Small RNA analysis
Practical for the smallRNA section of the NGS workshop.

While useful for those interested in small RNA analysis, this section will provide an example of what goes into the running of a typical bioinformatics pipeline.

____

## Before we get started...

+ \* -- wild card
+ # -- comment out
+ Relative paths versus absolute paths
+ On command line notation
+ On /proj/seq/data
+ On annotation files


## Setup

#### Log in

Login to Killdevil with your onyen

#### Load modules

The clusters (Killdevil and Kure) have many programs and applications already installed.  To call on them, you simply need to call on them.  Let's take a look at which applications are available to us, type:

```
$ module avail
```
miRQuant requires several applications.

The required applications are:  
cutadapt  
blast  
bedtools  
SHRiMP  

While modules can be loaded individually by typing:
```
$ module load blast
$ module load bedtools
ect...
```
this can be tedious when you have many modules required for a pipeline.  One alternative is to create a file containing the environmental variables (modules, paths, ect.) we want to use for the pipeline.  The uncENV.sh has all the necessary variables we need for the pipeline!  Lets take a look at this file.

`$ cat uncENV.sh`

uncENV.sh contains the location of programs used by the pipeline, the location of various files called on by the pipeline, and the commands to load the necessary modules.  We can execute the code within this file by typing `source uncENV.sh`.  If we type `module list` we should see the following output, confirming the necessary modules were loaded.

```
Currently Loaded Modulefiles:
  1) null              2) r/2.15.1          3) bowtie/1.1.0      4) samtools/1.3      5) bedtools/2.25.0   6) bamtools/1.0.2    7) perl/5.12.0
```

