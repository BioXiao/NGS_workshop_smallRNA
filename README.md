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
