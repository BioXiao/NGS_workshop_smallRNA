# Small RNA analysis
Practical for the smallRNA section of the NGS workshop.

While useful for those interested in small RNA analysis, this section will provide an example of what goes into the running of a typical bioinformatics pipeline.

____

## Before we get started...

+ \* -- wild card
+ # -- comment out
+ . -- current working directory
+ && -- chains commands together; do this command, then this one
+ Tab-completion
+ Relative paths versus absolute paths
+ On command line notation
+ On /proj/seq/data
+ On annotation files


## Setup

#### Log in

Login to Killdevil with your onyen.

#### Change to your scratch directory

Switch from the home directory to your scratch directory, by typing `$ cd /netscr/ONYEN`, where ONYEN is your UNC onyen.

#### Copy scripts and data to your scratch directory

The necessary miRquant scripts, resource files, and data have been compressed and can be copied from /netscr/mattkank/miRquant.tar.gz.  To copy the files over and un-compress them, we will submit the job to LSF by typing:

```
$ bsub 'cp /netscr/mattkank/miRquant.tar.gz . && tar -zxvf miRquant.tar.gz'

$ cd miRquant
$ pwd
/netscr/ONYEN/miRquant
```

#### Load pipeline requirements

The clusters (Killdevil and Kure) have many programs and applications already installed.  To call on them, you simply need to load the module.  Let's take a look at which applications are available to us, type:

```
$ module avail
```
miRquant requires several applications.

The required applications are:  
cutadapt  
blast  
bedtools  
SHRiMP  
R

While modules can be loaded individually by typing:
```
$ module load blast
$ module load bedtools
ect...
```
this can be tedious when you have many modules required for a pipeline.  One alternative is to create a file containing the environmental variables (modules, paths, ect.) we want to use for the pipeline.  The uncENV.sh has all the necessary variables we need for the pipeline!  Lets take a look at this file.
```
$ cat uncENV.sh

# Adds locations to system path, so the system knows where to look for programs
export PATH=/netscr/mattkank/miRquant/02_software/bin:./:$PATH
export PATH=$PATH:/netscr/mattkank/miRquant/02_software/SHRiMP_2_2_2/bin
export PYTHONPATH=/netscr/mattkank/miRquant/02_software/python/

# Sets location of required components to variables
export SHRIMP_FOLDER=/netscr/mattkank/miRquant/02_software/SHRiMP_2_2_2
export CHROMO_SIZE=/netscr/mattkank/miRquant/02_software/resources/mm9.chromSizes
export JB_GENOME_DIR=/proj/seq/data/mm9/bowtie_path/base/
export MMU_GENOME=/proj/seq/data/MM9_UCSC/Sequence/WholeGenomeFasta/genome.fa

# Loads the modules necessary for the pipeline
module load r/2.15.1
module load bowtie/1.1.0
module load samtools/1.3
module load bedtools/2.25.0
module load bamtools/1.0.2
module load perl/5.12.0
```
uncENV.sh contains the location of programs used by the pipeline, the location of various files called on by the pipeline, and the commands to load the necessary modules.  We can execute the code within this file by typing `source uncENV.sh`.  If we type `module list` we should see the following output, confirming the necessary modules were loaded.

```
Currently Loaded Modulefiles:
  1) null              2) r/2.15.1          3) bowtie/1.1.0      4) samtools/1.3      5) bedtools/2.25.0   6) bamtools/1.0.2    7) perl/5.12.0
```

OK, we should now be setup to run the pipeline on the non-descript SampleA.fastq.

## Running miRquant

Let's make sure first that we are in the correct location on the file system.

```
$ pwd
/netscr/ONYEN/miRquant
```

Take a look at what the directory structure is.

```
$ ls
01_logs      03_samples    05_scripts               chainSubmission.sh  process_all_summary2tab.pl  uncENV.sh
02_software  04_resources  06_final_output_scripts  post_runC.sh        runC.sh

where:
01_logs is where the pipeline logs will be stored as they are produced (currently empty)
02_software is where some of the programs necessary for the pipeline are stored
03_samples is where the sequencing data and pipeline outputs will be stored
04_resources is where miRNA resource files required for the pipeline are stored
05_scripts is where various scripts called on by the pipeline are stored
06_final_output_scripts is where scripts used for the final assembly of the data is stored
```

Next we'll take a look within our 03_samples folder.

```
$ cd 03_samples
$ ls
SampleA.adaptor  SampleA.fastq
$ head SampleA.adaptor
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
$ wc -l SampleA.fastq
400000 SampleA.fastq
$ head SampleA.fastq
@D1317JN1:309:C88CWACXX:1:1101:2446:2123 1:N:0:ATCACG
TGGAGTGTGACAATGGTGTTTTGGAATTCTCGGGTGCCAAGGAACTCCAGT
+
BB@FFDEDHHHFHIHIGHIJJJJJGHIJJJJJJJJIJJHHIJJJJJJJJHF
@D1317JN1:309:C88CWACXX:1:1101:2294:2127 1:N:0:ATCACG
GTCTACGGCCATACCACCCTGTGGAATTATCGGGTGCCAAGGAACTCCAGT
+
?=?DD=DD1C?ADE3CBEFI>AE3@??D*?DD?DDDDDDD9CDEDED7A@=
@D1317JN1:309:C88CWACXX:1:1101:2353:2170 1:N:0:ATCACG
AATACCGGGTGCTGTAGGCTTTTGGAATTCTCGGGTGCCAAGGAACTCCAG
```
For the fastq output, this is the information for our first two (and half of the third) reads.  For each read, the first line starts with a '@' and is followed by a read identifier containing various bits of information (instrument name, flowcell id, ect).  The second line contains the raw sequence letters.  The third line is a '+'.  The fourth line encodes the quality value for the sequence in the second line.  Each read contains a portion of the 3' adapter sequnce, which will have to be trimmed.

Change back to the main directory by typing `$ cd ..`.

#### Section 1: chainSubmission.sh

The first part of miRquant trims the adapter, aligns RNAs perfectly matching the genome with Bowtie, creates genomic 'windows' for the SHRiMP alignment using bedtools, and aligns small RNAs containing mismatches to these windows using SHRiMP.  This is something you'll come across often as you use various pipelines, the stringing together of well-written programs to achieve the goal (if a toilet is available, why dig a hole in the woods).  All these programs brought together and connected in a 'wrapper' script, in this case chainSubmission.sh.  To run chainSubmission.sh, enter:

```
$ bsub -o 01_logs/chainSubmission.log bash chainSubmission.sh mmu /netscr/ONYEN/miRquant/03_samples/SampleA.fastq
```

Pause for ~10 seconds, then confirm your job is running.

```
$ bjobs
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
186010  mattkan RUN   day        killdevil-l donor_pool2 *leA.fastq May 13 15:19
```
Continue to check until the job finishes, the type `$ more 01_logs/chainSubmission.log`.  This is the log file for the chainSubmission.sh run of the script.  In this file you'll see the output for the various programs which were 'wrapped' together in chainSubmission.sh.  Logs should be checked after every run, as this is where you are likely to catch errors which may have occured while running a pipeline.

#### Section 2 & 3: runC.sh & post_runC.sh

The second part of miRquant processes the results from the Bowtie alignment and SHRiMP alignment in the first stage and separates these based on length, internal nucleotide edits, and 3' additions.  This part also determines if a read corresponds to a known miRNA.  The third part combines the results from the second part into a single file.  To run this second part, be very mindful of tab-completion and enter:

```
$ bash runC.sh mmu /netscr/ONYEN/miRquant/03_samples/SampleA./IntermediateFiles/g1Results/CHR*.results
```
As you can see after entering this command, this script submits many jobs.  In the g1Results directory are the alignment results separated by chromosome.  The same analysis will be applied to each chromosome results, so instead of doing the analysis is sesquential order, each chromosome will be submitted as a individual job and the analysis will be done in parallel (which will speed up the process significantly).  Once each of these jobs is complete, enter:

```
$ bash post_runC.sh /netscr/ONYEN/miRquant/03_samples/SampleA./IntermediateFiles/g1Results/
```
In the previous step, runC.sh created results files on a chromosome by chromosome basis.  The above command combines the multiple chromosome results files into a single results files.  

#### Section 4: process_summary2tab.pl

The data has been processed at this point, but still isn't in a very human readable format.  To convert this data into a more usable format, enter:

```
$ bsub -o 01_logs/process_sum.log perl process_all_summary2tab.pl /netscr/ONYEN/miRquant/ mmu /netscr/ONYEN/miRquant/03_samples/SampleA./IntermediateFiles/g1Results/shift_summary.txt
```
This scripts combines the information into a more human readable format that can be imported into excel.  The files that are produced by this script can be seen by entering `$ ls /netscr/ONYEN/miRquant/03_samples/SampleA./`.  The files whose names begin with 'TAB' are the results files.  Lets look at the top of one, enter:

```
$ head 03_samples/SampleA./TAB_3p_summary_miR.txt
Name	tRNA	miRbaseOffset	Seed	Percent	Count	EM	E	T	AC	AT	G	AA	TT	TA	AG	CT	AAT	AAA	GT	TTT	TAT	CA	ATT	AAAT	GCT	GA	AC	TC	AGT	ACC	TTA	TGT	AGA	TAA	TTAT	TG	ATC	AAC	AATG	AAAA	ATG	CC	ACT	ATA	TTTT	AATT	TGA	ACCA	ACA	AAG	AAGT	GCGC	CG	TATT	GC	AATA	CCA	AACT	GCCC	CTC	AAGA	GGA	TAAT	GAC	CGT	ATTA	GAG	TCT	GGC	GAA	GCTA	ATCA	TCA
mmu-mir-122-5p	Mir122a	0	GGAGTGT	0.997220400488793	20451.544444444418030.8777777778	384.333333333333	513	1262.33333333333	21	79	26	52	2	11	2	7	10	18	10	1	0	0	4	0	1	2	3	1	00	2	2	1	0	1	0	1	2	0	20	1	2	0	0	0	0	0	0	1	01	0	0	1	0	0	0	0	1	0	00	0	0	1	1	0	0	0	0	0
mmu-mir-21-5p	Vmp1	0	AGCTTAT	1	2473.27777777778	2239.27777777778	17	9	170	24	0	9	1	1	00	0	0	0	0	0	0	0	0	0	10	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	1	0	0	0	0	10	0	0	0
mmu-mir-192-5p	Mir192	0	TGACCTA	0.83692628650904	2407	1639	246	385	24	17	3	1	11	40	15	1	42	0	3	0	2	1	1	0	0	1	12	0	0	3	0	0	1	2	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	10	0	1	0	0	0	0	0	0	0	00	0
mmu-mir-148a-3p	Mir148a	0	CAGTGCA	1	1489	1416	18	32	17	0	0	0	2	0	0	0	1	0	00	0	0	0	0	0	1	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	1	0	0	0	0	0	0	00	0	0	0	0	0	0	0	1	0	0
mmu-let-7b-5p	Mirlet7b	0	GAGGTAG	1	1432.36029411765	1131.61029411765	61.5833333333333	57.5	65.3333333333333	3	51.3333333333333	6	10.5	12.5	0	7	1	4	1.5	4.5	4	0	0	4.5	3	0	0	0	0	00	0	0	0	0	0.5	0	0	0	0	00	0	0	0	0.5	0	1	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	1	0	0	0	0	0	0.5	0
mmu-mir-22-3p	Mir22hg	0	AGCTGCC	1	1406	1338	19	5	33	0	5	2	2	0	0	0	0	1	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	1	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0
mmu-let-7f-2-5p	Huwe1	0	GAGGTAG	1	1240.26578947368	1200.81578947368	11.7	0.25	21	3	0	1	1.5	0	00	0	0	0	1	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0
mmu-let-7c-2-5p	Mirlet7c-2	0	GAGGTAG	1	1221.09362745098	1006.61029411765	131.2	25	38.5	1	0	6	3.5	4.2	02	0	0	0.5	0	0.5	0.5	0	0	0	00.75	0.333333333333333	0	0.5	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0
mmu-let-7g-5p	Wdr82	0	GAGGTAG	1	1162.53333333333	1093	11.5333333333333	3	40	2	7	2	2	0	0	00	0	0	1	0	0	0	0	0	1	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	00	0	0
```
This is the bulk of the pipeline, lets bring it all together and check how out sample sequencing analysis went.

## Bringing it all together - Final analyses scripts

These final scripts to get a complete output are in the '06_final_output_scripts', move there by entering `$ cd 06_final_output_scripts`.  We will also need to load a certain version of python, do this by entering `$ module load python/2.7.6`.  Before we start on the final steps, we will make sure our environment is ready to run the final scripts, enter:

```
$ ls
generate_adapter_files.py       genNormalRPM.py  lenDist.py
generate_mapping_info_table.pl  lenDistGraph.R   smRNAseq_correlation.R

$ module list
Currently Loaded Modulefiles:
  1) null              3) bowtie/1.1.0      5) bedtools/2.25.0   7) perl/5.12.0
  2) r/2.15.1          4) samtools/1.3      6) bamtools/1.0.2    8) python/2.7.6
```

#### Generate mapping statistics

COPY THE FINAL OUTPUTS
```cp -r /netscr/mattkank/miRquant/06_final_output_scripts/archived_results/ .```

Lets start by looking at some quality control statistics from the analyses.  The 'generate_mapping__info_table.pl' will bring together mapping statistics, including total mapping, miR mapping, and tRNA mapping.  Up to this point we have been using the absolute path for the files.  This script requires the relative path.  To run this script, enter:

```
$ perl generate_mapping_info_table.pl ../03_samples/SampleA./*.stats

$ cat MappingTableInfo.tsv

Sample name	SampleA
File name	/netscr/mattkank/miRquant/03_samples/test.fastq
Total reads	100000.00000000000000000000
Trimmed reads	90252.00000000000000000000
% Trimmed reads	90.252%
Short reads	7012.00000000000000000000
% Short	7.012%
Exact match to genome	57277
% EM	63.4634135531623%
No exact match to genome	29121
% NEM	32.2663209679564%
Total mapped reads	 61752.00247733
% Mapped	68.4217551714422%
Total mapped to miRs	 46831.3644750544
% of total mapped to miRs	75.8378070286009%
Total mapped to tRNAs	 353
% of total mapped to tRNAs	0.571641381394216%
```
For most samples we see the following: the % trimmed around 90%, % too short around 5-10%, % mapped around 80-85%, % mapped to miRs around 70-75%, and the % mapped to tRNAs under 10%.  If the statistics vary greatly from these expected numbers, it could be that the samples are of low quality or some other interesting biological process is occuring.

#### Read length distribution

Next we will look at the read length distribution for our sample.  As miRNAs are the major species of small RNA and the length of a mature miRNA is ~22 nucleotides long, we will expect a spike at 22 nucleotides.  To check the length distributions, enter:

```
$ python lenDist.py /netscr/mattkank/miRquant/03_samples/SampleA./IntermediateFiles/*1.fq --image
Loading required package: methods
null device 
          1 
$ cat lenDist.tsv
	SampleA._O10_E1.fq
	test._O10_E1.fq
14	0.0292957496787
15	0.0356778797146
16	0.0226698577317
17	0.0270686522182
18	0.0372512520498
19	0.0352014359793
20	0.039500509684
21	0.123565128751
22	0.264758675708
23	0.135996986216
24	0.0442981872978
25	0.0204095200106
26	0.0180494615078
27	0.016609050215
28	0.0175508575987
29	0.0131853033728
30	0.0151021584009
31	0.0216504897398
32	0.0128418206799
33	0.014758675708
34	0.0118556929486
35	0.010902805478
36	0.00699153481363
37	0.00591676638745
38	0.00717989629039
39	0.0055178832602
40	0.00278110180384
41	0.00341266675531
```
In addition to the text file, the script outputs an image, 'lenDistHistogram.png'.  To look at this, bring it to your local computer using filezilla, fugu, or scp.  On the graph there is a large peak at 22nt and a smaller peak around 30nt.  These are the miRNA and tRNA peaks, respectively.  If you see a high amount of other sequence lengths, this can be indicitive of degraded RNA.

#### miRNA expression level

We will have many expressed miRNAs, but likely those with only a couple reads aren't biologically relevant, so we will identify the highly expressed miRNAs.  To calculate the highly expressed miRNAs, enter:

```
$ python genNormalRPM.py mmu /netscr/mattkank/miRquant/03_samples/SampleA./TAB_lenDist_summary.txt

$ ls RPM*
RPM_all.tsv  RPM_miRs_only.tsv  RPM_miRs_over_100.tsv

$ head RPM_miRs_over_1000.tsv
	/netscr/mattkank/miRquant/03_samples/test./TAB_lenDist_summary.txt
mmu-mir-15a-5p	259.100909414
mmu-mir-191-5p	1829.90017274
mmu-let-7a-2-5p	13580.5630124
mmu-mir-10a-5p_+_1	485.814205151
mmu-mir-92a-1-3p	8064.5158055
mmu-mir-101b-3p	2671.97812833
mmu-let-7a-1-5p	14062.3287658
mmu-mir-30e-5p	1833.94862444
mmu-mir-22-3p	22768.4924147
```

These are our highly expressed miRNAs in our sample.  If we had a control, we could compare the two samples to see which miRNAs change between the treatment and the sample.  We could also look at the highly expressed miRNAs and see which mRNAs they target.

This is a beginning primer on smRNA-seq analysis according to the pipeline I use.  There are many other smRNA-seq analysis pipelines out there doing very similar things, and choosing the best bioinformatics tool can be half the challenge when approaching any analyses.  The best approach to making a choice (that I've found) is by reading the literature or googling extensively (SEQanswers and biostars especially).
