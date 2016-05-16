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
00_documentation  02_software  05_scripts        chainSubmission.sh  process_all_summary2tab.pl  runC.sh
01_logs           03_samples   06_final_scripts  post_runC.sh        resources                   uncENV.sh

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
$ wc -l *
4000000 SampleA.fastq
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
$ bsub -o 01_logs/chainSubmission.log bash chainSubmission.sh mmu /netscr/ONYEN/miRquant/03_Samples
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
$ bsub -o 01_logs/process_sum.log perl process_summary2tab.pl /netscr/ONYEN/miRquant/ mmu /netscr/ONYEN/miRquant/03_samples/SampleA./IntermediateResults/g1Results/shift_summary.txt
```
This scripts combines the information into a more human readable format that can be imported into excel.  The files that are produced by this script can be seen by entering `$ ls /netscr/ONYEN/miRquant/03_samples/SampleA./`.  The files whose names begin with 'TAB' are the results files.  Lets look at the top of one, enter:

```
$ head /netscr/ONYEN/miRquant/03_samples/SampleA./TAB_3p_summary_miR.txt
Name	tRNA	miRbaseOffset	Seed	Percent	Count	EM	E	T	A	C	G	AT	AA	TT	CT	TA	AAA	AG	AAT	GT	CCA	TTT	TAT	GCT	ATT	CA	AC	GA	TGT	AGT	CC	TG	TC	AAG	AAAT	TTA	AAAA	GC	ATA	TTAT	ACA	GG	TAA	AGA	TTTT	CG	AAGT	TCT	AAC	ATTT	ACC	TATT	AATT	GAC	ACT	GGA	ATAT	GCCA	TGA	ATG	CAA	ATGT	AAGA	ATC	CGC	CTT	TAAT	ACCA	CTA	TAGT	CCC	GGG	TTC	GAA	CTG	TCA	CGCT	AGCT	AATC	GAT	TGCT	CTAT	AACT	TGAC	CGT	ATTA	GTT	GCTT	GGC	CAAT	AACG	AGTA	GCCC	AGTT	CTGT	GTTT	CAT	AATG	AATA	CTC	TTG	TAG	AGAT	TTGT	AGG	GTA	TTCT	GCC	CAAG	ATCT	TTAG	CTTT	CTCT	CCAA	N	GAAA	GCTC	CACT	GAG	TGTA	CACC	TTGC	CCAT	GCAG	TAGA	TTTA	AGC	TCG	GCGC	GATC	TCTC	ACAC	AGAA	CGCA	ACG	CGTT	GTC	TAAA	GACC	GGT	TTAC	AAAG	CCGG	TCCA	GTAT	GCTA	NT	ATCA	TTAA	ATTC	GGCT	CCTG	GGGT
mmu-mir-122-5p	Mir122a	0	GGAGTGT	0.995552057919515	204080.277777778	180338.044444444	3814.03333333333	5218	12106.3333333333	213.7	227	757	549.333333333333	28	105	90	131	47.3333333333333	50	41	2	4	34	0	1	5	28.5	21.5	32	15	1	27	16	17	18	0	3	7	21	0	2	0	2	11.5	0	1	9	37	0	1	3	2	2	3	0	8	5	7	7	0	5	3	2	0	4	0	0	0	00	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	20	2	0	0	0	2	1	0	0	1	1	0	0	1	0	0	1	1	0	0	0	00	1	0	0	1	0	0	0	0	1	0	0	1	0	0	0	0	0	1	0	1	01	0	0	0	0	1	1	0	1	0	0	0	0	0	0	0	0
mmu-mir-21-5p	Vmp1	0	AGCTTAT	0.996357644610181	24892.8333333333	22479.8333333333	206	91	1684	249	82	3	22	15	1	1	0	0	1	0	0	0	32	0	4	1	1	4	0	0	4	0	0	0	00	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	9	0	0	0	10	0	0	0	0	0	0	0	0	2	0	0	0	0	0	1	0	0	0	0	0	00	0	0	3	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	00	0	0	0	0	0	0
mmu-mir-192-5p	Mir192	0	TGACCTA	0.834372869802318	24480.5	16580.5	2651.5	4079.5	230	126	31	10	84	312	28	178	0	13	11	15	0	0	16	2	2	7	3	11	0	3	0	3	7	2	0	19	0	0	018	0	0	5	1	0	0	0	2	0	0	0	1	1	0	3	0	0	0	0	0	20	0	0	0	0	2	0	3	2	0	0	1	0	0	0	0	0	0	0	0	3	00	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	1	2	0	0	0	10	0	0	0	0	1	0	1	0	0	0	0	1	0	0	0	1	0	0	1	0	00	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	00	0	0	0
mmu-mir-148a-3p	Mir148a	0	CAGTGCA	0.998540097548028	15047.5	14124	227	419.5	191	8	11	6	11	0	3	8	0	00	0	0	0	5	21	0	0	0	1	4	0	0	0	1	0	0	0	0	1	0	00	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	2	0	1	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	00	0	0
mmu-mir-22-3p	Mir22hg	0	AGCTGCC	0.997266610597141	14229	13511	155	83	388	6	6	46	22	1	0	2	0	04	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0
mmu-let-7b-5p	Mirlet7b	0	GAGGTAG	0.997871332488603	13912.6450980392	10874.0617647059	613.583333333333	580.833333333333	699.666666666667	25	49	493.333333333333	130.333333333333	129	2.5	1	12.8333333333333	53	55.5	28.8333333333333	032	0	1.5	41.5	1	5	5.5	0	23.5	0	0	0	2	8	0.999999999999999	0	0	0	3	00	0	1.5	4.5	0	6	1	0.5	9	0	0	4	0	0	0	0	0	1	0	0	0	30	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	2	0	0	0	0	01	0.666666666666666	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	1	0	00	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0.5	0.5	0	0	0
mmu-let-7f-2-5p	Huwe1	0	GAGGTAG	0.99504171826432	12392.1611111111	11955.6944444444	132.55	12	240.75	16.5	8.33333333333333	1.5	14.1666666666667	0	1	0.5	0.333333333333333	2	0	2	0.5	0	1	0	0.5	0.333333333333333	0.5	1	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
mmu-let-7c-2-5p	Mirlet7c-2	0	GAGGTAG	0.995775373609673	12335.349859944	10047.0617647059	1409.23333333333	279.333333333333	417.22619047619	5.5	44.7857142857143	2.14285714285714	38.8333333333333	39.65	2	2	5.5	18.5	0	0.5	0	7.5	3	0.833333333333333	0.5	0	2.5	1.75	0	1.5	0	0	2.5	1.5	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00
mmu-let-7g-5p	Wdr82	0	GAGGTAG	0.994159700569481	11376.6438596491	10782.0105263158	89.3	45	311.5	16	47.3333333333333	46.5	22	2	0	0	0	1	1	7	0	0	0	1	2	0	0	1	0	0	0	0	01	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0	0	0	0
```
This is the bulk of the pipeline, lets bring it all together and check how out sample sequencing analysis went.

## Bringing it all together - Final analyses scripts

These final scripts to get a complete output are in the '06_final_output_scripts', move there by entering `$ cd 06_final_output_scripts`.  We will also need to load a certain version of python, do this by entering `$ module load python/2.7.6`.

#### Generate mapping statistics

Lets start by looking at some quality control statistics from the analyses.  The 'generate_mapping__info_table.pl' will bring together mapping statistics, including total mapping, miR mapping, and tRNA mapping.  Up to this point we have been using the absolute path for the files.  This script requires the relative path.  To run this script, enter:

```
$ perl generate_mapping_info_table.pl ../03_samples/SampleA/IntermediateFiles/*O10_E1.fq

$ cat MappingTableInfo.tsv

Sample name	SampleA
File name	/netscr/mattkank/miRquant/03_samples/SampleA.fastq
Total reads	1000000.00000000000000000000
Trimmed reads	901782.00000000000000000000
% Trimmed reads	90.1782%
Short reads	70483.00000000000000000000
% Short	7.0483%
Exact match to genome	590247
% EM	65.4534022635182%
No exact match to genome	311535
% NEM	34.5465977364818%
Total mapped reads	 696249.747611942
% Mapped	77.2082108105886%
Total mapped to miRs	 474786.71503216
% of total mapped to miRs	68.1920125157136%
Total mapped to tRNAs	 16298.747979798
% of total mapped to tRNAs	2.34093413113626%
```
For most samples we see the following: the % trimmed around 10%, % too short around 5-10%, % mapped around 80-85%, % mapped to miRs around 70-75%, and the % mapped to tRNAs under 10%.  If the statistics vary greatly from these expected numbers, it could be that the samples are of low quality or some other interesting biological process is occuring.

#### Read length distribution

Next we will look at the read length distribution for our sample.  As miRNAs are the major species of small RNA and the length of a mature miRNA is ~22 nucleotides long, we will expect a spike at 22 nucleotides.  To check the length distributions, enter:

```
$ python lenDist.py /netscr/mattkank/miRquant/03_samples/SampleA./IntermediateFiles/*1.fq --image
Loading required package: methods
null device 
          1 
$ cat lenDist.tsv
	SampleA._O10_E1.fq
14	0.0284081962159
15	0.0354597896166
16	0.0226839746191
17	0.0265773767939
18	0.0366186062707
19	0.0349807381385
20	0.0400263034747
21	0.124606612241
22	0.264064929218
23	0.135104714887
24	0.0451927405958
25	0.0201234888255
26	0.018172906534
27	0.0174487847395
28	0.0179378164567
29	0.0132049652799
30	0.0155625195446
31	0.0211547802019
32	0.0132094009417
33	0.0152209735834
34	0.0116369588215
35	0.0102641214839
36	0.00735210948988
37	0.00621990680674
38	0.00705270231608
39	0.005678756063
40	0.00285323947473
41	0.00317815170407
42	4.43566183401e-06
```
In addition to the text file, the script outputs an image, 'lenDistGraph.png'.  To look at this, bring it to your local computer using filezilla, fugu, or scp.  On the graph there is a large peak at 22nt and a smaller peak around 30nt.  These are the miRNA and tRNA peaks, respectively.  If you see a high amount of other sequence lengths, this can be indicitive of degraded RNA.

#### miRNA expression level

We will have many expressed miRNAs, but likely the highly expressed miRNAs are the biologically relevant ones.  To calculate the highly expressed miRNAs, enter:

```
$ python genNormalRPM.py mmu /netscr/mattkank/miRquant/03_samples/SampleA./TAB_lenDist_summary.txt

$ ls RPM*
RPM_all.tsv  RPM_miRs_only.tsv  RPM_miRs_over_100.tsv

$ head RPM_miRs_over_1000.tsv
        /netscr/mattkank/miRquant/03_samples/SampleA./TAB_lenDist_summary.txt
mmu-mir-15a-5p	179.533278725
mmu-mir-191-5p	1612.92697606
mmu-let-7a-2-5p	11947.4368048
mmu-mir-10a-5p_+_1	393.536946964
mmu-mir-92a-1-3p	7371.63642444
mmu-mir-221-3p	1928.90554662
mmu-mir-101b-3p	2553.68135658
mmu-mir-126-3p	211.13113578
mmu-let-7a-1-5p	12303.2717632
```

These are our highly expressed miRNAs in our sample.  If we had a control, we could compare the two samples to see which miRNAs change between the treatment and the sample.  We could also look at the highly expressed miRNAs and see which mRNAs they target.

This is a beginning primer on smRNA-seq analysis according to the pipeline I use.  There are many other smRNA-seq analysis pipelines out there doing very similar things, and choosing the best bioinformatics tool can be half the challenge when approaching any analyses.  The best approach to making a choice (that I've found) is through reading the literature or googling extensively (SEQanswers and biostars especially).

Any questions?
