ChIP-seq 
========================================================
author:MRC LMS Bioinformatics Core
date:https://github.com/LMSBioinformatics/LMS_ChIPseq_short
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css


ChIP-seq introduction 
========================================================

Chromatin precipitation followed by deep sequencing (**ChIP-seq**) is a well established technique which allows for the genome wide identification of transcription factor binding sites and epigenetic marks. 

<div align="center">
<img src="imgs/igsss_1.png" alt="igv" height="600" width="1200">
</div>
ChIP-seq introduction (continued)
========================================================
In this course we will use a few of the packages from the comprehensive repository available from the [Bioconductor project](https://www.bioconductor.org/).

We will cover some of the basics of quality control, working with peaks, motif identification and functional analysis. 

For more details on alignment, working with ChIP-seq coverage and peak calling you can join us on [our extended course](https://github.com/LMSBioinformatics/LMS_chipseqcourse).

========================================================

* Where to find more information.
* [ChIP-seq file types covered](#/filetypes).
* [Story so far](#/background).
* [Materials](#/materials).
* [Assessing ChIP-seq quality](#/qc).
* [Working with peaks](#/peakpushing).
* [Functional annotation of peaks](#/functional).
* [Denovo motifs](#/motifs).
* [Getting hold of external data](#/external).
* [Exporting data for visualisation](#/viz).

Some extra work -
* [Complex Overlaps](#/complexOverlaps).
* [Differential ChIP-seq](#/diffchip).

Reminder of file types
========================================================
id: filetypes

In this session we will be dealing with two data types,

* [BED/BED6 files](http://mrccsc.github.io/genomicFormats.html#/18).

* [FASTA files](http://mrccsc.github.io/genomicFormats.html#/6).

Reminder of file types - BED files
========================================================

* BED or BED 6 files are used to store genomic locations. 
* A mimimum requirement is chromosome,start and end positions for intervals.
* BED6 can additionally store interval name, score and strand.


![alt text](imgs/bed6.png)

Reminder of file types - FASTA files
========================================================

* FASTA files store sequence information alongside names for stored sequences.
* Lines starting with ">" contains name and/or information on sequence.
* Lines following contain contig sequence

![alt text](imgs/fasta.png)

Story so far.
========================================================
id: background

In this course we will use some of the Encode data for Myc ChIP-seq in mouse aligned to the mm9 genome (GEO - GSM912934, GSM912906).

This data is composed of Myc chip for two cell lines, Mel and Ch12 cell lines, each with two replicates.

Due to the short time we have together the data has been processed from unaligned reads to called peaks. 

For full details on the analysis/processing of this data with all analysis steps covered in R/Bioconductor are available on [github](https://github.com/LMSBioinformatics/LMS_chipseqcourse/).

Materials.
========================================================
id: materials

All material for this course can be found on github.
* [ChIPseq_short](https://github.com/LMSBioinformatics/LMS_ChIPseq_short)

Or can be downloaded as a zip archive from here. 
* [Download zip](https://github.com/LMSBioinformatics/LMS_ChIPseq_short/archive/master.zip)

Materials. - Presentations, source code and practicals.
========================================================

Once the zip file in unarchived. All presentations as HTML slides and pages, their R code and HTML practical sheets will be available in the directories underneath.

* **presentations/slides/**
Presentations as an HTML slide show.
* **presentations/singlepage/** 
Presentations as an HTML single page.
* **presentations/rcode/**
R code in presentations.
* **presentations/practicals/**
Practicals as an HTML page. 

Materials. - Data for presentations, practicals.
========================================================

All data to run code in the presentations and in the practicals is available in the zip archive. This includes raw data (MACS peak calls) as well as R objects containing pre-compiled results.

**data/MacsPeaks/**
- MACS peak calls for the 4 Myc ChIP-seq datasets ending in *"_peaks.xls"*
+ 2 replicates of Myc ChIP from Mel cell-line (mycmelrep1,mycmelrep1) 
+ 2 replicates of Myc ChIP from Ch12 cell-line (mycch12rep1,mycch12rep1)

**data/robjects/**
- Robjects for further analysis or review

**data/peaksFromCourse/**
- Contains all GRanges generated in slides exported as a BED file for review in IGV or genome browser of choice.


Set the Working directory
========================================================

Before running any of the code in the practicals or slides we need to set the working directory to the folder we unarchived. 

You may navigate to the unarchived ChIPseq_1Day folder in the Rstudio menu

**Session -> Set Working Directory -> Choose Directory**

or in the console.


```r
setwd("/PathToMyDownload/ChIPseq_1Day/course")
# e.g. setwd("~/Downloads/ChIPseq_1Day/course")
```


Working With ChIP-seq data in R
========================================================
type:section

Quality Control.
========================================================
id: qc

ChIP-seq has many sources of potential noise including 
* Varying efficiency of antibodies
* Non-specific binding
* Library complexity
* ChIP artefacts and background.

Many of these sources of noise can be assessed using some now well-established methodology.

Quality Control. - Some references
========================================================

For some discussions:

* Encode quality metrics.

[Large-scale quality analysis of published ChIP-seq data. Marinov GK, Kundaje A, Park PJ, Wold BJ. G3 (Bethesda). 2014 Feb 19;4(2)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3931556/)

* Overestimation of artefact duplicates in ChIPseq.

[Systematic evaluation of factors influencing ChIP-seq fidelity.Nat Methods. Chen Y, Negre N, Li Q, Mieczkowska JO, Slattery M, Liu T, Zhang Y, Kim TK, He HH, Zieba J, Ruan Y, Bickel PJ, Myers RM, Wold BJ, White KP, Lieb JD, Liu XS. 2012 Jun;9(6)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3477507/)


* When and what QC is useful.

[Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data.Front Genet. 2014 Apr 10;5:75.Carroll TS, Liang Z, Salama R, Stark R, de Santiago I](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3989762/)

Quality Control - Always have an appropriate input.
========================================================

* Input samples are typically made from fragmented DNA prior to IP enrichment.

* Allows for control of artefact regions which occur across samples.

* NEVER run ChIP-seq without considering which input to use.

e.g. When using tumour samples for ChIP-seq, it is important to have matched input samples. 
Differing conditions of same tissue may share common input. 

Quality Control - Quality metrics for ChIP-seq.
========================================================

The ChIPQC package wraps some of the metrics into a Bioconductor package and takes care to measure these metrics under the appropriate condition. 

To run a single sample through ChIPQCsample function, we must provide the a set of peaks, the relevant unfiltered BAM file and we are recommended to supply a **blacklist** as a BED file or GRanges and Genome name.

You can find a Blacklist for most genomes at [Anshul Kundaje's site](https://sites.google.com/site/anshulkundaje/projects/blacklists) 


```r
QCresult <- ChIPQCsample(reads="/pathTo/myChIPreads.bam",
                         peaks="/pathTo/myChIPpeaks.bed",
                         genome="mm9",
                         blacklist = "/pathTo/mm9_Blacklist.bed")
```

Quality Control 
========================================================

Although we may not have time in this practical, we can look at the full course to see how to create ChIPQCexperiment objects containing multiple samples' QC metrics.

Here we can import the ChIPQCexperiment object from the course and take a look at some of the outputs. 

The first useful function is QCmetrics which will provide a table of QC scores.


```r
library(ChIPQC)
load("data/robjects/ChIPQCwithPeaks.RData")
QCmetrics(res)
```



|           |    Reads| Map%|   Filt%|  Dup%| ReadL| FragL| RelCC|  SSD| RiP%| RiBL%|
|:----------|--------:|----:|-------:|-----:|-----:|-----:|-----:|----:|----:|-----:|
|myc_ch12_1 | 10792905|  100| 0.0e+00| 10.20|    36|   176| 1.030| 5.22| 14.0|  13.9|
|myc_ch12_2 |  9880785|  100| 0.0e+00|  9.98|    36|   146| 1.370| 3.94| 19.5|  11.1|
|myc_Mel_1  |  9912979|  100| 0.0e+00|  9.70|    35|   169| 1.150| 4.57| 23.1|  13.0|
|myc_Mel_2  | 10475318|  100| 0.0e+00|  9.71|    35|   161| 0.973| 5.54| 21.7|  15.3|
|ch12       | 15907271|  100| 6.3e-06|  6.85|    36|   180| 0.744| 7.76|   NA|  16.2|
|MEL        | 18437914|  100| 0.0e+00|  5.65|    35|   173| 0.444| 8.61|   NA|  17.1|



Quality Control (Fraction of reads in peaks - FRIP/RIP)
========================================================

RIP/FRIP/PTIH/SPOT all are measurements of the fraction/percentage of reads landing in peaks. Variability in the proportion of reads in peaks for ChIP-seq samples can identify poorer quality replicates.


```r
frip(res)
```

```
myc_ch12_1 myc_ch12_2  myc_Mel_1  myc_Mel_2       ch12        MEL 
 0.1400971  0.1949768  0.2309460  0.2172935         NA         NA 
```

```r
plotFrip(res)
```

![plot of chunk unnamed-chunk-6](ChIPseq-figure/unnamed-chunk-6-1.png)


Quality Control (Assessing fragment length)
========================================================

The prediction of fragment length is an essential part of ChIP-seq affecting peaks calling, summit identification and coverage profiles. 

The use of cross-correlation or cross-coverage allows for an assessment of reads clustering by strand and so a measure of quality. 

Quality Control (Assessing fragment length)
========================================================
<div align="center">
<img src="imgs/ChIP-seq_biology_slides.png" alt="offset" height="900" width="1400">
</div>

Quality Control (Assessing fragment length)
========================================================
* In ChIP-seq typically short single end reads of dsDNA.

* **dsDNA single end sequencing means**
+ 5' will be sequenced on "+" strand
+ 3' end will be on "-" strand.

* **Although we only have partial sequence of strand, with predicted fragment length.**
+ "+" reads should extend only in positive direction 
+ "-" reads only in negative

Quality Control (Assessing fragment length)
========================================================
<div align="center">
<img src="imgs/pileup.png" alt="offset" height="900" width="900">
</div>

Quality Control (Assessing fragment length)
========================================================
<div align="center">
<img src="imgs/offset.jpg" alt="offset" height="900" width="900">
</div>


Quality Control (Assessing fragment length)
========================================================
<div align="center">
<img src="imgs/shifts.gif" alt="offset" height="300" width="1400">
</div>
<div align="center">
<img src="imgs/cor.gif" alt="offset" height="600" width="1400">
</div>


Quality Control (Assessing fragment length)
========================================================
<div align="center">
<img src="imgs/shifts.gif" alt="offset" height="900" width="1300">
</div>


Quality Control (Assessing fragment length)
========================================================
ChIPQC has already calculated change in coverage between successive shifts of the "+" strand reads towards "-" strand reads and we can plot these for inspection with the **plotCC** function. 


```r
ccplot <- plotCC(res)
ccplot$layers <- ccplot$layers[1]
ccplot
```

![plot of chunk unnamed-chunk-7](ChIPseq-figure/unnamed-chunk-7-1.png)


Quality Control - Blacklists and SSD.
========================================================

ChIP-seq will often show the presence of common artefacts such as ultra-high signal regions or **Blacklists**. Such regions can confound peak calling, fragment length estimation and QC metrics.

SSD is a measure of standard deviation of signal across the genome with higher scores reflecting significant pile-up of reads. SSD can therefore be used to assess both the extent of ultra high signals and the signal following removal of these blacklisted regions.

For a note on known Blacklisted regions and on associated resources.
* [Blacklisted Regions](http://mrccsc.github.io/analysisbeginings.html#/35)

For a note on SSD
* [SSD and Signal Pileup](http://mrccsc.github.io/analysisbeginings.html#/36).

Quality Control - Blacklist???
========================================================

<div align="center">
<img src="imgs/blacklist.png" alt="offset" height="900" width="1000">
</div>


Quality Control - Blacklist affects many metrics.
========================================================

<div align="center">
<img src="imgs/blacklistsAffects.jpg" alt="offset" height="900" width="1300">
</div>


Quality Control - Blacklist affects many metrics.
========================================================

<div align="center">
<img src="imgs/ssdAndBlacklist.png" alt="offset" height="900" width="600">
</div>



Quality Control - Standardised Standard Deviation.
========================================================

ChIPQC calculates SSD before and after removing signal coming from Blacklisted regions.

The plotSSD function plots samples's pre-blacklisting score in **red** and post-blacklisting score in **blue**.

Higher scores for pre-blacklisted SSD can suggest a strong background signal in blacklisted regions for that sample.


```r
plotSSD(res)+xlim(0,14)
```

![plot of chunk unnamed-chunk-8](ChIPseq-figure/unnamed-chunk-8-1.png)

Quality Control - Standardised Standard Deviation.
========================================================

Since SSD score is strongly affected by blacklisting it may be necessary to change the axis to see any differences between samples for post-blacklisting scores.

Higher post-blacklisted SSD scores reflect samples with stronger peak signal.


```r
plotSSD(res)+xlim(0.2,0.4)
```

![plot of chunk unnamed-chunk-9](ChIPseq-figure/unnamed-chunk-9-1.png)

Quality Control 
========================================================

For more details on assessing ChIP-seq quality you can visit the Bioconductor workshop we ran in Boston (Bioc 2014).

* Practical - https://www.bioconductor.org/help/course-materials/2014/BioC2014/Bioc2014_ChIPQC_Practical.pdf
* Practical Data - https://www.bioconductor.org/help/course-materials/2014/BioC2014/BCell_Examples.RData
* Theory part - https://www.bioconductor.org/help/course-materials/2014/BioC2014/ChIPQC_Presentation.pdf

Working with Peaks
========================================================
id: peakpushing

Macs2 is a frequently used peak caller and works well to identify both punctate and broad peaks.

For more details on peak calling steps for data in this course you can visit our [material](https://github.com/LMSBioinformatics/LMS_chipseqcourse/blob/master/Precticals_pres/peakCallingForMyc).

For more details on MACS2, see the github page for MACS2 software
* [MACS2 github page](https://github.com/taoliu/MACS).

Working with Peaks
========================================================
MACS peak calls can be found in the user specied output directory with the suffix and extension "_peaks.xls".

MACS peaks come as a tab seperated file thinly disguised as a ".xls".


|chr |   start|     end| length| abs_summit| pileup| X.log10.pvalue.| fold_enrichment| X.log10.qvalue.|name               |
|:---|-------:|-------:|------:|----------:|------:|---------------:|---------------:|---------------:|:------------------|
|1   | 4480665| 4480849|    185|    4480769|     12|         6.89435|         4.45956|         4.21774|mycch12rep1_peak_1 |
|1   | 4661593| 4661934|    342|    4661762|     36|        30.49564|        11.09288|        26.75324|mycch12rep1_peak_2 |
|1   | 4774202| 4774393|    192|    4774293|     10|         5.83769|         4.13574|         3.27058|mycch12rep1_peak_3 |

Working with Peaks
========================================================

In addition to the genomic coordinates of peaks, these files contain useful information on the samples, parameters and version used for peak calling at the top.


```
[1] "# Command line: callpeak -t wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam -c wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam -f BAM -n mycch12rep1 --nomodel --extsize 165"
[2] "# ARGUMENTS LIST:"                                                                                                                                                                        
[3] "# name = mycch12rep1"                                                                                                                                                                     
[4] "# format = BAM"                                                                                                                                                                           
[5] "# ChIP-seq file = ['wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam']"                                                                                                            
[6] "# control file = ['wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam']"                                                                                                                
```


Working with Peaks - Importing peaks
========================================================

We can import peak files therefore using read.delim function. Note we have set *comment.char* argument to **#** to exclude additional information on peak calling parameters stored within the MACS peak file.


```r
peakfile <- "data/MacsPeaks/mycch12rep1_peaks.xls"
macsPeaks_DF <- read.delim(peakfile,comment.char="#")
macsPeaks_DF[1:2,]
```

```
  chr   start     end length abs_summit pileup X.log10.pvalue.
1   1 4480665 4480849    185    4480769     12         6.89435
2   1 4661593 4661934    342    4661762     36        30.49564
  fold_enrichment X.log10.qvalue.               name
1         4.45956         4.21774 mycch12rep1_peak_1
2        11.09288        26.75324 mycch12rep1_peak_2
```

Working with Peaks - Importing peaks
========================================================

Now we have the information in a table we can create a GRanges object.

GRanges objects are made of chromosome names and intervals stored as IRanges.


```r
library(GenomicRanges)
macsPeaks_GR <- GRanges(
 seqnames=macsPeaks_DF[,"chr"],
 IRanges(macsPeaks_DF[,"start"],
         macsPeaks_DF[,"end"]
 )
)
macsPeaks_GR
```

```
GRanges object with 33498 ranges and 0 metadata columns:
          seqnames          ranges strand
             <Rle>       <IRanges>  <Rle>
      [1]        1 4480665-4480849      *
      [2]        1 4661593-4661934      *
      [3]        1 4774202-4774393      *
      [4]        1 4775399-4775792      *
      [5]        1 4775957-4776125      *
      ...      ...             ...    ...
  [33494]        Y   234019-234250      *
  [33495]        Y   307766-308084      *
  [33496]        Y   582005-582258      *
  [33497]        Y   622964-623320      *
  [33498]        Y 2721204-2721372      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Working with Peaks - Peaks as GRanges
========================================================

As we have seen before elements in GRanges can accessed and set using various GRanges functions.
Here we can deconstruct our object back to contig names and interval ranges.


```r
seqnames(macsPeaks_GR)
```

```
factor-Rle of length 33498 with 22 runs
  Lengths: 2257 1766 3226 1299 1509 1171 ... 1906 1571 1860   11  614    5
  Values :    1   10   11   12   13   14 ...    7    8    9   MT    X    Y
Levels(22): 1 10 11 12 13 14 15 16 17 18 19 2 3 4 5 6 7 8 9 MT X Y
```

```r
ranges(macsPeaks_GR)
```

```
IRanges object with 33498 ranges and 0 metadata columns:
              start       end     width
          <integer> <integer> <integer>
      [1]   4480665   4480849       185
      [2]   4661593   4661934       342
      [3]   4774202   4774393       192
      [4]   4775399   4775792       394
      [5]   4775957   4776125       169
      ...       ...       ...       ...
  [33494]    234019    234250       232
  [33495]    307766    308084       319
  [33496]    582005    582258       254
  [33497]    622964    623320       357
  [33498]   2721204   2721372       169
```

Working with Peaks - Peaks as GRanges
========================================================

GRanges objects may have metadata attached. Here we attach some useful information on our peaks including the summit position and the fold enrichment over input.


```r
mcols(macsPeaks_GR) <- macsPeaks_DF[,c("abs_summit", "fold_enrichment")]
macsPeaks_GR
```

```
GRanges object with 33498 ranges and 2 metadata columns:
          seqnames          ranges strand | abs_summit fold_enrichment
             <Rle>       <IRanges>  <Rle> |  <integer>       <numeric>
      [1]        1 4480665-4480849      * |    4480769         4.45956
      [2]        1 4661593-4661934      * |    4661762        11.09288
      [3]        1 4774202-4774393      * |    4774293         4.13574
      [4]        1 4775399-4775792      * |    4775544         7.31683
      [5]        1 4775957-4776125      * |    4776021         3.17407
      ...      ...             ...    ... .        ...             ...
  [33494]        Y   234019-234250      * |     234139         5.10991
  [33495]        Y   307766-308084      * |     307929        18.80093
  [33496]        Y   582005-582258      * |     582128         5.08412
  [33497]        Y   622964-623320      * |     623149         7.89867
  [33498]        Y 2721204-2721372      * |    2721341         6.14918
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```


Working with Peaks - GRanges manipulation.
========================================================

GRanges objects can be rows subset as vectors or data.frames/matrices using an index.

In this example, the first 3 genomic intervals are retrieved.


```r
macsPeaks_GR[1:3]
```

```
GRanges object with 3 ranges and 2 metadata columns:
      seqnames          ranges strand | abs_summit fold_enrichment
         <Rle>       <IRanges>  <Rle> |  <integer>       <numeric>
  [1]        1 4480665-4480849      * |    4480769         4.45956
  [2]        1 4661593-4661934      * |    4661762        11.09288
  [3]        1 4774202-4774393      * |    4774293         4.13574
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

```r
# or macsPeaks_GR[1:3,]
```

Working with Peaks - GRanges manipulation.
========================================================

GRanges objects can also be rows subset as vectors or data.frames/matrices using logical vectors.

Here we can retrieve all peaks on chromosome 1.


```r
macsPeaks_GR[seqnames(macsPeaks_GR) %in% "1"]
```

```
GRanges object with 2257 ranges and 2 metadata columns:
         seqnames              ranges strand | abs_summit fold_enrichment
            <Rle>           <IRanges>  <Rle> |  <integer>       <numeric>
     [1]        1     4480665-4480849      * |    4480769         4.45956
     [2]        1     4661593-4661934      * |    4661762        11.09288
     [3]        1     4774202-4774393      * |    4774293         4.13574
     [4]        1     4775399-4775792      * |    4775544         7.31683
     [5]        1     4775957-4776125      * |    4776021         3.17407
     ...      ...                 ...    ... .        ...             ...
  [2253]        1 196603937-196604107      * |  196604025          4.6006
  [2254]        1 196643096-196643377      * |  196643232        14.81234
  [2255]        1 196706013-196706236      * |  196706154          5.0807
  [2256]        1 196956950-196957403      * |  196957192         5.11378
  [2257]        1 196957597-196957944      * |  196957689          3.9841
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Working with Peaks - GRanges naming and indexing
========================================================

GRanges objects may have names attached to each genomic interval as seen with vector elements.

This name may also be used to retrieve an interval of interest.


```r
names(macsPeaks_GR) <- macsPeaks_DF[,"name"]
macsPeaks_GR["mycch12rep1_peak_33496"]
```

```
GRanges object with 1 range and 2 metadata columns:
                         seqnames        ranges strand | abs_summit
                            <Rle>     <IRanges>  <Rle> |  <integer>
  mycch12rep1_peak_33496        Y 582005-582258      * |     582128
                         fold_enrichment
                               <numeric>
  mycch12rep1_peak_33496         5.08412
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```






Working with peaks - Read in peaksets in one step.
========================================================

There are many helper functions for reading in peak sets.

The **rtracklayer** package has many tools for importing common file formats with the *import.bed* function being useful for peak sets in bed formats.

Since MACS has a slightly different format to BED or BED6 we can use a **ChIPQC** function *GetGRanges* to read here.


```r
peakfile <- "data/MacsPeaks/mycch12rep1_peaks.xls"
singlePeakSet <- ChIPQC:::GetGRanges(peakfile, sep="\t", simple=F)
singlePeakSet[1:2,]
```

```
GRanges object with 2 ranges and 7 metadata columns:
      seqnames          ranges strand |        ID     Score   Strand
         <Rle>       <IRanges>  <Rle> | <integer> <integer> <factor>
  [1]        1 4480665-4480849      * |       185   4480769        *
  [2]        1 4661593-4661934      * |       342   4661762        *
      X.log10.pvalue. fold_enrichment X.log10.qvalue.               name
            <numeric>       <numeric>       <numeric>           <factor>
  [1]         6.89435         4.45956         4.21774 mycch12rep1_peak_1
  [2]        30.49564        11.09288        26.75324 mycch12rep1_peak_2
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```


Manipulating Peak Sets - Finding Common peaks
=========================================================

Now we have our data in peak format we can start to do some simple but powerful analysis.

First lets read in the two replicate Myc ChIP replicates for Ch12 cell.


```r
firstPeakSet <- ChIPQC:::GetGRanges("data/MacsPeaks//mycch12rep1_peaks.xls", sep="\t", simple=F)
secondPeakSet <- ChIPQC:::GetGRanges("data/MacsPeaks//mycch12rep2_peaks.xls", sep="\t", simple=F)
```
Manipulating Peak Sets - Finding Common peaks
=========================================================

We have learnt how to identify overlapping intervals in the Bioconductor tutorial using the [%over%](?%over%) command.

Here we can apply this to identifying peaks in both replicates.


```r
OnlyfirstPeakSet <- firstPeakSet[!firstPeakSet %over% secondPeakSet]
firstANDsecondPeakSets <- firstPeakSet[firstPeakSet %over% secondPeakSet]
length(OnlyfirstPeakSet)
```

```
[1] 4697
```

```r
length(firstANDsecondPeakSets)
```

```
[1] 28801
```

Manipulating Peak Sets - Accessing data in GRanges metadata columns
=========================================================

Data in GRanges metadata columns may be accessed just as data.frame columns.

To access the fold_enrichment or abs_summit columns we can use the following syntax


```r
foldEnrichment <- firstPeakSet$fold_enrichment
# or foldEnrichment <- firstPeakSet[,"fold_enrichment"]
foldEnrichment[1:10]
```

```
 [1]  4.45956 11.09288  4.13574  7.31683  3.17407  3.73551  5.86629
 [8]  4.86149  6.01562  5.39653
```

Manipulating Peak Sets - Finding Common peaks
=========================================================

Now can plot the distribution of peaks' signal over input for those common to both replicates and those unique to replicate 1.

Here it is clear that the peaks with the highest fold enrichment are common to both replicates.


```r
FirstOnly_FE <- log2(OnlyfirstPeakSet$fold_enrichment)
FirstAndSecond_FE <- log2(firstANDsecondPeakSets$fold_enrichment)

boxplot(FirstOnly_FE,
        FirstAndSecond_FE,
        names=c("Only_in_First","Common_to_first_second"),
        ylab="log2 Fold_Enrichment")
```

![plot of chunk boxplotOfFE](ChIPseq-figure/boxplotOfFE-1.png)

Manipulating Peak Sets - Finding Common peaks
=========================================================

When looking at peaks which occur in both samples it is clear that the number of peaks in first replicate overlapping those in second is different from number of second replicate peaks overlapping first.

This is because 2 peaks from one replicate may overlap 1 peak in the other replicate.


```r
firstANDsecondPeakSets <- firstPeakSet[firstPeakSet %over% secondPeakSet]
secondANDfirstPeakSets <- secondPeakSet[secondPeakSet %over% firstPeakSet]

length(firstANDsecondPeakSets)
```

```
[1] 28801
```

```r
length(secondANDfirstPeakSets)
```

```
[1] 27858
```

=========================================================
![alt text](imgs/oneToMany.png)

Manipulating Peak Sets - Finding Common peaks
=========================================================

A common step with finding overlapping transcription factor peaks is to reduce peaksets to single non-overlapping peakset before interrogating whether a peak occurred in a sample.

This allows for a single peak set to be used as a consensus peakset between replicates.



```r
allPeaks <- c(firstPeakSet,secondPeakSet)
allPeaksReduced <- reduce(allPeaks)
length(allPeaks)
```

```
[1] 84856
```

```r
length(allPeaksReduced)
```

```
[1] 55427
```

=========================================================
![alt text](imgs/mel_Flattened.png)


=========================================================
id: makingcommonpeaks

Now we can use a logical expression to subset our reduced/flattened peak set to those overlapping peaks in both the first and second replicate.


```r
commonPeaks <- allPeaksReduced[allPeaksReduced %over% firstPeakSet 
                               & allPeaksReduced %over% secondPeakSet]
length(commonPeaks)
```

```
[1] 27232
```

=========================================================
![alt text](imgs/Ch12_highcon.png)

Time for a exercise.
=========================================================

Exercise on "Working with peaks" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/WorkingWithPeaks_Exercises.html).

Time for a solution.
=========================================================

Answers on "Working with peaks" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/WorkingWithPeaks_Solutions.html).

Rcode for "Working with peaks" solutions can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/WorkingWithPeaks.R).


<!-- Time to install a package. -->
<!-- ========================================================= -->

<!-- In this next section we will be using a package ,**rGREAT**, which has not been pre-installed on these machines. -->

<!-- We saw how to install packages in R and Bioconductor sessions. -->

<!-- You can go to the [rGREAT Bioconductor page](https://www.bioconductor.org/packages/release/bioc/html/rGREAT.html) and follow their instrutions or simply paste the below code into the console. -->

<!-- ```{r,eval=F} -->
<!-- source("https://bioconductor.org/biocLite.R") -->
<!-- biocLite("rGREAT") -->
<!-- ``` -->


Functional Annotation of peaks
=========================================================
id: functional

So far we have been working with ChIP-seq peaks corresponding to transcription factor binding. Transcription factors, as implied in the name, can affect the expression of their target genes.

The target of transcription factor is hard to assertain from ChIP-seq data alone and so often we will annotate peaks to genes by a simple set of rules:-

Peaks are typically annotated to a gene if
* They overlap the gene.
* The gene is the closest (and within a minimum distance.)




Peak annotation
=========================================================

A useful package for annotation of peaks to genes is **ChIPseeker**. 

By using pre-defined annotation in the from of a **TXDB** object for mouse (mm9 genome), ChIPseeker will provide us with an overview of where peaks land in the gene and distance to TSS sites.

First load the libraries we require for the next part.



```r
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)
```

Peak annotation
=========================================================

We use **GenomeInfoDb** package to fix chromosome naming style. We know we are using UCSC annotation (TxDb.Mmusculus.UCSC.mm9.knownGene) so we can style our previously defined [commonPeaks](#/makingcommonpeaks) to "UCSC" standard using the *seqlevelsStyle* function.


```r
commonPeaks[1:2, ]
```

```
GRanges object with 2 ranges and 0 metadata columns:
      seqnames          ranges strand
         <Rle>       <IRanges>  <Rle>
  [1]        1 4661367-4661934      *
  [2]        1 4775386-4776125      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

```r
seqlevelsStyle(commonPeaks) <- "UCSC"
commonPeaks[1:2, ]
```

```
GRanges object with 2 ranges and 0 metadata columns:
      seqnames          ranges strand
         <Rle>       <IRanges>  <Rle>
  [1]     chr1 4661367-4661934      *
  [2]     chr1 4775386-4776125      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Peak annotation
=========================================================

The annotatePeak function accepts a GRanges object of the regions to annotate, a TXDB object for gene locations and a database object name to retrieve gene names from.



```r
peakAnno <- annotatePeak(commonPeaks, tssRegion = c(-1000, 1000), TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene, 
    annoDb = "org.Mm.eg.db")
```

```
>> preparing features information...		 2019-09-24 12:41:12 
>> identifying nearest features...		 2019-09-24 12:41:12 
>> calculating distance from peak to TSS...	 2019-09-24 12:41:13 
>> assigning genomic annotation...		 2019-09-24 12:41:13 
>> adding gene annotation...			 2019-09-24 12:41:25 
>> assigning chromosome lengths			 2019-09-24 12:41:25 
>> done...					 2019-09-24 12:41:25 
```

Peak annotation
=========================================================

The result is a csAnno object containing annotation for peaks and overall annotation statistics.


```r
class(peakAnno)
```

```
[1] "csAnno"
attr(,"package")
[1] "ChIPseeker"
```

```r
peakAnno
```

```
Annotated peaks generated by ChIPseeker
27227/27232  peaks were annotated
Genomic Annotation Summary:
             Feature  Frequency
9           Promoter 29.4303449
4             5' UTR  0.5178683
3             3' UTR  1.7152092
1           1st Exon  0.5252139
7         Other Exon  3.1476108
2         1st Intron 12.5573879
8       Other Intron 19.7744886
6 Downstream (<=300)  1.5242223
5  Distal Intergenic 30.8076542
```


Peak annotation
=========================================================

The csAnno object contains the information on annotation of individual peaks to genes.

To extract this from the csAnno object the ChIPseeker functions *as.GRanges* or *as.data.frame* can be used to produce the respective object with peaks and their associated genes.


```r
peakAnno_GR <- as.GRanges(peakAnno)
peakAnno_GR[1:3, ]
```

```
GRanges object with 3 ranges and 12 metadata columns:
      seqnames          ranges strand |        annotation   geneChr
         <Rle>       <IRanges>  <Rle> |       <character> <integer>
  [1]     chr1 4661367-4661934      * | Distal Intergenic         1
  [2]     chr1 4775386-4776125      * |          Promoter         1
  [3]     chr1 4847417-4848050      * |          Promoter         1
      geneStart   geneEnd geneLength geneStrand      geneId transcriptId
      <integer> <integer>  <integer>  <integer> <character>  <character>
  [1]   4763279   4775807      12529          2       27395   uc007afd.2
  [2]   4763279   4775807      12529          2       27395   uc007afd.2
  [3]   4847775   4887990      40216          1       21399   uc007afi.2
      distanceToTSS            ENSEMBL      SYMBOL
          <numeric>        <character> <character>
  [1]        113873 ENSMUSG00000033845      Mrpl15
  [2]             0 ENSMUSG00000033845      Mrpl15
  [3]             0 ENSMUSG00000033813       Tcea1
                                       GENENAME
                                    <character>
  [1]       mitochondrial ribosomal protein L15
  [2]       mitochondrial ribosomal protein L15
  [3] transcription elongation factor A (SII) 1
  -------
  seqinfo: 22 sequences (1 circular) from mm9 genome
```



Visualising peak annotation
=========================================================
Now we have the annotated peaks from ChIPseeker we can use some of ChIPseeker's plotting functions to display distribution of peaks in gene features. Here we use the **plotAnnoBar** function to plot this as a bar chart but  **plotAnnoPie** would produce a similar plot as a pie chart.



```r
plotAnnoBar(peakAnno)
```

![plot of chunk unnamed-chunk-30](ChIPseq-figure/unnamed-chunk-30-1.png)

Visualising peak annotation
=========================================================
Similarly we can plot the distribution of peaks around TSS sites.



```r
plotDistToTSS(peakAnno)
```

![plot of chunk unnamed-chunk-31](ChIPseq-figure/unnamed-chunk-31-1.png)

Visualising peak annotation
=========================================================
ChIPseeker can also offer a succinct plot to describe the overlap between annotations.



```r
upsetplot(peakAnno, vennpie = F)
```

![plot of chunk unnamed-chunk-32](ChIPseq-figure/unnamed-chunk-32-1.png)


Gene Ontology and geneset testing.
=========================================================
id:gsa

Transcription factors or epigenetic marks may act on specific sets of genes grouped by a common biological feature (shared Biological function, common regulation in RNAseq experiment etc).

A frequent step in ChIP-seq analysis is to test whether common gene sets are enriched for transcription factor binding or epigenetic marks.

Sources of well curated genesets include [GO consortium](http://geneontology.org/) (gene's function, biological process and cellular localisation), [REACTOME](http://www.reactome.org/) (Biological Pathways) and [MsigDB](http://software.broadinstitute.org/gsea/msigdb/) (Computationally and Experimentally derived).

Geneset enrichment testing may be performed on the sets of genes with peaks associated to them. In this example we will consider genes with peaks within 1000bp of a gene's TSS. 


Gene ontology and geneset testing.
=========================================================

To perform geneset testing here, we will use the GOseq package and so must provide a named numeric vector of 1s or 0s to illustrate whether a gene had any peaks overlapping it's TSS.

In this example we use all TSS sites we found to be overlapped by our common peaks we [previously defined](#/makingcommonpeaks). So first lets find all peaks overlapping TSS regions.

The peaks landing in TSS regions will be marked as "Promoter" in the **annotation** column of our annotated GRanges object. We can extract the unique names of genes with peaks in their TSS by subsetting the annotated GRanges and retrieving gene names from the **geneId** column.


```r
genesWithPeakInTSS <- unique(peakAnno_GR[peakAnno_GR$annotation == "Promoter", 
    ]$geneId)

genesWithPeakInTSS[1:2]
```

```
[1] "27395" "21399"
```

Gene ontology and functional testing.
=========================================================

Next we can extract all genes which are included in the TxDb object to use as our universe of genes for pathway enrichment.



```r
allGenes <- unique(unlist(keys(TxDb.Mmusculus.UCSC.mm9.knownGene, "GENEID")))

length(allGenes)
```

```
[1] 21761
```

Gene ontology and functional testing.
=========================================================

Once we have a vector of all genes we can create a named vector of 1s or 0s representing whether a gene had peak in TSS or not.


```r
allGenesForGOseq <- as.integer(allGenes %in% genesWithPeakInTSS)
names(allGenesForGOseq) <- allGenes
allGenesForGOseq[1:3]
```

```
100009600 100009609 100009614 
        0         0         0 
```



Gene ontology and functional testing.
=========================================================

A little helper function you may recognise to add some useful KEGG names alongside the KEGG pathway IDs in GOseq results.


```r
library(KEGG.db)
library(goseq)
xx <- as.list(KEGGPATHID2NAME)
temp <- cbind(names(xx),unlist(xx))
addKeggTogoseq <- function(JX,temp){
  for(l in 1:nrow(JX)){
    if(JX[l,1] %in% temp[,1]){
      JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
      JX[l,"ontology"] <- "KEGG"
    }
    
  }
  return(JX)
}
```

Gene ontology and functional testing.
=========================================================

Now we have the the input for GOseq we can test against KEGG (or GO if we choose) using a standard hypergeometric test.


```r
library(goseq)

pwf = nullp(allGenesForGOseq, "mm9", "knownGene", plot.fit = FALSE)

Kegg_MycPeaks <- goseq(pwf, "mm9", "knownGene", test.cats = c("KEGG"), method = "Hypergeometric")

Kegg_MycPeaks <- addKeggTogoseq(Kegg_MycPeaks, temp)

Kegg_MycPeaks[1:4, ]
```

```
    category over_represented_pvalue under_represented_pvalue numDEInCat
90     03013            1.545710e-21                        1        111
96     03040            3.056583e-19                        1         91
89     03010            3.284752e-19                        1         70
112    04110            3.758630e-15                        1         84
    numInCat          term ontology
90       157 RNA transport     KEGG
96       125   Spliceosome     KEGG
89        87      Ribosome     KEGG
112      123    Cell cycle     KEGG
```


<!-- Gene ontology and functional testing. GREAT method. -->
<!-- ========================================================= -->

<!-- In addition to a standard enrichment tests, methods have been implemented specifically for ChIP-seq. Many of these tools aim to incorporate peaks distal to genes into their enrichment testing such as the popular [GREAT](http://bejerano.stanford.edu/great/public/html/splash.php) toolset. -->

<!-- Incorporating distal peaks by rules such as nearest gene results in some genes having a higher chance of being selected and hence some genesets as a whole having a higher chance of having its members selected. -->

<!-- [GREAT](http://bejerano.stanford.edu/great/public/html/splash.php) defines regulatory regions for each, individual gene and compares the proportion of peaks mapping to a geneset's regulatory regions to the proportion of the genome occupied by geneset's regulatory regions. -->

<!-- i.e. If a gene set's regulatory regions account for 1% of the genome then one might expect 1% of peaks to overlap these regions by chance. -->

<!-- rGREAT - R interface to GREAT server. -->
<!-- ========================================================= -->

<!-- We can use the GREAT Bioconductor interface available in the rGREAT package. Since GREAT uses UCSC annotation lets make sure our peaks our in UCSC style. -->

<!-- ```{r,eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T} -->
<!-- library(rGREAT) -->
<!-- seqlevelsStyle(commonPeaks) <- "UCSC" -->
<!-- ``` -->

<!-- Gene ontology and functional testing. GREAT method. -->
<!-- ========================================================= -->

<!-- To submit jobs we can use our GRanges of commonPeaks and specify a genome with the **submitGreatJob** function. -->

<!-- This function returns a GreatJob object containing a reference to our results on the GREAT server. To review the categories of results available we can use the availableCategories function on our GreatJob object. -->

<!-- ```{r,eval=T,echo=T, eval=T, echo=T,messages=F,message=F,warning=FALSE,tidy=T} -->
<!-- great_Job <- submitGreatJob(commonPeaks,species="mm9") -->
<!-- availableCategories(great_Job) -->
<!-- ``` -->

<!-- Gene ontology and functional testing. GREAT method. -->
<!-- ========================================================= -->

<!-- The results table can be retrieved using the getEnrichmentTables function and specifying which tables we wish to review. -->

<!-- Here we retrieve the results tables for the "MSigDB Predicted Promoter Motifs" genesets which contains 2 seperate database results. -->

<!-- ```{r,eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T} -->
<!-- great_ResultTable = getEnrichmentTables(great_Job,category= -->
<!--                           "Regulatory Motifs") -->
<!-- names(great_ResultTable) -->
<!-- great_ResultTable[["MSigDB Predicted Promoter Motifs"]][1:4,] -->

<!-- ``` -->

Time for an exercise
=========================================================

Exercise on "Functional Annotation of peaks" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/FunctionalAnnotationOfPeaks_Exercises.html).

Time for a solution.
=========================================================

Answers on "Functional Annotation of peaks" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/FunctionalAnnotationOfPeaks_Solutions.html).

Rcode for "Functional Annotation of peaks" solutions can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/FunctionalAnnotationOfPeaks.R).

Identifying Motifs
==========================================================
id: motifs

A common practice in transcription factor ChIP-seq is to investigate the motifs enriched under peaks. 

Denovo motif enrichment can be performed in R/Bioconductor but this can be very time consuming. Here we will use the Meme-ChIP suite available online to identify denovo motifs.

Meme-ChIP requires a FASTA file of sequences under peaks as input so we extract this using the **BSgenome** package.

Extracting sequences under regions
============================

First we need to load the BSgenome object for the genome we are working on, UCSC's mm9 build for the mouse genome. Again we ensure that the chromosome names and style for our common peaks [defined earlier](#/makingcommonpeaks) matches that seen in UCSC.


```r
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
genome <- BSgenome.Mmusculus.UCSC.mm9
seqlevelsStyle(commonPeaks) <- "UCSC"
```

Extracting sequences under regions
============================

The motif for the ChIP-ed transcription factor should in the centre of a peak. Meme-ChIP will trim our peaks to a common length internally if sequences are of different length.

It is best therefore to provide a peak set resized to a common length.


```r
commonPeaks <- resize(commonPeaks,300,fix="center")
commonPeaks[1:4,]
```

```
GRanges object with 4 ranges and 0 metadata columns:
      seqnames          ranges strand
         <Rle>       <IRanges>  <Rle>
  [1]     chr1 4661501-4661800      *
  [2]     chr1 4775606-4775905      *
  [3]     chr1 4847584-4847883      *
  [4]     chr1 5015913-5016212      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Extracting sequences under regions
============================

Once we have recentered our peaks we can use the **getSeq** function with our GRanges of resized common peaks and the BSgenome object for mm9.

The **getSeq** function returns a *DNAStringSet* object containing sequences under peaks. 

Here we provide names to the elements of the DNAStringSet as we did for GRanges objects earlier.


```r
commonPeaksSequences <- getSeq(genome,GRanges(commonPeaks))
names(commonPeaksSequences) <- paste0("peak_",seqnames(commonPeaks),"_",
                                         start(commonPeaks),
                                         "-",
                                         end(commonPeaks))

commonPeaksSequences[1:2,]
```

```
  A DNAStringSet instance of length 2
    width seq                                          names               
[1]   300 AGTAGTACACAGTTAAAGCAA...TGAGGCTTTGAAGTTGAGAC peak_chr1_4661501...
[2]   300 CAGAGTGACGCGGCCCCTGCA...ATCCGCGTCGGTAGGCTATG peak_chr1_4775606...
```

Writing to FASTA file
============================

The *writeXStringSet* function allows the user to write DNA/RNA/AA(aminoacid)StringSet objects out to file. 

By default the *writeXStringSet* function writes the sequence information in FASTA format (as required for Meme-ChIP).


```r
writeXStringSet(commonPeaksSequences,file="consensusPeaks.fa")
```

Now the file "consensusPeaks.fa" contains sequences around the geometric centre of peaks suitable for Motif analysis in Meme-ChIP. 

In your own work you will typically run this from a big machine (such as a cluster) with Meme installed locally but today we will upload our generated FASTA file to their [web portal](http://meme-suite.org/tools/meme-chip). 

Results files from Meme-ChIP can be found [here](http://mrccsc.github.io/myc_Meme_Example/meme-chip.html)

Time for an exercise
=========================================================
Exercise on "Identifying Motifs" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/IdentifingMotifs_Exercises.html).

Time for a solution.
=========================================================

Answers on "Identifying Motifs" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/IdentifingMotifs_Solutions.html).

Rcode for "Identifying Motifs" solutions can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/IdentifingMotifs.R).


Getting hold of Data
========================================================
id: external

Often external datasets for ChIP-seq would be useful to integrate with your own analysis.

AnnotationHub provides a nice interactive interface to retreive data from a range of repositories covering epigenetic or expression data into R.

Try the code below for yourself.



```r
library(AnnotationHub)
ah = AnnotationHub()
rowResults <- display(ah)
```

Getting hold of Data
========================================================

AnnotationHub may also be used in non-interactive modes.

To search AnnotationHub by keywords we can use the *query* function on the AnnotationHub object. This provides information on how to retrieve data too.

```r
query(ah, c("Myc","BED", "Mus Musculus"))
```

```
AnnotationHub with 1 record
# snapshotDate(): 2019-05-02 
# names(): AH28051
# $dataprovider: Haemcode
# $species: Mus musculus
# $rdataclass: GRanges
# $rdatadateadded: 2015-03-12
# $title: c-Myc_GSM912934_MEL.bed
# $description: peak file from Haemcode
# $taxonomyid: 10090
# $genome: mm10
# $sourcetype: BED
# $sourceurl: http://haemcode.stemcells.cam.ac.uk/blood/Peaks/mm10/c-My...
# $sourcesize: 303496
# $tags: c("c-Myc_GSM912934_MEL", "Mouse ErythroLeukaemic", "[CL]
#   MEL", "BTO:0004475", "Myc", "17869", "GSE36030", "GSM912934",
#   "tf") 
# retrieve record with 'object[["AH28051"]]' 
```

```r
# cmycAnnoHub <- ah[["AH28051"]]
# cmycAnnoHub[1:3,]
```
 
Exporting tracks for Visualisation
==========================================================
id: vis

Having produced our consensus sets or GRanges of any description it is useful to visualise this in a genome browser. 

One fast, locally installed browser is the Broad's Integrative Genome Browser (IGV).
IGV is available from [BROAD](http://www.broadinstitute.org/software/igv/download) and our quick course in IGV is available [here](https://github.com/LMSBioinformatics/LMS_IGV_course/).

To export GRanges from R into a ".bed" format acceptable to IGV (or other browser types) we can use the export.bed function from rtracklayer.


```r
library(rtracklayer)
export.bed(commonPeaks,con = "consensusPeaksForIGV.bed")
```

<!-- Time for an exercise -->
<!-- ========================================================= -->
<!-- Exercise on "External Data and Visualisation" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/External_Data_and_Visualisation_Exercises.html). -->

<!-- Installing a package to a personal directory. -->
<!-- ========================================================= -->

<!-- Sometimes you will be working on a machine you don't have full permissions for. -->

<!-- In these cases it is possible to install into and use a library from a local directory. -->

<!-- * Create a directory on your desktop called Rlibs. -->

<!-- * Download zipped file for [BSgenome.Mmusculus.UCSC.mm10  package](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Mmusculus.UCSC.mm10.html) required for the bonus question manually. [(Direct Download)](http://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Mmusculus.UCSC.mm10_1.4.0.tar.gz) -->

<!-- * Install using install.packages using the Rlibs directory as your library location. -->

<!-- ```{r,eval=F,echo=T,cache=T} -->
<!-- # Something like this -->
<!-- install.packages( "C:/Users/MYUSERNAME/Desktop/BSgenome.Mmusculus.UCSC.mm10_1.4.0.tar.gz", lib="C:/Users/MYUSERNAME/Desktop/Rlibs", repos = NULL, type = "source" -->
<!--                   ) -->

<!-- library(BSgenome.Mmusculus.UCSC.mm10, lib.loc="C:/Users/MYUSERNAME/Rlibs" -->
<!--         ) -->
<!-- ``` -->


<!-- Time for a solution. -->
<!-- ========================================================= -->

<!-- Answers on "External Data and Visualisation" can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/External_Data_and_Visualisation_Solutions.html). -->

<!-- Rcode for "External Data and Visualisation" solutions can be found [here](https://lmsbioinformatics.github.io/LMS_ChIPseq_short/course/presentations/practicals/External_Data_and_Visualisation.R). -->


Working with complex overlaps
=========================================================
id: complexOverlaps

When working with a larger number of transcription factor marks it can be useful to establish a common flattened peak set for all marks and to score the overlap for each peak set to flattened peaks.

So, first lets read in the data and flatten all peaksets into one set.


```r
listOfPeaks <- GRangesList(lapply(macsPeaksFiles,
                                  function(x)ChIPQC:::GetGRanges(x,sep="\t",simplify=T)
                                  )
                           )
flattenedPeaks <- reduce(unlist(listOfPeaks))
```

=========================================================
![alt text](imgs/flattened.png)

The next step would be to identify when samples shared peaks
========================================================

```r
matOfOverlaps <- sapply(listOfPeaks,function(x)
                          as.integer(flattenedPeaks %over% x)
                        )

colnames(matOfOverlaps) <- basename(gsub("_peaks\\.xls","",macsPeaksFiles))


mcols(flattenedPeaks) <- matOfOverlaps

flattenedPeaks[1:2,]
```

```
GRanges object with 2 ranges and 4 metadata columns:
      seqnames          ranges strand | mycch12rep1 mycch12rep2 mycmelrep1
         <Rle>       <IRanges>  <Rle> |   <integer>   <integer>  <integer>
  [1]        1 3049670-3049833      * |           0           0          1
  [2]        1 3435991-3436154      * |           0           0          1
      mycmelrep2
       <integer>
  [1]          0
  [2]          0
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

========================================================
We can get a quick idea about where overlaps occur using vennCounts


```r
library(limma)
vennCounts(as.data.frame(mcols(flattenedPeaks)))
```

```
   mycch12rep1 mycch12rep2 mycmelrep1 mycmelrep2 Counts
1            0           0          0          0      0
2            0           0          0          1   7155
3            0           0          1          0  18588
4            0           0          1          1  15314
5            0           1          0          0  16038
6            0           1          0          1   1053
7            0           1          1          0   2178
8            0           1          1          1   3327
9            1           0          0          0   3589
10           1           0          0          1    205
11           1           0          1          0    271
12           1           0          1          1    414
13           1           1          0          0  14409
14           1           1          0          1   1322
15           1           1          1          0   1238
16           1           1          1          1   9868
attr(,"class")
[1] "VennCounts"
```


========================================================
Or we can view as VennDiagram


```r
vennDiagram(as.data.frame(elementMetadata(flattenedPeaks)))
```

![plot of chunk unnamed-chunk-49](ChIPseq-figure/unnamed-chunk-49-1.png)

We can check the Venn to see our numbers add up

========================================================
Now we can identify common peaks


```r
mych12Peaks <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycch12rep1 + elementMetadata(flattenedPeaks)$mycch12rep2 == 2]
mycMelPeaks <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +  elementMetadata(flattenedPeaks)$mycmelrep2 == 2]

mych12Peaks[1:3,]
```

```
GRanges object with 3 ranges and 4 metadata columns:
      seqnames          ranges strand | mycch12rep1 mycch12rep2 mycmelrep1
         <Rle>       <IRanges>  <Rle> |   <integer>   <integer>  <integer>
  [1]        1 4661367-4661934      * |           1           1          0
  [2]        1 4775337-4776125      * |           1           1          1
  [3]        1 4847097-4848050      * |           1           1          1
      mycmelrep2
       <integer>
  [1]          0
  [2]          1
  [3]          1
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

========================================================
![alt text](imgs/highCon.png)


========================================================
And some unique peaks 

```r
mycMelPeaks_Only <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 + elementMetadata(flattenedPeaks)$mycmelrep2 == 2 &
elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 0]

mycMelPeaks_Only[1,]
```

```
GRanges object with 1 range and 4 metadata columns:
      seqnames          ranges strand | mycch12rep1 mycch12rep2 mycmelrep1
         <Rle>       <IRanges>  <Rle> |   <integer>   <integer>  <integer>
  [1]        1 7606348-7606524      * |           0           0          1
      mycmelrep2
       <integer>
  [1]          1
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

========================================================
![alt text](imgs/Mel_Only.png)





Simple Differential binding
========================================================
id: diffchip

Analysis of the differences in ChIP-seq data can often benefit from a more quantitative analysis. Many tools used in RNA-seq analysis studies can be applied to ChIP data including favorites such as DEseq2 and EdgeR.

Inorder to assess difference we first needed to identify peaks common within groups. Here we take our previously prepared set made of flattened peaks and identify those reproducible peaks in either group.



```r
highConfidence_Only <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +                  elementMetadata(flattenedPeaks)$mycmelrep2 == 2 |
elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 2]
```

========================================================
![alt text](imgs/forDB.png)

Simple Differential binding
========================================================
Now we can look to see if we need resizing.

```r
boxplot(width(highConfidence_Only))
abline(h=400,col="red")
```

![plot of chunk unnamed-chunk-53](ChIPseq-figure/unnamed-chunk-53-1.png)

The majority of peaks are around 400 so we will resize all peaks to this for ease here

Simple Differential binding
========================================================
Now we can resize to a sensible size

```r
PeaksToCount <- resize(highConfidence_Only,width = 400,fix = "center")
PeaksToCount[1:2,]
```

```
GRanges object with 2 ranges and 4 metadata columns:
      seqnames          ranges strand | mycch12rep1 mycch12rep2 mycmelrep1
         <Rle>       <IRanges>  <Rle> |   <integer>   <integer>  <integer>
  [1]        1 4661451-4661850      * |           1           1          0
  [2]        1 4775531-4775930      * |           1           1          1
      mycmelrep2
       <integer>
  [1]          0
  [2]          1
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Simple Differential binding
========================================================

Once we have our peakset we can can count the reads mapping to peaks from all samples.
Many options are available for counting including Rsubread package's FeatureCounts method and GenomicAlignments' summarizeOverlaps function.

For a more detailed description of counting, please take a look at some of our previous material.

For now, we can import the counts produced on this course into our present R session.


```r
load("data/robjects/MycCounts.Rdata")
countTable <- countTable[,-c(3,6)]

countTable[1:3,]
```

```
                      ch12myc ch12myc melmyc melmyc
ID1-1;4661451-4661850      45      81      0      1
ID2-1;4775531-4775930      56      68     59     47
ID3-1;4847374-4847773      51      46     64     70
```

Simple Differential binding - A simple DEseq2 DE analysis
========================================================
Here we set up a DEseq2 object much as you would for RNAseq.

We define the conditions in **colData** as celllines for Mel and ch12 as in RNAseq and provide a table of the counts in peak.

In contrast to the RNAseq DEseq2 setup we will provide additional information of the GRanges for peaks to the rowRanges argument.


```r
library("DESeq2")

colData <- data.frame(SampleName=paste0(colnames(countTable),c(1,2)), CellLine=c("ch12","ch12","mel","mel"))

colnames(countTable) <-  colData$SampleName
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = colData,
                              design = ~ CellLine,
                              rowRanges=PeaksToCount)

dds <- DESeq(dds)
```

Simple Differential binding - A simple DEseq2 DE analysis
========================================================

Now we have our dds object with rowRanges we can extract the results for differences between celllines as seen for RNA-seq.

In contrast to RNA-seq, we specify the "GRanges" format in our call to *results* function.

Finally we subset our GRanges to produce a GRanges of peaks with significantly higher signal in Mel cell line.


```r
test_cellline <- results(dds, contrast=c("CellLine","ch12","mel"),
                         format="GRanges")

UpinMel <- test_cellline[test_cellline$padj < 0.05 & !is.na(test_cellline$padj) 
                         & test_cellline$log2FoldChange > 0]

UpinMel
```

```
GRanges object with 14634 ranges and 6 metadata columns:
                                seqnames              ranges strand |
                                   <Rle>           <IRanges>  <Rle> |
          ID1-1;4661451-4661850        1     4661451-4661850      * |
          ID4-1;5015863-5016262        1     5015863-5016262      * |
          ID6-1;5210772-5211171        1     5210772-5211171      * |
          ID7-1;5273188-5273587        1     5273188-5273587      * |
          ID9-1;6252315-6252714        1     6252315-6252714      * |
                            ...      ...                 ...    ... .
  ID45884-X;166410701-166411100        X 166410701-166411100      * |
  ID45885-X;166417005-166417404        X 166417005-166417404      * |
  ID45886-X;166428102-166428501        X 166428102-166428501      * |
  ID45887-X;166434992-166435391        X 166434992-166435391      * |
        ID45889-Y;307700-308099        Y       307700-308099      * |
                                        baseMean   log2FoldChange
                                       <numeric>        <numeric>
          ID1-1;4661451-4661850 33.4987258065535 7.22030032288329
          ID4-1;5015863-5016262  11.492741487773 4.01128347553125
          ID6-1;5210772-5211171 12.5818260145118  4.7517271021384
          ID7-1;5273188-5273587 15.8972913963059 4.50076914167956
          ID9-1;6252315-6252714 24.4950641604032 5.73897203757724
                            ...              ...              ...
  ID45884-X;166410701-166411100 73.6663018288351 4.38898230258325
  ID45885-X;166417005-166417404 108.464233028677 3.15624745989976
  ID45886-X;166428102-166428501 89.0604007646183 2.84317416306789
  ID45887-X;166434992-166435391 48.2298270233863 6.15382638780086
        ID45889-Y;307700-308099 28.8586281179605 7.00076407045592
                                            lfcSE             stat
                                        <numeric>        <numeric>
          ID1-1;4661451-4661850  1.53570142126731 4.70163029277192
          ID4-1;5015863-5016262   1.2177153892343 3.29410592244675
          ID6-1;5210772-5211171  1.35143174162522 3.51606888885413
          ID7-1;5273188-5273587  1.12528610908257 3.99966648957301
          ID9-1;6252315-6252714  1.20609297125709 4.75831646012795
                            ...               ...              ...
  ID45884-X;166410701-166411100  0.53059084615728 8.27187716179002
  ID45885-X;166417005-166417404 0.389907618250926 8.09485968511788
  ID45886-X;166428102-166428501 0.423258322223824 6.71734969824028
  ID45887-X;166434992-166435391 0.955056503006767 6.44341603708996
        ID45889-Y;307700-308099  1.54835524166375 4.52141981508933
                                              pvalue                 padj
                                           <numeric>            <numeric>
          ID1-1;4661451-4661850 2.58092492556339e-06  1.6023242246206e-05
          ID4-1;5015863-5016262 0.000987352841511528  0.00246450359916256
          ID6-1;5210772-5211171 0.000437987324473888  0.00122018541217481
          ID7-1;5273188-5273587 6.34318107858661e-05 0.000233385124716184
          ID9-1;6252315-6252714 1.95214281072272e-06 1.26876841622556e-05
                            ...                  ...                  ...
  ID45884-X;166410701-166411100 1.31866129082079e-16 2.07958776489167e-14
  ID45885-X;166417005-166417404 5.73303977828751e-16 7.67057322172508e-14
  ID45886-X;166428102-166428501 1.85059486433119e-11 7.06521017636951e-10
  ID45887-X;166434992-166435391 1.16813906138065e-10 3.58104460954447e-09
        ID45889-Y;307700-308099  6.1426228115764e-06 3.29935915342772e-05
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Session Information
=========================================================


```r
sessionInfo()
```

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.2

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] limma_3.40.6                           
 [2] BSgenome.Mmusculus.UCSC.mm9_1.4.0      
 [3] BSgenome_1.52.0                        
 [4] rtracklayer_1.44.4                     
 [5] Biostrings_2.52.0                      
 [6] XVector_0.24.0                         
 [7] goseq_1.36.0                           
 [8] geneLenDataBase_1.20.0                 
 [9] BiasedUrn_1.07                         
[10] KEGG.db_3.2.3                          
[11] ChIPseeker_1.20.0                      
[12] org.Mm.eg.db_3.8.2                     
[13] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2
[14] GenomicFeatures_1.36.4                 
[15] AnnotationDbi_1.46.1                   
[16] DESeq2_1.24.0                          
[17] ChIPQC_1.20.0                          
[18] DiffBind_2.12.0                        
[19] SummarizedExperiment_1.14.1            
[20] DelayedArray_0.10.0                    
[21] BiocParallel_1.18.1                    
[22] matrixStats_0.55.0                     
[23] Biobase_2.44.0                         
[24] GenomicRanges_1.36.1                   
[25] GenomeInfoDb_1.20.0                    
[26] IRanges_2.18.2                         
[27] S4Vectors_0.22.1                       
[28] BiocGenerics_0.30.0                    
[29] ggplot2_3.2.1                          
[30] knitr_1.25                             

loaded via a namespace (and not attached):
  [1] tidyselect_0.2.5                         
  [2] RSQLite_2.1.2                            
  [3] htmlwidgets_1.3                          
  [4] grid_3.6.0                               
  [5] munsell_0.5.0                            
  [6] base64url_1.4                            
  [7] systemPipeR_1.18.2                       
  [8] withr_2.1.2                              
  [9] colorspace_1.4-1                         
 [10] GOSemSim_2.10.0                          
 [11] Category_2.50.0                          
 [12] highr_0.8                                
 [13] rstudioapi_0.10                          
 [14] DOSE_3.10.2                              
 [15] labeling_0.3                             
 [16] urltools_1.7.3                           
 [17] GenomeInfoDbData_1.2.1                   
 [18] hwriter_1.3.2                            
 [19] polyclip_1.10-0                          
 [20] bit64_0.9-7                              
 [21] farver_1.1.0                             
 [22] pheatmap_1.0.12                          
 [23] batchtools_0.9.11                        
 [24] vctrs_0.2.0                              
 [25] TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2  
 [26] xfun_0.9                                 
 [27] R6_2.4.0                                 
 [28] graphlayouts_0.5.0                       
 [29] locfit_1.5-9.1                           
 [30] gridGraphics_0.4-1                       
 [31] bitops_1.0-6                             
 [32] fgsea_1.10.1                             
 [33] assertthat_0.2.1                         
 [34] scales_1.0.0                             
 [35] ggraph_2.0.0                             
 [36] nnet_7.3-12                              
 [37] enrichplot_1.4.0                         
 [38] gtable_0.3.0                             
 [39] tidygraph_1.1.2                          
 [40] rlang_0.4.0                              
 [41] zeallot_0.1.0                            
 [42] genefilter_1.66.0                        
 [43] splines_3.6.0                            
 [44] lazyeval_0.2.2                           
 [45] acepack_1.4.1                            
 [46] europepmc_0.3                            
 [47] brew_1.0-6                               
 [48] checkmate_1.9.4                          
 [49] yaml_2.2.0                               
 [50] reshape2_1.4.3                           
 [51] TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2
 [52] backports_1.1.4                          
 [53] qvalue_2.16.0                            
 [54] Hmisc_4.2-0                              
 [55] RBGL_1.60.0                              
 [56] tools_3.6.0                              
 [57] gridBase_0.4-7                           
 [58] ggplotify_0.0.4                          
 [59] gplots_3.0.1.1                           
 [60] RColorBrewer_1.1-2                       
 [61] ggridges_0.5.1                           
 [62] Rcpp_1.0.2                               
 [63] plyr_1.8.4                               
 [64] base64enc_0.1-3                          
 [65] progress_1.2.2                           
 [66] zlibbioc_1.30.0                          
 [67] purrr_0.3.2                              
 [68] RCurl_1.95-4.12                          
 [69] prettyunits_1.0.2                        
 [70] rpart_4.1-15                             
 [71] viridis_0.5.1                            
 [72] cowplot_1.0.0                            
 [73] chipseq_1.34.0                           
 [74] ggrepel_0.8.1                            
 [75] cluster_2.1.0                            
 [76] magrittr_1.5                             
 [77] data.table_1.12.2                        
 [78] TxDb.Hsapiens.UCSC.hg18.knownGene_3.2.2  
 [79] DO.db_2.9                                
 [80] triebeard_0.3.0                          
 [81] amap_0.8-17                              
 [82] hms_0.5.1                                
 [83] evaluate_0.14                            
 [84] xtable_1.8-4                             
 [85] XML_3.98-1.20                            
 [86] gridExtra_2.3                            
 [87] compiler_3.6.0                           
 [88] biomaRt_2.40.4                           
 [89] tibble_2.1.3                             
 [90] KernSmooth_2.23-15                       
 [91] crayon_1.3.4                             
 [92] htmltools_0.3.6                          
 [93] GOstats_2.50.0                           
 [94] mgcv_1.8-28                              
 [95] Formula_1.2-3                            
 [96] tidyr_1.0.0                              
 [97] geneplotter_1.62.0                       
 [98] DBI_1.0.0                                
 [99] tweenr_1.0.1                             
[100] formatR_1.7                              
[101] MASS_7.3-51.4                            
[102] rappdirs_0.3.1                           
[103] boot_1.3-23                              
[104] ShortRead_1.42.0                         
[105] Matrix_1.2-17                            
[106] gdata_2.18.0                             
[107] igraph_1.2.4.1                           
[108] pkgconfig_2.0.2                          
[109] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2  
[110] rvcheck_0.1.3                            
[111] GenomicAlignments_1.20.1                 
[112] foreign_0.8-72                           
[113] TxDb.Celegans.UCSC.ce6.ensGene_3.2.2     
[114] xml2_1.2.2                               
[115] annotate_1.62.0                          
[116] AnnotationForge_1.26.0                   
[117] stringr_1.4.0                            
[118] VariantAnnotation_1.30.1                 
[119] digest_0.6.20                            
[120] graph_1.62.0                             
[121] fastmatch_1.1-0                          
[122] htmlTable_1.13.1                         
[123] edgeR_3.26.8                             
[124] GSEABase_1.46.0                          
[125] Rsamtools_2.0.0                          
[126] gtools_3.8.1                             
[127] rjson_0.2.20                             
[128] nlme_3.1-141                             
[129] lifecycle_0.1.0                          
[130] jsonlite_1.6                             
[131] viridisLite_0.3.0                        
[132] pillar_1.4.2                             
[133] lattice_0.20-38                          
[134] Nozzle.R1_1.1-1                          
[135] plotrix_3.7-6                            
[136] httr_1.4.1                               
[137] survival_2.44-1.1                        
[138] GO.db_3.8.2                              
[139] glue_1.3.1                               
[140] UpSetR_1.4.0                             
[141] bit_1.1-14                               
[142] Rgraphviz_2.28.0                         
[143] ggforce_0.3.1                            
[144] stringi_1.4.3                            
[145] blob_1.2.0                               
[146] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.7 
[147] latticeExtra_0.6-28                      
[148] caTools_1.17.1.2                         
[149] memoise_1.1.0                            
[150] dplyr_0.8.3                              
```

THE END!
=========================================================


