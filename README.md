# genomic.region.coverage
This folder contains two python scripts to compute sequencing depth for NGS data.
The input is a bam file from NA12878 WGS mapping to GRCH38 reference. If use this 
script to count other bam files, you will need to know which genomic reference 
is used and download its genome.size file, or use samtools to generate it 
from genome.fasta file. To find out which reference was used, use the command: 

        samtools view -H input.bam 

The goal here is to give a short coverage report and a pie chart to show how well 
the overall sequencing depth is. The coverage report tells how deep the sequencing
is and how well the genomic region is been sequenced. Also include the option to 
generate browsable track files (bedgraph/bigwig).

If you give a region, a bed, or a gene, I will only compute coverage in that region 
and draw the track picture to show how the sequencing depth is in that region. This is 
useful if you did sequencing and want to find out how well some gene or region is sequenced. 
Bedgraph/bigwig files can be loaded into a IGV browser or UCSC genome browser to see
the sequencing depth in whichever gene or region. See these two pages to learn more:
https://genome.ucsc.edu/goldenPath/help/bedgraph.html
https://genome.ucsc.edu/goldenPath/help/bigWig.html
Of course you if you are interested in certain gene or region but don't want use IGV,
follow the instruction of my script to generate the short report and pictures. 

Here are the steps to do this job:
1) I will generate a bunch of options to give input file, output file, if pictures will be
   drawn, if you want a (bigwig/bedgraph) track file, etc.
2) Read the bam file using pysam, it has a pileup option, similar to 'samtools pileup' which
   gives sequencing depth at each position. Write the depth at each position to a coverage file.
   store depth info in a list to calculate percentile etc.
   If a gene is given, a gtf/gff must be given to get the genomic region for that gene. It can
   be gene name or gene ID. 
4) The coverage report will be printed in a coverage.report.txt file.
5) matplotlib will be used to draw a piechart to show what's the percentage of each sequencing depth.
6) Read the coverage file and write it as bedgraph format.
7) Use bedGraphToBigWig tool to transform the bedgraph to bigwig file
8) The pyGenomeTracks tools will be introduced to draw a IGV like track file to show overall
   coverage for given window


   
Required packages or software (versions are not strickly required to be the same but 
I included what I used in case your versions conficts. For example, pyGenomeTracks may 
conflict with some of the installed packages, so better to create a new environment 
and install everything needed):

samtools 1.14 , bedGraphToBigWig v4,  Python 3.9.18, pyGenomeTracks 3.9, optparse 1.5.3, 
pybedtools 0.9.0, pysam 0.19.1, numpy 1.22.4, pandas 1.4.2, matplotlib 3.5.1, 
tempfile 3.9.18 (same as python), subprocess 3.9.18 (same as python)

Required files to run all functions:
1) bamfile, I only put a small bam file here for testing. It's generate by this command:
 
        samtools view -L test.2.bed -b -o NA12878.wgs.toy.bam NA12878.wgs.bam && samtools index NA12878.wgs.toy.bam

If you have a CRAM file and want to change it to bam file, use this command:

        samtools view -b -o NA12878.wgs.bam NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram && samtools index NA12878.wgs.bam

3) GRCh38_full_analysis_set_plus_decoy_hla.fa.fai (just my example, use yours when your bam file changes)
4) gencode.v35.annotation.gtf ( https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz)
(I only include a small portion of this file in the github repo to show case)

example usage: 

    python3 bam_coverage_report.py \
        -i NA12878.wgs.toy.bam \
        -o NA12878.wgs.ENSG00000227232.5 \
        -g ENSG00000227232.5 \
        -a gencode.v35.annotation.gtf \
        -f bigwig \
        -t True \
        -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
        -d 

    python3 bam_coverage_report.py \
        -i NA12878.wgs.toy.bam \
        -o NA12878.wgs \
        -f bigwig \
        -t True \
        -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
    
    python3 bam_coverage_report.py \
        -i NA12878.wgs.toy.bam \
        -o NA12878.wgs.bed \
        -b test.bed \
        -f bigwig \
        -d \
        -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai 
    
    python3 bam_coverage_report.py \
        -i NA12878.wgs.toy.bam \
        -o NA12878.wgs.region \
        -r chr1:29554-31109 \
        -f bigwig \
        -d \
        -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai 

 Options explained:

  -h, --help            show this help message and exit
  
  -i INPUT, --input=INPUT
                        Input BAM/CRAM file
                        
  -g GENE, --gene=GENE  Gene name or geneID to extract region from (requires
                        -a parameter, the annotation file needed to tell where
                        this gene is)
                        
  -r REGION, --region=REGION
                        Genomic region (chr:start-end)
                        
  -b BED, --bed=BED     BED file of regions
  
  -q QUALITY, --quality=QUALITY
                        min mapping quality to include
                        
  -o OUTPUT, --output=OUTPUT
                        Base output dir name, every thing will be generate here
                        
  -t, --track           output track file? [True or False] The file format can
                        be bedgraph or bigwig
                        
  -f TRACK_FILE_FORMAT, --track_file_format=TRACK_FILE_FORMAT
                        Output format for genomic track, default is bedgraph,
                        if you want a bigwig file for the output, please
                        provide -G parameter, the genome file, it can be
                        genome.fa.fai generated from samtools faidx. If this
                        parameter is on, --track will be automatically turned on
                        
  -d, --draw            Draw coverage track PDF, this option also need track
                        file, it automatically turn on the --track option
                        
  -p PICTURE_FORMAT, --picture_format=PICTURE_FORMAT
                        you want a PDF file or png file
                        
  -G GENOMEFILE, --Genomefile=GENOMEFILE
                        Genome file for BigWig conversion
                        
  -a ANNOTATION, --annotation=ANNOTATION
                        GTF or GFF file for gene region extraction


###################################################################################################

OUTPOUT files: 

output_prefix.coverage.txt: tab delimited coverage report for every non zero position, format: chr pos coverage

output_prefix.coverage.pdf: a picture showing real coverage within certain region

output_prefix.coverage.report.txt: stats about the coverage, mean, median, max, percentle, etc.

output_prefix.coverage.piechart.pdf: a pie chart of the sequencing depth

output_prefix.bedgraph: tab delimited bedgraph file for computed region, format: chr start end coverage

output_prefix.bw: bigwig file for the coverage, it's a indexed binary file for bedgraph, which is easy to load in to genome browser.


#############################################################################################
#############################################################################################

The second script depth.with.different.tools.py works slightly different , it uses many of the 
functions from bam_coverage_report.py, however it allows using other other tools (samtools or bamdst) to 
do the same thing. 

bamdst 1.0.9 can be found here: https://github.com/shiquan/bamdst

The benefit of bamdst is as it gives more comprehensive coverage report for a bed region. 
However, its check the whole bam file and will be slow if you only interested in a small
region. Check the bamdst.output folder for those report. This script has similar parameters 
as bam_coverage_report.py. 

Example usage:

    python3 depth.with.different.tools.py \
        -i input.bam \
        -o output \
        -t samools \
        -g ACTD \
        -a gencode.v35.annotation.gtf \
        -d
    python3 depth.with.different.tools.py  \
        -i input.bam \
        -o output \
        -t bamdst \
        -b a.bed
    python3 depth.with.different.tools.py \
        -i input.bam \
        -o output \
        -r chr:start-end \
        -d


Options:

  -h, --help            show this help message and exit
  
  -i INPUT, --input=INPUT
                        Input BAM/CRAM file
                        
  -g GENE, --gene=GENE  Gene name or geneID to extract region from (requires
                        -a parameter, the annotation file needed to tell where
                        this gene is)
                        
  -r REGION, --region=REGION
                        Genomic region (chr:start-end)
                        
  -b BED, --bed=BED     BED file of regions
  
  -q QUALITY, --quality=QUALITY
                        min mapping quality to include
                        
  -o OUTPUT, --output=OUTPUT
                        Base output name (no extension)
                        
  -t TOOLS, --tools=TOOLS
                        what tools to use? [samtools or bamdst] if bamdst is
                        used, the output will be a directory and bed file is
                        mandatory, you can gerate a large bed file for this, a
                        region will be transformed to a bed file, I wouldn't
                        use bamdst if the region I am interested is not large
                        as it's slow
                        
  -f OUTPUTFORMAT, --outputformat=OUTPUTFORMAT
                        Output format, default is bedgraph format, if you want
                        a bigwig file for the output, please provide -G
                        parameter, the genome file, it can be genome.fa.fai
                        generated from samtools faidx
                        
  -d, --draw            Draw coverage track PDF
  
  -G GENOMEFILE, --Genomefile=GENOMEFILE
                        Genome file for BigWig conversion
                        
  -a ANNOTATION, --annotation=ANNOTATION
                        GTF or GFF file for gene region extraction
                        



