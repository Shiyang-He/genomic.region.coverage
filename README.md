# genomic.region.coverage
  This folder contains the script to compute sequencing depth for NGS data. 
The input will be a bam file. You will need to know which genomic reference 
is used and download its genome.size file, or use samtools to generate it 
from genome.fasta file. To find out which reference was used, use the command: 
samtools view -H, you will see the genomic build info. 

Stat region can be given in 4 ways (a bed file, a region, a gene, or nothing 
given which will compute the whole genome region)

The idea is to give a short coverage report and a pie chart to show how well 
the overall sequencing depth is. If you give a region, a bed, or a gene, I will 
only show that region and draw the track picture to show how the sequencing depth 
is in that region. 

In addition, I also generate the browserable track files (bedgraph or bigwig file)
for future use. It can be loaded into a IGV browser or UCSC genome browser to see
the sequencing depth in whichever gene or region. Of course you if you are interested
in certain gene or region, follow the instruction of my script to generate the 
report and pictures. 

Required packages or software:

(versions are not strickly required to be the same but I included what I used in 
case your versions conficts. For example, pyGenomeTracks may conflict with some of 
the installed packages, so better to create a new environment and install everything needed)

samtools 1.14

bedGraphToBigWig v4

Python 3.9.18

pyGenomeTracks 3.9

optparse 1.5.3

pybedtools 0.9.0

pysam 0.19.1

numpy 1.22.4

pandas 1.4.2

matplotlib 3.5.1

tempfile 3.9.18 (same as python)

subprocess 3.9.18 (same as python)


Required files to run all functions:
bamfile, I only put a small bam file here for testing. It's generate by this command: 
samtools view -L test.2.bed -b -O NA12878.wgs.toy.bam NA12878.wgs.bam && samtools index NA12878.wgs.toy.bam
GRCh38_full_analysis_set_plus_decoy_hla.fa.fai (just my example, use yours when your bam file changes)
gencode.v35.annotation.gtf ( https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz)

example usage: 
python3 bam_coverage_report.py \
    -i NA12878.wgs.bam \
    -o NA12878.wgs.ENSG00000227232.5 \
    -g ENSG00000227232.5 \
    -a gencode.v35.annotation.gtf \
    -f bigwig \
    -t True \
    -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
    -d 

python3 bam_coverage_report.py \
    -i NA12878.wgs.bam \
    -o NA12878.wgs \
    -f bigwig \
    -t True \
    -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
    
python3 bam_coverage_report.py \
    -i NA12878.wgs.bam \
    -o NA12878.wgs.bed \
    -b test.bed \
    -f bigwig \
    -d \
    -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai 
    
python3 bam_coverage_report.py \
    -i NA12878.wgs.bam \
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
                        Base output name (no extension)
                        
  -t, --track           output track file? [True or False] The file format can
                        be bedgraph or bigwig
                        
  -f OUTPUTFORMAT, --outputformat=OUTPUTFORMAT
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

The script depth.with.different.tools.py works as well, it uses many of the functions from bam_coverage_report.py, 
just to use other tools (samtools or bamdst) to do the samething. 
bamdst 1.0.9 can be found here: https://github.com/shiquan/bamdst

The benefit of bamdst is as it gives more comprehensive coverage report for a bed region. Check the bamdst.output folder for those report.
This script has similar parameters as bam_coverage_report.py. 
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
                        




