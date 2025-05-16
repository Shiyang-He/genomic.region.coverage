#!usr/bin/python3
#coding:utf-8
import sys,os
from optparse import OptionParser
import tempfile
from pybedtools import BedTool
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import numpy as np
from bam_coverage_report import compute_coverage_stats
from bam_coverage_report import write_coverage_report
from bam_coverage_report import transform_depth_file_to_bedgraph
from bam_coverage_report import parse_gtf_for_gene
from bam_coverage_report import convert_to_bigwig 
from bam_coverage_report import draw_piechart 
from bam_coverage_report import draw_pyGenomeTracks

def parseOptions():
    desc=''' this script compute the coverage for a region or the whole genome for a bam file, using samtools, you can give a bed file, a region, or a gene to compute the coverage'''
    epilog='''usages:
python3 depth.with.different.tools.py -i input.bam -o output -t samools -g ACTD -d -a gencode.v35.annotation.gtf
python3 depth.with.different.tools.py  -i input.bam -o output -t bamdst -b a.bed -f bigwig -G GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
python3 depth.with.different.tools.py -i input.bam -o output -r chr:start-end -d
'''
    parser=OptionParser(description=desc,epilog=epilog)
    parser.add_option("-i", "--input", help="Input BAM/CRAM file")
    parser.add_option("-g", "--gene", help="Gene name or geneID to extract region from (requires -a parameter, the annotation file needed to tell where this gene is)")
    parser.add_option("-p", "--picture_format", default="pdf", help="you want a PDF file or png file")
    parser.add_option("-r", "--region", help="Genomic region (chr:start-end)")
    parser.add_option("-b", "--bed", help="BED file of regions")
    parser.add_option("-q", "--quality", default = 0, help="min mapping quality to include")
    parser.add_option("-o", "--output", default="output", help="Base output name (no extension)")
    parser.add_option("-t", "--tools", default="samtools",help="what tools to use? [samtools or bamdst] if bamdst is used, the output will be a directory and bed file is mandatory, you can gerate a large bed file for this, a region will be transformed to a bed file, I wouldn't use bamdst if the region I am interested is not large as it's slow")
    parser.add_option("-f", "--outputformat", choices=["bedgraph", "bigwig"], default="bedgraph", help="Output format, default is bedgraph format, if you want a bigwig file for the output, please provide -G parameter, the genome file, it can be genome.fa.fai generated from samtools faidx")
    parser.add_option("-d", "--draw", action="store_true", help="Draw coverage track PDF")
    parser.add_option("-G", "--Genomefile", help="Genome file for BigWig conversion")
    parser.add_option("-a", "--annotation", help="GTF or GFF file for gene region extraction")
    return parser.parse_args()

def compute_depth_with_samtools(input_bam,depth_file,quality=0,bed=None,region=None):
    cmd=[]
    if bed:
        cmd=['samtools','depth','-b',bed,'-o',depth_file,'-Q',quality,input_bam]
    elif region:
        cmd=['samtools','depth','-r',region,'-o',depth_file,'-Q',quality,input_bam]
    else:
        cmd=['samtools','depth','-o',depth_file,'-Q',quality,input_bam]
    subprocess.run(cmd, check=True)
def compute_depth_with_bamdst(input_bam,output_dir=None,quality=0,bed=None,region=None):
    cmd=[]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    coverage_file=f"{output_dir}/coverage.txt"
    if bed:
        cmd=['bamdst','-p',bed,'-o',output_dir,'-q',quality,input_bam]
    if region:
        #gerate a bed file to read
        bed_region=region.replace(":","\t").replace("-","\t")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
            tmp.write(bed_region+"\n")
            tmp_path=tmp.name
        cmd=['bamdst','-p',tmp_path,'-o',output_dir,'-q',quality,input_bam]
    subprocess.run(cmd, check=True)
    os.remove(tmp_path)
    import gzip
    gzipped_depth=f"{output_dir}/depth.tsv.gz"
    coverage_file_handle=open(coverage_file,'w')
    with gzip.open(gzipped_depth,'rb') as f:
        for line in f:
            region_line=line.decode('utf-8').split()
            if not region_line[0].startswith("#"):
                coverage_file_handle.write("\t".join(map(str,region_line[0:3]))+"\n")

def read_coverage_file(coverage_file):
    df=pd.read_csv(coverage_file,header=None,sep="\t")
    df.columns=["chr","pos","depth"]
    depths=np.array(df['depth'])
    return depths
        
def main():
    if len(sys.argv)<2:
        print('''too few options, please type: 
python3 depth.with.different.tools.py -h
to see more info''')
        sys.exit(1)
    else:
        (args,options)=parseOptions()
        region = args.region
        quality = str(args.quality)
        bed = args.bed
        annotation= args.annotation
        tools = args.tools
        gene = args.gene
        draw = args.draw
        genome_file = args.Genomefile
        input_bam = args.input
        picture_format=args.picture_format
        output_prefix = args.output

        if gene and not annotation:
            print("A gene is given, -a (the annotation) must be provided")
            sys.exit(1)

        # Step 1 calculate coverage and draw a piechart
        coverage_file= f"{output_prefix}.coverage.txt"
        report_file = f"{output_prefix}.coverage.report.txt"
        if gene:
            # read gene region from a gtf/gff file
            region=parse_gtf_for_gene(annotation,gene)
            print(f"a gene is provided, use its region {region} to calculate")

        depths=[]
        depth_stat={} 
        if tools == "samtools":
            compute_depth_with_samtools(input_bam,coverage_file,quality=quality,bed=bed,region=region)
            depths=read_coverage_file(f"{output_prefix}.coverage.txt")
            piechart_value, coverage_report = compute_coverage_stats(depths,len(depths))

        if tools == "bamdst":
            output_dir=output_prefix
            compute_depth_with_bamdst(input_bam,output_dir=output_prefix,quality=quality,bed=bed,region=region)
            coverage_file = f"{output_dir}/coverage.txt"
            depths=read_coverage_file(f"{output_dir}/coverage.txt")
            piechart_value, coverage_report = compute_coverage_stats(depths,len(depths))
        write_coverage_report(coverage_report,report_file)

        # draw a piechart
        piechart_file = f"{output_prefix}.coverage.piechart.pdf"
        if picture_format == "PDF":
            picture_format = "pdf"
        if picture_format == "PNG" or picture_format =="png":
            piechart_file = f"{output_prefix}.coverage.piechart.png"
        labels=piechart_value.keys()
        sizes=piechart_value.values()
        draw_piechart(labels,sizes,piechart_file)

        # Step 2 Generate coverage track
        # convert the depth file to bedgraph format
        bedgraph_file = f"{output_prefix}.bedgraph"
        transform_depth_file_to_bedgraph(coverage_file,bedgraph_file)
        if args.outputformat == "bigwig":
            # covert to bigwig file, the genome size is required
            if not genome_file:
                print("please give me a genome size file with argument: -G")
                sys.exit(1)
            bigwig_file = f"{output_prefix}.bw"
            convert_to_bigwig(bedgraph_file, genome_file, bigwig_file)
#            os.remove(bedgraph_file)

        # Step 3: draw a track picture for the given region or gene
        if draw :
            bigwig_input = bigwig_file if args.outputformat == "bigwig" else bedgraph_file
            output_pdf = f"{output_prefix}.coverage.pdf"
            if picture_format == "PNG" or picture_format == "png": 
                output_pdf = f"{output_prefix}.coverage.png"
            if region:
                draw_pyGenomeTracks(bigwig_input, region, output_pdf)
            if bed:
                output_pdf_dir=f"{output_prefix}.track.files"
                if not os.path.exists(output_pdf_dir):
                    os.mkdir(output_pdf_dir)
                bed_regions = BedTool(bed)
                for interval in bed_regions:
                    tmp_region=interval.chrom+":"+str(interval.start)+"-"+str(interval.end)
                    tmp_region_file=f"{output_prefix}.track.files/{interval.chrom}_{interval.start}_{interval.end}.pdf"
                    if picture_format == "PNG" or picture_format == "png":
                        tmp_region_file=f"{output_prefix}.track.files/{interval.chrom}_{interval.start}_{interval.end}.png"
                    draw_pyGenomeTracks(bigwig_input, tmp_region, tmp_region_file)
    
            if gene:
                if not annotation: # quit if annotation file is not provided
                    print("Annotation file required for gene plotting. Use -a to specify GTF/GFF. or you can just provide the region of that gene to draw this plot")
                    sys.exit(1)
        
                region_from_gene = parse_gtf_for_gene(annotation, gene)
                if not region_from_gene:
                    print(f"Gene {gene} not found in annotation.")
                    sys.exit(1)
        
                bigwig_input = bigwig_file if args.outputformat == "bigwig" else bedgraph_file
                output_pdf = f"{output_prefix}.coverage.pdf"
                if picture_format == "PNG" or picture_format == "png": 
                    output_pdf = f"{output_prefix}.coverage.png"
                draw_pyGenomeTracks(bigwig_input, region_from_gene, output_pdf)


if __name__ == "__main__":
    main()
