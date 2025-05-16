#!usr/bin/python3
#coding:utf-8
import sys,os
from optparse import OptionParser
#import pybedtools
from pybedtools import BedTool
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tempfile
import subprocess

def parseOptions():
    desc='''This script computes the coverage for a region or the whole genome for a bam file, you can give a bed, a region (chr:start-end), or a gene to compute coverage, it generate a report include the depth level (how many necleotides are covered by how many reads) and depth percentle info, also can draw a overview if you add --draw parameter for a given region or a gene'''
    epilog='''usages:
python3 bam_coverage_report.py -i input.bam -o output -g ACTD -d
python3 bam_coverage_report.py -i input.bam -o output -b a.bed
python3 bam_coverage_report.py -i input.bam -o output -r chr:start-end -d
'''
    parser=OptionParser(description=desc,epilog=epilog)
    parser.add_option("-i", "--input", help="Input BAM/CRAM file")
    parser.add_option("-g", "--gene", help="Gene name or geneID to extract region from (requires -a parameter, the annotation file needed to tell where this gene is)")
    parser.add_option("-r", "--region", help="Genomic region (chr:start-end)")
    parser.add_option("-b", "--bed", help="BED file of regions")
    parser.add_option("-q", "--quality", default = 0, help="min mapping quality to include")
    parser.add_option("-o", "--output", default="output", help="Base output dir name, every thing will be generate here")
    parser.add_option("-t", "--track", action="store_true",help="output track file? [True or False] The file format can be bedgraph or bigwig")
    parser.add_option("-f", "--track_file_format", choices=["bedgraph", "bigwig"], default="bedgraph", help="Output format for genomic track, default is bedgraph, if you want a bigwig file for the output, please provide -G parameter, the genome file, it can be genome.fa.fai generated from samtools faidx. If this parameter is on, --track will be automatically on")
    parser.add_option("-d", "--draw", action="store_true", help="Draw coverage track PDF, this option also need track file, it automatically turn on the --track option")
    parser.add_option("-p", "--picture_format", default="pdf", help="you want a PDF file or png file")
    parser.add_option("-G", "--Genomefile", help="Genome file for BigWig conversion")
    parser.add_option("-a", "--annotation", help="GTF or GFF file for gene region extraction")
    return parser.parse_args()

def compute_coverage_stats(depths, total_bases):
    '''Compute basic and percentile coverage statistics. store the results in two dictionary,
in reality, this function is not good for large region compute, need to refine: best to store these values in a dict when generate them, do not use depth in a array, which comsume a lot of memory'''
    depths = np.array(depths)
    avg_depth = np.mean(depths)
    med_depth = np.median(depths)
    zero_cov = np.sum(depths == 0)
    over_1 = np.sum(depths >= 1)
    over_1_less_5 = np.sum((depths >=1 ) & (depths < 5))
    over_5 = np.sum(depths >= 5)
    over_5_less_10 = np.sum((depths >=5 ) & (depths < 10))
    over_10 = np.sum(depths >= 10)
    over_10_less_30 = np.sum((depths >=10 ) & (depths < 30))
    over_30 = np.sum(depths >= 30)
    over_30_less_50 = np.sum((depths >=30 ) & (depths < 50))
    over_50 = np.sum(depths >= 50)
    over_50_less_100 = np.sum((depths >=50 ) & (depths < 100))
    over_100 = np.sum(depths >= 100)
    over_100_less_200 = np.sum((depths >=100 ) & (depths < 200))
    over_200 = np.sum(depths >= 200)
    max_depth = np.max(depths)
    p10 = np.percentile(depths, 10)
    p25 = np.percentile(depths, 25)
    p50 = np.percentile(depths, 50)
    p75 = np.percentile(depths, 75)
    p90 = np.percentile(depths, 90)
    piechart_value = { 
    "0X": zero_cov,
    "[1,5X)": over_1_less_5,
    "[5,10X)": over_5_less_10,
    "[10,30X)": over_10_less_30,
    "[30,50X)": over_30_less_50,
    "[50,100X)": over_50_less_100,
    "[100,200X)": over_100_less_200,
    ">=200X": over_200
    }
    coverage_report = {
        "Average Depth": avg_depth,
        "Median Depth": med_depth,
        "Max Depth": max_depth,
        ">=1X": str(over_1) + "\t" + str(round(over_1/total_bases*100,2)) + "%",
        ">=5X": str(over_5) + "\t" + str(round(over_5/total_bases*100,2)) + "%",
        ">=10X": str(over_10) + "\t" + str(round(over_10/total_bases*100,2)) + "%",
        ">=30X": str(over_30) + "\t" + str(round(over_30/total_bases*100,2)) + "%",
        ">=50X": str(over_50) + "\t" + str(round(over_50/total_bases*100,2)) + "%",
        ">=100X": str(over_100) + "\t" + str(round(over_100/total_bases*100,2)) + "%",
        ">=200X": str(over_200) + "\t" + str(round(over_200/total_bases*100,2)) + "%",
        "0X": str(zero_cov) + "\t" + str(round(zero_cov/total_bases*100,2)) + "%",
        "[1,5X)": str(over_1_less_5) + "\t" + str(round(over_1_less_5/total_bases*100,2)) + "%",
        "[5,10X)": str(over_5_less_10) + "\t" + str(round(over_5_less_10/total_bases*100,2)) + "%",
        "[10,30X)": str(over_10_less_30) + "\t" + str(round(over_10_less_30/total_bases*100,2)) + "%",
        "[30,50X)": str(over_30_less_50) + "\t" + str(round(over_30_less_50/total_bases*100,2)) + "%",
        "[50,100X)": str(over_50_less_100) + "\t" + str(round(over_50_less_100/total_bases*100,2)) + "%",
        "[100,200X)": str(over_100_less_200) + "\t" + str(round(over_100_less_200/total_bases*100,2)) + "%",
        "10th Percentile": p10,
        "25th Percentile": p25,
        "50th Percentile": p50,
        "75th Percentile": p75,
        "90th Percentile": p90,
        "Total Bases computed": total_bases
    }
    return piechart_value,coverage_report 

def write_coverage_report(stats, output_file):
    '''write coverage report to output file'''
    with open(output_file, "w") as f:
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")

def calculate_depth(bam_file, coverage_file, annotation=None, region=None, bed=None, gene=None, quality=0):
    bam = pysam.AlignmentFile(bam_file, "rb") # chagne rc to rb if it bam file
    out = open(coverage_file,'w')
    depths = []
    total_bases = 0
    if region or gene:
        tmp_region=None
        if gene:
            tmp_region=parse_gtf_for_gene(annotation,gene)
        elif region:
            tmp_region=region
        chrom, coords = tmp_region.split(":")
        start, end = map(int, coords.split("-"))
        for pileupcolumn in bam.pileup(chrom, start, end, truncate=True,min_mapping_quality = quality):
            depths.append(pileupcolumn.nsegments)
            out.write("\t".join(map(str,[pileupcolumn.reference_name,pileupcolumn.reference_pos,pileupcolumn.nsegments]))+"\n")
        total_bases = end - start
    elif bed:
        bed_regions = BedTool(bed)
        for interval in bed_regions:
            for pileupcolumn in bam.pileup(interval.chrom, interval.start, interval.end, truncate=True, min_mapping_quality = quality):
                depths.append(pileupcolumn.nsegments)
                out.write("\t".join(map(str,[pileupcolumn.reference_name,pileupcolumn.reference_pos,pileupcolumn.nsegments]))+"\n")
        total_bases = sum([i.length for i in bed_regions])
    else:
        for pileupcolumn in bam.pileup(min_mapping_quality = quality):
            depths.append(pileupcolumn.nsegments)
            out.write("\t".join(map(str,[pileupcolumn.reference_name,pileupcolumn.reference_pos,pileupcolumn.nsegments]))+"\n")
        total_bases = sum(bam.lengths)
    bam.close()
    depths.extend([0] * (total_bases - len(depths)))  ## other position will be zero coverage, add them to depths
    piechart_value, coverage_report = compute_coverage_stats(depths, total_bases)
    return piechart_value, coverage_report

def transform_depth_file_to_bedgraph(coverage_file,bedgraph_file):
    '''convert the depth format (chr pos depth) to bedgraph format (chr start end depth)'''
    inputfile=open(coverage_file,'r')
    outfile=open(bedgraph_file,'w')
    prev_chr, prev_pos, prev_depth = None, None, None
    for line in iter(inputfile):
        chrom, pos_str, depth_str = line.strip().split()
        pos=int(pos_str)
        if chrom == prev_chr and depth_str == prev_depth and pos == prev_pos + 1:
            prev_pos = pos # extend the current region
        else:
            #write privious region
            if prev_chr is not None:
                outfile.write("\t".join([prev_chr,str(start_pos),str(prev_pos),prev_depth])+"\n")
            # start new a region 
            prev_chr, prev_pos, prev_depth,start_pos = chrom, pos, depth_str, pos - 1
    # process the last regioin
    if prev_chr is not None:
        outfile.write("\t".join([prev_chr,str(start_pos),str(prev_pos + 1),prev_depth])+"\n")

def convert_to_bigwig(bedgraph_file, genome_file, output_file):
    cmd = ["bedGraphToBigWig", bedgraph_file, genome_file, output_file]
    subprocess.run(cmd, check=True)

def parse_gtf_for_gene(gtf_file, gene_name):
    ''' if a gene is provided, read the gtf/gff file and get the region like chr:TSS-TES to calculate''' 
    tss, tes, chrom = None, None, None
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attr_field = fields[8]
            if f'gene_name "{gene_name}"' in attr_field or f'gene_id "{gene_name}"' in attr_field or f'gene_id={gene_name}' or f'gene_name={gene_name}':
                chrom = fields[0]
                start = int(fields[3]) - 1
                end = int(fields[4])
                strand = fields[6]
                return f"{chrom}:{start}-{end}"
                break
    return None

def draw_pyGenomeTracks(bigwig_file, region, output_pdf):
    config = f"""
[track bigwig]
file = {bigwig_file}
title = Coverage
height = 4
color = black
[x-axis]
fontsize=10
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".ini", delete=False) as tmp:
        tmp.write(config)
        tmp_path = tmp.name
    cmd = [
        "pyGenomeTracks",
        "--tracks", tmp_path,
        "--region", region,
        "--outFileName", output_pdf]
    subprocess.run(cmd, check=True)
    os.remove(tmp_path)

def draw_piechart(labels,sizes,piechart_file):
    total = sum(sizes)
    fig, ax = plt.subplots()
    def autopct_func(pct):
        return f'{pct:.1f}%' if pct > 5 else ''
    wedges, texts, autotexts = ax.pie(
        sizes, 
        labels=[label if size / total > 0.05 else '' for label, size in zip(labels, sizes)], 
        startangle=140, 
        autopct=autopct_func)
    ax.legend(wedges, 
        labels,
        title="Coverage",
        loc="best",
        bbox_to_anchor=(1, 0.5))
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(piechart_file)

def main():
    if len(sys.argv)<2:
        print('''too few options, type: 
        python3 bam_coverage_report.py -h 
        to see help info''')
        sys.exit(-1)
    else:
        (args,options)=parseOptions()

        region = args.region
        quality = str(args.quality)
        bed = args.bed
        annotation= args.annotation
        gene = args.gene
        draw = args.draw
        track = args.track
        genome_file = args.Genomefile
        input_bam = args.input
        output_prefix = args.output
        picture_format=args.picture_format
        if picture_format == "PDF":
            picture_format = "pdf"
        if gene and not annotation:
            print("A gene is given, -a (the annotation) must be provided")
            sys.exit(1)
    
        # Step 1: Calculate coverage and draw a simple pie chart
        if not os.path.exists(output_prefix):
            os.mkdir(output_prefix)
        coverage_file= f"{output_prefix}/{output_prefix}.coverage.txt"
        piechart_value,stats = calculate_depth(input_bam, coverage_file, annotation=annotation, region=region, bed=bed, gene=gene)
        report_file = f"{output_prefix}/{output_prefix}.coverage.report.txt"
        write_coverage_report(stats, report_file)

        # draw a piechart to show how deep the sequencing is
        piechart_file = f"{output_prefix}/{output_prefix}.coverage.piechart.pdf"
        if picture_format == "PNG" or picture_format =="png":
            piechart_file = f"{output_prefix}/{output_prefix}.coverage.piechart.png"
        labels=piechart_value.keys()
        sizes=piechart_value.values()
        draw_piechart(labels,sizes,piechart_file)
    
        # Step 2: Generate coverage track
        # convert the depth file to bedgraph format
        if args.track_file_format:
            # if a file format given, automatically give bigwig/bedgraph file
            track = True
        if track:
            bedgraph_file = f"{output_prefix}/{output_prefix}.bedgraph"
            transform_depth_file_to_bedgraph(coverage_file,bedgraph_file)

#        # alternative way to generate bedgraph file using pybedtools
#        bam_bedtool = BedTool(input_bam)
#        if track:
#            # if you only want to get the bed graph file, use this argument
#            coverage = bam_bedtool.genome_coverage(bg=True)
#            bedgraph_file = f"{output_prefix}/{output_prefix}.bedgraph"
#            coverage.saveas(bedgraph_file)
#    
            if args.track_file_format == "bigwig":
                if not genome_file:
                    print("to generate a bigwig file, genome size file is required, use -G paramete")
                    sys.exit(1)
                bigwig_file = f"{output_prefix}/{output_prefix}.bw"
                convert_to_bigwig(bedgraph_file, genome_file, bigwig_file)
            #os.remove(bedgraph_file)
    
        # Step 3: draw a track picture for the given region or gene
        if draw :
            bigwig_input = bigwig_file if args.track_file_format == "bigwig" else bedgraph_file
            output_pdf = f"{output_prefix}/{output_prefix}.coverage.pdf"
            if picture_format == "PNG" or picture_format == "png": 
                output_pdf = f"{output_prefix}/{output_prefix}.coverage.png"
            if region:
                draw_pyGenomeTracks(bigwig_input, region, output_pdf)
            if bed:
                bed_regions = BedTool(bed)
                for interval in bed_regions:
                    tmp_region=interval.chrom+":"+str(interval.start)+"-"+str(interval.end)
                    tmp_region_file=f"{output_prefix}/{output_prefix}.{interval.chrom}_{interval.start}_{interval.end}.pdf"
                    if picture_format == "PNG" or picture_format == "png":
                        tmp_region_file=f"{output_prefix}/{output_prefix}.{interval.chrom}_{interval.start}_{interval.end}.png"
                    draw_pyGenomeTracks(bigwig_input, tmp_region, tmp_region_file)
    
            if gene:
                if not annotation: # quit if annotation file is not provided
                    print("Annotation file required for gene plotting. Use -a to specify GTF/GFF. or you can just provide the region of that gene to draw this plot")
                    sys.exit(1)
        
                region_from_gene = parse_gtf_for_gene(annotation, gene)
                if not region_from_gene:
                    print(f"Gene {gene} not found in annotation.")
                    sys.exit(1)
        
                bigwig_input = bigwig_file if args.track_file_format == "bigwig" else bedgraph_file
                output_pdf = f"{output_prefix}/{output_prefix}.coverage.pdf"
                if picture_format == "PNG" or picture_format == "png": 
                    output_pdf = f"{output_prefix}/{output_prefix}.coverage.png"
                draw_pyGenomeTracks(bigwig_input, region_from_gene, output_pdf)

if __name__ == "__main__":
    main()

