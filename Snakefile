import os
###############################################################################
 # SETTINGS                                                                    #
 #                                                                             #
 # Set the sample name                                                         #
 #                                                                             #
SAMPLE_NAME = "testSample"
 #                                                                             #
 #                                                                             #
 # By default, the pipeline will expect file to be in a subfolder called       #
 # 'fastq' and to be names *_1_sequence.txt.gz *_2_sequence.txt.gz             #
 #                                                                             #
FASTQ, = glob_wildcards("fastq/{cell}_1_sequence.txt.gz")
 #                                                                             #
 #                                                                             #
abbreviate_names = False
 #                                                                             #
 ###############################################################################

configfile: "Snake.config.json"

 ###############################################################################
 # PRE-PROCESSING                                                              #
 #                                                                             #
 # Remove common substrs from front and end of fastq names to determine NAMES  #
 #                                                                             #
def longest_pref(strs):
    if len(strs) == 1:
        return strs[0]
    ref = strs[0]
    lpos = len(ref)
    for i in range(1,len(strs)):
        pos = 0
        while pos<len(ref) and pos<len(strs[i]) and strs[i][pos]==ref[pos]:
            pos+=1
        if pos < lpos:
            lpos = pos
        if lpos == 0:
            break
    return lpos
pref = longest_pref(FASTQ)
suff = longest_pref([x[::-1] for x in FASTQ])
NAMES = dict([(f,f) for f in FASTQ])
if abbreviate_names:
    NAMES = dict([(fq, fq[pref:(len(fq)-suff)]) for fq in FASTQ])
    print("NOTE")
    print("I will try to uniquely abbreviate cell names. Here is an example:")
    for cell, name in NAMES.items():
        print("  ", cell, " -> ", name)
        break
 #                                                                             #
 # Finished pre-processing                                                     #
 #                                                                             #
 ###############################################################################




rule all:
    input:
        expand("plot_overview/{s}.{window}.pdf", s = SAMPLE_NAME, window = [20000, 50000, 100000, 200000, 500000]),
        expand("plot_chromosomes/{s}/window_{window}.chr1.pdf", s = SAMPLE_NAME, window = [20000, 50000, 100000, 200000, 500000])
		# dynamic(expand("plot_chromosomes/{s}/window_{window}.{{chrom}}.pdf", s = SAMPLE_NAME, window = [20000, 50000, 100000, 200000, 500000]))
        # Note: dynamic rule does not work !


 #
 # PART I
 # Read mapping
 #

def get_time_for_map_reads(wc, input):
    fsize = os.stat(input[1]).st_size
    return max(1,int(round(fsize/400000000)))

rule map_reads:
    input:
        "fastq/{cell}_1_sequence.txt.gz", "fastq/{cell}_2_sequence.txt.gz"
    output:
        temp("bam/{cell}.bam")
    params:
        name = lambda wc: NAMES[wc.cell],
        genome = config["genome"]
    threads:
        4
    resources:
        time = get_time_for_map_reads
    message:
        "Start mapping cell {input} using {threads} threads and giving {resources.time} hours"
    shell:
        """
        module load BWA/0.7.15-foss-2016b SAMtools/1.3.1-foss-2016b
        bwa mem -t {threads} \
            -R '@RG\tID:{params.name}\tSM:{SAMPLE_NAME}\tLB:{wildcards.cell}' \
            {params.genome} \
            {input} \
        | samtools view -bT {params.genome} - \
        > {output}
        """

rule sort_bam:
    input:
        "bam/{cell}.bam"
    output:
        temp("bam/{cell}.sort.bam")
    threads:
        2
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools sort -@ {threads} -O BAM -o {output} {input}
        """

rule markdups:
    input:
        "bam/{cell}.sort.bam"
    output:
        "bam/{cell}.sort.mdup.bam"
    threads:
        2
    shell:
        """
        module load biobambam2/2.0.76-foss-2016b
        bammarkduplicates markthreads={threads} I={input} O={output} index=1 rmdup=0
        """

rule index:
    input:
        "bam/{cell}.sort.mdup.bam"
    output:
        "bam/{cell}.sort.mdup.bam.bai"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools index {input}
        """





 #
 # PART II
 # MosaiCatcher counts & plots
 #

rule mosaic_count:
    input:
        bam = expand("bam/{cell}.sort.mdup.bam", cell = FASTQ),
        bai = expand("bam/{cell}.sort.mdup.bam.bai", cell = FASTQ)
    output:
        counts = "counts/" + SAMPLE_NAME + ".{window}.txt.gz",
        info   = "counts/" + SAMPLE_NAME + ".{window}.info"
    params:
        mosaic = config["mosaic"],
        exclude = config["exclude"]
    shell:
        """
        {params.mosaic} count \
            -o {output.counts} \
            -i {output.info} \
            -x {params.exclude} \
            -w {wildcards.window} \
            {input.bam}
        """

rule plot_mosaic_counts:
    input:
        counts = "counts/" + SAMPLE_NAME + ".{window}.txt.gz",
        info   = "counts/" + SAMPLE_NAME + ".{window}.info"
    output:
        "plot_overview/" + SAMPLE_NAME + ".{window}.pdf"
    params:
        qc_plot = config["qc_plot"]
    shell:
        """
        module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0
        Rscript {params.qc_plot} {input.counts} {input.info} {output}
        """

rule plot_chromosomes:
    input:
        counts = "counts/{SAMPLE_NAME}.{window}.txt.gz"
    output:
        dynamic("plot_chromosomes/{SAMPLE_NAME}/window_{window}.{chrom}.pdf")
    params:
        sv_plot = config["sv_plot"]
    shell:
        """
        module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0
        Rscript {params.sv_plot} {input} plot_chromosomes/{SAMPLE_NAME}/window_{wildcards.window}
        """

