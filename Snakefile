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
 # 'bam' and to be names *.bam and *.bai                                       #
 #                                                                             #
BAMFILE, = glob_wildcards("bam/{cell}.bam")
 #                                                                             #
 #                                                                             #
abbreviate_names = False
 #                                                                             #
 ###############################################################################

configfile: "Snake.config.json"

rule all:
    input:
        expand("result/{s}.tab", s = SAMPLE_NAME),
        expand("result/{s}.npz", s = SAMPLE_NAME)


 #
 # PART I
 # Read preprocessing
 #
	
rule remove_low_quality_reads:
    input:
        bam = "bam/{cell}.bam"
    output:
        bam_pre = "bam/{cell}.sc_pre_mono.bam",
        bam_header = "bam/{cell}.header_test.sam"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
		samtools view -H {input} > {output.bam_header} 
		samtools view -F 2304 {input.bam} | awk -f utils/awk_1st.awk | cat {output.bam_header} - | samtools view -Sb - > {output.bam_pre}	
        """

rule sort_bam:
    input:
        "bam/{cell}.sc_pre_mono.bam"
    output:
        "bam/{cell}.sc_pre_mono_sort_for_mark.bam"
    threads:
        2
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools sort -@ {threads} -O BAM -o {output} {input}
        """

rule index_num1:
    input:
        "bam/{cell}.sc_pre_mono_sort_for_mark.bam"
    output:
        "bam/{cell}.sc_pre_mono_sort_for_mark.bam.bai"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools index {input}
        """	
	
rule remove_dup:
    input:
        bam="bam/{cell}.sc_pre_mono_sort_for_mark.bam"
    output:
        bam_uniq="bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bam_metrix="bam/{cell}.sc_pre_mono.metrix_dup.txt"
    shell:
        """
        module load biobambam2/2.0.76-foss-2016b
		bammarkduplicates markthreads=2 I={input.bam} O={output.bam_uniq} M={output.bam_metrix} index=1 rmdup=1
        """

rule index_num2:
    input:
        "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam"
    output:
        "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools index {input}
        """

 #
 # PART II
 # Read counting
 #

rule count_reads:
    input:
        bam = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        bai = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai", cell=BAMFILE)
    output:
        tab = "result/" + SAMPLE_NAME + ".tab",
        npz = "result/" + SAMPLE_NAME + ".npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/regions_all_hg38_v2_resize_2kb_sort.bed --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """

