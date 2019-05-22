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
BAMFILE, = glob_wildcards("bam/{cell}.bam")
 #                                                                             #
 #                                                                             #
abbreviate_names = False
 #                                                                             #
 ###############################################################################

configfile: "Snake.config.json"

 #
 # PART I
 # Read preprocessing
 #
	
rule remove_low_quality_reads:
    input:
        bam = "bam/{cell}.bam"
	awk = awk_1st.awk
    output:
        temp("bam/{cell}.sc_pre_mono.bam")
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
	samtools view -H {input} > header_test.sam
	samtools view -F 2304 {input.bam} | awk -f {input.awk} | cat header_test.sam - | samtools view -Sb - > {output}	
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
        bam = "bam/{cell}.sc_pre_mono_sort_for_mark.bam"
	bai = "bam/{cell}.sc_pre_mono_sort_for_mark.bam.bai"
    output:
        bam_uniq = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam"
	metrix = "bam/{cell}.sc_pre_mono.metrix_dup.txt"
    shell:
        """
        module load biobambam2/2.0.76-foss-2016b
        bammarkduplicates markthreads=2 I={input.bam} O={output.bam_uniq} M={output.metrix} index=1 rmdup=1
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
 # In this step, users can choose the set of regulatory elements they want to analyze
 
#rule count_reads:
#   input:
#       counts = "counts/{SAMPLE_NAME}.{window}.txt.gz"
#   output:
#       dynamic("plot_chromosomes/{SAMPLE_NAME}/window_{window}.{chrom}.pdf")
#   params:
#       sv_plot = config["sv_plot"]
#   shell:
#       """
#       module load deeptools/2.5.1-foss-2016b-Python-2.7.12
#	perl Strand_seq_deeptool_DHS_chromVAR.pl > Strand_seq_deeptool_bin_auto.pl
#        perl Strand_seq_deeptool_bin_auto.pl
#        """

