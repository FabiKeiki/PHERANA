configfile: "../config/config.yaml"

rule all:
    input:
        "../results/data_validation/motus_output/Summary/motus_combined.txt"


#################### Main Functions #######################################
import os
def get_direction_r1(wildcards) :
#This  function returns a list of R1 reads when the same sample was sequenced on multiple lanes
    samname=wildcards.sample
    l=config["samples"][samname]

    r_list=[]
    r1_l1=[s for s in l if "L1_R1" in s]
    r1_l2=[s for s in l if "L2_R1" in s]

    r_list.extend(r1_l1)
    r_list.extend(r1_l2)
    return(r_list)


def get_direction_r2(wildcards) :
#This  function returns a list of R2 reads when the same sample was sequenced on multiple lanes
    samname=wildcards.sample
    l=config["samples"][samname]

    r_list=[]
    r2_l1=[s for s in l if "L1_R2" in s]
    r2_l2=[s for s in l if "L2_R2" in s]

    r_list.extend(r2_l1)
    r_list.extend(r2_l2)
    return(r_list)

def get_files_commas(path, sep=","):
    file_l=os.listdir(path)
    file_l2=[]
    for f in file_l:
        file_l2.append(os.path.join(path,f))
    out=sep.join(file_l2)
    return(out)

def get_wildcard_commas(paths, sep=","):
    file_l=[paths]
    out=sep.join(file_l)
    return(out)

def get_type(wildcards) :
    type_w=wildcards.type
    return(type_w)

def join_list(li, car):
    fin_l=[]
    for i in li:
        add=i+car
        fin_l.append(add)
    return(fin_l)


################################################################################
########################### DATA VALIDATION ####################################
################################################################################
# The following pipeline is used to ascertain that the metagenomics samples
# are composed of what we expect (viruses in the virome fraction, bacteria in the
# bacteriome one)
################################################################################
################################################################################


################################### Concat Lanes ###################################


#Concatenate raw reads from the same sample that were sequenced on different lanes
rule concat_lanes:
    input:
        R1s=get_direction_r1,
        R2s=get_direction_r2
    output:
        R1_concat="../scratch_link/concat_raw_reads/{sample}_R1_concat.fastq.gz", # these files go in scratch because they can be quickly recreated if needed
        R2_concat="../scratch_link/concat_raw_reads/{sample}_R2_concat.fastq.gz"
    threads: 2
    log:
        "logs/data_validation/lane_concatenation/{sample}_concat.log"
    resources:
        account = "pengel_beemicrophage",
        runtime= "02:00:00"
    shell:
        "cat {input.R1s} > {output.R1_concat}; cat {input.R2s} > {output.R2_concat}; "


################################### Kraken2 ###################################

# This rule takes the concatenated raw reads files and runs them against the krakend db
rule run_Krakern2:
    input:
        R1="../scratch_link/concat_raw_reads/{sample}_R1_concat.fastq.gz",
        R2="../scratch_link/concat_raw_reads/{sample}_R2_concat.fastq.gz",
        db="../../../mndiaye1/PHOSTER/workflow/resources/databases/220131_costum_kraken2db"
    output:
        tab=temp("../results/data_validation/kraken2_output/{sample}_kraken2_report.kraken"),
        rep="../results/data_validation/kraken2_output/Reports/{sample}_kraken2_report"
    conda:
        "envs/Kraken2.yaml"
    threads: 24
    log:
        "logs/data_validation/kraken2/run/{sample}_kraken2.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "03:00:00"
    shell:
        "kraken2 --use-names --threads {threads} \
         --db {input.db} \
         --fastq-input --report {output.rep}  --gzip-compressed \
         --paired {input.R1} {input.R2} \
         > {output.tab}"
    
# This rule parses the kraken output for further analyes
rule parse_kraken_report:
    input:
        "scripts/data_validation/parse_kraken_report.py", # if you change the script, the rule runs again
        expand("../results/data_validation/kraken2_output/Reports/{sample}_kraken2_report", sample=config["samples"])
    output:
        "../results/data_validation/kraken2_output/Summary/all_samples_report.txt"
    threads: 2
    params:
        "../results/data_validation/kraken2_output/"
    log:
        "logs/data_validation/kraken2/parsing/report_parser_kraken2.log"
    resources:
        account = "pengel_beemicrophage",
        runtime= "00:15:00"
    script:
        "scripts/data_validation/parse_kraken_report.py"



############################# QC and Trimming ##################################

# This rule does a fastQC on the raw reads
rule fastQC_PreTrimming:
    input:
        R1=get_direction_r1,
        R2=get_direction_r2
    output:
        temp(directory("../results/data_validation/QC/preTrimming/QC_{sample}/"))
    threads: 2
    log:
        "logs/data_validation/QC/{sample}_QC.log"
    conda:
        "envs/fastqc.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 6000,
        runtime= "01:00:00"
    shell:
        "mkdir -p {output}; "
        "fastqc -o {output} {input.R1};"
        "fastqc -o {output} {input.R2};"


# This rule runs the trimming of the raw reads
rule rawreads_trimming:
    input:
        R1="../scratch_link/concat_raw_reads/{sample}_R1_concat.fastq.gz",
        R2="../scratch_link/concat_raw_reads/{sample}_R2_concat.fastq.gz"
    output:
        R1_paired="../data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        R2_paired="../data/trimmed_reads/{sample}_R2_paired.fastq.gz",
        R1_unpaired="../data/trimmed_reads/{sample}_R1_unpaired.fastq.gz",
        R2_unpaired="../data/trimmed_reads/{sample}_R2_unpaired.fastq.gz"
    threads: 8
    params:
        nextera="../data/reference_assemblies/short_RefSeqs/NexteraPE-PE.fa",
        q=28,
        min_length=40
    log:
        "logs/data_validation/trimming/{sample}_trimming.log"
    conda:
        "envs/trimmomatic.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 6000,
        runtime= "03:00:00"
    shell:
        "trimmomatic PE -phred33 -threads {threads} {input.R1} {input.R2} \
         {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
         ILLUMINACLIP:{params.nextera}:2:30:10 \
         LEADING:{params.q} TRAILING:{params.q} MINLEN:{params.min_length}"

# this rule does a post trimming fast QC
rule fastQC_PostTrimming:
    input:
        R1="../data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        R2="../data/trimmed_reads/{sample}_R2_paired.fastq.gz"
    output:
        temp(directory("../results/data_validation/QC/postTrimming/QC_{sample}/"))
    threads: 2
    log:
        "logs/data_validation/QC/{sample}_postQC.log"
    conda:
        "envs/fastqc.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 6000,
        runtime= "01:00:00"
    shell:
        "mkdir -p {output}; "
        "fastqc -o {output} {input.R1} {input.R2}"

rule parse_fastQC:
    input:
        preT=expand("../results/data_validation/QC/preTrimming/QC_{sample}/", sample=config["samples"]),
        postT=expand("../results/data_validation/QC/postTrimming/QC_{sample}/", sample=config["samples"])
    output:
        "../results/data_validation/QC/Summary/fastQC_summary.txt"
    threads: 2
    log:
        "logs/data_validation/QC/summarize_fastQC.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 6000,
        runtime= "00:20:00"
    script:
        "scripts/data_validation/parse_fastqc_output.py"


############################### Host Filtering #################################

# This rule uses bbsplit to map remove reads from the honeybee genome and/or the human genomes
rule host_filtering:
    input:
        R1="../data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        R2="../data/trimmed_reads/{sample}_R2_paired.fastq.gz"
    output:
        dir=temp(directory("../results/data_validation/host_filtering/discarded/{sample}_discarded/")), # I don't need the sam of the mapping so I delete them immediatly
        unmapped_R1="../data/host_filtered_reads/{sample}_R1_HF.fastq.gz", # these are the filtered reads
        unmapped_R2="../data/host_filtered_reads/{sample}_R2_HF.fastq.gz",
        refstats="../results/data_validation/host_filtering/HF_mappings_stats/{sample}_refstats.out"
    conda:
        "envs/bwa_mapping.yaml"
    threads: 25
    params:
        ref_Acer="../data/reference_assemblies/A_cerana/GCF_001442555.1_ACSNU-2.0_genomic.fna.gz",
        ref_Hsap="../data/reference_assemblies/H_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        xmx="50g"
    log:
        "logs/data_validation/HF/{sample}_HF.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "08:00:00"
    shell:
        "bbsplit.sh in1={input.R1} in2={input.R2} ref={params.ref_Acer},{params.ref_Hsap} \
        basename={output.dir}/{wildcards.sample}_HF_discarded_%.sam \
        refstats={output.refstats} rebuild=t nodisk=t \
        outu1={output.unmapped_R1} outu2={output.unmapped_R2} nzo=f -Xmx{params.xmx} threads={threads}"


# This rule parses the refstats output of the filering
rule parse_filtering_refstats:
    input:
        files=expand("../results/data_validation/host_filtering/HF_mappings_stats/{sample}_refstats.out", sample=config["samples"])
    output:
        "../results/data_validation/host_filtering/HF_mappings_stats/HF_refstats.txt"
    threads: 2
    log:
        "logs/data_validation/HF/HF_refstats_parsing.log"
    params:
        file1="../results/data_validation/host_filtering/HF_mappings_stats/file1.txt",
        file2="../results/data_validation/host_filtering/HF_mappings_stats/file2.txt",
        tmp="../results/data_validation/host_filtering/HF_mappings_stats/tmp.txt",
        sams=expand("{sample}", sample=config["samples"])
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500,
        runtime= "00:30:00"
    shell:
        "echo -e 'sample\tname\tperc_unambiguousReads\tunambiguousMB\tperc_ambiguousReads\tambiguousMB\tunambiguousReads\tambiguousReads\tassignedReads\tassignedBases' > {output}; "
        "tail -n +2 -q {input.files} > {params.file1}; "
        "printf '%s\n' {params.sams} > {params.file2}; "
        "awk '{{for(i=0;i<2;i++)print}}' {params.file2} > {params.tmp}; "
        "paste -d '\t' {params.tmp} {params.file1} >> {output}"


############################### count_reads ##################################
# After trimming and host filtering the reads, I wanna know how much I lost in
# Terms of reads and bases


# this rules returns a table of read count before and after trimming
rule count_reads_qc:
    input:
        preT="../scratch_link/concat_raw_reads/{sample}_R1_concat.fastq.gz",
        postT="../data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        postF="../data/host_filtered_reads/{sample}_R1_HF.fastq.gz"
    output:
        temp("../results/data_validation/QC/{sample}_read_count.txt")
    log:
        "logs/data_validation/QC/read_count/{sample}_read_count.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 4000,
        runtime= "00:30:00"
    shell:
        "touch {output}; "
        "(./scripts/data_validation/count_reads.sh {input.preT} {input.postT} {input.postF} {output})2>{log}"


rule count_reads_summary:
    input:
        sams=expand("../results/data_validation/QC/{sample}_read_count.txt", sample=config["samples"])
    output:
        "../results/data_validation/QC/Summary/read_count.txt"
    log:
        "logs/data_validation/QC/read_count/summary_read_count.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 4000,
        runtime= "00:30:00"
    shell:
        "(awk 'FNR!=NR && FNR==1 {{next}} 1' {input.sams} > {output})2>{log}"



################################################################################
###############################  mOTUS  ########################################
################################################################################
# mOTUs is more has more taxonomical resolution than kraken, so once the reads
# are all cleaned and filtered, I run mOTUs to obtain a genus-level composition
# of the remaining reads (for the bacterial samples)
################################################################################
################################################################################

# this rule runs the motu profiling
# it is better to launch this rule with one sample first and then all the others,
# mOTUs db download doesn't handle weell multiple files trying to download it and
# access it at the same time and the jobs might fail.
rule run_motus:
    input:
        reads1 = "../data/host_filtered_reads/{sample}_R1_HF.fastq.gz",
        reads2 = "../data/host_filtered_reads/{sample}_R2_HF.fastq.gz"
    output:
        motus_temp = temp("../results/data_validation/motus_output/map/{sample}_map.motus")
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "02:30:00"
    threads: 24
    log:
        "logs/data_validation/motus/{sample}_motus.log"
    conda:
        "envs/motus-env.yaml"
    shell:
        "motus downloadDB; " # just leave it here to be safe - it warns and continues
        "motus profile -f {input.reads1} -r {input.reads2} -n {wildcards.sample} -o {output.motus_temp} -t {threads}"

# this roule run the mouts counting of the marker genes that map to the datbase
rule run_motus_count:
    input:
        reads1 = "../data/host_filtered_reads/{sample}_R1_HF.fastq.gz",
        reads2 = "../data/host_filtered_reads/{sample}_R2_HF.fastq.gz"
    output:
        motus_temp = temp("../results/data_validation/motus_output/count/{sample}_count.motus")
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "02:30:00"
    threads: 24
    log:
        "logs/data_validation/motus/{sample}_motus.log"
    conda:
        "envs/motus-env.yaml"
    shell:
        "motus downloadDB; " # just leave it here to be safe - it warns and continues
        "motus profile -f {input.reads1} -r {input.reads2} -c -n {wildcards.sample} -o {output.motus_temp} -t {threads}"

rule merge_motus:
    input:
        motus_temp = expand("../results/data_validation/motus_output/map/{sample}_map.motus", sample=config["samples"]),
        motus_count= expand("../results/data_validation/motus_output/count/{sample}_count.motus", sample=config["samples"])
    output:
        motus_merged = "../results/data_validation/motus_output/Summary/samples_merged_map.motus",
        mouts_merged_count = "../results/data_validation/motus_output/Summary/samples_merged_count.motus"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 500000,
        runtime= "02:00:00"
    threads: 12
    log:
        "logs/data_validation/motus/merge_motus.log"
    conda:
        "envs/motus-env.yaml"
    shell:
        "motus merge -a bee -i $(echo \"{input.motus_temp}\" | sed -e 's/ /,/g' ) > {output.motus_merged}; "
        "motus merge -a bee -c -i $(echo \"{input.motus_count}\" | sed -e 's/ /,/g' ) > {output.mouts_merged_count}"

rule parse_motus:
    input:
        motus_tab = "../results/data_validation/motus_output/Summary/samples_merged_map.motus",
        motus_count = "../results/data_validation/motus_output/Summary/samples_merged_count.motus"
    output:
        "../results/data_validation/motus_output/Summary/motus_combined.txt"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 50000,
        runtime= "03:30:00"
    log:
        "logs/data_validation/motus/parse_motus.log"
    conda:
        "envs/base_R_env.yaml"
    script:
        "scripts/data_validation/parse_motus.R"



