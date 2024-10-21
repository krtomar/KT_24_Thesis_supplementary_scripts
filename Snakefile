configfile: "config.yaml"

import pandas as pd
import os, glob


sample_table =  pd.read_table("SampleTable.txt")
samples = list(sample_table['LibraryID'])

sample_table['Pairs'] = sample_table['Fraction'].str.replace(' |fraction', '') + "_" + sample_table['Time'].astype(str).str.zfill(4) + "_" + sample_table['Dataset']
pairs = list(sample_table['Pairs'])

sample_table['TimePoints'] = sample_table['Fraction'].str.replace(' |fraction', '') + "_" + sample_table['Genotype'] + "_" + sample_table['Dataset']
timepoints = list(sample_table['TimePoints'])


# everything that need network connection for
localrules: get_index_files


rule all:
    input:
        expand("fastqc/{ID}_fastqc.txt",              ID=samples),
        expand("Output/BAM/{ID}.ReadsPerGene.out.tab", ID=samples),
        # expand("Output/rsem/{ID}_rsem.genes.results",  ID=samples),
        # expand("Output/bedgraph_frags/{ID}.norm.bedgraph",  ID=samples),
        expand("Output/frag_size/{ID}.pdf",  ID=samples),
        # expand("Output/composite/plots/composite.{ID}.pdf", ID=samples),
        # expand("Output/composite/pairs/composite.{PID}.pdf", PID=pairs),
        # expand("Output/composite2/pairs/composite.{PID}.pdf", PID=pairs),
        # expand("Output/composite/timepoints/composite.{TID}.pdf", TID=timepoints),
        "Output/deseq2_SpomNorm/plots/sessionInfo.txt",
        expand("Output/deseq2_SpomNorm/tables/paramtests.{LT}.{DT}.{FT}.txt",
        LT = ['log2', 'cnt'], DT = ['norm_pullout','lennorm_pullout'], FT = ['expo','logistic']),
        expand("Output/scaledmat/{MID}/{ID}_GB_matrix.rds", ID=samples, MID = ['wmulti','nomulti']),
        expand("Output/scaledcomp/{MID}/composite.{PID}.pdf", PID=pairs, MID = ['wmulti','nomulti']),
        # "Output/deseq2_ScerNorm/plots/sessionInfo.txt",


rule fit_tests:
    input:
        anno="genome/combined.gtf",
        table="SampleTable.txt",
        lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
    params:
        norm="Spom"
    output:
        "Output/deseq2_SpomNorm/tables/paramtests.{LT}.{DT}.{FT}.txt",
    shell:
        """
        Rscript --vanilla scripts/fit_test.R {input} {params.norm} {output}
        """

rule deseq2_SpomNorm:
    input:
        anno="genome/combined.gtf",
        table="SampleTable.txt",
        counts=expand("Output/BAM/{SN}.ReadsPerGene.out.tab", SN=sample_table["LibraryID"]),
    params:
        norm="Spom"
    output:
        lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
        pca="Output/deseq2_SpomNorm/plots/PCA.pdf",
        session="Output/deseq2_SpomNorm/plots/sessionInfo.txt",
    shell:
        """
        Rscript --vanilla scripts/deseq2.R {input} {params.norm} {output}
        """
#
# rule deseq2_ScerNorm:
#     input:
#         anno="genome/combined.gtf",
#         table="SampleTable.txt",
#         counts=expand("Output/BAM/{SN}.ReadsPerGene.out.tab", SN=sample_table["LibraryID"]),
#     params:
#         norm="Scer"
#     output:
#         lnc="Output/deseq2_ScerNorm/tables/log2_norm_counts.txt",
#         pca="Output/deseq2_ScerNorm/plots/PCA.pdf",
#         session="Output/deseq2_ScerNorm/sessionInfo.txt",
#     shell:
#         """
#         Rscript --vanilla scripts/deseq2.R {input} {params.norm} {output}
#         """



# rule RSEM:
#     input:
#         "Output/BAM/{ID}.Aligned.toTranscriptome.out.bam",
#         index="star_index/genome"
#     output:
#         "Output/rsem/{ID}_rsem.genes.results"
#     threads: 12
#     shell:
#         """
#         rsem-calculate-expression --bam \
#         --paired-end --strandedness none \
#         -p {threads} \
#         Output/BAM/{wildcards.ID}.Aligned.toTranscriptome.out.bam {input.index} Output/rsem/{wildcards.ID}_rsem
#         """


rule frag_size:
    input:
        "Output/BAM/{ID}.Aligned.sortedByCoord.out.bam",
    output:
        "Output/frag_size/{ID}.pdf",
    threads: 12
    shell:
        """
        Rscript --vanilla scripts/frag_size.R {input} {output}
        """


### Note: this does not work when reads don´t cover junction ###
# rule bedtools_frags:
#     input:
#         "Output/BAM/{ID}.Aligned.sortedByCoord.out.bam",
#     output:
#         pc="Output/bedgraph_frags/{ID}.pc.bedgraph",
#         unsplit="Output/bedgraph_frags/{ID}.unsplit.bedgraph",
#         split="Output/bedgraph_frags/{ID}.split.bedgraph",
#         substr="Output/bedgraph_frags/{ID}.substr.bedgraph",
#         final="Output/bedgraph_frags/{ID}.norm.bedgraph",
#     threads: 12
#     shell:
#         """
#         samtools index {input}
#
#         TOTAL=`samtools view -c -q 255 {input}`
#         SCALER=`echo "scale=10; 1/(${{TOTAL}}/1000000)" | bc`
#
#         bedtools genomecov -ibam {input} -bga -pc    > {output.pc}
#         bedtools genomecov -ibam {input} -bga        > {output.unsplit}
#         bedtools genomecov -ibam {input} -bga -split > {output.split}
#
#         bedtools unionbedg -i {output.pc} {output.unsplit} {output.split} | awk '{{print $1,$2,$3,$4-$5+$6}}' > {output.substr}
#
#         cat {output.substr} | awk -v SCALER="$SCALER" '{{print $1,$2,$3,$4*SCALER}}' > {output.final}
#         """


rule scaled_comp_pairs:
    input:
        anno="genome/combined.gtf",
        table="SampleTable.txt",
        lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
        mats=lambda wildcards: expand("Output/scaledmat/{MID}/{LID}_GB_matrix.rds", LID = list(sample_table['LibraryID'][sample_table['Pairs'] == wildcards.PID]), MID = wildcards.MID),
    output:
        plot="Output/scaledcomp/{MID}/composite.{PID}.pdf",
    shell:
        """
        set +o pipefail;
        Rscript --vanilla scripts/compscaled_pairs.R {input} {output}
        """


rule bedgraph2scaledmat:
    input:
        anno="genome/combined.gtf",
        ncov="Output/bedgraph/{MID}/{ID}.norm.bedgraph",
    output:
        mat="Output/scaledmat/{MID}/{ID}_GB_matrix.rds",
    shell:
        """
        set +o pipefail;
        Rscript --vanilla scripts/bedgraph2scaledmat.R {input.anno} {input.ncov} {output.mat}
        """


# rule composite_timepoints:
#     input:
#         anno="genome/combined.gtf",
#         table="SampleTable.txt",
#         lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
#         mats=lambda wildcards: expand("Output/matrix/{LID}_TSS_matrix.rds", LID = list(sample_table['LibraryID'][sample_table['TimePoints'] == wildcards.TID])),
#     output:
#         plot="Output/composite/timepoints/composite.{TID}.pdf",
#     shell:
#         """
#         Rscript --vanilla scripts/composite_timepoints.R {input} {output}
#         """


# rule composite_pairs2:
#     input:
#         anno="genome/combined.gtf",
#         table="SampleTable.txt",
#         lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
#         mats1=lambda wildcards: expand("Output/matrix/{LID}_TSS_matrix.rds", LID = list(sample_table['LibraryID'][sample_table['Pairs'] == wildcards.PID])),
#         mats2=lambda wildcards: expand("Output/matrix/{LID}_TTS_matrix.rds", LID = list(sample_table['LibraryID'][sample_table['Pairs'] == wildcards.PID])),
#     output:
#         plot="Output/composite2/pairs/composite.{PID}.pdf",
#     shell:
#         """
#         Rscript --vanilla scripts/composite_pairs2.R {input} {output}
#         """


# rule composite_pairs:
#     input:
#         anno="genome/combined.gtf",
#         table="SampleTable.txt",
#         lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
#         mats=lambda wildcards: expand("Output/matrix/{LID}_TSS_matrix.rds", LID = list(sample_table['LibraryID'][sample_table['Pairs'] == wildcards.PID])),
#     output:
#         plot="Output/composite/pairs/composite.{PID}.pdf",
#     shell:
#         """
#         Rscript --vanilla scripts/composite_pairs.R {input} {output}
#         """


# rule composite:
#     input:
#         anno="genome/combined.gtf",
#         table="SampleTable.txt",
#         lnc="Output/deseq2_SpomNorm/tables/log2_norm_counts.txt",
#         TSS="Output/matrix/{ID}_TSS_matrix.rds",
#         TTS="Output/matrix/{ID}_TTS_matrix.rds",
#     output:
#         plot="Output/composite/plots/composite.{ID}.pdf",
#     shell:
#         """
#         Rscript --vanilla scripts/composite.R {input} {output}
#         """

# rule bedgraph2matrix:
#     input:
#         anno="genome/combined.gtf",
#         ncov="Output/bedgraph/{ID}.norm.bedgraph",
#     output:
#         matTSS="Output/matrix/{ID}_TSS_matrix.rds",
#         matTTS="Output/matrix/{ID}_TTS_matrix.rds",
#     shell:
#         """
#         Rscript --vanilla scripts/bedgraph2matrix.R {input.anno} {input.ncov} {output.matTSS} {output.matTTS}
#         """


rule bedtools_wmulti:
    input:
        "Output/BAM/{ID}.Aligned.sortedByCoord.out.bam",
    output:
        "Output/bedgraph/wmulti/{ID}.norm.bedgraph",
    threads: 12
    shell:
        """
        samtools index {input}

        set +o pipefail;

        TOTAL=`samtools view -c -q 255 {input} spom_I spom_II spom_III`
        SCALER=`echo "scale=10; 1/(${{TOTAL}}/1000000)" | bc`

        genomeCoverageBed -ibam {input} -bg -pc -scale ${{SCALER}} > {output}
        """


rule bedtools_nomulti:
    input:
        "Output/BAM/{ID}.Aligned.sortedByCoord.out.bam",
    output:
        "Output/bedgraph/nomulti/{ID}.norm.bedgraph",
    threads: 12
    shell:
        """
        samtools index {input}

        set +o pipefail;

        TOTAL=`samtools view -c -q 255 {input} spom_I spom_II spom_III`
        SCALER=`echo "scale=10; 1/(${{TOTAL}}/1000000)" | bc`

        samtools view -bS -@ {threads} -q 255 {input} | samtools sort -@ {threads} - | tee Output/BAM/{wildcards.ID}_unique.bam | samtools index - Output/BAM/{wildcards.ID}_unique.bam.bai

        genomeCoverageBed -ibam Output/BAM/{wildcards.ID}_unique.bam -bg -pc -scale ${{SCALER}} > {output}
        """


rule star_align:
    input:
        mate1=lambda wildcards: glob.glob('FastQ/' + wildcards.ID + '*_R1.fastq.gz'),
        mate2=lambda wildcards: glob.glob('FastQ/' + wildcards.ID + '*_R2.fastq.gz'),
        index="star_index/Genome"
    output:
        # "Output/BAM/{ID}.Aligned.toTranscriptome.out.bam",
        "Output/BAM/{ID}.ReadsPerGene.out.tab",
        "Output/BAM/{ID}.Aligned.sortedByCoord.out.bam",
    threads: 12
    params:
        quantMode="GeneCounts",
        outSAMtype="BAM SortedByCoordinate",
        limitBAMsortRAM="40000000000",
        readFilesCommand="gunzip -c"
    shell:
        """
        STAR --genomeDir star_index \
        --sjdbGTFfile genome/combined.gtf \
        --runThreadN {threads} \
        --readFilesCommand {params.readFilesCommand} \
        --alignIntronMax 10000 \
        --quantMode {params.quantMode} --outSAMtype {params.outSAMtype} \
        --limitBAMsortRAM {params.limitBAMsortRAM} \
        --readFilesIn {input.mate1} {input.mate2} \
        --outFileNamePrefix Output/BAM/{wildcards.ID}.
        """


rule fastqc:
    input:
        fastq1=lambda wildcards: glob.glob('FastQ/' + wildcards.ID + '*_R1.fastq.gz'),
        fastq2=lambda wildcards: glob.glob('FastQ/' + wildcards.ID + '*_R2.fastq.gz'),
    output:
        txt="fastqc/{ID}_fastqc.txt",
    threads: 12
    shell:
        """
        fastqc {input.fastq1} -t {threads} -o fastqc
        fastqc {input.fastq2} -t {threads} -o fastqc
        echo 'FastQC Finished!' > {output.txt}
        """


rule star_index:
    input:
        fasta="genome/combined_masked.fa",
        gtf="genome/combined.gtf"
    output:
        "star_index/genome",
        "star_index/Genome",
        "star_index/genome.transcripts.fa",
        "star_index/chrNameLength.txt"
    threads: 12
    shell:
        """
        rsem-prepare-reference --gtf {input.gtf} --star -p {threads} {input.fasta} star_index/genome
        touch star_index/genome #hack for RSEM input
        """


rule combine_genomes:
    input:
        genomeSc="genome/genomeSc.fa",
        genomeSp="genome/genomeSp.fa",
        genomeEc="genome/genomeEc.fa",
        gtfSc="genome/genomeSc.gtf",
        gtfSp="genome/genomeSp.gtf",
        gtfEc="genome/genomeEc.gtf",
    output:
        "genome/combined_masked.fa",
        "genome/combined.gtf",
    shell:
        """
        cat {input.genomeSc} | sed 's/^>/>scer_/' > genome/tmp_scer.fa
        cat {input.genomeSp} | sed 's/^>/>spom_/' > genome/tmp_spom.fa
        cat {input.genomeEc} | sed 's/^>/>ecol_/' > genome/tmp_ecol.fa

        cat genome/tmp_scer.fa genome/tmp_spom.fa genome/tmp_ecol.fa > genome/combined.fa

        rm genome/tmp_scer.fa
        rm genome/tmp_spom.fa
        rm genome/tmp_ecol.fa

        cat {input.gtfSc} | grep -v '#!gen' | sed 's/^/scer_/' > genome/tmp_scer.gtf
        cat {input.gtfSp} | grep -v '#!gen' | sed 's/^/spom_/' > genome/tmp_spom.gtf
        cat {input.gtfEc} | grep -v '#!gen' | sed 's/^/ecol_/' > genome/tmp_ecol.gtf


        cat genome/tmp_scer.gtf genome/tmp_spom.gtf genome/tmp_ecol.gtf | grep -v "RDN37" > genome/combined.gtf

        rm genome/tmp_scer.gtf
        rm genome/tmp_spom.gtf
        rm genome/tmp_ecol.gtf

        bedtools maskfasta -fi genome/combined.fa -bed external/rRNA_masked.bed -fo genome/combined_masked.fa
        """


rule get_index_files:
    params:
        genomeSc=config['genomeScFTP'],
        genomeSp=config['genomeSpFTP'],
        genomeEc=config['genomeEcFTP'],
        gtfSc=config['gtfScFTP'],
        gtfSp=config['gtfSpFTP'],
    output:
        "genome/genomeSc.fa",
        "genome/genomeSp.fa",
        "genome/genomeEc.fa",
        "genome/genomeSc.gtf",
        "genome/genomeSp.gtf",
    shell:
        """
        wget -O genome/genomeSc.fa.gz {params.genomeSc}
        gunzip genome/genomeSc.fa.gz

        wget -O genome/genomeSc.gtf.gz {params.gtfSc}
        gunzip genome/genomeSc.gtf.gz

        wget -O genome/genomeSp.fa.gz {params.genomeSp}
        gunzip genome/genomeSp.fa.gz

        wget -O genome/genomeSp.gtf.gz {params.gtfSp}
        gunzip genome/genomeSp.gtf.gz

        wget -O genome/genomeEc.fa.gz {params.genomeEc}
        gunzip genome/genomeEc.fa.gz
        """

onsuccess:
        print("Finished!")
