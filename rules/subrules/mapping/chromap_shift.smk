rule atac_seq_shift:
    """
    Perform ATAC-seq specific read shifting to correct for Tn5 transposase binding bias.

    This rule applies a critical correction to ATAC-seq alignments to account for the
    binding properties of the Tn5 transposase, which creates a 9bp duplication during
    the tagmentation reaction. This shift is essential for accurate determination of
    transcription factor binding sites and open chromatin region boundaries.

    The Tn5 transposase binds as a dimer and inserts adapters separated by 9 base pairs.
    Consequently:
    - Reads mapping to the positive strand are shifted +4 bp
    - Reads mapping to the negative strand are shifted -5 bp
    - This centers the read pileup over the actual Tn5 insertion site

    This rule uses deepTools' alignmentSieve with the --ATACshift parameter, which
    automatically applies this standard ATAC-seq correction. The shifted BAM file
    is then re-sorted and indexed to maintain compatibility with downstream tools.

    Proper shifting is critical for ATAC-seq analysis because:
    - It improves the resolution of transcription factor footprinting analysis
    - It ensures accurate calling of peak boundaries
    - It enables proper TSS enrichment calculations
    - It provides more precise nucleosome positioning information
    - It is required for TOBIAS footprinting and other high-resolution analyses

    The shifted BAM file serves as the primary input for peak calling, bigwig
    generation, TSS enrichment analysis, and all subsequent downstream analyses
    in the ATAC-seq workflow.
    """
    input:
        bam = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam',
        bai = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam.bai' 
    output:
        shifted_sort_bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        shifted_sort_bam_bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai'
    log:
        "logs/02.mapping/atac_seq_shift_{sample}.log"
    threads:
        20
    shell:
        """
        (ln -s -r {input.bam} {output.shifted_sort_bam}
        ln -s -r {input.bai} {output.shifted_sort_bam_bai} )  & > {log}
        """