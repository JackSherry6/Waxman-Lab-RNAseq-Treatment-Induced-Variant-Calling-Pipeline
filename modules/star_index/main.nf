process STAR_INDEX {
    label 'process_high'
    conda 'envs/star_env.yml'

    input:
    path(gtf)
    path(ref_genome)

    output:
    path "star", emit: index

    script:
    """
    mkdir -p star
    STAR --runMode genomeGenerate \
         --runThreadN $task.cpus \
         --genomeDir star \
         --genomeFastaFiles $ref_genome \
         --sjdbGTFfile $gtf \
         --sjdbOverhang 150
    """

    stub:
    """
    mkdir -p star
    """
}