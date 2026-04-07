process CSV_TO_XLSX {
    label 'process_single'
    conda 'envs/biopython_env.yml'
    publishDir "${params.outdir}/individual_excels", mode: 'copy'

    input:
    path(variants_csv)

    output:
    path("${variants_csv.simpleName}.xlsx")

    script:
    """
    csv_to_excel_with_key.py $variants_csv -o ${variants_csv.simpleName}.xlsx
    """
}

// ask claude to modify the python script to check what the column names are and build the index accordingly (but first I have to actually make the index)
