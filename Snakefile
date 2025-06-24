CHROMOSOMES = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


rule all:
    input:
        "output/all_chroms_combined.tsv.gz"
        # expand("output/merged_dist_match_valid/{chromosome}.tsv.gz", chromosome=CHROMOSOMES),
        # expand("output/valid/{chromosome}.tsv.gz", chromosome=CHROMOSOMES),
        # "output/pairwise_distances.tsv"


rule extract_header:
    input:
        "data/phASER_GTEx_v8_matrix.txt.gz"
    output:
        "output/phased/header.txt"
    shell:
        """
        gzip -dc {input} | sed -n 1p > {output}
        """

rule split_by_chromosome:
    input:
        phased="data/phASER_GTEx_v8_matrix.txt.gz",
        header="output/phased/header.txt"
    output:
        "output/phased/{chromosome}.tsv.gz"
    shell:
        """
        # mkdir -p $(dirname {output})
        temp_file=$(mktemp)
        zcat {input.phased} | grep -P "^{wildcards.chromosome}\t" > $temp_file
        cat {input.header} $temp_file | gzip > {output}
        rm $temp_file
           """


rule h0_up_analysis:
    input:
        "output/phased/{chromosome}.tsv.gz"
    output:
        "output/h0_up/{chromosome}.tsv.gz"
    log:
        "logs/h0_up/{chromosome}.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/is_h0_up.py {input} {output} 2> {log}
        """


rule compute_match_valid:
    input:
        "output/h0_up/{chromosome}.tsv.gz"
    output:
        match = "output/match/{chromosome}.tsv.gz",
        valid = "output/valid/{chromosome}.tsv.gz"
    shell:
        """
        mkdir -p $(dirname {output.match})
        python scripts/compute_match_valid.py {input} {output.match} {output.valid}
        """


rule compute_pairwise_distances:
    input:
        "data/phASER_GTEx_v8_matrix.txt.gz"
    output:
        "output/pairwise_distances.tsv"
    shell:
        "python scripts/compute_pairwise_distance.py {input} {output}"


rule merge_scores:
    input:
        distances = "output/pairwise_distances.tsv",
        match = "output/match/{chrom}.tsv.gz",
        valid = "output/valid/{chrom}.tsv.gz"
    output:
        "output/merged_dist_match_valid/{chrom}.tsv.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/merge_scores.py \
            --chrom {wildcards.chrom} \
            --distances {input.distances} \
            --match {input.match} \
            --valid {input.valid} \
            --output {output}
        """

rule concat_scores:
    input:
        expand("output/merged_dist_match_valid/{chrom}.tsv.gz", chrom=CHROMOSOMES)
    output:
        "output/all_chroms_combined.tsv.gz"
    shell:
        """
        # Extract header from first file
        gzip -dc {input[0]} | sed -n 1p > output/combined.tmp

        # Append data from all files (skipping headers)
        for f in {input}; do
            gzip -dc "$f" | tail -n +2
        done >> combined.tmp

        # Compress the result
        gzip -c combined.tmp > {output}
        rm combined.tmp
        """