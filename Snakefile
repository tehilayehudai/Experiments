import os.path
import os
import subprocess

def read_replicates(in_path):
    with open(in_path) as in_file:
        return [k.rstrip() for k in in_file]


EXPERIMENTS = {
    "Infected_vs_Naive": ("Infected", "Naive"),
}


rule all:
    input:
        "data/model/06_best_model/best_model_Infected_vs_Naive/svg/"



rule remove_cysteine_loops:
    input:
        "input/peptides/{replicate}.faa",
    output:
        clashes="data/peptides/nocys/{replicate}.csv",
        out="data/peptides/nocys/{replicate}.faa",
    shell:
        "motifier legacy remove-cysteine-loop {input} {output.clashes} {output.out}"


rule cluster_peptides:
    input:
        "data/peptides/nocys/{replicate}.faa",
    output:
        "data/clusters/00_raw/{replicate}.csv",
    params:
        threshold=0.6,
        word_length=2,
        discard=4,
        mode=1,
    shell:
        "motifier legacy cluster {input} --threshold {params.threshold} "
        "--word_length {params.word_length} --discard {params.discard} "
        "--cluster_algorithm_mode {params.mode} > {output}"


rule reduce_clusters:
    input:
        "data/clusters/00_raw/{replicate}.csv",
    output:
        "data/clusters/01_top_1000/{replicate}.csv",
    params:
        top_clusters=100,
        top_peptides=1000,
    shell:
        "motifier legacy filter-top-peptides {input} {params.top_peptides} | "
        "motifier legacy filter-top-clusters - {params.top_clusters} > {output}"


rule restore_sequences:
    input:
        cluster="data/clusters/01_top_1000/{replicate}.csv",
        faa="input/peptides/{replicate}.faa",
        nocys="data/peptides/nocys/{replicate}.csv",
    output:
        "data/clusters/02_top_1000_orig/{replicate}.csv",
    shell:
        "motifier legacy replace-sequences {input.cluster} {input.faa} {input.nocys} > {output}"


rule align_clusters:
    input:
        "data/clusters/02_top_1000_orig/{replicate}.csv",
    output:
        "data/clusters/03_aligned/{replicate}.csv",
    threads: workflow.cores,
    benchmark:
        "benchmark/align_clusters/{replicate}.txt"
    shell:
        "motifier legacy align-sequences {input} "
        "-p {threads} --ignore-version > {output}"


rule remove_gaps:
    input:
        "data/clusters/03_aligned/{replicate}.csv",
    output:
        "data/clusters/04_gaps_removed/{replicate}.csv",
    params:
        threshold=0.5,
    shell:
        "motifier legacy remove-gappy-columns "
        "{input} {params.threshold} > {output}"


rule concat_meme_clusters:
    input:
        lambda wc: expand(
            "data/clusters/04_gaps_removed/{replicate}.csv",
            replicate=read_replicates(f"input/groups/{wc.group}.txt"),
        ),
    output:
        "data/clusters/05_concatenated_meme/{group}.meme",
    shell:
        "motifier legacy concat-clusters --prefix hiv_neg_cor_pos {input} | "
        "motifier legacy create-meme - > {output}"


rule run_unitepssms:
    input:
        "data/clusters/05_concatenated_meme/{group}.meme",
    output:
        "data/clusters/06_unite_pssms/{group}",
    params:
        aln_cutoff=24,
        pcc_cutoff=0.7,
    shell:
        "unite_pssms --pssm {input} --out {output} "
        "--aln_cutoff {params.aln_cutoff} --pcc_cutoff {params.pcc_cutoff}"


rule sort_unitepssms:
    input:
        "data/clusters/06_unite_pssms/{group}",
    output:
        "data/clusters/07_sorted_unite_pssms/{group}",
    shell:
        "motifier legacy sort-unitepssm-groups {input} > {output}"


rule concat_reduced_clusters:
    input:
        lambda wc: expand(
            "data/clusters/02_top_1000_orig/{replicate}.csv",
            replicate=read_replicates(f"input/groups/{wc.group}.txt"),
        ),
    output:
        "data/clusters/08_concatenated_top_1000/{group}.csv",
    shell:
        "motifier legacy concat-clusters --prefix hiv_neg_cor_pos {input} > {output}"


rule merge_clusters:
    input:
        csv="data/clusters/08_concatenated_top_1000/{group}.csv",
        merging_list="data/clusters/07_sorted_unite_pssms/{group}",
    output:
        "data/clusters/09_merged/{group}.csv",
    shell:
        "motifier legacy merge-clusters --prefix {wildcards.group} "
        "{input.csv} {input.merging_list} > {output}"


rule align_merged_clusters:
    input:
        "data/clusters/09_merged/{group}.csv",
    output:
        "data/clusters/10_merged_aligned/{group}.csv",
    threads: workflow.cores,
    shell:
        "motifier legacy align-sequences {input} "
        "-p {threads} --ignore-version > {output}"


rule remove_gaps_from_merged_clusters:
    input:
        "data/clusters/10_merged_aligned/{group}.csv",
    output:
        "data/clusters/11_merged_gaps_removed/{group}.csv",
    params:
        threshold=0.5,
    shell:
        "motifier legacy remove-gappy-columns "
        "{input} {params.threshold} > {output}"


rule create_merged_meme:
    input:
        "data/clusters/11_merged_gaps_removed/{group}.csv",
    output:
        "data/clusters/12_merged_meme/{group}.txt",
    shell:
        "motifier legacy create-meme {input} > {output}"


rule get_top_meme:
    input:
        "data/clusters/12_merged_meme/{group}.txt",
    output:
        "data/clusters/13_top_merged_meme/{group}.txt",
    params:
        num_motifs=400,
    shell:
        "motifier legacy head-meme {input} {params.num_motifs} > {output}"


rule calculate_cutoffs:
    input:
        "data/clusters/13_top_merged_meme/{group}.txt",
    output:
        "data/clusters/14_cutoffs/{group}.txt",
    params:
        perc=0.05,
        min_length=4,
        max_length=14,
    log:
        "logs/PSSM_score_Peptide/{group}.log",
    benchmark:
        "benchmark/PSSM_score_Peptide/{group}.txt"
    threads: workflow.cores
    shell:
        "pssm_score_peptide CalcPSSM_Cutoff "
        "--pssm {input} --pssm_cutoffs {output} --total_memes 0 "
        "--cutoff_random_peptitdes_percentile {params.perc} "
        "--min_library_length_cutoff {params.min_length} "
        "--max_library_length_cutoff {params.max_length} > {log}"


rule assign_peptides_to_motifs:
    input:
        motifs="data/clusters/13_top_merged_meme/{group}.txt",
        cutoffs="data/clusters/14_cutoffs/{group}.txt",
        peptides="input/peptides/{replicate}.faa",
    output:
        "data/model/00_hits/{replicate}_peptides_vs_{group}_motifs.txt",
    log:
        "logs/PSSM_score_Peptide/hits/{replicate}_peptides_vs_{group}_motifs.log",
    benchmark:
        "benchmark/PSSM_score_Peptide/hits/{replicate}_peptides_vs_{group}_motifs.txt"
    threads: workflow.cores
    shell:
        "pssm_score_peptide CalcPSSM_Hits "
        "--pssm {input.motifs} --pssm_cutoffs {input.cutoffs} --seq {input.peptides} "
        "--out {output} > {log}"



rule merge_hits:
    input:
        lambda wc: expand(
            f"data/model/00_hits/{{replicate}}_peptides_vs_{EXPERIMENTS[wc.label][0]}_motifs.txt",
            replicate=[
                *read_replicates(f"input/groups/{EXPERIMENTS[wc.label][0]}.txt"),
                *read_replicates(f"input/groups/{EXPERIMENTS[wc.label][1]}.txt"),
            ],
        ),
    output:
        "data/model/01_merged_hits/{label}.csv",
    params:
        condition=lambda wc: " ".join(
            f"-c {k}"
            for k in read_replicates(f"input/groups/{EXPERIMENTS[wc.label][0]}.txt")
        ),
        label=lambda wc: EXPERIMENTS[wc.label][0],
    shell:
        "motifier legacy merge-hits {input} {params.condition} "
        "-l {params.label} > {output}"

rule single_feature:
    input:
        "data/model/01_merged_hits/{label}.csv"
    output:
        "data/model/02_single_feature/{label}_single_feature.txt"
    shell: 
        "evaluate-single-features {input} > {output}"

rule hyperparameter:
    output:
        "data/model/03_hyperparameter/hyperparameter.json"
    shell: 
        "generate-hyperparameters > {output}"

rule train_random_forest:
    input:
        hits = "data/model/01_merged_hits/{label}.csv", 
        conf = "data/model/03_hyperparameter/hyperparameter.json"
    output:
        directory("data/model/04_train_rf/model_{label}_{number}/")
    shell: 
        "mkdir -p {output} && train-rf {input.hits} {output} {input.conf} {wildcards.number}"

rule select_best_model:
    input:
        expand("data/model/04_train_rf/model_{{label}}_{number}/", 
        number = range(0,100))
    output:
        "data/model/05_select_best_model/best_model_{label}.txt"
    params:
        error_rate = lambda w, input: " ".join(os.path.join(k,"error_rate.txt") for k in input)
    shell:
        "select-best-model {params.error_rate} > {output}"

def read_best_model(best_model_path):
    with open (str(best_model_path)) as f:
        best_model = os.path.split(f.readline())[0]
    return best_model

rule copy_best_model:
    input:
        "data/model/05_select_best_model/best_model_{label}.txt"
    output:
        directory("data/model/06_best_model/best_model_{label}")
    params:
        dir_name = lambda w, input: read_best_model(input)
    shell:
        "cp -R {params.dir_name} {output}"

rule create_csv_files:
    input:
        hits = "data/model/01_merged_hits/{label}.csv",
        error_rate = "data/model/06_best_model/best_model_{label}/error_rate.txt",
        sorted_feature_path = "data/model/06_best_model/best_model_{label}/sorted_feature_importance.txt"
    output:
        directory("data/model/06_best_model/best_model_{label}/csv/")
    shell:
        "mkdir -p {output} && create-csv-files {input.hits} {input.error_rate} {input.sorted_feature_path} {output}"

rule generate_heatmap:
    input:
        "data/model/06_best_model/best_model_{label}/csv/"
    output:
        "data/model/06_best_model/best_model_{label}/svg/"
    run:
        os.mkdir(str(output))
        for file in os.listdir(str(input)):
            if file.endswith(".csv"):
                subprocess.run(["generate-heatmap", os.path.join(str(input), file), os.path.join(str(output), os.path.splitext(file)[0])])


