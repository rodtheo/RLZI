import pandas as pd

configfile: "configs/config-test.yaml"

samples = pd.read_table(config['sample_tsv']).set_index("sample", drop=False)

species = config['species']

ms = config['parameters']['m']
ks = config['parameters']['k']


def get_sample_names(wildcards):
        return samples.loc[wildcards.samplesamp, "id"]

def get_ref_name(wildcards):
        return samples.loc[wildcards.sampler, "id"]

if config['plot_distribution']:
    rule all:
        input:
            expand("{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/synteny_to_draw.txt", species=config['species'], sampler=samples['sample'][0], samplesamp=samples['sample'][1:], k=ks, m=ms), expand("{species}/final_stats.txt", species=config['species']), expand("{species}/{sampler}_{samplesamp}/distribution_factors.png", species=config['species'], sampler=samples['sample'][0], samplesamp=samples['sample'][1])
else:
    rule all:
        input:
            expand("{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/chromosome.svg", species=config['species'], sampler=samples['sample'][0], samplesamp=samples['sample'][1:], k=ks, m=ms), expand("{species}/final_stats.txt", species=config['species'])


if config['download_fasta']:
    rule get_sample_genomes:
        params:
            genome=get_sample_names
        output:
            "{species}/{samplesamp}.fa"
        conda:
            "envs/downloadEnv.yaml"
        shell:
            "ncbi-acc-download -m nucleotide --out /dev/stdout -F fasta -v {params.genome} > {output}"

# rule get_ref_genome:
#     output:
#         "{species}/{sampler}.fa"
#     conda:
#         "envs/downloadEnv.yaml"
#     shell:
#         "ncbi-acc-download -m nucleotide --out /dev/stdout -F fasta -v {input} > {output}"

rule RLZ_parsing:
    input:
        ref="{species}/{sampler}.fa",
        samp="{species}/{samplesamp}.fa"
    output:
        "{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/RLZ_k{k}_m{m}.log"
    benchmark:
        "{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/benchmark.txt"
    shell:
        "RLZI.x --reverse-complement -k {wildcards.k} -m {wildcards.m} -R {input.ref} -S {input.samp} -o {wildcards.species}/{wildcards.sampler}_{wildcards.samplesamp}/RLZ_k{wildcards.k}_m{wildcards.m}/RLZ_k{wildcards.k}_m{wildcards.m} > {output}"

rule prepare_to_draw_synteny:
    input:
        flag_finished_rlzi="{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/RLZ_k{k}_m{m}.log",
        ref="{species}/{sampler}.fa",
        samp="{species}/{samplesamp}.fa"
    output:
        chr_info="{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/chr_infos.txt",
        synt="{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/synteny_to_draw.txt",
        syntall="{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/synteny_to_draw_all.txt"
    conda:
        "envs/flexidotEnv.yaml"
    shell:
        "python scripts/parse_to_draw_synteny.py {input.ref} {input.samp} {wildcards.species}/{wildcards.sampler}_{wildcards.samplesamp}/RLZ_k{wildcards.k}_m{wildcards.m}/RLZ_k{wildcards.k}_m{wildcards.m}_sdsl.txt {output.chr_info} {output.synt} {output.syntall}"

rule draw_synteny:
    input:
        chr_info="{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/chr_infos.txt",
        synt="{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/synteny_to_draw.txt"
    output:
        "{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/chromosome.svg"
    conda:
        "envs/ideogramEnv.yaml"
    script:
        "scripts/draw_ideogram.R"

if config['plot_distribution']:
    rule plot_distribution:
        input:
            expand("{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/synteny_to_draw.txt", species=["synthC"], sampler=["RefL"], samplesamp=["TarL"], k=[10], m=[2,5,7,9]+[i for i in range(1,101,10)])
        output:
            "{species}/{sampler}_{samplesamp}/distribution_factors.png"
        shell:
            "python scripts/plot_distributions.py {output} {input}"

rule get_stats:
    input:
        expand("{species}/{sampler}_{samplesamp}/RLZ_k{k}_m{m}/synteny_to_draw_all.txt", species=config['species'], sampler=samples['sample'][0], samplesamp=samples['sample'][1:], k=config['parameters']['k'][0], m=config['parameters']['m'][0])
    output:
        "{species}/final_stats.txt"
    shell:
        "python scripts/rank_translocations_ecoli.py {input} > {output}"
