configfile: "config.yaml"

WEB_PATH = "/var/www/html/genotyping_figures/"

include:
    "analysis.smk"

include:
    "figures.smk"

include:
    "read_simulation.smk"

include:
    "paragraph.smk"

include:
    "vg.smk"

include:
    "glimpse.smk"

include:
    "low_coverage_genotyping.smk"

include:
    "kage_with_mapped_reads.smk"

ruleorder:
    downsample_real_reads15x > convert_fa_to_fq

ruleorder:
    downsample_real_reads30x > convert_fa_to_fq


rule all:
    input:
        "data/dataset1/happy-hg002-kageNoHelperModelN250all_hg002_simulated_reads_15x.extended.csv",
        #"data/dataset1/debugging-usN20all-hg002-hg002_simulated_reads_15x.txt"
        #"data/dataset1/happy-hg002-usN250all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset2/happy-hg002-usN25all_hg002_real_reads_15x.extended.csv"
        #"data/dataset2/happy-hg002-pangenie_hg002_real_reads_15x.extended.csv"
        #"data/dataset1/debugging-naivekageN2548all-hg002-hg002_simulated_reads_15x.txt"
        #"data/dataset2/haplotype_to_nodes_2548all.npz"
        #"data/dataset1/tricky_variants_2548all.npy"
        #"figure1.html",
        #"data/dataset1/happy-hg002-naivekageN2548all_hg002_simulated_reads_15x.extended.csv"
        #"glimpse/GLIMPSE_phase_static"
        #"data/dataset1/happy-hg002-kageNoHelperModelN250all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset4/happy-hg002-usN250all_hg002_simulated_reads_15x.extended.csv"
        #"data/svdataset1/truth_NA12889_regions.bed"
        #"data/dataset1/debugging-usN250all-hg002-hg002_simulated_reads_15x.txt"
        #"data/dataset1/haplotype_to_nodes_250all.npz"
        #"data/dataset1/happy-hg002-usN250all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset2/happy-hg002-usN2548all_hg002_real_reads_15x.extended.csv"
        #"data/dataset2/debugging-usN250all-hg002-hg002_simulated_reads_15x.txt"
        #"data/svdataset1/debugging-usN154all-NA12889-NA12889_simulated_reads_15x.txt"
        #"data/svdataset2/happy-NA12889-usN154all_NA12889_simulated_reads_15x.extended.csv"
        #"data/svdataset2/happy-NA12889-paragraph_NA12889_simulated_reads_15x.extended.csv"
        #"data/svdataset1/happy-NA12889-usN154all_NA12889_simulated_reads_15x.extended.csv"
        #"data/dataset4/happy-hg002-usN250all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset1/debugging-usN2548all-hg002-hg002_simulated_reads_15x.txt"
        #"data/simulated_dataset1/debugging-usN100all-seed1-seed1_simulated_reads_15x.txt"
        #"data/simulated_dataset1/happy-seed1-usN100all_seed1_simulated_reads_15x.extended.csv"
        #"data/dataset1/happy-hg002-usN500all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset1/sampling_count_model_100all.npz"
        #"data/dataset6/happy-hg002-kageMappedReadsN2548all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset6/hg002_simulated_reads_15x.node_counts_giraffe.npy"
        #"data/dataset6/obgraph.npz"
        #"data/dataset1/debugging-kageMappedReadsN2548all-hg002-hg002_simulated_reads_15x.txt"
        #"data/dataset1/model_for_read_mapping_with_sampled_reads.npz"
        #"data/dataset1/happy-hg002-kageMappedReadsN2548all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset1/hg002_simulated_reads_15x.node_counts_giraffe.npy"
        #"data/dataset1/obgraph.npz"
        #"data/dataset1/happy-hg002-usN2548all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset3/usN2548all_hg002_simulated_reads_15x.vcf.gz"
        #"data/dataset3/obgraph.npz"
        #"data/dataset5/gatkRefined_hg002_simulated_reads_4x.vcf.gz"
        #"data/dataset5/variants.map"
        #"data/dataset4/happy-hg002-usN2058all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset1/happy-hg002-usN2548all_hg002_simulated_reads_15x.extended.csv"
        #"data/simulated_dataset1/obgraph.npz",
        #"data/simulated_dataset1/happy-seed1-usN800all_seed1_simulated_reads_15x.extended.csv",
        #"data/dataset1/haplotype_to_nodes_2548all.npz"
        #"data/dataset1/haplotype_matrix_2548all.npy"
        #"data/svdataset1/debugging-usN150all-NA12889-NA12889_simulated_reads_15x.txt"
        #WEB_PATH + "svtable1.html"
        #"data/svdataset2/happy-NA12889-paragraph_NA12889_simulated_reads_15x.extended.csv"
        #"data/svdataset1/paragraph_NA12889_simulated_reads_15x/genotypes.vcf.gz"
        #"data/svdataset1/paragraph_samples_NA12889_simulated_reads_15x.txt"
        #"data/svdataset2/happy-NA12889-usN150all_NA12889_simulated_reads_15x.extended.csv"
        #"data/svdataset1/NA12889_simulated_reads_15x.fa"
        #WEB_PATH + "table2.html"
        #"data/svdataset1/happy-NA12889-usN150all_NA12889_simulated_reads_15x.extended.csv"
        #"data/svdataset1/variant_to_nodes.npz"
        #"data/svdataset1/variants.vcf.gz"
        #"data/sv-variants.vcf.gz"
        #"data/dataset2/happy-hg002-usN2548all_hg002_real_reads_15x.extended.csv"
        #"data/dataset4/happy-hg002-usN2058all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset4/index_bundle.npz"
        #"data/dataset2/flat_kmers_for_model.npz"
        #"data/dataset3/happy-hg002-usN2548all_hg002_simulated_reads_15x.extended.csv"
        #"data/dataset1/linear_kmers_with_reverse_complements.npz",
        #"data/dataset1/linear_kmers_counter.npz",
        #"data/dataset1/happy-hg002-usN2548all_hg002_simulated_reads_15x.extended.csv",
        #"data/dataset2/benchmarks/mapI1000_hg002_real_reads_15x.tsv"
        #WEB_PATH  + "table2.html"
        #"data/dataset2/usN2548all_hg002_real_reads_15x.vcf.gz"
        #WEB_PATH + "supplementary_table2.html",
        #WEB_PATH + "table20.html",
        #WEB_PATH + "table2.html",
        #"data/dataset1/happy-hg002-usN2548all_hg002_simulated_reads_15x.extended.csv",
        #"data/dataset2/happy-hg002-usN2548all_hg002_real_reads_15x.extended.csv",
        #"data/dataset2/happy-hg006-usN85chinese_hg006_real_reads_15x.extended.csv",
        #"data/dataset1/genotype_matrix_50chinese.npy"
        #"data/dataset1/variants_50chinese.vcf.gz"
        #WEB_PATH + "supplementary_table2.html"
        #WEB_PATH + "table20.html",
        #"data/dataset2/happy-hg006-usN2548_hg006_real_reads_15x.extended.csv",
        #"WEB_PATH + "table1.html",
        #"data/simulated_dataset1/happy-seed1-usN800_seed1_simulated_reads_15x.extended.csv",
        #WEB_PATH + "figure2.html"
        #"data/simulated_dataset1/happy-seed1-usN800_seed1_simulated_reads_15x.extended.csv",
        #WEB_PATH + "supplementary_table.html"
        #"data/simulated_dataset0/index_2058individuals.npz"
        #WEB_PATH + "supplementary_table.html"
        #"data/dataset1/happy-hg002-usN2548_hg002_simulated_reads_15x.extended.csv",
        #"data/dataset2/happy-hg002-malva_hg002_real_reads_15x.extended.csv",
        #"data/dataset1/happy-hg002-usN2548_hg002_simulated_reads_15x.extended.csv",
        #WEB_PATH + "f1_figure_whole_genome.html"
        #"data/dataset3/happy-hg002-usN2058_hg002_simulated_reads_15x.extended.csv",
        #"data/simulated_dataset1/happy-seed1-usN800_seed1_simulated_reads_15x.extended.csv",
        #"data/dataset1/happy-hg002-gatk_hg002_simulated_reads_15x.extended.csv",
        #WEB_PATH + "figure1.html",
        #"data/dataset1/debugging-usN2058-hg002-hg002_simulated_reads_15x.txt",
        #"data/dataset2/debugging-usN2058-hg002-hg002_real_reads_15x.txt"
        #"data/dataset2/happy-hg002-usN2058_hg002_real_reads_15x.extended.csv",
        #"WEB_PATH + "table1.html",
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.extended.csv",
        #"data/dataset1/happy-hg002-usN30_hg002_simulated_reads_15x.summary.csv",
        #WEB_PATH + "table2.html",
        #"data/dataset2/usN2058_hg002_real_reads_15x.vcf.gz"
        #WEB_PATH + "table10.html",
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv",
        #"data/dataset1/debugging-usN2058-hg002-hg002_simulated_reads_15x.txt",
        #"data/simulated_dataset2/debugging-usN1000-seed1-seed1_simulated_reads_15x.txt",
        #WEB_PATH + "table2.html"
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv",
        #"data/dataset1/debugging-usN2058-hg002-hg002_simulated_reads_15x.txt",
        #"table_simulated_data.html",
        #"data/simulated_dataset1/happy-seed1-usN800_seed1_simulated_reads_15x.summary.csv",
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv"
        #"data/dataset1/truth_hg002.only_callable.vcf.gz"
        #WEB_PATH + "figure2.html",
        #WEB_PATH + "table2.html"
        #"data/dataset3/debugging-usN2058-hg002-hg002_simulated_reads_15x.txt"
        #WEB_PATH + "table9.html"
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv"
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv"
        #"data/dataset2/debugging-usN15-seed1-seed1_simulated_reads_15x.txt"
        #"data/dataset1/debugging-usN2058-hg002-hg002_simulated_reads_15x.txt"
        #"data/dataset2/debugging-usN2058-hg002-hg002_real_reads_15x.txt"
        #WEB_PATH + "table2.html"
        #"data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv"
        #"data/simulated_dataset1/happy-seed1-usN800_seed1_simulated_reads_15x.summary.csv",
        #"data/simulated_dataset2/happy-seed1-usN50_seed1_simulated_reads_15x.summary.csv",
        #"data/simulated_dataset2/happy-seed1-pangenieN15_seed1_simulated_reads_15x.summary.csv",
        #"data/simulated_dataset1/happy-seed1-usN500_seed1_simulated_reads_15x.summary.csv",
        #"table_simulated_data.html",
        #"table_simulated_dataset1-seed1_simulated_reads_15x-seed1.html"
        #"data/simulated_dataset1/variants.vcf.gz",
        #"data/simulated_dataset1/usN100_seed1_simulated_reads_10x.vcf.gz",
        #"data/simulated_dataset1/happy-seed1-usN500_seed1_simulated_reads_15x.summary.csv",
        #"data/simulated_dataset1/happy-seed1-pangenie_seed1_simulated_reads_15x.summary.csv"
        ##"data/dataset2/pangenie_hg002_simulated_reads_15x.15individuals.vcf"
        #WEB_PATH + "table2.html"
        #WEB_PATH + "table4.html",
        #WEB_PATH + "table3.html",
        #WEB_PATH + "figure4.html",
        #WEB_PATH + "table2.html"
        #"/var/www/html/genotyping_figures/table1.html"
        #"/var/www/html/genotyping_figures/figure3.html"
        #"data/dataset1/hg002_real_reads_15x.fa"
        #"data/dataset1/gatk_hg002_simulated_reads_15x.vcf.gz"
        #"figure3.html"
        #"data/dataset1/bayestyper_hg002_simulated_reads_15x.vcf.gz"
        #"/var/www/html/genotyping_figures/figure3.html"
        #"data/dataset1/hg002_simulated_reads_1x.fa"
        #"data/dataset1/hg002_coordinate_map_chromosome1_haplotype0.npz"
        #"figure1.html"
        #"data/dataset1/happy-hg002-malva_simulated_reads15x.summary.csv"
        #"data/dataset1/malva_genotypes_simulated_reads15x.fa"
        #"data/dataset1/happy-hg002-genotypes_simulated_reads15x.summary.csv"

