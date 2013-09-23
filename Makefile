 # Some steps just require a lot of number-crunching, but no plotting
 # or other requirements which require standard cpython.  For these
 # steps, use pypy if it's available.

.SECONDARY: # Don't delete intermediate files

PP = $(shell which pypy || which python)

git: 
	git submodule init
	git submodule update

data/NC_000913.fna:
	curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna -o data/NC_000913.fna

# nomenclature: tf_{genome|control}_{exact|approximation}

results/binding_landscapes/%_genome_binding_landscape.dat: src/generate_binding_landscapes.py data/NC_000913.fna
	echo using: $(PP) for fast numerics
	$(PP) src/generate_binding_landscapes.py $* genome

results/binding_landscapes/%_control_binding_landscape.dat: src/generate_binding_landscapes.py
	echo using: $(PP) for fast numerics
	$(PP) src/generate_binding_landscapes.py $* control

results/mu_k_tables/%_exact_mu_k_table.csv: src/generate_mu_k_relation.py results/binding_landscapes/%_binding_landscape.dat # this matches both tf_genome and tf_control
	$(PP) src/generate_mu_k_relation.py \
	results/binding_landscapes/$*_binding_landscape.dat \
	results/mu_k_tables/$*_exact_mu_k_table.csv

results/mu_k_tables/%_approximation_mu_k_table.csv: src/generate_mu_k_relation.py results/binding_landscapes/%_binding_landscape.dat # this matches both tf_genomme and tf_control
	$(PP) src/generate_mu_k_relation.py \
	results/binding_landscapes/$*_binding_landscape.dat \
	results/mu_k_tables/$*_approximation_mu_k_table.csv \
	approximation

results/fig/mu_k_figs/%_mu_k_fig.png: src/generate_mu_k_fig.py results/mu_k_tables/%_genome_exact_mu_k_table.csv results/mu_k_tables/%_control_exact_mu_k_table.csv results/mu_k_tables/%_genome_approximation_mu_k_table.csv results/mu_k_tables/%_control_approximation_mu_k_table.csv
	python src/generate_mu_k_fig.py results/mu_k_tables/$*_genome_exact_mu_k_table.csv results/mu_k_tables/$*_control_exact_mu_k_table.csv results/mu_k_tables/$*_genome_approximation_mu_k_table.csv results/mu_k_tables/$*_control_approximation_mu_k_table.csv results/fig/mu_k_figs/$*_mu_k_fig.png # use python for matplotlib...

results/fig/mu_mu_hat_figs/%_mu_mu_hat_fig.png: src/generate_mu_mu_hat_fig.py results/mu_k_tables/%_genome_exact_mu_k_table.csv results/mu_k_tables/%_genome_approximation_mu_k_table.csv 
	python src/generate_mu_mu_hat_fig.py results/mu_k_tables/$*_genome_exact_mu_k_table.csv results/mu_k_tables/$*_genome_approximation_mu_k_table.csv results/fig/mu_mu_hat_figs/$*_mu_mu_hat_fig.png # use python for matplotlib...

results/fig/k_occupancy_figs/%_k_occupancy_fig.png: src/generate_k_occupancy_fig.py results/mu_k_tables/%_genome_exact_mu_k_table.csv results/mu_k_tables/%_genome_approximation_mu_k_table.csv 
	python src/generate_k_occupancy_fig.py results/mu_k_tables/$*_genome_exact_mu_k_table.csv results/mu_k_tables/$*_genome_approximation_mu_k_table.csv results/fig/k_occupancy_figs/$*_k_occupancy_fig.png # use python for matplotlib...
tables: $(all_exact_mu_k_tables)
	echo $^
	echo $(all_exact_mu_k_tables)

all_mu_mu_hats = results/fig/mu_mu_hat_figs/Fis_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/OxyR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/IciA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/GlpR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/ArgR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/DnaA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/ArcA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/YhiX_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/TyrR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/MarA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Lrp_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/FruR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/NarL_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/LexA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/SoxS_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Crp_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Fur_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Fnr_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/FadR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/GlnG_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/MetJ_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/IHF_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/MalT_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/CpxR_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/PhoB_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/FliA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Rob_mu_mu_hat_fig.png
all_mu_mu_hats : $(all_mu_mu_hats)
# shameful variable lists go here...
all_exact_mu_k_tables = results/mu_k_tables/Fis_genome_exact_mu_k_table.csv results/mu_k_tables/Fis_genome_approximation_mu_k_table.csv results/mu_k_tables/OxyR_genome_exact_mu_k_table.csv results/mu_k_tables/OxyR_genome_approximation_mu_k_table.csv results/mu_k_tables/IciA_genome_exact_mu_k_table.csv results/mu_k_tables/IciA_genome_approximation_mu_k_table.csv results/mu_k_tables/GlpR_genome_exact_mu_k_table.csv results/mu_k_tables/GlpR_genome_approximation_mu_k_table.csv results/mu_k_tables/ArgR_genome_exact_mu_k_table.csv results/mu_k_tables/ArgR_genome_approximation_mu_k_table.csv results/mu_k_tables/DnaA_genome_exact_mu_k_table.csv results/mu_k_tables/DnaA_genome_approximation_mu_k_table.csv results/mu_k_tables/ArcA_genome_exact_mu_k_table.csv results/mu_k_tables/ArcA_genome_approximation_mu_k_table.csv results/mu_k_tables/YhiX_genome_exact_mu_k_table.csv results/mu_k_tables/YhiX_genome_approximation_mu_k_table.csv results/mu_k_tables/TyrR_genome_exact_mu_k_table.csv results/mu_k_tables/TyrR_genome_approximation_mu_k_table.csv results/mu_k_tables/MarA_genome_exact_mu_k_table.csv results/mu_k_tables/MarA_genome_approximation_mu_k_table.csv results/mu_k_tables/Lrp_genome_exact_mu_k_table.csv results/mu_k_tables/Lrp_genome_approximation_mu_k_table.csv results/mu_k_tables/FruR_genome_exact_mu_k_table.csv results/mu_k_tables/FruR_genome_approximation_mu_k_table.csv results/mu_k_tables/NarL_genome_exact_mu_k_table.csv results/mu_k_tables/NarL_genome_approximation_mu_k_table.csv results/mu_k_tables/LexA_genome_exact_mu_k_table.csv results/mu_k_tables/LexA_genome_approximation_mu_k_table.csv results/mu_k_tables/SoxS_genome_exact_mu_k_table.csv results/mu_k_tables/SoxS_genome_approximation_mu_k_table.csv results/mu_k_tables/Crp_genome_exact_mu_k_table.csv results/mu_k_tables/Crp_genome_approximation_mu_k_table.csv results/mu_k_tables/Fur_genome_exact_mu_k_table.csv results/mu_k_tables/Fur_genome_approximation_mu_k_table.csv results/mu_k_tables/Fnr_genome_exact_mu_k_table.csv results/mu_k_tables/Fnr_genome_approximation_mu_k_table.csv results/mu_k_tables/FadR_genome_exact_mu_k_table.csv results/mu_k_tables/FadR_genome_approximation_mu_k_table.csv results/mu_k_tables/GlnG_genome_exact_mu_k_table.csv results/mu_k_tables/GlnG_genome_approximation_mu_k_table.csv results/mu_k_tables/MetJ_genome_exact_mu_k_table.csv results/mu_k_tables/MetJ_genome_approximation_mu_k_table.csv results/mu_k_tables/IHF_genome_exact_mu_k_table.csv results/mu_k_tables/IHF_genome_approximation_mu_k_table.csv results/mu_k_tables/MalT_genome_exact_mu_k_table.csv results/mu_k_tables/MalT_genome_approximation_mu_k_table.csv results/mu_k_tables/CpxR_genome_exact_mu_k_table.csv results/mu_k_tables/CpxR_genome_approximation_mu_k_table.csv results/mu_k_tables/PhoB_genome_exact_mu_k_table.csv results/mu_k_tables/PhoB_genome_approximation_mu_k_table.csv results/mu_k_tables/FliA_genome_exact_mu_k_table.csv results/mu_k_tables/FliA_genome_approximation_mu_k_table.csv results/mu_k_tables/Rob_genome_exact_mu_k_table.csv results/mu_k_tables/Rob_genome_approximation_mu_k_table.csv


results/fig/misclassification_rates.png: src/generate_misclassification_plot.py $(all_exact_mu_k_tables)
	echo $^
	echo $(all_exact_mu_k_tables)
	python src/generate_misclassification_plot.py results/fig/misclassification_rates.png

paper_mu_k_figs = results/fig/mu_k_figs/LexA_mu_k_fig.png results/fig/mu_k_figs/Crp_mu_k_fig.png results/fig/mu_k_figs/FliA_mu_k_fig.png results/fig/mu_k_figs/Fur_mu_k_fig.png results/fig/mu_k_figs/OxyR_mu_k_fig.png
paper_mu_mu_hat_figs = results/fig/mu_mu_hat_figs/LexA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Crp_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/FliA_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/Fur_mu_mu_hat_fig.png results/fig/mu_mu_hat_figs/OxyR_mu_mu_hat_fig.png
paper_k_occupancy_figs = results/fig/k_occupancy_figs/LexA_k_occupancy_fig.png results/fig/k_occupancy_figs/Crp_k_occupancy_fig.png results/fig/k_occupancy_figs/FliA_k_occupancy_fig.png results/fig/k_occupancy_figs/Fur_k_occupancy_fig.png results/fig/k_occupancy_figs/OxyR_k_occupancy_fig.png
pub/paper/paper.pdf: pub/paper/chem_pot_paper.tex pub/paper/bib/  $(paper_mu_k_figs)  $(paper_mu_mu_hat_figs)  $(paper_k_occupancy_figs)
	cd pub/paper; pdflatex chem_pot_paper.tex