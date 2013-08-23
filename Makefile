 # Some steps just require a lot of number-crunching, but no plotting
 # or other requirements which require standard cpython.  For these
 # steps, use pypy if it's available.

.SECONDARY: # Don't delete intermediate files

PP = $(shell which pypy || which python)

git: 
	git submodule init
	git submodule update

pub/paper/paper.pdf: pub/paper/chem_pot_paper.tex pub/paper/bib/
	pdflatex pub/paper/chem_pot_paper.tex pub/paper/chem_pot_paper.pdf

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

results/fig/mu_k_diagrams/%_mu_k_diagram.png: src/generate_mu_k_diagram.py results/mu_k_tables/%_genome_exact_mu_k_table.csv results/mu_k_tables/%_control_exact_mu_k_table.csv results/mu_k_tables/%_genome_approximation_mu_k_table.csv results/mu_k_tables/%_control_approximation_mu_k_table.csv
	python src/generate_mu_k_diagram.py results/mu_k_tables/$*_genome_exact_mu_k_table.csv results/mu_k_tables/$*_control_exact_mu_k_table.csv results/mu_k_tables/$*_genome_approximation_mu_k_table.csv results/mu_k_tables/$*_control_approximation_mu_k_table.csv results/fig/mu_k_diagrams/$*_mu_k_diagram.png # use python for matplotlib...

results/fig/mu_mu_hat_diagrams/%_mu_mu_hat_diagram.png: src/generate_mu_mu_hat_diagram.py results/mu_k_tables/%_genome_exact_mu_k_table.csv results/mu_k_tables/%_genome_approximation_mu_k_table.csv 
	python src/generate_mu_mu_hat_diagram.py results/mu_k_tables/$*_genome_exact_mu_k_table.csv results/mu_k_tables/$*_genome_approximation_mu_k_table.csv results/fig/mu_mu_hat_diagrams/$*_mu_mu_hat_diagram.png # use python for matplotlib...

tables: $(all_exact_mu_k_tables)
	echo $^
	echo $(all_exact_mu_k_tables)

all_mu_mu_hats = results/fig/mu_mu_hat_diagrams/Fis_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/OxyR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/IciA_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/GlpR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/ArgR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/DnaA_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/ArcA_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/YhiX_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/TyrR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/MarA_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/Lrp_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/FruR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/NarL_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/LexA_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/SoxS_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/Crp_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/Fur_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/Fnr_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/FadR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/GlnG_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/MetJ_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/IHF_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/MalT_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/CpxR_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/PhoB_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/FliA_mu_mu_hat_diagram.png results/fig/mu_mu_hat_diagrams/Rob_mu_mu_hat_diagram.png
all_mu_mu_hats : $(all_mu_mu_hats)
# shameful variable lists go here...
all_exact_mu_k_tables = results/mu_k_tables/Fis_genome_exact_mu_k_table.csv results/mu_k_tables/Fis_genome_approximation_mu_k_table.csv results/mu_k_tables/OxyR_genome_exact_mu_k_table.csv results/mu_k_tables/OxyR_genome_approximation_mu_k_table.csv results/mu_k_tables/IciA_genome_exact_mu_k_table.csv results/mu_k_tables/IciA_genome_approximation_mu_k_table.csv results/mu_k_tables/GlpR_genome_exact_mu_k_table.csv results/mu_k_tables/GlpR_genome_approximation_mu_k_table.csv results/mu_k_tables/ArgR_genome_exact_mu_k_table.csv results/mu_k_tables/ArgR_genome_approximation_mu_k_table.csv results/mu_k_tables/DnaA_genome_exact_mu_k_table.csv results/mu_k_tables/DnaA_genome_approximation_mu_k_table.csv results/mu_k_tables/ArcA_genome_exact_mu_k_table.csv results/mu_k_tables/ArcA_genome_approximation_mu_k_table.csv results/mu_k_tables/YhiX_genome_exact_mu_k_table.csv results/mu_k_tables/YhiX_genome_approximation_mu_k_table.csv results/mu_k_tables/TyrR_genome_exact_mu_k_table.csv results/mu_k_tables/TyrR_genome_approximation_mu_k_table.csv results/mu_k_tables/MarA_genome_exact_mu_k_table.csv results/mu_k_tables/MarA_genome_approximation_mu_k_table.csv results/mu_k_tables/Lrp_genome_exact_mu_k_table.csv results/mu_k_tables/Lrp_genome_approximation_mu_k_table.csv results/mu_k_tables/FruR_genome_exact_mu_k_table.csv results/mu_k_tables/FruR_genome_approximation_mu_k_table.csv results/mu_k_tables/NarL_genome_exact_mu_k_table.csv results/mu_k_tables/NarL_genome_approximation_mu_k_table.csv results/mu_k_tables/LexA_genome_exact_mu_k_table.csv results/mu_k_tables/LexA_genome_approximation_mu_k_table.csv results/mu_k_tables/SoxS_genome_exact_mu_k_table.csv results/mu_k_tables/SoxS_genome_approximation_mu_k_table.csv results/mu_k_tables/Crp_genome_exact_mu_k_table.csv results/mu_k_tables/Crp_genome_approximation_mu_k_table.csv results/mu_k_tables/Fur_genome_exact_mu_k_table.csv results/mu_k_tables/Fur_genome_approximation_mu_k_table.csv results/mu_k_tables/Fnr_genome_exact_mu_k_table.csv results/mu_k_tables/Fnr_genome_approximation_mu_k_table.csv results/mu_k_tables/FadR_genome_exact_mu_k_table.csv results/mu_k_tables/FadR_genome_approximation_mu_k_table.csv results/mu_k_tables/GlnG_genome_exact_mu_k_table.csv results/mu_k_tables/GlnG_genome_approximation_mu_k_table.csv results/mu_k_tables/MetJ_genome_exact_mu_k_table.csv results/mu_k_tables/MetJ_genome_approximation_mu_k_table.csv results/mu_k_tables/IHF_genome_exact_mu_k_table.csv results/mu_k_tables/IHF_genome_approximation_mu_k_table.csv results/mu_k_tables/MalT_genome_exact_mu_k_table.csv results/mu_k_tables/MalT_genome_approximation_mu_k_table.csv results/mu_k_tables/CpxR_genome_exact_mu_k_table.csv results/mu_k_tables/CpxR_genome_approximation_mu_k_table.csv results/mu_k_tables/PhoB_genome_exact_mu_k_table.csv results/mu_k_tables/PhoB_genome_approximation_mu_k_table.csv results/mu_k_tables/FliA_genome_exact_mu_k_table.csv results/mu_k_tables/FliA_genome_approximation_mu_k_table.csv results/mu_k_tables/Rob_genome_exact_mu_k_table.csv results/mu_k_tables/Rob_genome_approximation_mu_k_table.csv


results/fig/misclassification_rates.png: src/generate_misclassification_plot.py $(all_exact_mu_k_tables)
	echo $^
	echo $(all_exact_mu_k_tables)
	python src/generate_misclassification_plot.py results/fig/misclassification_rates.png
