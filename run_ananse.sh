# ANANSE BINDING

# Source cell type
ananse binding -C enhancers_source.tsv -o source_binding -g hg38

# Target cell type
ananse binding -C enhancers_target.tsv -o target_binding -g hg38


# ANANSE NETWORK

# Source cell type
ananse network source_binding/binding.h5 -e TPM_source_*.txt -o network_source.tsv -g hg38 -n 6

# Target cell type
ananse network target_binding/binding.h5 -e TPM_target_*.txt -o network_target.tsv -g hg38 -n 6


# ANANSE INFLUENCE

ananse influence -s network_source.tsv -t network_target.tsv -d DE_genes.txt -o influence.tsv -i 500000 -n 8

ananse plot influence.tsv -d influence_diffnetwork.tsv -o influence_results
