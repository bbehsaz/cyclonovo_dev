# CycloNovo: Algorithm for finding and de novo sequencing cyclospectra

###Prerequisites

1. Copy (or symlink) `print_score` and `print_spectrum` under `./scripts/` subdirectory
2. Copy (or symlink) `./Fragmentation_rule/`
3. Copy (or symlink) `./configs/common/` under `./configs/` subdirectory

###Examples
`python cyclonovo.py -s ../cycloquest_minimal/test_data/varquest/surugamide_769.mgf -o outdirectory --pname prefix_for_generated_files --monomers standard  --fragment_ion_thresh 0.02 --precursor_ion_thresh 0.02 --kmer_score 4 --cyclointensity 60 --num_frequent_clusters 2 --alpha 0.0067 --beta -1 --preprocess`