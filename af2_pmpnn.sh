#!/bin/bash
# @jderoo

##################
### EDIT BELOW ###
##################

folder_with_pdbs=input_dir           # where are our input PDBs to start the cycle? This is also where data will end up

seqs_per_run=15                         # how many sequences should we make per input structure? Remember a structure will 
                                       # be predicted for each seq here! Time intensive variable

determine_CDRs=martin                 # this variable is ONLY allowed to be 'structure' OR an antibody numbering scheme
                                       # I have only tested kabat, martin, and chothia for these values. 

ss_near_CDRs=3                         # locate all residues within X Angstroms of the CDRs to group together with as the secondary shell. 

simple_grab=true                       # When doing the distance calculations for the CDRs, do a simple grab or do an intelligent grab of the neighbors

output_dir=output_dir                # temp output directory where less important (mediary ProteinMPNN) data ends up

to_design=framework                  # what are we designing? ONLY 3 options: "framework", "loops", and "lss". lss = loops and secondary shell,
                                       # or the Vernier zone region. The part of the scFv thats 'in between' the CDR loops and the framework.

linker_seq=GGGGSGGGGSGGGGS                          # specify linker sequence: if none is specified and left blank, automatically try to detect
                                       # linker sequences from the following list: 
				       # ['GGGGSGGGGSGGGGS', 'GGSGGSGGSGGSGGS', 'GSGSGSGSGSGSGS', 'GGGSGGGSGGGSGGS', 'GGGSGGGSGGGS']

PMPNN=home/c21094846/Proteinmpnn # the global path to where the ProteinMPNN github was cloned to locally

# How was this used to generate data?   nohup bash af2_pmpnn.sh &

##################
### EDIT ABOVE ###
##################


if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

# make this input sequence and not input structure
pdb_count=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.pdb" | wc -l)
fasta_count=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.fasta" | wc -l)

pdbs=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.pdb")
fastas=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.fasta")

echo "will process $pdbs"
echo "will process $fastas"


# ensure we're in right env
# this is my env - ignored for prod
# if [ "$CONDA_DEFAULT_ENV" != "pmpnn" ]; then
#     source /home/tagteam/anaconda3/etc/profile.d/conda.sh
#     conda activate pmpnn
# fi


if [ "$pdb_count" -gt 0 ]; then
    for pdb_file in $pdbs; do
        # Extract the base name without the extension
        base_name=$(basename "$pdb_file" .pdb)
    
        if [[ "$base_name" == *"seq"* ]]; then
            echo "Error: 'seq' pattern found not allowed in input structures! Found 'seq' in: $base_name"
            exit 1
        fi
    
        # Create a directory for this base name if it doesn't already exist
        mkdir -p "$folder_with_pdbs/$base_name"
    
        # Move the .pdb file into its corresponding directory
        mv "$pdb_file" "$folder_with_pdbs/$base_name"
    done
fi

echo "finished moving $pdbs"

if [ "$fasta_count" -gt 0 ]; then
    for fasta_file in $fastas; do
        base_name=$(basename "$fasta_file" .fasta)
    
        if [[ "$base_name" == *"seq"* ]]; then
            echo "Error: 'seq' pattern found not allowed in input structures! Found 'seq' in: $base_name"
            exit 1
        fi
    
        # Create a directory for this base name if it doesn't already exist
        mkdir -p "$folder_with_pdbs/$base_name"
        
        echo "$base_name has no structure - creating one..."
        scripts/scfv_anarci.py $fasta_file --output ${base_name}_scfv.fasta --polyG
        echo $PWD
	colabfold_batch --templates --num-recycle 6 --model-type alphafold2_multimer_v2 ${base_name}_scfv.fasta ${base_name}_initial_structures 
        cp ${base_name}_initial_structures/*rank_001*pdb ./${base_name}.pdb  
        mv ./${base_name}.pdb "$folder_with_pdbs/$base_name" 
        mv ./${base_name}_scfv.fasta "$folder_with_pdbs/$base_name" 
        mv ./${fasta_file} "$folder_with_pdbs/$base_name" 
	rm -rf ${base_name}_initial_structures 
    done
fi
echo "finished moving $fastas"

for dir in "$folder_with_pdbs"/*; do
    dir=${dir%/}
    echo "Processing directory $dir"

    pdb_file=$(find "$dir" -type f -name "*.pdb")
    if [ -n "$pdb_file" ]; then
        echo "Found PDB file: $pdb_file"
    else
        echo "No PDB file found in $dir"
    fi

    echo "The protein that is being fed into proteinmpnn is: $dir or $pdb_file"

    path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
    path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
    path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
    chains_to_design=$(scripts/longest_chain.py $pdb_file)


    #The first amino acid in the chain corresponds to 1 and not PDB residues index for now.

    if [ "$determine_CDRs" = "structure" ]; then 
        IFS=" " read design_only_positions <<< $(scripts/find_loops.py $pdb_file --output $to_design)
    fi

    #if [ "$determine_CDRs" = "kabat" ] || [ "$determine_CDRs" = "chothia" ] || [ "$determine_CDRs" = "martin" ]; then
#	    IFS=" " read design_only_positions <<< $(scripts/loops_from_sequence.py $pdb_file --scheme $determine_CDRs --output $to_design --dist $ss_near_CDRs --linker-seq $linker_seq --simple_grab $simple_grab)
#    fi

    if [ "$determine_CDRs" = "kabat" ] || [ "$determine_CDRs" = "chothia" ] || [ "$determine_CDRs" = "martin" ]; then
        cmd="scripts/loops_from_sequence.py $pdb_file --scheme $determine_CDRs --output $to_design --dist $ss_near_CDRs --simple_grab $simple_grab"

        if [ -n "$linker_seq" ]; then
            cmd+=" --linker-seq $linker_seq"
        fi

        IFS=" " read -r design_only_positions <<< $($cmd)
    fi


    python $PMPNN/helper_scripts/parse_multiple_chains.py --input_path=$dir --output_path=$path_for_parsed_chains
    
    python $PMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
    
    python $PMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed
    
    python $PMPNN/protein_mpnn_run.py \
            --jsonl_path $path_for_parsed_chains \
            --chain_id_jsonl $path_for_assigned_chains \
            --fixed_positions_jsonl $path_for_fixed_positions \
            --out_folder $output_dir \
            --num_seq_per_target $seqs_per_run \
            --sampling_temp "0.1" \ 
            --seed 37 \
            --batch_size 1 \
	    --use_soluble_model 

    # make fasta files of just the WT and the ProteinMPNN sequences for folding purposes.
    trimmed=$(basename "$dir")
    scripts/simplify_fasta.py $output_dir/seqs/$trimmed.fa 
    tail -n +3 $trimmed.fa > $dir/${trimmed}_noWT.fa
    head -n +2 $trimmed.fa > $dir/${trimmed}_WT.fa


    colabfold_batch --templates --num-recycle 1 --cache-mmseq-results $dir/${trimmed}.pkl $dir/${trimmed}_WT.fa $dir/structures
    colabfold_batch --templates --num-recycle 6 --model-type alphafold2_multimer_v2 --use-cached-mmseq-results $dir/${trimmed}.pkl $dir/${trimmed}_noWT.fa $dir/structures

    mv $trimmed.fa $dir
    scripts/make_logo_cli.py $dir/$trimmed.fa --output $dir/${trimmed}_sequence_logo.png
    scripts/analyze_structures.py $dir --prefix $trimmed
    mv $dir $output_dir/$trimmed

done

