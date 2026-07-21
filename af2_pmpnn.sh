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
    
done

