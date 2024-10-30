#!/bin/bash

pdbpath='ring_file.pdb'
original_pdbpath='../../pdbs//2dsf.pdb'
sugar_ring_ResID='2.C'
bfmp_path='/home/rajan/apps/BFMP/detect_shape'

# if [[ $# -ne 2 ]]; then
#   echo "Error !!!"
#   echo "USAGE: $0 <pdbpath> <sugar_ring_ResID>"
#   exit 1
# fi

echo "$original_pdbpath $sugar_ring_ResID" >>monosaccharide_ring_conformations.txt
echo "$original_pdbpath $sugar_ring_ResID" >>monosaccharide_ring_conformations.log

sed -i "9s/.*/$sugar_ring_ResID/" pdb_input.txt
$bfmp_path $pdbpath pdb_input.txt

cat ring_conformations.txt >>monosaccharide_ring_conformations.txt
cat ring_conformations.txt >>monosaccharide_ring_conformations.log
