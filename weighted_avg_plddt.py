#!/usr/bin/env python

import sys
from Bio.PDB import PDBParser, NeighborSearch, Selection
from collections import defaultdict

def calculate_bfactor_for_contact_residues_with_zeros(pdb_file, chain_a_id='A', chain_b_id='B', distance_threshold=5.0):
    # Initialize the PDB parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    # Extract chains A and B
    chain_a = None
    chain_b = None
    for model in structure:
        chain_a = model[chain_a_id]
        chain_b = model[chain_b_id]
        break  # Assuming only one model is needed
    
    # Identify contact residues in chain B (within distance threshold to any residue in chain A)
    atoms_in_chain_a = Selection.unfold_entities(chain_a, 'A')  # List of all atoms in chain A
    atoms_in_chain_b = Selection.unfold_entities(chain_b, 'A')  # List of all atoms in chain B

    # Using NeighborSearch to find contacts
    ns = NeighborSearch(atoms_in_chain_a)
    contact_residues_in_chain_b = set()
    
    for atom in atoms_in_chain_b:
        neighbors = ns.search(atom.coord, distance_threshold, 'R')  # 'R' for residues
        if neighbors:
            contact_residues_in_chain_b.add(atom.get_parent().get_id())

    # To store B-factors for each residue in chain B (including zeros for non-contact residues)
    residue_bfactors = defaultdict(list)

    # Loop over all residues in chain B
    for residue in chain_b:
        residue_id = residue.get_id()
        if residue_id in contact_residues_in_chain_b:
            # Get B-factors for all atoms in the residue
            for atom in residue:
                residue_bfactors[residue_id].append(atom.get_bfactor())
        else:
            # If not in contact, add a zero for the residue
            residue_bfactors[residue_id].append(0)

    # Calculate average B-factor per residue (averaging over all atoms within a residue)
    avg_bfactors_per_residue = {res_id: sum(bfactors)/len(bfactors) 
                                for res_id, bfactors in residue_bfactors.items()}
    
    # Calculate the overall average (Direct average of all atoms in contact and non-contact residues)
    all_bfactors = [bfactor for bfactors in residue_bfactors.values() for bfactor in bfactors]
    overall_avg_bfactor = sum(all_bfactors) / len(all_bfactors) if all_bfactors else 0
    
    # Calculate the average of residue averages
    avg_bfactor_from_residues = sum(avg_bfactors_per_residue.values()) / len(avg_bfactors_per_residue) if avg_bfactors_per_residue else 0
    
    return overall_avg_bfactor, avg_bfactor_from_residues

# Example usage
if __name__=='__main__':
 pdb_file = sys.argv[1]
 chain_a_id = sys.argv[2]
 chain_b_id = sys.argv[3]
 overall_avg, avg_residue = calculate_bfactor_for_contact_residues_with_zeros(pdb_file, chain_a_id, chain_b_id)
 print(f"Overall weighted average B-factor for chain B: {overall_avg}")
 print(f"Weighted average of residue average B-factors for chain B: {avg_residue}")
