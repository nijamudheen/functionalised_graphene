"""_summary_

script to add halogen atoms to a selected graphene surface.

returns:
    coordinate file of halographene  
    with halogen atoms added on either sides of a graphene surface randomly. 
    script takes care of some local geometry changes
    and considers vdW repulsion as defined below.
"""

import itertools

import numpy as np
import pandas as pd

from ase import Atom
from ase import neighborlist
# from ase.io import read, write
# from ase.build import graphene_nanoribbon



# plausible buckling of carbon due to halogene addition
buckling_halogen_dict = {
    'F' : np.arange(0.25,0.36,0.01),
    'Cl' : np.arange(0.15,0.26,0.01),
    'Br' : np.arange(0.10,0.16,0.01),
    'I' : np.arange(0.05,0.11,0.01)
}


# dictionary of approx. atomic radii
atomic_radii_dict = {
    'H' : 1.20 / 2,
    'C' : 1.85 / 2,
    'O' : 1.52 / 2,
    'F' : 1.60 / 2,
    'Cl' : 1.90 / 2,
    'Br' : 2.0 / 2,
    'I' : 2.2 / 2
}


# Generate all possible pairwise combinations of atomic radii 
vdW_parameters = list(atomic_radii_dict.keys())

# a threshold is defined to avoid steric hindrance when adding atoms
pairwise_distance_threshold_dict = {}

# create new keys for each pairwise combinations of atomic radii
# assign the sum of the values of the parameters in the pair to the new keys
for pair in itertools.product(vdW_parameters, repeat=2):
    key = pair[0] + pair[1]
    value = atomic_radii_dict[pair[0]] + atomic_radii_dict[pair[1]]
    pairwise_distance_threshold_dict[key] = value

pairwise_distance_threshold_df = pd.DataFrame(list(pairwise_distance_threshold_dict.items()), columns=['Parameter', 'Value'])


def updated_nearest_atoms_list(graphene, cutoff=2.85):
    i, j = neighborlist.neighbor_list(
        "i" "j", graphene, cutoff=cutoff, self_interaction=False
    )
    nearest_atoms_list = np.array((i, j)).T
    
    return nearest_atoms_list



def add_halogen(
    graphene, C_atoms, functionalised_carbons, halogen,
    maximum_iterations, close_contact=False
): 
    """_summary_
    
    pick a carbon (non-edge for 0D/1D) atom randomly for halogenation and add requested halogen
    to either sides of it as perpendicular to the basal plane. 
    
    Any added atom should pass a collision threshold to avoid extreme steric interactions.

        Args:
        graphene : ase.Atom
            graphene layer 
        C_atoms : int
            number of C atoms in graphene layer
        functionalised_carbons : list
            list of the indices of the carbon atoms which have been functionalised
        halogen : str
            halogen to add (F, Cl, Br, I)
        maximum iterations : int
            attempts to avoid close contacts during functionalisation
        close_contact : boolean
            checks close contacts

    """
    
    
    iterations = 0
    close_contact=False
    
    # C_choices = list(range(C_atoms))
    C_choices = [c for c in C_atoms if c not in functionalised_carbons]
    
    if len(C_choices) == 0:
        print(f"All carbon atoms have been saturated!")
        return graphene, functionalised_carbons
    
    else:
        while iterations < maximum_iterations:
            
            carbon_atom = np.random.choice(C_choices)
            
            rand = np.random.random()
            
            if rand < 0.5:
                graphene[carbon_atom].position += [0, np.random.choice(buckling_halogen_dict[halogen]), 0]
                halogen_pos = graphene[carbon_atom].position + [0, 2 * atomic_radii_dict[halogen], 0]
                
                graphene.append(Atom(halogen, halogen_pos))
                
                halogen_index = len(graphene)-1
                
                nearest_atoms_list = updated_nearest_atoms_list(graphene)
                halogen_neighbors_list = nearest_atoms_list[np.where(nearest_atoms_list[:, 0] == halogen_index)]
                
                halogen_neighbors_list = [n for n in halogen_neighbors_list[:, 1] if graphene[n].symbol != 'C']
                
                if len(halogen_neighbors_list) == []:
                    pass
                
                else:
                    for i in halogen_neighbors_list:
                        if graphene.get_distance(
                            i, halogen_index, mic=True
                            ) > pairwise_distance_threshold_dict[graphene[i].symbol + graphene[halogen_index].symbol]:
                            functionalised_carbons.append(carbon_atom)
                            pass
                        
                        else:
                            del graphene[halogen_index]
                            graphene[carbon_atom].position += [0, -1 * np.random.choice(buckling_halogen_dict[halogen]), 0]
                            # print("close contacts found!")
                            close_contact=True
                            break
                        
                    
            elif rand >= 0.5:
                graphene[carbon_atom].position += [0, -1 * np.random.choice(buckling_halogen_dict[halogen]), 0]
                halogen_pos = graphene[carbon_atom].position + [0, -1 * 2 * atomic_radii_dict[halogen], 0]
                
                graphene.append(Atom(halogen, halogen_pos))
                
                halogen_index = len(graphene)-1
                
                nearest_atoms_list = updated_nearest_atoms_list(graphene)
                halogen_neighbors_list = nearest_atoms_list[np.where(nearest_atoms_list[:, 0] == halogen_index)]

                halogen_neighbors_list = [n for n in halogen_neighbors_list[:, 1] if graphene[n].symbol != 'C']
                
                
                if len(halogen_neighbors_list) == []:
                    pass
                
                
                else:
                    for i in halogen_neighbors_list:
                        if graphene.get_distance(
                            i, halogen_index, mic=True
                            ) > pairwise_distance_threshold_dict[graphene[i].symbol + graphene[halogen_index].symbol]:
                            functionalised_carbons.append(carbon_atom)
                            pass
                        
                        else:
                            del graphene[halogen_index]
                            graphene[carbon_atom].position += [0, np.random.choice(buckling_halogen_dict[halogen]), 0]
                            # print("close contacts found!")
                            close_contact=True
                            break                            
            
            
            if close_contact:
                iterations += 1
            
            else:
                return graphene, functionalised_carbons
    
    
    print(f"max iterations reached")
    return graphene, functionalised_carbons

