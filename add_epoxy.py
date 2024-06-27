# from math import sin, cos, pi

import numpy as np
import pandas as pd
import itertools

from ase import Atom
from ase import neighborlist
# from ase import Atoms
# from ase.io import write
# from ase.build import graphene_nanoribbon


# list of values for possible buckling of carbon due to epoxidation
buckling_epoxy = np.arange(0.10,0.25,0.01)

atomic_radii_dict = {
    'H' : 1.20 / 2.0,
    'C' : 1.85 / 2.0,
    'O' : 1.85 / 2.0,
    'F' : 1.50 / 2.0,
    'Cl' : 1.90 / 2.0,
    'Br' : 2.0 / 2.0,
    'I' : 2.2 / 2.0
}

# Generate all possible pairwise combinations 
vdW_parameters = list(atomic_radii_dict.keys())

pairwise_distance_threshold_dict = {}

# Create new keys for each pairwise combinations 
# assign the sum of the values of the parameters in the pair to the new keys
for pair in itertools.product(vdW_parameters, repeat=2):
    key = pair[0] + pair[1]
    value = atomic_radii_dict[pair[0]] + atomic_radii_dict[pair[1]]
    pairwise_distance_threshold_dict[key] = value

pairwise_distance_threshold_df = pd.DataFrame(list(pairwise_distance_threshold_dict.items()), columns=['Parameter', 'Value'])


def updated_nearest_atoms_list(graphene, cutoff_distance = 2.85):
    i, j = neighborlist.neighbor_list(
        "i" "j", graphene, cutoff=cutoff_distance, self_interaction=False
    )
    nearest_atoms_list = np.array((i, j)).T
    
    return nearest_atoms_list



def select_carbon_epoxy(
    graphene, C_atoms, functionalised_carbons,
    maximum_iterations=100
):
    """_summary_
    
    select a carbon atom and its neighbor randomly for epoxidation.
    if the carbon atom or its neighbor is already functionalised, iterate.
    if no pair of unfunctionalised carbon atoms are found after
    maximum iterations, return none. 

    Args:
        graphene : ase Atoms
            graphene system to functionalise
        C_atoms : int
            number of arbon atoms in graphene layer
        functionalised_carbons : list
            indices of carbon atoms already functionalised
        maximum_iterations : int, default = 100
            maximum iterations to find a pair of unfunctionalised carbon atoms
    """
    
    
    iterations = 0
    
    while iterations < maximum_iterations:
        
        # available_carbon_atoms = list(range(C_atoms))
        available_carbon_atoms = [c for c in C_atoms if
                                  c not in functionalised_carbons]
        
        
        if len(available_carbon_atoms) == 0:
            print(
                f"all carbon atoms already functionalised. can not add epoxy groups"
            )
            
            return graphene, functionalised_carbons, iterations
        
        else:
            carbon_atom = np.random.choice(available_carbon_atoms)
            
            nearest_atoms_list = updated_nearest_atoms_list(graphene, cutoff_distance = 1.85)
            
            carbon_neighbors = nearest_atoms_list[np.where(nearest_atoms_list[:, 0] == carbon_atom)]
            neighbors = [x for x in carbon_neighbors[:, -1] if x not in functionalised_carbons]
            
            if len(neighbors) == 0:
                print(f"All neighbors saturated, start new iteration")
                
                iterations += 1
                
            else:
                neighbor = np.random.choice(neighbors)
                
                functionalised_carbons.append(carbon_atom)
                functionalised_carbons.append(neighbor)
                
                return graphene, functionalised_carbons, carbon_atom, neighbor
            
    
    print(
        f"maximum iterations reached. could not find a pair of carbon atoms for epoxidation"
    )
    
    return graphene, functionalised_carbons, None, None




def add_epoxy(
    graphene, C_atoms, functionalised_carbons,
    maximum_iterations
):
    """_summary_
    
    Add an oxygen in epoxy fashion above or below the plane to selected
    pair of carbon atoms

    Args:
        graphene : ase Atoms
            graphene layer
        C_atoms : int
            number of carbon atoms in graphene
        functionalised_carbons : list
            indices of carbon atoms already functionalised
        buckling_epoxy : list
            list of values for possible buckling of carbon due to epoxidation
        maximum_iterations : int
            maximum attempts to add an epoxide without close contacts. 
            
            (close contacts check not implemented for epoxides, but can be 
            done following the same approach as used for hydroxyls and halogens)
    """
    
    graphene, functionalised_carbons, carbon_atom, neighbor = select_carbon_epoxy(
        graphene, C_atoms, functionalised_carbons,
        maximum_iterations=100
        )
    
    box_size = graphene.get_cell().diagonal()
    
    if carbon_atom == None:
        print(
            f"all carbon atoms already functionalised. can not add epoxide"
        )
        return graphene, functionalised_carbons
    
    elif neighbor == None:
        print(
            "all neighbors are already functionalised. can not add epoxide"
        )
        return graphene, functionalised_carbons
        
    else:
        bond_vector = graphene[neighbor].position - graphene[carbon_atom].position
        

        for i in range(3):
            if abs(bond_vector[i]) > box_size[i] / 2:
                bond_vector[i] -= box_size[i] * round(bond_vector[i] / box_size[i])
        
        midpoint = graphene[carbon_atom].position + bond_vector / 2

        
        rand = np.random.random()
        
        if rand < 0.5:
            
            graphene[carbon_atom].position += [
                0,
                np.random.choice(buckling_epoxy),
                0
            ]
            
            graphene[neighbor].position += [
                0,
                -1 * np.random.choice(buckling_epoxy),
                0
            ]
            
            graphene.append(
                Atom("O", midpoint + [0, 1.26, 0])
            )
            
        elif rand >= 0.5:
            graphene[carbon_atom].position += [
                0,
                -1 * np.random.choice(buckling_epoxy),
                0
            ]
            
            graphene[neighbor].position += [
                0,
                np.random.choice(buckling_epoxy),
                0
            ]
            
            graphene.append(
                Atom("O", midpoint + [0, -1.26, 0])
            )
            
        return graphene, functionalised_carbons


