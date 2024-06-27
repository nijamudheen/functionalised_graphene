import itertools

import numpy as np
import pandas as pd

from math import sin, cos, pi

from ase import Atom
from ase import neighborlist
# from ase import Atoms
# from ase.io import write
# from ase.build import graphene_nanoribbon


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


def add_hydroxyl(
    graphene, C_atoms, functionalised_carbons, maximum_iterations
):
    """_summary_
    
    add a hydroxyl functional group on graphene surface in a random fashion.
    consider local structural changes and chemical information, avoid close contacts.

    Args:
        graphene : ase Atoms
            graphene layer
        C_atoms : list
            carbon atoms in graphene plane
        functionalised_carbons : list
            indices of carbon atoms already functionalised
        maximum_iterations : int
            maximum attempts to add an epoxide without close contacts
    """
    
    box_size = graphene.get_cell().diagonal()
    
    # buckling of carbon due to hydroxyl addition
    buckling_hydroxyl = np.arange(0.10,0.31,0.01)
    
    # available_carbon_atoms = list(range(C_atoms))
    available_carbon_atoms = [c for c in C_atoms if
                                c not in functionalised_carbons]


    if len(available_carbon_atoms) == 0:
        print(
            f"all carbon atoms already functionalised. can not add hydroxy groups"
        )
        
        return graphene, functionalised_carbons
    
    
    else:
        
        rand = np.random.random()
        selected_buckling_oh = np.random.choice(buckling_hydroxyl)
        
        
        iterations_o = 0 # iterations to attempt O addition
        
        while iterations_o < maximum_iterations:
            
            close_contact = False
            
            carbon_atom = np.random.choice(available_carbon_atoms)

            if rand < 0.5:
                graphene[carbon_atom].position += [0, selected_buckling_oh, 0]
                oxygen_position = graphene[carbon_atom].position + [0, 1.47, 0]
            
            elif rand >= 0.5:
                graphene[carbon_atom].position += [0, -1.0 * selected_buckling_oh, 0]
                oxygen_position = graphene[carbon_atom].position + [0, -1.47, 0]
            
            graphene.append(Atom("O", oxygen_position))
            
            oxygen_index = len(graphene)-1    
            
            
            nearest_atoms_list = updated_nearest_atoms_list(graphene, cutoff_distance=1.85).copy()
            oxygen_neighbors_list = nearest_atoms_list[np.where(nearest_atoms_list[:, 0] == oxygen_index)]
            
            oxygen_neighbors_list = [n for n in oxygen_neighbors_list[:, 1] if graphene[n].symbol != 'C']

            if len(oxygen_neighbors_list) == []:
                pass
            
            else:
                for i in oxygen_neighbors_list:
                    if graphene.get_distance(
                        i, oxygen_index, mic=True
                        ) > pairwise_distance_threshold_dict[graphene[i].symbol + graphene[oxygen_index].symbol]:
                        
                        pass
                    
                    else:
                        
                        # delete added oxygen atom
                        del graphene[oxygen_index]
                        
                        # undo buckling of carbon_atom
                        if rand < 0.5:
                            graphene[carbon_atom].position += [0, -1 * selected_buckling_oh, 0]
                            
                        elif rand >= 0.5:
                            graphene[carbon_atom].position += [0, selected_buckling_oh, 0]
                        
                        print("close contacts found!")
                        
                        close_contact = True
                        
                        break
                    
                if close_contact:
                    pass
                
                else:
                    break
                
            if close_contact:
                iterations_o += 1
                
            else:
                break
            
        functionalised_carbons.append(carbon_atom)
        
        print(f"iterations = {iterations_o}")
            
        if iterations_o < maximum_iterations:
            pass
        
        else:
            print(f"can not add hydroxyl groups")
            
            return graphene, functionalised_carbons
        
        
        iterations_oh = 0
        
        while iterations_oh < maximum_iterations:
        
            # add H of OH group with a random orientation
            oh_bond = 0.98
            
            theta_deg = 100.0 + np.random.random() * 15.0
            
            if rand < 0.5:
                theta = np.deg2rad(180-theta_deg)
                
            elif rand >= 0.5:
                theta = np.deg2rad(theta_deg)
            
            # add a random orientation to H of OH group
            # apply spherical to cartesian coordinate system conversion
            phi = np.deg2rad(np.random.random() * 360.0) 
            x = oh_bond * np.sin(theta) * np.cos(phi)            
            y = oh_bond * np.cos(theta)
            z = oh_bond * np.sin(theta) * np.sin(phi)                       
            
            h_position = oxygen_position + [x,y,z]
            
            graphene.append(Atom("H", h_position))
            
            hydrogen_index = len(graphene) - 1
            #nearest_atoms_list = updated_nearest_atoms_list(graphene, cutoff_distance=1.85)
            
            hydrogen_neighbors_list = nearest_atoms_list[np.where(nearest_atoms_list[:, 0] == hydrogen_index)]
            
            # hydrogen_neighbors_list = [n for n in hydrogen_neighbors_list if n != hydrogen_index - 1]
            
            if len(hydrogen_neighbors_list) == []:
                pass
            
            else:
                for i in hydrogen_neighbors_list:
                    if graphene.get_distance(
                        i, hydrogen_index, mic=True
                        ) > pairwise_distance_threshold_dict[graphene[i].symbol + graphene[hydrogen_index].symbol] ** 2.0:
                        
                        pass
                    
                    else:
                        
                        # delete added oxygen atom
                        del graphene[hydrogen_index]
                        del graphene[hydrogen_index - 1]
                        
                        # undo buckling of carbon_atom
                        if rand < 0.5:
                            graphene[carbon_atom].position += [0, -1 * selected_buckling_oh, 0]
                            
                        elif rand >= 0.5:
                            graphene[carbon_atom].position += [0, selected_buckling_oh, 0]
                        
                        print("close contacts found!")
                        
                        close_contact = True
                        
                        break
                    

                
            if close_contact:
                print(f"iterations = {iterations_oh}")
                iterations_oh += 1
                
            else:
                functionalised_carbons.append(carbon_atom)
                return graphene, functionalised_carbons
                
            
        
        
    print(f"can not add hydroxyl groups")
    return graphene, functionalised_carbons
    

    