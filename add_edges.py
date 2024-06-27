"""_summary_

script to add random functional groups to the edges of a 
selected graphene nanoribbon/wire/quantum dot (1D/0D).

returns:
    coordinate file of edge functionalised graphene structure  
    with randomly selected functional groups added on random positions. 
    script takes care of close-contacts.
    
    functional groups are either defined locally or may be imported from an external file.

    some of the functional groups defined locally are: CO2H, CHO and OH

    one may be interested in generating functional groups using simple python functions,
    RDkit or import from existing databases for further applications.

    another interesting usecase could be extending the workflow to PAHs from chemistry databases.

    some applications: bioconjugation, amino acids and polymer attachments.

"""

from math import sin, cos, pi
import itertools

import numpy as np
import pandas as pd

from ase import Atoms
from ase import neighborlist
# from ase import Atom
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


def functional_groups():
    """_summary_
    Generate some of the general functional groups

    Returns:
        _type_: _description_
    """

    # carboxyl group centered at (0,0,0)
    carboxyl = Atoms('COOH', positions=[(0,0,0), (0,1,0), (0,0,1), (1,0,0)])
    carboxyl.set_positions([(0,0,0), (sin(5*pi/6)*1.21,0 ,cos(5*pi/6)*1.21), (sin(pi/6)*1.30,0,cos(pi/6)*1.30), (sin(pi/6)*1.30 + sin(pi/2)*0.96,0,cos(pi/6)*1.30 + cos(pi/2)*0.96)])
    
    # aldehyde group centered at (0,0,0) 
    aldehyde = Atoms('CHO', positions=[(0,0,0), (0,1,0), (0,0,1)])
    aldehyde.set_positions([(0,0,0), (sin(5*pi/6)*1.09,0 ,cos(5*pi/6)*1.09), (sin(pi/6)*1.20,0,cos(pi/6)*1.20)])
    
    # hydroxyl group centered at (0,0,0)
    hydroxyl = Atoms('OH', positions=[(0,0,0), (0,1,0)])
    hydroxyl.set_positions([(0,0,0), (sin(pi/6)*0.96,0 ,cos(pi/6)*0.96)])

    funct_groups = (carboxyl, aldehyde, hydroxyl)
    
    return funct_groups


functional_groups = functional_groups()


def functional_group_position(functional_group, functional_group_pos_init, angle, H_position):
    """_Summary_
    
    Add rotated functional group to ribbon


    Args:
    
        functional_group: ase.Atoms
            functional group added to graphene nanoriibon
        functional_group_pos_init : numpy array
            Positions of functional group if origin was (0,0,0)
            Indices of atoms that can collide with functional group
        angle:
            Angle by which to rotate group
        h_pos:
            Position of removed H atom
    """
    
    # Generate group centered at (0,0,0) with a random orientation
    functional_group.set_positions(functional_group_pos_init)
    functional_group.rotate(angle, 'x') # along 'x' axis

    # Connect group to graphene network
    functional_group_pos = functional_group.get_positions()
    functional_group_pos +=  np.broadcast_to(H_position, (len(functional_group), 3))
    
    functional_group.set_positions(functional_group_pos)


    return functional_group



def updated_nearest_atoms_list(graphene, cutoff_distance = 2.85):
    i, j = neighborlist.neighbor_list(
        "i" "j", graphene, cutoff=cutoff_distance, self_interaction=False
    )
    nearest_atoms_list = np.array((i, j)).T
    
    return nearest_atoms_list



def edge_functionalisation(graphene, edge_FG_atoms, maximum_iterations = 50):
    """_summary_
    
    Substitute randomly selected H atoms of graphene nanoribbon 
    by randomly selected functional group from a list of options.
     

    Args:
        graphene : ASE atoms object
            graphene nanoribbon with H on edges
        edge_FG_atoms : list
            list of indices of atoms of functional groups already added to the edge
        maximum_iterations : int
            maximum attempts to add functional group without any close contacts
            as defined by the pairwise_distance_threshold_dict
    """
    
    box_size = graphene.get_cell().diagonal()
        
    iterations = 0
    
    close_contact = False
    
    original_graphene_atoms = len(graphene)
    
    nearest_atoms_list = updated_nearest_atoms_list(graphene, cutoff_distance = 1.85)
    CH_carbon_atoms = updated_nearest_atoms_list(graphene, cutoff_distance = 1.20)
    
    
    functional_group = functional_groups[np.random.randint(0, len(functional_groups))].copy()
    
    if functional_group[0].symbol == 'C':
        distance = 1.60
    elif functional_group[0].symbol == 'O':
        distance = 1.50
        
    C_H_bond_length = 1.09

    indices_H = [i for i, x in enumerate(graphene.numbers) if x == 1 and i not in edge_FG_atoms]
    
    
    if len(indices_H) == 0:
        print(f"no more edge H atoms available for functionalisation!!")
        return graphene, edge_FG_atoms
    
    else:    
        selected_H = np.random.choice(indices_H)
    
    H_position = np.array(graphene[selected_H].position)
    
    if H_position[0] <= box_size[0]/2:
        orientation = -1
    
    else:
        orientation = 1      
    
    functional_group_pos_init = functional_group.get_positions()
    functional_group_pos_init *= orientation
    
    # Move the selected functional group along X axis 
    # to ensure realistic bond length than original C-H bond. 
    # Bond angle of 120 degree corresponding to sp2 hybridisation is used.
    H_position[0] += orientation * cos((1/6)*pi) * (distance - C_H_bond_length)
    
    close_CH_carbon = CH_carbon_atoms[np.where(CH_carbon_atoms[:, 0] == selected_H)].copy()
    
    H_neighbors_list = nearest_atoms_list[np.where(nearest_atoms_list[:, 0] == selected_H)]
    H_neighbors_updated = [x for x in H_neighbors_list[:,1]].copy()
    
    # Assign a random orientation to add the functional group.
    angle = np.random.uniform(0, 180)
    
    while iterations < maximum_iterations:
        
        print("iterations = ", iterations)    
        
        
        # neighbors of selected edge H and previously added functional group atoms
        # remove alpha-carbon from the list as it is already bonded within close contact distance
        collision_area = H_neighbors_updated + [x for x in edge_FG_atoms]
        collision_area = [x for x in collision_area if x not in close_CH_carbon]
        
        # during each iteration, functional group position is updated with a new angle to avoid close contacts
        functional_group = functional_group_position(functional_group, functional_group_pos_init, angle, H_position)
        
        # append selected functional group with the initial atom positions and orientation
        graphene.extend(functional_group)
        
        # newly added atoms
        functional_group_atoms = np.arange(original_graphene_atoms,len(graphene))
        
        for atom in functional_group_atoms:

            if len(H_neighbors_updated) == 0:
                print("empty neighbors!!")
        
            else:
                
                # check if any newly added atom is in close contact with any atoms in collision area.
                for i in collision_area:
                    
                    if  graphene.get_distance(
                        i, atom, mic=True
                        ) < (pairwise_distance_threshold_dict[graphene[i].symbol + graphene[atom].symbol]):
                        
                        close_contact = True
                        
                        break
                    
                    else:
                        close_contact = False
                        
                        pass
                        

            if close_contact:
                
                # initiate a new iteration with updated angle 
                # break iterations
                angle += np.random.uniform(5, 15)
                
                del graphene[original_graphene_atoms:len(graphene)]
                
                iterations += 1

                break
                
            else:
                # continue iterations for all atoms
                close_contact = False
                
                pass
            
        
        if close_contact:
            # enter new iterations.
            
            pass
            
        else:
            # if no close contacts for all atoms, break iterations.
            
            break
        
    
                    
    if close_contact:
        # all iterations finished
        print("can't add functional group!!")
        
        return graphene, edge_FG_atoms

    else:
        
        del graphene[selected_H]
        
        edge_FG_atoms = [i if i < selected_H else i - 1 for i in edge_FG_atoms]
        new_fg_atoms = np.arange(original_graphene_atoms - 1, len(graphene)).tolist()
        edge_FG_atoms.extend(new_fg_atoms)
                

        print("added functional group")
        
        
        return graphene, edge_FG_atoms
    



    
    
    
    