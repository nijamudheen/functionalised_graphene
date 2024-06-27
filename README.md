# functionalised_graphene
this repository contains some python scripts to generate input files 
for running simulations of functionalised graphene materials. 

# dependencies
(tested with python '3.12.1' and ASE - versoin = '3.22.1'.)

# tutorials
notebooks starting with tut* are tutorials.

# sample coordinate files
these are provided in the directory './structure_files'.

# description of codes
chemically functionalised systems such as fluorographene, chlorographene, graphane and 
graphene oxide are some of the classes of highly disordered materials with varying 
atomic structures depending on how they were made, stored and other factors. 

functionalised_graphene code can be used to generate relatively simple 
modified graphene structures and some complex structures. 

these scripts are written in a modular fashion, so they can be easily adapted and modified 
to generate more complex structures, for instance extending to 3D layered structures 
or add some specific domains of pristince, defect, specially designed / functionalised structures.

mainly, the code can be used for high throughput generation of thousands of random 
chemically / physically reasonable structures that are simulations ready. 
these were written with the aim of using them for high-throughput DFT calculations, 
training general purpose MLIPs and/or for running specific molecular dynamics simulations.

code can be used to functionalise 0D, 1D and 2D graphene structures.
both pristine and amorphous graphene structures either generated using 
ase functions or imported from a database can be used.

(Some utility scripts for MLIP training are also included in the repository (to be added!!)).

### Functionalisation strategy
both edges and bulk carbon atoms can be functionalised separately. 

### Functional groups
bulk carbon atoms can be functionalised with 
OH (hydroxyl), O (epoxides, bridge fashion) and / or halogens.
edge carbon atoms can be functionalised from a list of stored functional groups.
some of the groups stored are -OH, CO2H and CHO.
RDkit can be used to generate specific organic groups and used;
for instance, for bioconjugation or polymer attachments to graphene
and to use for subsequent molecular dynamics simulations. 


