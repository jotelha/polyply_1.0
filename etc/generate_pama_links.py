# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Backbone links

# %% [markdown]
# We generate links along the polymer backbone.

# %%
from vermouth.pdb.pdb import read_pdb
import networkx as nx
import matplotlib.pyplot as plt

import itertools


# %%
# --- 1. Find the atom index by atom name (optionally with residue filter) ---
def find_atom_index(mol, atomname, resname=None, resid=None):
    for idx, atom in mol.nodes(data=True):
        if atom['atomname'] == atomname:
            if (resname is None or atom['resname'] == resname) and \
               (resid   is None or atom['resid'] == resid):
                return idx
    raise ValueError("Atom not found.")


# %%
monomer_set = ('EH', 'LA', 'MMA', 'OC', 'ST')   # fixed the unmatched quote

# %%
monomer_pdb_dict = {
    'EH': './PAMA-monomers/eh.pdb',
    'LA': './PAMA-monomers/la.pdb',
    'MMA': './PAMA-monomers/mma.pdb',
    'OC': './PAMA-monomers/oc.pdb',
    'ST': './PAMA-monomers/st.pdb'
}

# %%
monomer_molecule_dict = {monomer: read_pdb(filename)[0] for monomer, filename in monomer_pdb_dict.items()}

# %%
Calpha = {
    'EH': 'C04',
    'LA': 'C04',
    'MMA': 'C07',
    'OC': 'C03',
    'ST': 'C04'
}

# %%
Cbeta = {
    'EH': 'C05',
    'LA': 'C05',
    'MMA': 'C08',
    'OC': 'C01',
    'ST': 'C05'
}

# %%
backbone_hydrocarbons = {
    'alpha': Calpha,
    'beta': Cbeta
}

# %%
terminal_set = ('CH3','CH3n')

# %%
terminal_pdb_dict = {
    'CH3': './terminals/CH3.pdb',
    'CH3n': './terminals/CH3.pdb'
}

# %%
terminal_molecule_dict = {terminal: read_pdb(filename)[0] for terminal, filename in terminal_pdb_dict.items()}

# %%
Cterminal = {
    'CH3': 'C1',
    'CH3n': 'C1',
}

# %%
lines = []

# %%
output_ff = 'PAMA.oplsaa.LigParGen_links.ff'

# %% [markdown]
# ### Harmonic bond potential:
#
# $V_a(\theta) = \frac{1}{2} k_{\theta}(\theta-\theta_0)^2$

# %%
bond_type_parameter_dict = {}

# %%
# harmonic bond potential: V_b(r) = 1/2*k_b(r-b_0)^2
# default CT-CT bond from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 192
default_bond_func_type = 1  # harmonic bond
default_bond_b_0 = 0.15290  # b_0 (nm)
default_bond_k_b = 224262.4 # k_b (kJ mol^-1 nm^-2)

bond_type_parameter_dict[('C','C')] = {
    'func_type': default_bond_func_type,
    'b_0': default_bond_b_0,
    'k_0': default_bond_k_b,
    'comment': "CT-CT bond from 'oplsaa.ff/ffbonded.itp', line 192"
}

# %% [markdown]
# ### Harmonic angle potential:
#
# $V_b(r) = \frac{1}{2} k_{b}(r-b_0)^2$

# %%
angle_type_parameter_dict = {}

# %%
# default CT-CT-CT angle from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 873
# extracted with
#    cat -n ffbonded.itp | grep -E 'CT[[:space:]]+CT[[:space:]]+CT'
default_angle_func_type = 1  # harmonic angle
default_angle_theta_0 = 112.700 # theta_0 (deg)
default_angle_k_theta = 488.273 # k_theta (kJ mol^-1 rad^-2)

angle_type_parameter_dict[('C','C','C')] = {
    'func_type': default_angle_func_type,
    'theta_0': default_angle_theta_0,
    'k_theta': default_angle_k_theta,
    'comment': "CT-CT-CT angle from 'oplsaa.ff/ffbonded.itp', line 873"
}

# %%
# default CT-CT-HC angle from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 881
#    881	  CT     CT     HC      1   110.700    313.800   ; CHARMM 22 parameter file
angle_type_parameter_dict[('C','C','H')] = {
    'func_type': 1,
    'theta_0': 110.700,
    'k_theta': 313.800,
    'comment': "CT-CT-HC angle from 'oplsaa.ff/ffbonded.itp', line 881"
}

# %% [markdown]
# ### Ryckaert-Bellmans dihedral potential:
#
# $ V_{rb}(\phi) = \sum_{n=0}^5 C_n (\cos(\psi))^n$ with $\psi = \phi - 180\degree$
#
# OPLS-AA uses only four coefficients $n=0\dots3$

# %%
dihedral_type_parameter_dict = {}

# %%
# default CT-CT-CT-CT dihedral from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 1596
# extracted with
#    cat -n ffbonded.itp | grep -E 'CT[[:space:]]+CT[[:space:]]+CT[[:space:]]+CT'
# 1596	  CT     CT     CT     CT      3      2.92880  -1.46440   0.20920  -1.67360   0.00000   0.00000 ; hydrocarbon all-atom
default_dihedral_func_type = 3  # Ryckaert-Bellemans dihedral
default_dihedral_C = [2.92880, -1.46440, 0.20920, -1.67360, 0.00000, 0.00000] # C_n (kJ mol^-1)

dihedral_type_parameter_dict[('C','C','C','C')] = {
    'func_type': default_dihedral_func_type,
    'C': default_dihedral_C,
    'comment': "CT-CT-CT-CT dihedral from 'oplsaa.ff/ffbonded.itp', line 1596"
}

# %%
# CT-CT-CT-OS dihedral from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 1608
#  1608	  CT     CT     CT     OS      3      2.87441   0.58158   2.09200  -5.54799   0.00000   0.00000 ; alcohols, ethers AA
dihedral_type_parameter_dict[('C','C','C','O')] = {
    'func_type': 3,
    'C': [2.87441, 0.58158, 2.09200, -5.54799, 0.00000, 0.00000],
    'comment': "CT-CT-CT-OS dihedral from 'oplsaa.ff/ffbonded.itp', line 1608"
}

# %%
# CT-CT-CT-HC dihedral from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 1599
#  1599	  CT     CT     CT     HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ; hydrocarbon all-atom
dihedral_type_parameter_dict[('C','C','C','H')] = {
    'func_type': 3,
    'C': [0.62760, 1.88280, 0.00000, -2.51040, 0.00000, 0.00000],
    'comment': "CT-CT-CT-HC dihedral from 'oplsaa.ff/ffbonded.itp', line 1599"
}

# %%
# HC-CT-CT-HC dihedral from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 1824
#  1824	  HC     CT     CT     HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ; hydrocarbon *new* 11/99
dihedral_type_parameter_dict[('H','C','C','H')] = {
    'func_type': 3,
    'C': [0.62760, 1.88280, 0.00000, -2.51040, 0.00000, 0.00000],
    'comment': "HC-CT-CT-HC dihedral from 'oplsaa.ff/ffbonded.itp', line 1824"
}

# %%
for dihedral_type, parameter_dict in dihedral_type_parameter_dict.items():
    print("{}: {}".format(dihedral_type, sum(parameter_dict["C"])))

# %% [markdown]
# ### visualize connectivity graphs

# %%
for monomer_name, mol in monomer_molecule_dict.items():
    # ---------------------------------------------------------
    # CONFIG: atom names you want to emphasize
    # ---------------------------------------------------------
    highlight_atomnames = {Calpha[monomer_name], Cbeta[monomer_name]}
    
    # ---------------------------------------------------------
    # 1. Create node labels (atom name or atomname+resid)
    # ---------------------------------------------------------
    labels = {
        idx: atom["atomname"]
        for idx, atom in mol.nodes(data=True)
    }
    
    # ---------------------------------------------------------
    # 2. Determine which nodes should be emphasized
    # ---------------------------------------------------------
    highlight_nodes = [
        idx for idx, atom in mol.nodes(data=True)
        if atom["atomname"] in highlight_atomnames
    ]
    
    # ---------------------------------------------------------
    # 3. Visual styling
    # ---------------------------------------------------------
    node_colors = [
        "red" if idx in highlight_nodes else "lightgray"
        for idx in mol.nodes()
    ]
    
    node_sizes = [
        500 if idx in highlight_nodes else 200
        for idx in mol.nodes()
    ]
    
    # ---------------------------------------------------------
    # 4. Layout (graphviz-like spring layout)
    # ---------------------------------------------------------
    pos = nx.spring_layout(mol, seed=42)   # deterministic layout
    
    # ---------------------------------------------------------
    # 5. Draw the graph
    # ---------------------------------------------------------
    plt.figure(figsize=(10, 10))
    
    nx.draw_networkx_edges(mol, pos, alpha=0.4)
    
    nx.draw_networkx_nodes(
        mol, pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="black"
    )
    
    nx.draw_networkx_labels(
        mol, pos,
        labels=labels,
        font_size=8,
        font_color="black"
    )

    plt.title(monomer_name)
    plt.axis("off")
    plt.tight_layout()
    plt.show()


# %% [markdown]
# ### evaluate first and second neighbors of backbone hydrocarbons

# %%
# iterate over Calpha and Cbeta
connectivity_dict = {}
for backbone_hydrocarbon_type, backbone_hydrocarbon_dict in backbone_hydrocarbons.items():
    print(f"Treating C{backbone_hydrocarbon_type}...")

    connectivity_dict[backbone_hydrocarbon_type] = {}
    
    #iterate over all monomers
    for monomer_name, backbone_hydrocarbon_name in backbone_hydrocarbon_dict.items():
        print(f"    Treating {monomer_name}:{backbone_hydrocarbon_name}...")

        connectivity_dict[backbone_hydrocarbon_type][monomer_name] = {}

        mol = monomer_molecule_dict[monomer_name]

        start = find_atom_index(mol, atomname=backbone_hydrocarbon_name)
        
        first_neighbors = list(mol.neighbors(start))
        
        neighbor_index_map = {}

        for nbr in first_neighbors:
            # Second neighbors: neighbors of this first neighbor,
            # excluding the start atom (so only true 2-step nodes)
            second = [x for x in mol.neighbors(nbr) if x != start]
            neighbor_index_map[nbr] = second

        neighbor_name_map = {mol.nodes[idx1]['atomname']: [mol.nodes[idx2]['atomname'] for idx2 in second_neighbors] for idx1, second_neighbors in neighbor_index_map.items()}

        print(f"        Identified first and second neighbors of {monomer_name}:{backbone_hydrocarbon_name}:{neighbor_name_map}.")

        connectivity_dict[backbone_hydrocarbon_type][monomer_name] = neighbor_name_map


# %%
connectivity_dict

# %%
# iterate over Cterminal
terminal_connectivity_dict = {}

#iterate over all monomers
for terminal_name, backbone_hydrocarbon_name in Cterminal.items():
    print(f"    Treating {terminal_name}:{backbone_hydrocarbon_name}...")

    terminal_connectivity_dict[terminal_name] = {}

    mol = terminal_molecule_dict[terminal_name]

    start = find_atom_index(mol, atomname=backbone_hydrocarbon_name)
    
    first_neighbors = list(mol.neighbors(start))
    
    neighbor_index_map = {}

    for nbr in first_neighbors:
        # Second neighbors: neighbors of this first neighbor,
        # excluding the start atom (so only true 2-step nodes)
        second = [x for x in mol.neighbors(nbr) if x != start]
        neighbor_index_map[nbr] = second

    neighbor_name_map = {mol.nodes[idx1]['atomname']: [mol.nodes[idx2]['atomname'] for idx2 in second_neighbors] for idx1, second_neighbors in neighbor_index_map.items()}

    print(f"        Identified first and second neighbors of {terminal_name}:{backbone_hydrocarbon_name}:{neighbor_name_map}.")

    terminal_connectivity_dict[terminal_name] = neighbor_name_map


# %%
terminal_connectivity_dict

# %%
monomer_2_tuples = list(itertools.product(monomer_set, repeat=2))

# %%
monomer_3_tuples = list(itertools.product(monomer_set, repeat=3))

# %%
terminal_monomer_2_tuples = list(itertools.product(terminal_set, monomer_set))

# %%
monomer_3_tuples = list(itertools.product(monomer_set, repeat=3))

# %%
terminal_monomer_monomer_3_tuples = list(itertools.product(terminal_set, monomer_set, monomer_set))

# %% [markdown]
# # Bonds

# %% [markdown]
# ## bonds between residues $i$ and $i+1$

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta)_{i}$

# %% [markdown]
# For bonds, we write entries like this:
#
# ```
# ; EH-headed links
# [ link ]
# ; EH-EH C^beta - C^alpha link
# [ bonds ]
# +C05 {"resname": "EH"} C04 {"resname": "EH"} 1 0.15290 224262.4
# ```

# %%
bond_tuple_list = []

# %%
for tail_monomer, head_monomer in monomer_2_tuples:
    
    # atom name, res name, res increment:
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)
    
    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, 0)
    
    bond_tuple_list.append((left_atom_tuple, right_atom_tuple))

# %%
len(bond_tuple_list)

# %% [markdown]
# ## bonds between residue $i$ and terminal $i+1$

# %% [markdown]
# ### $(T)_{i+1} - (C^\beta)_{i}$

# %%
for terminal, tail_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    left_atom_name = Cterminal[terminal]
    left_atom_tuple = (left_atom_name, terminal, 1)
    
    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, 0)
    
    bond_tuple_list.append((left_atom_tuple, right_atom_tuple))

# %%
len(bond_tuple_list)

# %% [markdown]
# ## bonds between terminal $i$ and residue $i+1$

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (T)_{i}$

# %%
for terminal, head_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)
    
    right_atom_name = Cterminal[terminal]
    right_atom_tuple = (right_atom_name, terminal, 0)
    
    bond_tuple_list.append((left_atom_tuple, right_atom_tuple))

# %%
len(bond_tuple_list)

# %% [markdown]
# ## Map bond types

# %%
bond_type_dict = {}
bond_type_set = set()
bond_type_count = {}

for bond_tuple in bond_tuple_list:
    bond_type = (bond_tuple[0][0][0], bond_tuple[1][0][0])
    
    # canonical ordering, reverse if necessary
    if bond_type[0] > bond_type[1]:
        bond_type = (bond_type[1], bond_type[0])

    bond_type_dict[bond_tuple] = bond_type

    if bond_type not in bond_type_count:
        bond_type_count[bond_type] = 1
    else:
        bond_type_count[bond_type] += 1

# %%
bond_type_count

# %%
bond_parameter_dict = {bond_tuple: bond_type_parameter_dict[bond_type] for bond_tuple, bond_type in bond_type_dict.items()}

# %% [markdown]
# ## Generate bond links

# %%
lines.append('')
lines.append("; backbone bonds")
for bond_tuple, parameter_dict in bond_parameter_dict.items():
    atom_tokens = []
    for atom_tuple in bond_tuple:
        if atom_tuple[2] > 0:
            prefix = '+'
        elif atom_tuple[2] < 0:
            prefix = '-'
        else:
            prefix = ''
            
        atom_tokens.append(f'{prefix}{atom_tuple[0]} {{"resname": "{atom_tuple[1]}"}}')
    atoms_token = ' '.join(atom_tokens)
    
    parameter_token = f'{parameter_dict["func_type"]:d} {parameter_dict["b_0"]:f} {parameter_dict["k_0"]:f}'

    comment_token = ' ; {}'.format(parameter_dict["comment"]) if 'comment' in parameter_dict else ''

    lines.append('')
    lines.append('[ link ]')
    lines.append('[ bonds ]')
    lines.append(f'{atoms_token} {parameter_token}{comment_token}')

# %% [markdown]
# # Angles

# %% [markdown]
# ## angles between residues $i$ and $i+1$

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta - A)_{i}$

# %% [markdown]
# A is any neighbor of $C^\beta$ in residue i

# %%
angle_tuple_list = []

# %%
for tail_monomer, head_monomer in monomer_2_tuples:
    
    # atom name, res name, res increment:
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)
    
    center_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    center_atom_tuple = (center_atom_name, tail_monomer, 0)
    
    for right_atom_name in connectivity_dict["beta"][tail_monomer].keys():
        right_atom_tuple = (right_atom_name, tail_monomer, 0)
        
        angle_tuple_list.append((left_atom_tuple, center_atom_tuple, right_atom_tuple))

# %% [markdown]
# ### $ (A - C^\alpha)_{i+1} - (C_\beta)_i $

# %% [markdown]
# A is any neighbor of $C^\alpha$ in residue $i+1$

# %%
for tail_monomer, head_monomer in monomer_2_tuples:
    
    # atom name, res name, res increment:
    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, 0)
    
    center_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    center_atom_tuple = (center_atom_name, head_monomer, 1)
    
    for left_atom_name in connectivity_dict["alpha"][head_monomer].keys():
        left_atom_tuple = (left_atom_name, head_monomer, 1)
        angle_tuple_list.append((left_atom_tuple, center_atom_tuple, right_atom_tuple))

# %%
len(angle_tuple_list)

# %% [markdown]
# ## angles between terminal $i$ and residue $i+1$

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (T - A)_{i}$

# %% [markdown]
# A is any neighbor of $T$ in terminal i

# %%
for terminal, head_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)
    
    center_atom_name = Cterminal[terminal]
    center_atom_tuple = (center_atom_name, terminal, 0)
    
    for right_atom_name in terminal_connectivity_dict[terminal].keys():
        right_atom_tuple = (right_atom_name, terminal, 0)
        
        angle_tuple_list.append((left_atom_tuple, center_atom_tuple, right_atom_tuple))

# %%
len(angle_tuple_list)

# %% [markdown]
# ### $ (A - C^\alpha)_{i+1} - (T)_i $

# %% [markdown]
# A is any neighbor of $C^\alpha$ in residue $i+1$

# %%
for terminal, head_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    right_atom_name = Cterminal[terminal]
    right_atom_tuple = (right_atom_name, terminal, 0)
    
    center_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    center_atom_tuple = (center_atom_name, head_monomer, 1)
    
    for left_atom_name in connectivity_dict["alpha"][head_monomer].keys():
        left_atom_tuple = (left_atom_name, head_monomer, 1)
        angle_tuple_list.append((left_atom_tuple, center_atom_tuple, right_atom_tuple))

# %%
len(angle_tuple_list)

# %% [markdown]
# ## angles between residue $i$ and terminal $i+1$

# %% [markdown]
# ### $(T)_{i+1} - (C^\beta - A)_{i}$

# %% [markdown]
# A is any neighbor of $C^\beta$ in residue i

# %%
for terminal, tail_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    left_atom_name = Cterminal[terminal]
    left_atom_tuple = (left_atom_name, terminal, 1)
    
    center_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    center_atom_tuple = (center_atom_name, tail_monomer, 0)
    
    for right_atom_name in connectivity_dict["beta"][tail_monomer].keys():
        right_atom_tuple = (right_atom_name, tail_monomer, 0)
        
        angle_tuple_list.append((left_atom_tuple, center_atom_tuple, right_atom_tuple))

# %%
len(angle_tuple_list)

# %% [markdown]
# ### $ (A - T)_{i+1} - (C_\beta)_i $

# %% [markdown]
# A is any neighbor of $T$ in terminal $i+1$

# %%
for terminal, tail_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, 0)
    
    center_atom_name = Cterminal[terminal]
    center_atom_tuple = (center_atom_name, terminal, 1)
    
    for left_atom_name in terminal_connectivity_dict[terminal].keys():
        left_atom_tuple = (left_atom_name, terminal, 1)
        angle_tuple_list.append((left_atom_tuple, center_atom_tuple, right_atom_tuple))

# %%
len(angle_tuple_list)

# %% [markdown]
# ## Map angle types

# %%
angle_type_dict = {}
angle_type_set = set()
angle_type_count = {}

for angle_tuple in angle_tuple_list:
    angle_type = (angle_tuple[0][0][0], angle_tuple[1][0][0], angle_tuple[2][0][0])
    
    # canonical ordering, reverse if necessary
    if angle_type[0] > angle_type[2]:
        angle_type = (angle_type[2], angle_type[1], angle_type[0])

    angle_type_dict[angle_tuple] = angle_type

    if angle_type not in angle_type_count:
        angle_type_count[angle_type] = 1
    else:
        angle_type_count[angle_type] += 1

# %%
angle_type_count

# %%
angle_parameter_dict = {angle_tuple: angle_type_parameter_dict[angle_type] for angle_tuple, angle_type in angle_type_dict.items()}

# %% [markdown]
# ## Generate angle links

# %%
lines.append('')
lines.append('')
lines.append("; backbone angles")
for angle_tuple, parameter_dict in angle_parameter_dict.items():
    atom_tokens = []
    for atom_tuple in angle_tuple:
        if atom_tuple[2] > 0:
            prefix = '+'
        elif atom_tuple[2] < 0:
            prefix = '-'
        else:
            prefix = ''
            
        atom_tokens.append(f'{prefix}{atom_tuple[0]} {{"resname": "{atom_tuple[1]}"}}')
    atoms_token = ' '.join(atom_tokens)
    
    parameter_token = f'{parameter_dict["func_type"]:d} {parameter_dict["theta_0"]:f} {parameter_dict["k_theta"]:f}'

    comment_token = ' ; {}'.format(parameter_dict["comment"]) if 'comment' in parameter_dict else ''

    lines.append('')
    lines.append('[ link ]')
    lines.append('[ angles ]')
    lines.append(f'{atoms_token} {parameter_token}{comment_token}')

# %% [markdown]
# # Dihedrals

# %% [markdown]
# ## dihedrals between residues $i$ and $i+1$

# %% [markdown]
# ### $(A - C^\alpha)_{i+1} - (C^\beta - B)_{i}$

# %% [markdown]
# A is any neighbor of $C^\alpha$ in residue i+1, B is any neighbor of $C^\beta$ in residue i

# %%
dihedral_tuple_list = []

# %%
for tail_monomer, head_monomer in monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    center_left_atom_tuple = (center_left_atom_name, head_monomer, 1)
    
    center_right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    center_right_atom_tuple = (center_right_atom_name, tail_monomer, 0)

    for left_atom_name in connectivity_dict["alpha"][head_monomer].keys():
        left_atom_tuple = (left_atom_name, head_monomer, 1)
        
        for right_atom_name in connectivity_dict["beta"][tail_monomer].keys():
            right_atom_tuple = (right_atom_name, tail_monomer, 0)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(B - A - C^\alpha)_{i+1} - (C^\beta)_{i}$

# %% [markdown]
# A is any neighbor of $C^\alpha$ in residue i+1, B is any neighbor of A in residue i+1

# %%
for tail_monomer, head_monomer in monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_right_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    center_right_atom_tuple = (center_right_atom_name, head_monomer, 1)
    
    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, 0)

    for center_left_atom_name, left_atom_names in connectivity_dict["alpha"][head_monomer].items():
        center_left_atom_tuple = (center_left_atom_name, head_monomer, 1)
        
        for left_atom_name in left_atom_names:
            left_atom_tuple = (left_atom_name, head_monomer, 1)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta - A - B)_{i}$

# %% [markdown]
# A is any neighbor of $C^\beta$ in residue i, B is any neighbor of A in residue i

# %%
for tail_monomer, head_monomer in monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_left_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    center_left_atom_tuple = (center_left_atom_name, tail_monomer, 0)
    
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)

    for center_right_atom_name, right_atom_names in connectivity_dict["beta"][tail_monomer].items():
        center_right_atom_tuple = (center_right_atom_name, tail_monomer, 0)
        
        for right_atom_name in right_atom_names:
            right_atom_tuple = (right_atom_name, tail_monomer, 0)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta - C^\alpha)_{i} - (C^\beta)_{i-1}$

# %%
len(monomer_3_tuples)

# %%
for tail_monomer, center_monomer, head_monomer in monomer_3_tuples:
    
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)
    
    center_left_atom_name = backbone_hydrocarbons["beta"][center_monomer]
    center_left_atom_tuple = (center_left_atom_name, center_monomer, 0)
    
    center_right_atom_name = backbone_hydrocarbons["alpha"][center_monomer]
    center_right_atom_tuple = (center_right_atom_name, center_monomer, 0)

    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, -1)
    
    dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ## dihedrals between terminal $i$ and residue $i+1$

# %% [markdown]
# ### $(A - C^\alpha)_{i+1} - (T - B)_{i}$

# %% [markdown]
# A is any neighbor of $C^\alpha$ in residue i+1, B is any neighbor of $T$ in terminal i

# %%
for terminal, head_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    center_left_atom_tuple = (center_left_atom_name, head_monomer, 1)
    
    center_right_atom_name = Cterminal[terminal]
    center_right_atom_tuple = (center_right_atom_name, terminal, 0)

    for left_atom_name in connectivity_dict["alpha"][head_monomer].keys():
        left_atom_tuple = (left_atom_name, head_monomer, 1)
        
        for right_atom_name in terminal_connectivity_dict[terminal].keys():
            right_atom_tuple = (right_atom_name, terminal, 0)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(B - A - C^\alpha)_{i+1} - (T)_{i}$

# %% [markdown]
# A is any neighbor of $C^\alpha$ in residue i+1, B is any neighbor of A in residue i+1

# %%
for terminal, head_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_right_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    center_right_atom_tuple = (center_right_atom_name, head_monomer, 1)
    
    right_atom_name = Cterminal[terminal]
    right_atom_tuple = (right_atom_name, terminal, 0)

    for center_left_atom_name, left_atom_names in connectivity_dict["alpha"][head_monomer].items():
        center_left_atom_tuple = (center_left_atom_name, head_monomer, 1)
        
        for left_atom_name in left_atom_names:
            left_atom_tuple = (left_atom_name, head_monomer, 1)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (T - A - B)_{i}$

# %% [markdown]
# A is any neighbor of $T$ in terminal i, B is any neighbor of A in terminal i

# %%
for terminal, head_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_left_atom_name = Cterminal[terminal]
    center_left_atom_tuple = (center_left_atom_name, terminal, 0)
    
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)

    for center_right_atom_name, right_atom_names in terminal_connectivity_dict[terminal].items():
        center_right_atom_tuple = (center_right_atom_name, terminal, 0)
        
        for right_atom_name in right_atom_names:
            right_atom_tuple = (right_atom_name, terminal, 0)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta - C^\alpha)_{i} - (T)_{i-1}$

# %%
len(terminal_monomer_monomer_3_tuples)

# %%
for terminal, center_monomer, head_monomer in terminal_monomer_monomer_3_tuples:
    
    left_atom_name = backbone_hydrocarbons["alpha"][head_monomer]
    left_atom_tuple = (left_atom_name, head_monomer, 1)
    
    center_left_atom_name = backbone_hydrocarbons["beta"][center_monomer]
    center_left_atom_tuple = (center_left_atom_name, center_monomer, 0)
    
    center_right_atom_name = backbone_hydrocarbons["alpha"][center_monomer]
    center_right_atom_tuple = (center_right_atom_name, center_monomer, 0)

    right_atom_name = Cterminal[terminal]
    right_atom_tuple = (right_atom_name, terminal, -1)
    
    dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ## dihedrals between residue $i$ and terminal $i+1$

# %% [markdown]
# ### $(A - T)_{i+1} - (C^\beta - B)_{i}$

# %% [markdown]
# A is any neighbor of $T$ in terminal i+1, B is any neighbor of $C^\beta$ in residue i

# %%
for terminal, tail_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_left_atom_name = Cterminal[terminal]
    center_left_atom_tuple = (center_left_atom_name, terminal, 1)
    
    center_right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    center_right_atom_tuple = (center_right_atom_name, tail_monomer, 0)

    for left_atom_name in terminal_connectivity_dict[terminal].keys():
        left_atom_tuple = (left_atom_name, terminal, 1)
        
        for right_atom_name in connectivity_dict["beta"][tail_monomer].keys():
            right_atom_tuple = (right_atom_name, tail_monomer, 0)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(B - A - T)_{i+1} - (C^\beta)_{i}$

# %% [markdown]
# A is any neighbor of $T$ in terminal i+1, B is any neighbor of A in residue i+1

# %%
for terminal, tail_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_right_atom_name = Cterminal[terminal]
    center_right_atom_tuple = (center_right_atom_name, terminal, 1)
    
    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, 0)

    for center_left_atom_name, left_atom_names in terminal_connectivity_dict[terminal].items():
        center_left_atom_tuple = (center_left_atom_name, terminal, 1)
        
        for left_atom_name in left_atom_names:
            left_atom_tuple = (left_atom_name, terminal, 1)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(T)_{i+1} - (C^\beta - A - B)_{i}$

# %% [markdown]
# A is any neighbor of $C^\beta$ in residue i, B is any neighbor of A in residue i

# %%
for terminal, tail_monomer in terminal_monomer_2_tuples:
    
    # atom name, res name, res increment:
    center_left_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    center_left_atom_tuple = (center_left_atom_name, tail_monomer, 0)
    
    left_atom_name = Cterminal[terminal]
    left_atom_tuple = (left_atom_name, terminal, 1)

    for center_right_atom_name, right_atom_names in connectivity_dict["beta"][tail_monomer].items():
        center_right_atom_tuple = (center_right_atom_name, tail_monomer, 0)
        
        for right_atom_name in right_atom_names:
            right_atom_tuple = (right_atom_name, tail_monomer, 0)
            dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ### $(T)_{i+1} - (C^\beta - C^\alpha)_{i} - (C^\beta)_{i-1}$

# %%
len(terminal_monomer_monomer_3_tuples)

# %%
for terminal, tail_monomer, center_monomer in terminal_monomer_monomer_3_tuples:
    
    left_atom_name = Cterminal[terminal]
    left_atom_tuple = (left_atom_name, terminal, 1)
    
    center_left_atom_name = backbone_hydrocarbons["beta"][center_monomer]
    center_left_atom_tuple = (center_left_atom_name, center_monomer, 0)
    
    center_right_atom_name = backbone_hydrocarbons["alpha"][center_monomer]
    center_right_atom_tuple = (center_right_atom_name, center_monomer, 0)

    right_atom_name = backbone_hydrocarbons["beta"][tail_monomer]
    right_atom_tuple = (right_atom_name, tail_monomer, -1)
    
    dihedral_tuple_list.append((left_atom_tuple, center_left_atom_tuple, center_right_atom_tuple, right_atom_tuple))

# %%
len(dihedral_tuple_list)

# %%
len(set(dihedral_tuple_list))

# %% [markdown]
# ## Map dihedral types

# %%
dihedral_type_dict = {}
dihedral_type_set = set()
dihedral_type_count = {}

for dihedral_tuple in dihedral_tuple_list:
    dihedral_type = (dihedral_tuple[0][0][0], dihedral_tuple[1][0][0], dihedral_tuple[2][0][0], dihedral_tuple[3][0][0])
    
    # canonical ordering, reverse if necessary
    if dihedral_type[0] > dihedral_type[3]:
        dihedral_type = (dihedral_type[3], dihedral_type[2], dihedral_type[1], dihedral_type[0])

    dihedral_type_dict[dihedral_tuple] = dihedral_type

    if dihedral_type not in dihedral_type_count:
        dihedral_type_count[dihedral_type] = 1
    else:
        dihedral_type_count[dihedral_type] += 1

# %%
dihedral_type_count

# %%
dihedral_parameter_dict = {dihedral_tuple: dihedral_type_parameter_dict[dihedral_type] for dihedral_tuple, dihedral_type in dihedral_type_dict.items()}

# %% [markdown]
# ## Generate dihedral links

# %%
lines.append('')
lines.append('')
lines.append("; backbone dihedrals")
for dihedral_tuple, parameter_dict in dihedral_parameter_dict.items():
    atom_tokens = []
    for atom_tuple in dihedral_tuple:
        if atom_tuple[2] > 0:
            prefix = '+'
        elif atom_tuple[2] < 0:
            prefix = '-'
        else:
            prefix = ''
            
        atom_tokens.append(f'{prefix}{atom_tuple[0]} {{"resname": "{atom_tuple[1]}"}}')
    atoms_token = ' '.join(atom_tokens)
    
    parameter_token = f'{parameter_dict["func_type"]:d} {parameter_dict["C"][0]:f} {parameter_dict["C"][1]:f} {parameter_dict["C"][2]:f} {parameter_dict["C"][3]:f} {parameter_dict["C"][4]:f} {parameter_dict["C"][5]:f}'

    comment_token = ' ; {}'.format(parameter_dict["comment"]) if 'comment' in parameter_dict else ''

    lines.append('')
    lines.append('[ link ]')
    lines.append('[ dihedrals ]')
    lines.append(f'{atoms_token} {parameter_token}{comment_token}')

# %%
lines[-1]

# %% [markdown]
# ## write links file

# %%
with open(output_ff, "w") as f:
    for line in lines:
        f.write(line + "\n")

# %%
