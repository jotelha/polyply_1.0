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
lines = []

# %%
output_ff = 'PAMA.oplsaa.LigParGen_links.ff'

# %% [markdown]
# Harmonic bond potential:
#
# $V_a(\theta) = \frac{1}{2} k_{\theta}(\theta-\theta_0)^2$

# %%
# harmonic bond potential: V_b(r) = 1/2*k_b(r-b_0)^2
# default CT-CT bond from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 192
default_bond_func_type = 1  # harmonic bond
default_bond_b_0 = 0.15290  # b_0 (nm)
default_bond_k_b = 224262.4 # k_b (kJ mol^-1 nm^-2)

# %% [markdown]
# Harmonic angle potential:
#
# $V_b(r) = \frac{1}{2} k_{b}(r-b_0)^2$

# %%
# default CT-CT-CT angle from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 873
# extracted with
#    cat -n ffbonded.itp | grep -E 'CT[[:space:]]+CT[[:space:]]+CT'
default_angle_func_type = 1  # harmonic angle
default_angle_theta_0 = 112.700 # theta_0 (deg)
default_angle_k_theta = 488.273 # k_theta (kJ mol^-1 rad^-2)

# %% [markdown]
# Ryckaert-Bellmans dihedral potential:
#
# $ V_{rb}(\phi) = \sum_{n=0}^5 C_n (\cos(\psi))^n$ with $\psi = \phi - 180\degree$
#
# OPLS-AA uses only four coefficients $n=0\dots3$

# %%
# default CT-CT-CT-CT dihedral from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 1596
# extracted with
#    cat -n ffbonded.itp | grep -E 'CT[[:space:]]+CT[[:space:]]+CT[[:space:]]+CT'
default_dihedral_func_type = 3  # Ryckaert-Bellemans dihedral
default_dihedral_C = [12.92880, -1.46440, 0.20920, -1.67360, 0.00000, 0.00000] # C_n (kJ mol^-1)

# %% [markdown]
# ## bonds between residues $i$ and $i+1$

# %% [markdown]
# ### $(C^\beta)_{i+1} - (C^\alpha)_{i}$

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
backbone_bond_tuples = []
for tail_monomer, tail_backbone_atom in Calpha.items():
    for head_monomer, head_backbone_atom in Cbeta.items():
        backbone_bond_tuples.append(((tail_monomer, tail_backbone_atom), (head_monomer, head_backbone_atom)))

# %%
len(backbone_bond_tuples)

# %%
backbone_bond_parameter_dict = {
    backbone_bond_tuple: {
        "func_type": default_bond_func_type,
        "b_0": default_bond_b_0,
        "k_b": default_bond_k_b
    } for backbone_bond_tuple in backbone_bond_tuples}

# %%
lines.append("; default CT-CT bond from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 192")
for backbone_bond_tuple, parameter_dict in backbone_bond_parameter_dict.items():
    lines.append('')
    lines.append('[ link ]')
    lines.append(f'; {backbone_bond_tuple[1][0]}-{backbone_bond_tuple[0][0]} (C^beta)_i+1 - (C^alpha)_i link')
    lines.append('[ bonds ]')
    lines.append(f'+{backbone_bond_tuple[1][1]} {{"resname": "{backbone_bond_tuple[1][0]}"}} {backbone_bond_tuple[0][1]} {{"resname": "{backbone_bond_tuple[0][0]}"}} {parameter_dict["func_type"]:d} {parameter_dict["b_0"]:f} {parameter_dict["k_b"]:f}')


# %% [markdown]
# ## angles between residues $i$ and $i+1$

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta - C^\alpha)_{i}$

# %%
backbone_angle_tuples = []

# alpha-beta-alpha
for (tail_monomer, tail_backbone_atom), (center_monomer, center_backbone_atom) in zip(Calpha.items(), Cbeta.items()):
    for head_monomer, head_backbone_atom in Calpha.items():
        backbone_angle_tuples.append(((tail_monomer, tail_backbone_atom), (center_monomer, center_backbone_atom), (head_monomer, head_backbone_atom)))

# %%
len(backbone_angle_tuples)

# %%
backbone_angle_parameter_dict = {
    backbone_angle_tuple: {
        "func_type": default_angle_func_type,
        "theta_0": default_angle_theta_0,
        "k_theta": default_angle_k_theta
    } for backbone_angle_tuple in backbone_angle_tuples}

# %%
lines.append('')
lines.append('')
lines.append("; default CT-CT-CT angle from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 873")
for backbone_angle_tuple, parameter_dict in backbone_angle_parameter_dict.items():
    lines.append('')
    lines.append('[ link ]')
    lines.append(f'; {backbone_angle_tuple[2][0]}-{backbone_angle_tuple[0][0]} (C^alpha)_i+1 - (C^beta - C^alpha)_i link')
    lines.append('[ angles ]')
    lines.append(f'+{backbone_angle_tuple[2][1]} {{"resname": "{backbone_angle_tuple[2][0]}"}} {backbone_angle_tuple[1][1]} {{"resname": "{backbone_angle_tuple[1][0]}"}} {backbone_angle_tuple[0][1]} {{"resname": "{backbone_angle_tuple[0][0]}"}} {parameter_dict["func_type"]:d} {parameter_dict["theta_0"]:f} {parameter_dict["k_theta"]:f}')

# %% [markdown]
# ### $ (C^\beta - C^\alpha)_{i+1} - (C_\beta)_i $

# %%
backbone_angle_tuples = []

# beta-alpha-beta
for tail_monomer, tail_backbone_atom in Cbeta.items():
    for (center_monomer, center_backbone_atom), (head_monomer, head_backbone_atom) in zip(Calpha.items(), Cbeta.items()) :
        backbone_angle_tuples.append(((tail_monomer, tail_backbone_atom), (center_monomer, center_backbone_atom), (head_monomer, head_backbone_atom)))

# %%
len(backbone_angle_tuples)

# %%
backbone_angle_parameter_dict = {
    backbone_angle_tuple: {
        "func_type": default_angle_func_type,
        "theta_0": default_angle_theta_0,
        "k_theta": default_angle_k_theta
    } for backbone_angle_tuple in backbone_angle_tuples}

# %%
for backbone_angle_tuple, parameter_dict in backbone_angle_parameter_dict.items():
    lines.append('')
    lines.append('[ link ]')
    lines.append(f'; {backbone_angle_tuple[2][0]}-{backbone_angle_tuple[0][0]} (C^beta - C^alpha)_i+1 - (C^beta)_i link')
    lines.append('[ angles ]')
    lines.append(f'+{backbone_angle_tuple[2][1]} {{"resname": "{backbone_angle_tuple[2][0]}"}} +{backbone_angle_tuple[1][1]} {{"resname": "{backbone_angle_tuple[1][0]}"}} {backbone_angle_tuple[0][1]} {{"resname": "{backbone_angle_tuple[0][0]}"}} {parameter_dict["func_type"]:d} {parameter_dict["theta_0"]:f} {parameter_dict["k_theta"]:f}')

# %% [markdown]
# ## dihedrals

# %% [markdown]
# ### $(C^\beta - C^\alpha)_{i+1} - (C^\beta - C^\alpha)_{i}$

# %%
backbone_dihedral_tuples = []

# alpha-beta-alpha-alpha
for (tail_monomer, tail_backbone_atom), (center_left_monomer, center_left_backbone_atom) in zip(Calpha.items(), Cbeta.items()):
    for (center_right_monomer, center_right_backbone_atom), (head_monomer, head_backbone_atom) in zip(Calpha.items(), Cbeta.items()) :
        backbone_dihedral_tuples.append(((tail_monomer, tail_backbone_atom), (center_left_monomer, center_left_backbone_atom), (center_right_monomer, center_right_backbone_atom), (head_monomer, head_backbone_atom)))


# %%
len(backbone_dihedral_tuples)

# %%
backbone_dihedral_parameter_dict = {
    backbone_dihedral_tuple: {
        "func_type": default_dihedral_func_type,
        "C": default_dihedral_C,
    } for backbone_dihedral_tuple in backbone_dihedral_tuples}

# %%
lines.append('')
lines.append('')
lines.append("; default CT-CT-CT-CT dihedral from CHARMM 22 parameter file in 'oplsaa.ff/ffbonded.itp', line 1596")
for backbone_dihedral_tuple, parameter_dict in backbone_dihedral_parameter_dict.items():
    lines.append('')
    lines.append('[ link ]')
    lines.append(f'; {backbone_dihedral_tuple[3][0]}-{backbone_dihedral_tuple[0][0]} (C^beta - C^alpha)_i+1 - (C^beta - C^alpha)_i link')
    lines.append('[ dihedrals ]')
    lines.append(f'+{backbone_dihedral_tuple[3][1]} {{"resname": "{backbone_dihedral_tuple[3][0]}"}} +{backbone_dihedral_tuple[2][1]} {{"resname": "{backbone_dihedral_tuple[2][0]}"}} {backbone_angle_tuple[1][1]} {{"resname": "{backbone_angle_tuple[1][0]}"}} {backbone_angle_tuple[0][1]} {{"resname": "{backbone_angle_tuple[0][0]}"}} {parameter_dict["func_type"]:d} {parameter_dict["C"][0]:f} {parameter_dict["C"][1]:f} {parameter_dict["C"][2]:f} {parameter_dict["C"][3]:f} {parameter_dict["C"][4]:f} {parameter_dict["C"][5]:f}')

# %% [markdown]
# ### $(C^\alpha)_{i+1} - (C^\beta - C^\alpha)_{i} - (C^\beta)_{i-1}$

# %%
backbone_dihedral_tuples = []

# beta-alpha-beta-alpha
for tail_monomer, tail_backbone_atom in Cbeta.items():
    for (center_left_monomer, center_left_backbone_atom), (center_right_monomer, center_right_backbone_atom) in zip(Calpha.items(), Cbeta.items()):
        for head_monomer, head_backbone_atom in Calpha.items() :
            backbone_dihedral_tuples.append(((tail_monomer, tail_backbone_atom), (center_left_monomer, center_left_backbone_atom), (center_right_monomer, center_right_backbone_atom), (head_monomer, head_backbone_atom)))

# %%
backbone_dihedral_tuples

# %%
len(backbone_dihedral_tuples)

# %%
backbone_dihedral_parameter_dict = {
    backbone_dihedral_tuple: {
        "func_type": default_dihedral_func_type,
        "C": default_dihedral_C,
    } for backbone_dihedral_tuple in backbone_dihedral_tuples}

# %%
for backbone_dihedral_tuple, parameter_dict in backbone_dihedral_parameter_dict.items():
    lines.append('')
    lines.append('[ link ]')
    lines.append(f'; {backbone_dihedral_tuple[3][0]}-{backbone_dihedral_tuple[2][0]}-{backbone_dihedral_tuple[0][0]} (C^alpha)_i+1 - (C^beta - C^alpha)_i (C^beta)_i-1 link')
    lines.append('[ dihedrals ]')
    lines.append(f'+{backbone_dihedral_tuple[3][1]} {{"resname": "{backbone_dihedral_tuple[3][0]}"}} {backbone_dihedral_tuple[2][1]} {{"resname": "{backbone_dihedral_tuple[2][0]}"}} {backbone_dihedral_tuple[1][1]} {{"resname": "{backbone_dihedral_tuple[1][0]}"}} -{backbone_dihedral_tuple[0][1]} {{"resname": "{backbone_dihedral_tuple[0][0]}"}} {parameter_dict["func_type"]:d} {parameter_dict["C"][0]:f} {parameter_dict["C"][1]:f} {parameter_dict["C"][2]:f} {parameter_dict["C"][3]:f} {parameter_dict["C"][4]:f} {parameter_dict["C"][5]:f}')

# %% [markdown]
# ## write links file

# %%
with open(output_ff, "w") as f:
    for line in lines:
        f.write(line + "\n")
