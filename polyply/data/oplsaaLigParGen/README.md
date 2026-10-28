Provenance of PAMA monomer topologies and parameters.

We use renumberd and charge-normalized itp,

    cp ~/ws/2025-12-01-pama-aa/monomers/EH/normalized.itp PAMA_EH.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/LA/normalized.itp PAMA_LA.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/MMA/normalized.itp PAMA_MMA.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/OC/normalized.itp PAMA_OC.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/ST/normalized.itp PAMA_ST.oplsaa.LigParGen.ff

    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_EH.oplsaa.LigParGen.ff
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_LA.oplsaa.LigParGen.ff 
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_MMA.oplsaa.LigParGen.ff 
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_OC.oplsaa.LigParGen.ff 
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_ST.oplsaa.LigParGen.ff 

## 2026-10-28 backbone hybridization fix + OPLS-AA/L-OPLS re-parametrization

All six monomer `.ff` files (and `etc/custom-ffnonbonded.itp`) were updated
by `etc/fix_backbone_hybridization.py`, `etc/assign_library_charges.py`,
and `etc/retype_alkyl_tail_lopls.py`, followed by `etc/renormalize-charges.py`.
Full rationale, the OPLS-AA atom-type identification trail, and the
deviation found during implementation (ester-junction torsions use OPLS-AA's
`C`/`OS` bond types, not the `C_2`/`O_2` variant L-OPLS extends) are in
`etc/report.md` §6.

Reference data for both the OPLS-AA library charges and the L-OPLS deltas
was vendored from `etc/lopls-reference/` (`wesbarnett/lopls`,
commit `59ffe9819f53598eb96576c75a6417fef52c1755`). Cite, in addition to
the LigParGen/OPLS-AA citations already covering these files:

- S. W. I. Siu, K. Pluhackova, and R. A. Böckmann, *J. Chem. Theory
  Comput.* **8**, 1459–1470 (2012). https://doi.org/10.1021/ct200908r
- K. Pluhackova et al., *J. Phys. Chem. B* (2015).
  https://doi.org/10.1021/acs.jpcb.5b08569

