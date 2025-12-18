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

