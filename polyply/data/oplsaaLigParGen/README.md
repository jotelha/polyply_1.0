Provenance of PAMA monomer topologies and parameters.

    cp ~/ws/2025-12-01-pama-aa/monomers/EH/EH_34067F/tmp/EH_34067F.itp PAMA_EH.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/LA/LA_7F9BE9/tmp/LA_7F9BE9.itp PAMA_LA.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/MMA/MMA_7EA60E/tmp/MMA_7EA60E.itp PAMA_MMA.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/OC/OC_76D421/tmp/OC_76D421.itp PAMA_OC.oplsaa.LigParGen.ff
    cp ~/ws/2025-12-01-pama-aa/monomers/ST/ST_C8C92E/tmp/ST_C8C92E.itp PAMA_ST.oplsaa.LigParGen.ff

    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_EH.oplsaa.LigParGen.ff
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_LA.oplsaa.LigParGen.ff 
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_MMA.oplsaa.LigParGen.ff 
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_OC.oplsaa.LigParGen.ff 
    sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' PAMA_ST.oplsaa.LigParGen.ff 
