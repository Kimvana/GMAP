Chlorophyll a (Chl a) map
-------------------------

Chlorophyll a is a central pigment in photosynthesis and features a porphyrin ring with a methyl group (-CH₃) at the C7 position, distinguishing it from Chlorophyll b. It absorbs light primarily in the blue-violet and red regions, contributing significantly to the light-harvesting process.

sections:
- intended use
- how to use
- available parameters
- warnings that can be raised by this map


intended use
------------

This map is intended for treating Chlorophyll a (Chl a) molecules in protein environments. The parameters included are derived from literature, specifically from DFT calculations of the S₀ and S₁ (Qᵧ) transitions of Chl a. Parameters are provided for the electrostatic solvent shift, transition dipole coupling, and TrEsp coupling.


how to use
----------

This map has no special prerequisites aside from the presence of the relevant dynamic structure. It has been validated in a LHCII monomer system that includes both Chl a and Chl b pigments. The model follows a united atom representation, meaning only acidic protons are explicitly included, while hydrogens bound to carbon atoms are omitted.


Available parameters
--------------------

Currently, the frequency map does not provide user-adjustable options. However, you can toggle between B3LYP and HF-CIS difference charges by manually editing the core.txt file to select the preferred charge set.


warnings that can be raised by this map
---------------------------------------

None at this time.
