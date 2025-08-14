## Chlorophyll b (CLB) map

Like Chlorophyll a, Chlorophyll b contains a porphyrin ring, but it differs from Chlorophyll a in having a formyl group (-CHO) at the C7
position instead of a methyl group (-CH₃). This molecule plays a vital role in photosynthesis, acting as a light-harvesting pigment.

sections:
- intended use
- how to use
- available parameters
- warnings that can be raised by this map

Intended use
------------
This map is intended for treating Chlorophyll B (CLB) molecules in proteins. The parameters in this map are taken from literature, specifically from
DFT calculations of the S0 and S1 (QY) transitions of Chl b. Parameters are available both for the electrostatic solvent shift, transition
dipole coupling, and TrEsp coupling.  


How to use
----------
This map has no specific requirements beyond the availability of the corresponding dynamic structure. This map was tested in a LHCII
monomer (not trimer) model system containing both Chl a and Chl b. The model is a united atom model, so only acidic protons are present,
no carbon bound hydrogens.


Available parameters
--------------------
At this point the frequency map has no options or parameters. It is possible to change between B3LYP and HF-CIS difference charges,
by manually selecting the desired choice in the core.txt file.


Warnings that can be raised by this map
---------------------------------------
At this point none.