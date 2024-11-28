.. _UserGuide_page_references_file:

##############
references.bib
##############

(This applies both to singles and pairs maps)

This is an optional file containing all references that might need to be cited when a paper uses a GMAP calculation. Each reference is present in the default .bib format, with the exeption of two extra (mandatory!) fields being present.

When a calculation is performed, GMAP checks if a map is involved in the calculation at all. Only if the map was both requested (either by the user or by another map) and used (i.e. its molecules are present), its references are considered.

Figuring out what references to print for a map is a two-step process:
- The first step is to ask the map if it is okay with the current .bib file contents. This step is important for maps which have multiple models with separate references (like the amide maps that come with GMAP), as only the model used should be referenced, and the others shouldn't. A map influences this step through the ``GM_report_references`` function in the main.py file of the map. Map makers will probably want to use the 'mapkey' field for this purpose.
- The second step is for GMAP, to see if a reference was actually used for that specific calculation.

