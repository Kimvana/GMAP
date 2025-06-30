.. _DevGuide_page_home:

################
Developers guide
################

These pages will contain information specifically for developers of the program.

.. note::
    Map != map()! Two main programs in the GMAP package (GEM and AIM) serve a specific purpose: converting an md trajectory into a Hamiltonian trajectory. In the spectroscopic community, this conversion is traditionally done using so-called maps. In this package, :class:`~GMAP.src.tools.MapReader.Map` objects contain all information on such a spectroscopic map. Do not confuse these objects with python's built-in map() function - it doesn't occur yet in the program as of writing this text (january 2024), and is not expected to, either.


.. grid:: 1 2 2 3

    .. grid-item-card::
        :margin: 0 3 0 0
        :link: useful_resources
        :link-type: doc

        **useful resources**
        ^^^^^^^^^^^^^^^^
        All kinds of links and information 

    .. grid-item-card::
        :margin: 0 3 0 0
        :link: Program_flow/index
        :link-type: doc

        **Program flow**
        ^^^^^^^^^^^^^^
        A rough outline of how the program is structured.
    
    .. grid-item-card::
        :margin: 0 3 0 0
        :link: code_style
        :link-type: doc

        **Code Style**
        ^^^^^^^^^^^^^^^^^^^^^^^^^
        Guidelines on how to style your code to keep the entire codebase consistent.
    
    .. grid-item-card::
        :margin: 0 3 0 0
        :link: print_colors
        :link-type: doc

        **print_colors**
        ^^^^^^^^^^^^^^^^^^^^^^^^^
        A chart of the available colors and how they'll render.
    
    .. grid-item-card::
        :margin: 0 3 0 0
        :link: VScode_setup
        :link-type: doc

        **VScode setup**
        ^^^^^^^^^^^^^^^^^^^^^^^^^
        How to recreate my exact coding environment.


.. toctree::
    :hidden:
    
    useful_resources
    Program_flow/index
    code_Style
    print_colors
    VScode_setup