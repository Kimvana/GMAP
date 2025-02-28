
################################################################
Calculating the amide-I spectrum of a Protein
################################################################

This tutorial will demonstrate how to use GMAP's GEM to calculate the spectrum of the amide-I vibration of a protein. Or, more specifically, how to convert an MD trajectory into a :ref:`Hamiltionian <Theory_page_Hamiltonians>` trajectory, one of the important steps in the spectrum-calculating process.

Before this tutorial, make sure you have/did the following:

- You know how to do simple things on the command line. Things like navigating between directories, creating a file/directory, or listing all files present in the current directory.
- You have MD trajectory files. Not .pdb files directly from the protein data bank. Trajectories from GROMACS, CHARMM and Amber have been confirmed to work. This tutorial does not cover creating these files or how to perform an MD calculation.
- You have done the GMAP installation, and have the enviroment activated. This means that typing 'GMAP' on the command line shows the GMAP help.
- You have NISE installed. NISE will convert our hamiltonian trajectory into a spectrum. You can use the output files from GEM to follow the NISE tutorial.


*******************************************************************************
Defining the filetree
*******************************************************************************

This tutorial will be pointing to file locations quite a bit, so lets establish first what our filetree looks like:

.. code-block:: text
    
    D:\spectra\my_protein_spectra
    ├── MD_files
    │   ├── my_protein.tpr
    │   ├── my_protein.xtc
    │   └── (optionally other files, too)
    ├── GEM_files
    │   └── (this directory is empty for now)
    └── (optionally other files, too)

In this case, we are working on a system called 'my_protein'. A special directory for calculating spectra has been established, containing a directory for each software package used. This structure is not the only way of doing this, but we'll use it here. Of course, you should pick a name more indicative of your actual system, make sure to keep replacing 'my_protein' with the same term of your choice throughout this tutorial.

This main directory (my_protein_spectra) will contain the scripts, parameter files and perhaps slurm files (when running on a slurm-style cluster). There might be an additional folder with protein data bank files, or maybe you placed these in the MD_files directory. The MD_files directory might contain other files generated while computing the trajectory.

This example uses GROMACS-style MD files, and windows-style file paths. This matters little, however, just use your own MD package with the corresponding output files, and replace the windows-style paths with those of your OS.



*******************************************************************************
Creating our first input parameter file.
*******************************************************************************

The input parameter file should be created in the my_protein_spectra directory, and we'll call it GEM_input.txt:

.. code-block:: text
    
    D:\spectra\my_protein_spectra
    ├── MD_files
    │   ├── my_protein.tpr
    │   └── my_protein.xtc
    ├── GEM_files
    │   └── (this directory is empty for now)
    └── GEM_input.txt

The contents of the file should (for now) look like this:

.. code-block:: text
    
    topology_file        D:\spectra\my_protein_spectra\MD_files\my_protein.tpr
    trajectory_file      D:\spectra\my_protein_spectra\MD_files\my_protein.xtc

    output_directory     D:\spectra\my_protein_spectra\GEM_files
    log_directory        D:\spectra\my_protein_spectra\GEM_files

    maps_to_use          AmideBB AmideSC
    couplings_to_use     ProteinAmide   :All

    number_frames        1

Note that a file like this is great to get an initial feel for the program, but doesn't generate files that are usable for NISE. More information on how to specify parameters in a file like this :ref:`can be found here. <UserGuide_page_specifying_parameters>`

So, to explain what's happening:

The first two parameters/lines are to tell GEM what input files it should read from. In the example, the files have been specified using absolute paths (this is also how GMAP returns paths when printing/logging), but they can also be specified relative to the input parameter files. More information on these parameters can be found at :ref:`topology_file <UserGuide_page_parameter_overview_topfile>` and :ref:`trajectory_file <UserGuide_page_parameter_overview_trjfile>`.

The second pair of parameters/lines instruct GEM where it should place its output. Using a different directory for this helps keep our main directory clean and organized. More information on these parameters can be foud at :ref:`output_directory <UserGuide_page_parameter_overview_outdir>` and :ref:`log_directory <UserGuide_page_parameter_overview_logdir>`.

Then, we get to the meat of the instructions. We've defined the system we'd like GEM to look at, now we need to tell GEM what to look _for_ when analyzing the system. GEM will look for so-called singles, which are dyes, choromophores, oscillators, anything that gives a spectroscopic signal. For each single, GEM will calculate an absorption/resonance frequency (the diagonal of the :ref:`Hamiltionian <Theory_page_Hamiltonians>`). For each pair of singles, a coupling between the two can be calculated.

The first step is to tell GMAP which kind of single we're interested in. The options available depend on the maps present. GMAP ships with a few options, and more can be downloaded. The maps required for calculating the Amide-I mode of proteins are included with the GMAP download. So, we select them using the parameter :ref:`maps_to_use <UserGuide_page_parameter_overview_usemaps>`. We select two, as AmideBB is meant for amide groups in the protein backbone, while AmideSC works for amide groups in the protein side chains (asparagine and glutamine residues).

Next, we would like to calculate couplings, so we request a coupling map using the parameter :ref:`couplings_to_use <UserGuide_page_parameter_overview_usecoup>`. In this case, we opt for the map designed to couple protein amide groups, dubbed 'ProteinAmide', and we tell it that all singles should be computed with it. 

As the very last thing, we tell GEM to only treat a single frame. The reason for doing this is that a single frame will be computed fast, so we can quickly see if the program works the way we think it does. We can also look at the file size of the output files - the size scales linearly with the amount of frames. This can give an estimate of the amount of storage space needed to fit all frames.




*******************************************************************************
Running our first calculation
*******************************************************************************

To use the input file, open a command prompt (windows) or terminal (unix). Navigate to the directory my_protein_spectra, so that typing dir (windows) or ls (unix) lists the GEM_input.txt file we created during the following step. Let's run the following command:

``GMAP GEM run GEM_input.txt``

This should start the program. If not, check out the troubleshooting section at the bottom of this page. 






*******************************************************************************
Troubleshooting
*******************************************************************************

There's always something that can go wrong. If you get stuck on any of the steps above, you can find some tips on getting unstuck here.


Running our first calculation
===============================================================================
If you get stuck here, try the following steps. Start at the top, work down the list until you find something that solves your issue.

- very carefully check all input file paths used in your file. Are they correct? You could try using absolute paths, and triple-checking the path to use. Navigate to the directory containing the files - such that typing ``dir`` (windows) or ``ls`` (unix) lists the file. Then, by typing ``cd`` (windows) or ``pwd`` (unix), you are shown the full absolute path. This path, plus the filename, should match the one in your input file.
- While the output files don't exist, the output directory must exist. GMAP will not create a directory to put your output files in, so if it's missing, make sure to create it.
- Try running the ``GMAP`` command again in your command line. Does it give the GMAP logo, and some information? There should be no errors raised here. If there are, make sure to resolve those first.
- Try using the default MD files shipped with GMAP. If the following input file works, but the one before doesn't, there might be an issue with your MD files:

  .. code-block:: text
    
    # topology_file        D:\spectra\my_protein_spectra\MD_files\my_protein.tpr
    # trajectory_file      D:\spectra\my_protein_spectra\MD_files\my_protein.xtc

    # output_directory     D:\spectra\my_protein_spectra\GEM_files
    # log_directory        D:\spectra\my_protein_spectra\GEM_files

    maps_to_use          AmideBB AmideSC
    couplings_to_use     ProteinAmide   :All

    number_frames        1

  These '#' means that GMAP should ignore the lines. Then, it'll use the default files instead.





