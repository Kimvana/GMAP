.. _UserGuide_page_Hamiltonians:

=============
Hamiltonians
=============

Use Case
--------
The main reason to use GMAP is to obtain the strange creatures known as Hamiltionians.
They are an important input for programs like NISE that compute absorption spectra.
The Hamiltonian contains the frequency at which a an oscillator absorbs light and
how strongly each of these oscillators talk to each other and thus influence each
others absorption frequency.

Physical description
--------------------

The Hamiltonian is given by[1]:

.. math::
    H(t) = \sum_n \epsilon_n(t)B_n^\dagger B_n - \frac{1}{2}\sum_n \Delta_n(t)B_n^\dagger B_n^\dagger B_nB_n+
    \sum_{nm} J_{nm}(t)B_n^\dagger B_m - \sum_n \vec{E}(t)\cdot \vec{\mu}_n(t)[B_n^\dagger + B_n] ,

where :math:`B_n^\dagger` and :math:`B_n` are the creation and annihilation operators respectively.
:math:`\epsilon_n(t)` represents the energy gap between the ground state and the first excited states on
amide group :math:`n` . The coupling strength of amide group :math:`n` with the external electric field :math:`\vec{E}(t)`
is controlled by the transition dipole :math:`\mu_n`. The anharmonicity :math:`\Delta_n(t)` defines the
difference in energy energy gap between the ground and first excited state and the first excited and second
excited state, where a positive anharmonicity indicates that the latter is smaller than the former. 
:math:`J_{nm}` represents the transition dipole coupling. The reason why :math:`\vec{E}(t)`, :math:`J_{nm}(t)`
and :math:`\vec{\mu}_n(t)` are dependent on time is because the investigated structure is not necessarily
completely static and constantly shifts ever so slightly, and is thus examined at discrete frames.

This formula is impressive, however, it is not used directly. The program works by using maps which use physical 
characteristics of the system to compute the values inside the Hamiltionians.

File layout
-----------

A Hamiltonian is computed for every single frame. In the txt format, the entire Hamiltonian lives on
a single line. Therefore, the file contains a line for each frame. On a single line, the first number is the
number of the frame, followed by a space. Then, the Hamiltonian itself follows. As the Hamiltonian is a
diagonal matrix, only roughly half of the entries need to be written. The Hamiltonian is written row by
row, starting at the diagonal entry. Therefore, after the frame number, the entry at 0,0 is written to the
output file, followed by the one at 0,1, etc. For the next line, the first item is 1,1 (1,0 is skipped as it is
identical to the value at 0,1), the second is 1,2, etc. The number of entries total in the Hamiltonian stored
this way is the triangular number of the width of the Hamiltonian. The values at the diagonal of the Hamiltonian
are the expected vibrational frequencies in wavenumbers. These include electrostatic effects from the 
environment, but exclude any influences from coupling with neighbouring amide groups. The values are usually 
around 1600. The off-diagonal elements of the Hamiltonian list the coupling between amide groups. The coupling 
can be both positive and negative, and occasionally reach values of 30. However, virtually all values are smaller
than 5, with most of them being smaller than 0.1.

.. math::

   \begin{pmatrix}
   C_1 & 1-2 & 1-3 & 1-4 & ... & 1-n \\
   2-1 & C_2 & 2-3 & 2-4 & ... & 2-n \\
   3-1 & 3-2 & C_3 & 3-4 & ... & 3-n \\
   4-1 & 4-2 & 4-3 & C_4 & ... & 4-n \\
   ... & ... & ... & ... & ... & ... \\
   n-1 & n-2 & n-3 & n-4 & ... & C_n
   \end{pmatrix}

The binary file is very similar. It contains no spaces, only 32-bit floats. The first one indicates the frame number, and is followed by the hamiltonian entries
in the same order as the text file. After one hamiltonian has been written, no spacer is used: the next
float is the number of the next frame.

Below is an example Hamiltonian:

[1]: Thomas L.C. Jansen. “Computational Spectroscopy of Complex Systems”. In: The Journal of Chemical
Physics 155 (170901 2021)