.. _Theory_page_Hamiltonians:

############
Hamiltonians
############

********
Use Case
********
The main reason to use GMAP is to obtain the Hamiltonians and transition
dipole moments describing the energy landscape of chromophores and their interaction with light.
They are an important input for programs like NISE that compute absorption spectra.
The Hamiltonian contains the frequency at which each chromophore absorbs light, and
the couplings between these chromophores describe how strongly each of 
these chromophores talk to each other and thus influence each others absorption frequency.

********************
Physical description
********************

The Hamiltonian is given by the following equation [1]:

.. math::
    H(t) = \sum_n \epsilon_n(t)B_n^\dagger B_n - \frac{1}{2}\sum_n \Delta_n(t)B_n^\dagger B_n^\dagger B_nB_n+
    \sum_{n,m} J_{nm}(t)B_n^\dagger B_m - \sum_n \vec{E}(t)\cdot \vec{\mu}_n(t)[B_n^\dagger + B_n] ,

where :math:`B_n^\dagger` and :math:`B_n` are the creation and annihilation operators respectively.
:math:`\epsilon_n(t)` represents the energy gap between the ground state and the first excited state
in a chromophore :math:`n`. The coupling strength of amide group :math:`n` with the external electric field :math:`\vec{E}(t)`
is affected by the transition dipole :math:`\mu_n`. The anharmonicity :math:`\Delta_n(t)` is the
difference in energy energy gaps between the ground and first excited states and the first excited and second
excited states, when considering the chromophores as three level systems. A positive anharmonicity indicates that
the latter is smaller than the former. :math:`J_{nm}(t)` represents the coupling between chromophores :math:`n` and :math:`m`.
The time-dependent electric field, :math:`\vec{E}(t)`, is not a system property and not considered by GEM.
The coupling and transition dipoles dependent on time as the physical structure that we derive the Hamiltonian
from is not necessarily static and can change over time. 

GEM obtains some of the terms in this equation, which are then used by programs
like NISE to compute spectral simulations. The program works by using maps which use physical 
characteristics of the system to compute the values that make up the Hamiltionian trajectory. 
The way these maps work and how to add them to the program is explained :ref:`here<UserGuide_page_adding_map>`

An example Hamiltonian with :math:`n` chromophores is given below. :math:`\epsilon_m`
represents the frequency of the :math:`m` th chromophore. :math:`J_{kl}`` is 
the coupling between the :math:`k` th and the  :math:`l` th chromophores.
Note that :math:`J_{kl}` and :math:`J_{lk}` are identical. In the following the number of chromophores is :math:`(N+1)`.

.. math::

   \begin{pmatrix}
   \epsilon_0 & J_{01} & J_{02} & J_{03} & ... & J_{0N} \\
   J_{10} & \epsilon_1 & J_{11} & J_{13} & ... & J_{1N} \\
   J_{20} & J_{21} & \epsilon_2 & J_{23} & ... & J_{2N} \\
   J_{30} & J_{31} & J_{32} & \epsilon_3 & ... & J_{3N} \\
   ... & ... & ... & ... & ... & ... \\
   J_{N0} & J_{N1} & J_{N2} & J_{N3} & ... & \epsilon_N
   \end{pmatrix}

***********
File layout
***********

Each frame has its own Hamiltonian. In the txt format, each frame has its
own line. Thus, after each Hamiltonian a new line is started. The first number on
each line is the number of the frame. Then, the values that make up the Hamiltonian follow.
As the Hamiltonian is a matrix which is symmetrical across the diagonal, only half of
the couplings need to be written saving almost a factor two in storage space.
The Hamiltonian is written starting at the diagonal entry followed by the couplings
beween the first unit and all the other units. The amount of entries total in the Hamiltonian stored
this way is the triangular number of the size of the Hamiltonian (:math:`(N+2)(N+1)/2`). The values at the diagonal of the Hamiltonian
are the expected frequencies in wavenumbers [cm :math:`^{-1}` ]. These frequencies include electrostatic effects from the 
environment, but exclude any influences from coupling with neighbouring amide groups. 
The off-diagonal elements of the Hamiltonian denote the coupling between amide groups,
which can be both positive and negative.

Thus the Hamiltonian trajectory in text format looks like:

    | 0  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0N}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1N}` ... :math:`\epsilon_N`
    | 1  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0N}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1N}` ... :math:`\epsilon_N`
    | 2  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0N}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1N}` ... :math:`\epsilon_N`
    | ...  ...
    | m  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0N}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1N}` ... :math:`\epsilon_N`

Where :math:`m` is the number of Hamiltionians. :math:`J_{kl}` and :math:`\epsilon_n` keep the same meaning 
as in the Hamiltonian above, except now they are per frame in the trajectory.

The binary file contains only 32-bit floats. The first one indicates the 
frame number, which is still a 32-bit float even though it represents an integer,
and is followed by the Hamiltonian entries that belong to that frame number.
After one Hamiltonian has been written, no spacer is used: the next
float is the number of the next frame. Thus, in binary format the Hamiltonian trajectory looks like:

    | 0  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0N}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1n}` ... :math:`\epsilon_n` 1  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0n}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1n}` ... :math:`\epsilon_n` 2  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0n}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1n}` ... :math:`\epsilon_n` ...  ... m  :math:`\epsilon_0`  :math:`J_{01}` :math:`J_{02}` :math:`J_{03}` ...  :math:`J_{0N}` :math:`\epsilon_1`  :math:`J_{12}` :math:`J_{13}` ...  :math:`J_{1N}` ... :math:`\epsilon_N`

Where :math:`m`, :math:`J_{kl}` and :math:`\epsilon_N` keep the same meaning as in the
representation of the text trajectory above.

The transition dipole moments are stored in a separate file. The values are
expected the units of Debye. As for the Hamiltonian in the text format The
for the transition dipoles every frame is stored at a separate line. Note that
first the x-components of the vectors are stored together, then the y-components
and finally the z-components.

The transition dipole in the text format looks like:

    | 0  :math:`\mu_{0,x}`  :math:`\mu_{1,x}`  :math:`\mu_{2,x}`  ...  :math:`\mu_{N,x}`  :math:`\mu_{0,x}`  :math:`\mu_{1,y}`  :math:`\mu_{2,y}`  ...  :math:`\mu_{N,y}`  :math:`\mu_{0,z}`  :math:`\mu_{1,z}`  :math:`\mu_{2,z}`  ...  :math:`\mu_{N,z}`
    | 1  :math:`\mu_{0,x}`  :math:`\mu_{1,x}`  :math:`\mu_{2,x}`  ...  :math:`\mu_{N,x}`  :math:`\mu_{0,x}`  :math:`\mu_{1,y}`  :math:`\mu_{2,y}`  ...  :math:`\mu_{N,y}`  :math:`\mu_{0,z}`  :math:`\mu_{1,z}`  :math:`\mu_{2,z}`  ...  :math:`\mu_{N,z}`
    | 2  :math:`\mu_{0,x}`  :math:`\mu_{1,x}`  :math:`\mu_{2,x}`  ...  :math:`\mu_{N,x}`  :math:`\mu_{0,x}`  :math:`\mu_{1,y}`  :math:`\mu_{2,y}`  ...  :math:`\mu_{N,y}`  :math:`\mu_{0,z}`  :math:`\mu_{1,z}`  :math:`\mu_{2,z}`  ...  :math:`\mu_{N,z}`
    | ... ...
    | m  :math:`\mu_{0,x}`  :math:`\mu_{1,x}`  :math:`\mu_{2,x}`  ...  :math:`\mu_{N,x}`  :math:`\mu_{0,x}`  :math:`\mu_{1,y}`  :math:`\mu_{2,y}`  ...  :math:`\mu_{N,y}`  :math:`\mu_{0,z}`  :math:`\mu_{1,z}`  :math:`\mu_{2,z}`  ...  :math:`\mu_{N,z}`

As for the Hamiltonian the binary file contains only 32-bit floats. The first one indicates the 
frame number, which is still a 32-bit float even though it represents an integer,
and is followed by the transition dipole entries that belong to that frame number.
After one frame has been written, no spacer is used: the next
float is the number of the next frame.

[1]: T.L.C. Jansen, J. Chem. Phys. 155:170901 (2021) https://doi.org/10.1063/5.0064092