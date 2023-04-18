"""
aimsChain module for usage with FHI-aims
========================================

A package for finding minimum energy paths and transition states.
For usage please refer to the FHI-aims manual.


==== Rough Update Log ====

- Summer 2013 - 
Original development by Yinguy Yao under the direction of Luca M. Ghiringhelli

- Next 2.5 years - 
Minor updates including additional comments and optimizers

- February 2018 - 
Updates from Yair Litman and Julian Gebhardt to make compatible with
newest version of FHI-aims.

- March 2018 - 
Justin C Smith merged the February 2018 updates with the updates not included 
in the "Next 2.5 years" category. Created this __init__ file, cleaned up the 
commenting and style throughout the code, and 
implemented outputs to inform the user on what is going.



==== From the FHI-aims manual ====
This project aims to provide an all-in-one package for various flavours of 
the chain of states methods for finding the minimum energy path(MEP). 
Currently the nudged elastic band method (NEB), the string method, 
and the growing string method are included.

The aimsChain code is distributed under the Lesser General Public License as 
published in Version 3, 29 June 2007, at http://www.gnu.org/licenses/lgpl.html.
Some of the optimizer routines in the code originated within the 
Atomic Simulation Environment (ASE). We wish to give full credit to the 
developers of these routines. The aimsChain code can also be found in a 
separate distribution maintained by Yingyu Yao. The reason we distribute them 
directly within the FHI-aims repository is for the convenience of all FHI-aims 
users, but again: We wish to give full credit to the work and the license of 
the original authors within ASE.
"""

