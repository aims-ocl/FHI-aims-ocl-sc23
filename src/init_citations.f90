!****s* FHI-aims/init_citations
!  NAME
!   init_citations
!  SYNOPSIS
      subroutine init_citations()
!  PURPOSE
!     This functions helps you to define a citation for your feature at the end
!     of the output of an FHI-aims run. Just define the text for your citation
!     using the citetext array, choose a tag for your feature and call
!
!        register_citation("your tag", citetext, "priority")
!
!     where priority can be either "core", "feature" or "technical". "core" is
!     reserved for special references, you should use either "feature" or
!     "technical", where "technical" aims at stuff like integration routines or
!     grids.
!     After registering your citation, you can place activators for it in your
!     code. When the run enters your codeblock, call
!
!        cite_reference("your tag")
!
!     to signal the citation manager that your reference should be printed at
!     the of the run. That's already everything you have to do to add your cite!
!  USES
      use applicable_citations, only: init_citationlist, register_citation
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2014).
!  SOURCE
      implicit none
      character(len=130) :: citetext(15)

      call init_citationlist()
      citetext(:) = "____"

      citetext(1) = "  For any use of FHI-aims, please cite:"
      citetext(2) = ""
      citetext(3) = "    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,"
      citetext(4) = "    Xinguo Ren, Karsten Reuter, and Matthias Scheffler"
      citetext(5) = "    'Ab initio molecular simulations with numeric atom-centered orbitals'"
      citetext(6) = "    Computer Physics Communications 180, 2175-2196 (2009)"
      citetext(7) = "    http://dx.doi.org/10.1016/j.cpc.2009.06.022"
      call register_citation("FHI-aims-CPC", citetext, "core")

      citetext(1) = "  For Hartree-Fock, hybrid functionals, or many-body perturbation theory used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    Xinguo Ren, Patrick Rinke, Volker Blum, Juergen Wieferink, Alex Tkatchenko,"
      citetext(4) = "    Andrea Sanfilippo, Karsten Reuter, and Matthias Scheffler,"
      citetext(5) = "    'Resolution-of-identity approach to Hartree-Fock, hybrid density functionals,"
      citetext(6) = "    RPA, MP2, and GW with numeric atom-centered orbital basis functions'"
      citetext(7) = "    New Journal of Physics 14, 053020 (2012)."
      citetext(8) = "    http://dx.doi.org/10.1088/1367-2630/14/5/053020"
      call register_citation("Ren_NJP_2012", citetext, "feature")

      citetext(1) = "  For the NAO-VCC-nZ basis sets used in this calculation, please cite:"
      citetext(2) = ""
      citetext(3) = "    Igor Ying Zhang, Xinguo Ren, Patrick Rinke, Volker Blum, and Matthias Scheffler"
      citetext(4) = "    'Numeric atom-centered-orbital basis sets with valence-correlation consistency from H to Ar'"
      citetext(5) = "    New Journal of Physics 15, 123033 (2013)."
      citetext(6) = "    http://dx.doi.org/10.1088/1367-2630/15/12/123033"
      call register_citation("NAO-VCC-2013", citetext, "feature")

      citetext(1) = "  For pseudopotentials or QM/MM embedding infrastructure used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    D. Berger, A. J. Logsdail, H. Oberhofer, M. R. Farrow, C. R. A. Catlow, P. Sherwood,"
      citetext(4) = "    A. A. Sokol, V. Blum, and K. Reuter,"
      citetext(5) = "    'Embedded-cluster calculations in a numeric atomic orbital'"
      citetext(6) = "    'density-functional theory framework'"
      citetext(7) = "    The Journal of Chemical Physics 141, (2014)."
      citetext(8) = "    http://dx.doi.org/10.1063/1.4885816"
      call register_citation("Berger_JCP_2014", citetext, "feature")

      citetext(1) = "  For the real-space grid partitioning and parallelization used in this calculation, please cite:"
      citetext(2) = ""
      citetext(3) = "    Ville Havu, Volker Blum, Paula Havu, and Matthias Scheffler,"
      citetext(4) = "    'Efficient O(N) integration for all-electron electronic structure calculation'"
      citetext(5) = "    'using numerically tabulated basis functions'"
      citetext(6) = "    Journal of Computational Physics 228, 8367-8379 (2009)."
      citetext(7) = "    http://dx.doi.org/10.1016/j.jcp.2009.08.008"
      call register_citation("FHI-aims-grids", citetext, "technical")

      citetext(1) = "  The scalable eigensolver library ELPA was used in your run."
      citetext(2) = "  ELPA is essential especially for large systems on hundreds or thousands of CPUs."
      citetext(3) = "  For ELPA, please cite:"
      citetext(4) = ""
      citetext(5) = "    Andreas Marek, Volker Blum, Rainer Johanni, Ville Havu, Bruno Lang,"
      citetext(6) = "    Thomas Auckenthaler, Alexander Heinecke, Hans-Joachim Bungartz, and Hermann Lederer,"
      citetext(7) = "    'The ELPA Library - Scalable Parallel Eigenvalue Solutions"
      citetext(8) = "    for Electronic Structure Theory and Computational Science'"
      citetext(9) = "    The Journal of Physics: Condensed Matter 26, 213201 (2014)."
      citetext(10) = "    http://dx.doi.org/10.1088/0953-8984/26/21/213201"
      call register_citation("ELPA", citetext, "technical")

      citetext(1) = "  For ELPA's efficient two-stage tridiagonalization used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    T. Auckenthaler, V. Blum, H.-J. Bungartz, T. Huckle, R. Johanni, L. Kraemer,"
      citetext(4) = "    B. Lang, H. Lederer, and P. R. Willems,"
      citetext(5) = "    'Parallel solution of partial symmetric eigenvalue problems'"
      citetext(6) = "    'from  electronic structure calculations'"
      citetext(7) = "    Parallel Computing 37, 783-794 (2011)."
      citetext(8) = "    http://dx.doi.org/10.1016/j.parco.2011.05.002"
      call register_citation("ELPA-2stage", citetext, "technical")

      citetext(1) = "  For second-variational spin-orbit coupling used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    William P. Huhn and Volker Blum," 
      citetext(4) = "    'One-hundred-three compound band-structure benchmark of post-self-consistent'"
      citetext(5) = "    'spin-orbit coupling treatments in density functional theory'"
      citetext(6) = "    Phys. Rev. Materials 1, 033803 (2017)."
      citetext(8) = "    https://dx.doi.org/10.1103/PhysRevMaterials.1.033803"
      call register_citation("Spin_Orbit_Coupling", citetext, "feature")

      citetext(1) = "  For the SMPB implicit solvation architecture used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    S. Ringe, H. Oberhofer, C. Hille, S. Matera, K. Reuter,"
      citetext(4) = "    'Function-Space-Based Solution Scheme for the Size-Modified '"
      citetext(5) = "    'Poisson-Boltzmann Equation in Full-Potential DFT'"
      citetext(6) = "    Journal of Chemical Theory and Computation 12, 4052-4066 (2016)."
      citetext(7) = "    http://dx.doi.org/10.1021/acs.jctc.6b00435"
      citetext(8) = ""
      citetext(9) = "  For the SMPB ionic parameters used in your run, please cite:"
      citetext(10)= "    S. Ringe, H. Oberhofer, K. Reuter,"
      citetext(11)= "    'Transferable ionic parameters for first-principles Poisson-Boltzmann'"
      citetext(12)= "    'solvation calculations: Neutral solutes in aqueous monovalent'"
      citetext(13)= "    'salt solutions'"
      citetext(14)= "    Journal of Chemical Physics 146, 134103 (2017)."

      
      call register_citation("Ringe_JCTC_2016", citetext, "feature")

      citetext(1) = "  For the MPE implicit solvation feature used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    M. Sinstein, C. Scheurer, S. Matera, V. Blum, K. Reuter, and H. Oberhofer"
      citetext(4) = "    'An efficient implicit solvation method for full potential DFT'"
      citetext(5) = "    Journal of Chemical Theory and Computation (2017)."
      citetext(6) = "    http://dx.doi.org/10.1021/acs.jctc.7b00297"
      call register_citation("MPE_implicit_solvation", citetext, "feature")

      citetext(1) = "  For the analytical stress tensor used in your run, please cite:"
      citetext(2) = ""
      citetext(3) = "    Franz Knuth, Christian Carbogno, Viktor Atalla, Volker Blum, Matthias Scheffler"
      citetext(4) = "    'All-electron formalism for total energy strain derivatives and"
      citetext(5) = "    stress tensor components for numeric atom-centered orbitals'"
      citetext(6) = "    Computer Physics Communications 190, 33-50 (2015)."
      citetext(7) = "    http://dx.doi.org/10.1016/j.cpc.2015.01.003"
      call register_citation("FHI-aims-analytical-stress-tensor", citetext, "feature")

      citetext(1) = "  The high-level electronic structure method used in your calculation"
      citetext(2) = "  employed our efficient localized resolution of identity for the Coulomb operator."
      citetext(3) = "  Some calculations, especially large-scale hybrid density functional theory,"
      citetext(4) = "  are only possible thanks to this development. For this feature, please cite:"
      citetext(5) = ""
      citetext(6) = "    Arvid Ihrig, Juergen Wieferink, Igor Ying Zhang, Matti Ropo, Xinguo Ren,"
      citetext(7) = "    Patrick Rinke, Matthias Scheffler, and Volker Blum,"
      citetext(8) = "    'Accurate localized resolution of identity approach for linear-scaling"
      citetext(9) = "    hybrid density functionals and for many-body perturbation theory'"
      citetext(10)= "    New Journal of Physics 117, 093020 (2015)."
      citetext(11)= "    http://dx.doi.org/10.1088/1367-2630/17/9/093020"
      call register_citation("FHI-aims-lvl", citetext, "technical")

      citetext(1) = "  The provided symmetry information was generated with SPGlib:"
      citetext(2) = ""
      citetext(3) = "    Atsushi Togo, Yusuke Seto, Dimitar Pashov"
      citetext(4) = "    SPGlib 1.7.3 obtained from http://spglib.sourceforge.net"
      citetext(5) = "    Licence: New BSD license"
      citetext(5) = "    Copyright (C) 2008 Atsushi Togo"
      call register_citation("SPGlib", citetext, "feature")

      citetext(1) = "  For calculations of electronic coupling values using the "
      citetext(2) = "  Fragment molecular orbital scheme (FO-DFT), please cite:"
      citetext(3) = ""
      citetext(4) = "    Christoph Schober, Karsten Reuter, Harald Oberhofer "
      citetext(5) = "    'Critical analysis of fragment-orbital DFT schemes"
      citetext(6) = "    for the calculation of electronic coupling values'"
      citetext(7) = "    Journal of Chemical Physics, 144, 054103 (2016)."
      call register_citation("FODFT", citetext, "feature")

      citetext(1) = "  For the 'Strongly Constrained and Appropriately Normed' (SCAN)"
      citetext(2) = "  meta-GGA functional, please cite:"
      citetext(3) = ""
      citetext(4) = "    Jianwei Sun and Adrienn Ruzsinszky and John P. Perdew,"
      citetext(5) = "    'Strongly Constrained and Appropriately Normed Semilocal Density Functional'"
      citetext(6) = "     Physical Review Letters 115, 036402 (2015)."
      call register_citation("SCAN", citetext, "feature")

      citetext(1) = "  For the use of the dfauto autogenerated XC functionals, please cite:"
      citetext(2) = ""
      citetext(3) = "    R. Strange, F. R. Manby and P. J. Knowles,"
      citetext(4) = "    'Automatic code generation in density functional theory',"
      citetext(5) = "    Comp. Phys. Comm. 136, 310-318 (2001)."
      citetext(6) = "    http://dx.doi.org/10.1016/S0010-4655(01)00148-5"
      call register_citation("dfauto", citetext, "feature")


      citetext(1) = " Please note that your FHI-aims calculation has made use of functionals "
      citetext(2) = " obtained from the Density Functional Repository:"
      citetext(3) = " http://www.cse.scitech.ac.uk/ccg/dft/ "
      citetext(4) = " In addition to the original scientific reference" 
      citetext(5) = "of the functionals, users of this repository are"
      citetext(6) = "required to include the following citation in any "
      citetext(7) = "publications resulting from its use:"
      citetext(8) = "Functionals were obtained from the Density Functional Repository"
      citetext(9) = "as developed and distributed by the Quantum Chemistry Group,"
      citetext(10) = "CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD United Kingdom."
      citetext(11) = "Contact Huub van Dam (hvandam@bnl.gov) or Paul Sherwood for further information."
      call register_citation("DFPT_Density_Functional_Repository", citetext, "feature")


      citetext(1) = "  The atom_sphere code by Stefan Goedecker, Santanu Saha, and others was used "
      citetext(2) = "  to create some of the basis functions used in your FHI-aims run. This code is an"
      citetext(3) = "  evolution of a program first used in"
      citetext(4) = ""
      citetext(5) = "    TRANSFERABILITY OF PSEUDOPOTENTIALS"
      citetext(6) = "    By: GOEDECKER, S; MASCHKE, K"
      citetext(7) = "    PHYSICAL REVIEW A   Volume: 45   Issue: 1   Pages: 88-93   Published: JAN 1 1992"
      citetext(8) = ""
      citetext(9) = "  Please reference this publication as part of the FHI-aims basis functions used to "
      citetext(10)= "  obtain your results."
      call register_citation("atom_sphere", citetext, "technical")

      citetext(1) = "  For calculations of tensorial nonadiabatic friction using "
      citetext(2) = "  time-dependent perturbation theory, please cite:"
      citetext(3) = ""
      citetext(4) = "    Reinhard J. Maurer, Mikhail Askerka, Victor S.Batista, John C. Tully"
      citetext(5) = "    'Ab initio tensorial electronic friction for molecules"
      citetext(6) = "    on metal surfaces: Nonadiabatic vibrational relaxation' "
      citetext(7) = "    Physical Review B, 94, 115432 (2016)"
      call register_citation("friction", citetext, "feature")

      citetext(1) = "  The DFPT (perturbation by an atomic displacement) part of the code was used."
      citetext(2) = "  Here is the corresponding reference:     "
      citetext(3) = "                                                                             "
      citetext(4) = "    Honghui Shang, Christian Carbogno, Patrick Rinke and Matthias Scheffler "
      citetext(5) = "    'Lattice dynamics calculations based on density-functional"
      citetext(6) = "    perturbation theory in real space' "
      citetext(7) = "    Computer Physics Communications 215, 26-46 (2017)"
      citetext(8) = "    https://doi.org/10.1016/j.cpc.2017.02.001"
      call register_citation("dfpt", citetext, "feature")
      
      citetext(1) = "  The DFPT (perturbation by an electric field) part of the code was used."
      citetext(2) = "  Here is the corresponding reference:                                  "
      citetext(3) = "                                                                             "
      citetext(4) = "    Honghui Shang, Nathaniel Raimbault, Mariana Rossi,                       "
      citetext(5) = "    Christian Carbogno, Patrick Rinke and Matthias Scheffler,                "
      citetext(6) = "    'All-Electron, Real-Space Perturbation Theory for Homogeneous            "
      citetext(7) = "    Electric Fields: Theory, Implementation, and Application within dft'     "
      citetext(8) = "    New Journal of Physics, 20(7):073040, 2018                               "
      call register_citation("dfpt_electricfield", citetext, "feature")

      citetext(1) = "  The ELSI infrastructure was used in your run to solve the Kohn-Sham electronic structure."
      citetext(2) = "  Please check out http://elsi-interchange.org to learn more."
      citetext(3) = "  If scalability is important for your project, please acknowledge ELSI by citing:"
      citetext(4) = ""
      citetext(5) = "    V. W-z. Yu, F. Corsetti, A. Garcia, W. P. Huhn, M. Jacquelin, W. Jia,"
      citetext(6) = "    B. Lange, L. Lin, J. Lu, W. Mi, A. Seifitokaldani, A. Vazquez-Mayagoitia,"
      citetext(7) = "    C. Yang, H. Yang, and V. Blum"
      citetext(8) = "    'ELSI: A unified software interface for Kohn-Sham electronic structure solvers'"
      citetext(9) = "    Computer Physics Communications 222, 267-285 (2018)."
      citetext(10) = "    http://dx.doi.org/10.1016/j.cpc.2017.09.007"
      call register_citation("ELSI",citetext,"technical")

      citetext(1) = "  The Conductive Heat Flux was calculated using the formalism introduced in"
      citetext(2) = "    C. Carbogno, R. Ramprasad, and M. Scheffler"
      citetext(3) = "    'Ab Initio Green-Kubo Approach for the Thermal Conductivity of Solids'"
      citetext(4) = "    Phys. Rev. Lett. 118, 175901 (2017)."
      citetext(5) = "     http://dx.doi.org/10.1103/PhysRevLett.118.175901"
      call register_citation("HeatFluxGK",citetext,"feature")

      citetext(1) = " DFT+U infrastructure was used in your run "
      citetext(2) = " please cite: "
      citetext(3) = ""
      citetext(4) = " Matthias Kick, Karsten Reuter and Harald Oberhofer,"
      citetext(5) = " Intricacies Of DFT+U, not only in a Numeric Atom Centered Orbital Framework "
      citetext(6) = ""
      citetext(7) = " Journal of Chemical Theory and Computation "
      citetext(8) = " https:// "
      call register_citation("dftpU", citetext, "feature")


      end subroutine init_citations
