! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine read_options( na, ns, nspin )

  ! Subroutine to read the options for the SIESTA program
  !
  !     It uses the FDF (Flexible Data Format) package 
  !     of J.M.Soler and A.Garcia
  !
  ! Writen by P.Ordejon, December'96
  ! Modified for introduction of dynamic memory in SIESTA by JDG Sept 99
  ! Wrapping of most fdf and broadcast calls: A. Garcia, June 2005
  !

  use siesta_options
  use precision, only : dp, grid_p
  use parallel,  only : IOnode, Nodes
  use fdf
  use files,     only : slabel
  use files,     only : filesOut_t   ! derived type for output file names
  use sys
  use units,     only : eV, Ang, Kelvin
  use siesta_cml
  use m_target_stress, only: set_target_stress
  use m_spin, only: print_spin_options

  use m_charge_add, only : read_charge_add
  use m_hartree_add, only : read_hartree_add
  
  use m_mixing_scf, only: mixers_scf_init
  use m_mixing_scf, only: mixers_scf_print, mixers_scf_print_block

  use m_cite, only: add_citation
  
  implicit none
  !----------------------------------------------------------- Input Variables
  ! integer na               : Number of atoms
  ! integer ns               : Number of species
  ! integer nspin            : Number of spin-components
  !                            1=non-polarized, 2=polarized, 4=non-collinear,
  !                            8=spin-orbit

  integer, intent(in)  :: na, ns, nspin

  ! This routine sets variables in the 'siesta_options' module

  ! The following are comment lines that should be merged into 'siesta_options'.

  ! real*8 charnet           : Net charge (in units of |e|)
  ! logical outlng           : Long (true) or short (false) output
  ! real*8 g2cut             : PW cutoff energy (Ry)
  ! logical negl             : True = Neglect interactions between
  !                            non-overlaping orbitals (coming from
  !                            KB projectors)
  ! integer nscf             : Maximum number of SCF cycles per time step
  ! real*8 dDtol             : Maximum Density Matrix tolerance in SCF
  ! real*8 Energy_tolerance  : Maximum Total energy tolerance in SCF
  ! real*8 Harris_tolerance  : Maximum Harris energy tolerance in SCF
  ! logical mix              : Perform mix in first SCF step
  ! real*8 wmix              : Amount of output DM for new DM
  ! integer isolve           : Method of solution.  0   = Diagonalization
  !                                                 1   = Order-N
  !                                                 2   = Transiesta
  !                                                 3   = OMM
  !                                                 4   = PEXSI
  !                                                 5   = (Matrix write)
  ! real*8 temp              : Temperature for Fermi smearing (Ry)
  ! logical fixspin          : Fix the spin of the system?
  ! real*8  total_spin       : Total spin of the system
  ! integer ncgmax           : Maximum number of CG steps for 
  !                            band structure energy minimization
  ! real*8 etol              : Relative tolerance in CG minimization
  !                            of band structure energy
  ! real*8 eta(2)            : Fermi level parameter of Kim functional
  ! real*8 rcoor             : Cutoff radius of LWF's (Bohr)
  ! integer ioptlwf          : Option to build LWF's according to:
  !                             0 = Read blindly from disk
  !                             1 = Functional of Kim et al.
  !                             2 = Functional of Ordejon-Mauri
  ! logical chebef          : Compute the chemical potential 
  ! logical noeta            : Use computed Chem.pot. instead of eta
  ! real*8 rcoorcp           : Cutoff (Bohr) to compute the chem.pot.
  ! real*8 beta              : Inverse temperature to compute chem.pot.
  ! integer pmax             : Order of Chebi expansion for chem.pot.
  ! integer idyn             : Atomic dynamics option:
  !                             0 = Geometry optimization
  !                             1 = Standard MD run (Verlet)
  !                             2 = Nose thermostat MD
  !                             3 = Parrinello-Rahman MD
  !                             4 = Nose thermostat + Parrinello-Rahman MD
  !                             5 = Annealing MD
  !                             6 = Force constants
  !                             7 = Deprecated (Forces for PHONON program)
  !                             8 = Force evaluation
  !                             9 = Explicit set of coordinates
  !                            10 = Lua controlled dynamics
  ! integer istart           : Initial time step for MD
  ! integer ifinal           : Final time step for MD
  ! integer nmove            : Number of steps in *any* MD/optimization
  ! real*8 ftol              : Maximum force for structural optimization
  ! real*8 strtol            : Maximum stress for structural optimization
  ! integer ianneal          : Annealing option for idyn = 5
  !                             1 = Temperature 
  !                             2 = Pressure
  !                             3 = Temperature and Pressure
  ! integer iquench          : Quench option: 0 = No; 1 = Yes; 2 = Fire
  ! real*8 dt                : Length of time step (fs)
  ! real*8 dx                : Atomic displacement for Force Constants
  !                             calculation
  ! integer ia1              : First atom to displace for force constants
  ! integer ia2              : Last atom to displace for force constants
  ! real*8 dxmax             : Maximum atomic displacement in one atomic move
  ! real*8 tt                : Target temperature (Kelvin)
  ! real*8 tp                : Target Pressure (Ry/Bohr**3)
  ! real*8 mn                : Mass of Nose variable (Ry/fs**2)
  ! real*8 mpr               : Mass of Parrinello-R. variable (Ry/fs**2)
  ! real*8 bulkm             : Estimate of bulk modulus (Ry/Bohr**3)
  ! real*8 taurelax          : Annealing time to reach targer T and P (fs)
  ! logical usesavelwf       : True = try to use continuation LWF files 
  !                              from disk
  ! logical usesavedm        : True = try to use continuation DM files 
  !                              from disk
  ! logical usesavecg        : True = try to use continuation CG files
  !                              from disk
  ! integer mullipop         : Option for Mulliken Pop. analysis
  ! logical init_anti_ferro  : Spin initialization for spin-polarized
  !                              .true.  -> Antiferro
  !                              .false. -> Ferro
  ! integer maxsav           : Number of density-matrices stored for Pulay
  !                            mixing. .lt.2 => linear mixing only
  !                                    .ge.2 => pulay mixing
  ! integer nkick            : Perform a linear mixing eack nkick scf cycles
  ! real*8 wmixkick          : Mixing parameter for linear mixing each nkick scf
  !                            cycles
  ! logical pulfile          : Use file (.true.) or memory (.false.)
  !                            to store Pulay miximg intermediate vectors
  !                            Default: .false.
  ! real*8 tempinit          : Initial temperature (Kelvin) of the MD simulation
  ! logical dumpcharge       : True: Dump information to plot charge contours
  !                            by the external DENCHAR application program.
  !     (This is now obsolete: info will appear in .RHO file)
  ! logical varcel           : variable shape for optimization or dynamics
  ! logical harrisfun        : swith that indicates if harris functional will
  !                            be used or not
  ! real*8  occtol           : Occupancy threshold for DM build
  ! integer broyden_maxit    : Number of histories saved in Broyden SCF mixing
  ! logical require_energy_convergence  : Impose E. conv. criterion?
  ! logical broyden_optim    : Use Broyden method for optimization
  ! logical want_domain_decomposition:  Use domain decomposition for orbitals in O(N)
  ! logical want_spatial_decomposition:  Use spatial decomposition for orbitals in O(N)


  !----------------------------------------------------------- Local Variables
  real(dp) :: tcp

  character(len=22) :: annop, dyntyp
  character(len=13) :: lwfopt
  character(len=30) :: ctmp
  character(len=6) :: method

  logical :: qnch, qnch2
  logical :: tBool

  !--------------------------------------------------------------------- BEGIN
  ! New template, using fdf
  !
  !      param = fdf_get('ParamName', param_default)
  !      if (ionode)  write(6,'(a,i)'),
  !     .    'redata: ParamName           = ',param
  !      if (cml_p) call cmlAddParameter(xf=mainXML, name='ParamName',
  !     .                 value=param, dictref='siesta:param')

  !
  !      cml_p is only true in the master node
  !
  ! Start of code

  if (cml_p) then
     call cmlStartParameterList(mainXML, title='Input Parameters')
  endif

  ! for cml output, find the system name & label
  if (cml_p) then
     call cmlAddParameter(xf=mainXML, name='SystemName',             &
          value=trim(sname), dictref='siesta:sname')
     call cmlAddParameter(xf=mainXML, name='SystemLabel',            &
          value=trim(slabel), dictref='siesta:slabel')
  endif

  ! Start by printing out spin-configuration
  call print_spin_options()

  ! H setup only
  h_setup_only = fdf_get('HSetupOnly', .false.)
  if (ionode .and. h_setup_only) then
     write(6,1) 'redata: H Setup Only', h_setup_only
  endif

  ! Type of output
  outlng = fdf_get('LongOutput', .false.)
  if (ionode) then
     write(6,1) 'redata: Long output', outlng
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='LongOutput',            &
          value=outlng, dictRef='siesta:verbosity' )
  endif

  ! Write about Number of species, as before
  if (ionode) then
     write(6,4) 'redata: Number of Atomic Species', ns
  endif

  if (ns .le. 0) then
     call die( 'redata: ERROR: Number of species must be larger than zero.' )
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, title='NumberOfSpecies', &
          value=ns, dictRef='siesta:ns', units="cmlUnits:countable" )
  endif

  ! Dump information to plot charge contours
  ! by the external DENCHAR application program.
  dumpcharge = fdf_get('WriteDenchar',.false.)
  if (ionode) then
     write(6,2) 'redata: Charge density info will appear in .RHO file'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='WriteDenChar', value=dumpcharge)
  endif

  ! Perform Mulliken Population Analysis
  mullipop = fdf_get('WriteMullikenPop', 0)
  if (mullipop == 0 .and. outlng) then
     mullipop = 1
  endif
  ! <L> output
  orbmoms                = fdf_get( 'WriteOrbMom'   , .false. )

  if (ionode) then
     select case (mullipop)
     case(0)
        write(6,3) 'redata: Write Mulliken Pop.','NO'
     case(1)
        write(6,3) 'redata: Write Mulliken Pop.','Atomic and Orbital charges'
     case(2)
        write(6,3)'redata: Write Mulliken Pop.','Atomic and Orbital charges'
        write(6,10) 'plus Atomic Overlap Populations'
     case(3)
        write(6,3)'redata: Write Mulliken Pop.','Atomic and Orbital charges'
        write(6,10) 'plus Atomic Overlap Populations'
        write(6,10) 'plus Orbital Overlap Populations'
     case default
        call die( 'redata: Invalid value for WriteMullikenPop' )
     end select
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='WriteMullikenPop', value=mullipop, &
          units="cmlUnits:dimensionless" )
  endif

  ! Perform Hirshfeld and/or Voronoi Population Analysis
  hirshpop= fdf_get('WriteHirshfeldPop',.false.)
  voropop=  fdf_get('WriteVoronoiPop',.false.)
  partial_charges_at_every_geometry =  &
       fdf_get('PartialChargesAtEveryGeometry',.false.)
  partial_charges_at_every_scf_step =  &
       fdf_get('PartialChargesAtEveryScfStep',.false.)


  if ( fdf_get('Compat.Matel.NRTAB', .false.) ) then
    matel_NRTAB = 128
  else
    matel_NRTAB = 1024
  end if
  if ( IONode ) then
    write(6,4) 'redata: Matel table size (NRTAB)', matel_NRTAB
  end if
  if (cml_p) then
    call cmlAddParameter( xf=mainXML, name='MatelNRTAB',value=matel_NRTAB, &
        dictRef='siesta:matel_nrtab', units="cmlUnits:countable")
  end if

  ! Planewave cutoff of the real space mesh ...
  g2cut = fdf_get('MeshCutoff',300._dp,'Ry')
  if (ionode) then
     write(6,6) 'redata: Mesh Cutoff', g2cut,' Ry'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='MeshCutOff', value=g2cut,     &
          dictRef='siesta:g2max', units='siestaUnits:Ry' )
  endif

  ! Net charge in the cell ...
  charnet = fdf_get('NetCharge',0.0_dp)
  if (ionode) then
     write(6,6) 'redata: Net charge of the system',charnet,' |e|'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='NetCharge', value=charnet, &
          dictRef='siesta:NetCharge', units='siestaUnits:e__')
  endif

  ! SCF Loop parameters ...
  !     Minimum/Maximum number of SCF iterations
  min_nscf = fdf_get('MinSCFIterations',0)
  nscf     = fdf_get('MaxSCFIterations',1000)
  SCFMustConverge = fdf_get('SCF.MustConverge', .true.)
  if (ionode) then
     write(6,4) 'redata: Min. number of SCF Iter',min_nscf
     write(6,4) 'redata: Max. number of SCF Iter',nscf
     if (SCFMustConverge) then
        write(6,2) 'redata: SCF convergence failure will abort job'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='MaxSCFIterations',  &
          value=nscf, dictRef='siesta:maxscf',  &
          units="cmlUnits:countable")
     call cmlAddParameter( xf=mainXML, name='MinSCFIterations',  &
          value=min_nscf, dictRef='siesta:minscf',  &
          units="cmlUnits:countable")
  endif

  call fdf_deprecated('TS.MixH','SCF.Mix')
  call fdf_deprecated('MixHamiltonian','SCF.Mix')
  call fdf_deprecated('MixCharge','SCF.Mix')

  ! Note, since 4.1 mixing the Hamiltonian is the default option!
  mixH = fdf_get('TS.MixH',.true.) ! Catch old-style keyword (prefer new key)
  mixH = fdf_get('MixHamiltonian',mixH)
  mix_charge = fdf_get('MixCharge',.false.)

  if ( mix_charge ) then
     ctmp = 'charge'
  else if ( mixH ) then
     ctmp = 'Hamiltonian'
  else
     ctmp = 'density'
  end if
  
  ctmp = fdf_get('SCF.Mix', trim(ctmp))
  if ( leqi(ctmp, 'charge') .or. &
       leqi(ctmp, 'rho') ) then
     mix_charge = .true.
     mixH = .false.
  else if ( leqi(ctmp, 'Hamiltonian') &
       .or. leqi(ctmp, 'H') ) then
     mix_charge = .false.
     mixH = .true.
  else if ( leqi(ctmp, 'density') &
       .or. leqi(ctmp, 'density-matrix') &
       .or. leqi(ctmp, 'DM') ) then
     mix_charge = .false.
     mixH = .false.
  else
     call die('Unrecognized option for: SCF.Mix. Please see the manual.')
  end if
  
  if ( IONode ) then
  if ( mix_charge ) then
     write(6,3) 'redata: SCF mix quantity', 'charge'
  else if ( mixH ) then
     write(6,3) 'redata: SCF mix quantity', 'Hamiltonian'
  else
     write(6,3) 'redata: SCF mix quantity', 'density-matrix'
  end if
  end if

  ! Options for pre-4.0 compatibility
  compat_pre_v4_DM_H  = fdf_get('Compat-pre-v4-DM-H',.false.)
  mix_after_convergence = fdf_get('SCF.MixAfterConvergence',compat_pre_v4_DM_H)
  recompute_H_after_scf = fdf_get('SCF.Recompute-H-After-Scf',compat_pre_v4_DM_H)

  if (ionode) then
     if (compat_pre_v4_DM_H) then
        write(6,2) ':!:Next two options activated by pre-4.0 compat. switch'
     endif
     write(6,1) 'redata: Mix DM or H after convergence',mix_after_convergence
     write(6,1) 'redata: Recompute H after scf cycle', recompute_H_after_scf
  endif

  ! Pulay mixing, number of iterations for one Pulay mixing (maxsav)
  maxsav = fdf_get('DM.NumberPulay', 0)

  ! Broyden SCF mixing, number of iterations 
  broyden_maxit = fdf_get('DM.NumberBroyden',0)
  ! FIRE SCF mixing, no parameters
  fire_mix = fdf_get('DM.FIRE.Mixing',.false.)

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.NumberPulay',     &
          value=maxsav, dictRef='siesta:maxsav', &
          units="cmlUnits:countable" )
     call cmlAddParameter( xf=mainXML, name='DM.NumberBroyden',    &
          value=broyden_maxit, dictRef='siesta:broyden_maxit', &
          units="cmlUnits:countable" )
  endif

  ! Mix density matrix on first SCF step (mix)
  mix_scf_first = fdf_get('DM.MixSCF1', &
       .not. compat_pre_v4_DM_H)
  mix_scf_first = fdf_get('SCF.Mix.First', mix_scf_first)
  mix_scf_first_force = fdf_get('SCF.Mix.First.Force', .false.)
  if ( mix_scf_first_force ) then
    ! Also set this, to note the user of mixing first SCF regardless
    ! of flag.
    mix_scf_first = .true.
  end if
  if (ionode) then
    write(6,1) 'redata: Mix DM in first SCF step',mix_scf_first
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.MixSCF1',   &
          value=mix_scf_first, dictRef='siesta:mix' )
  endif

  ! Use disk or memory to store intermediate Pulay mixing vectors
  ! (pulfile)
  pulfile = fdf_get('DM.PulayOnFile',.false.)
  if (ionode) then
     if (pulfile) then
        call die( 'redata: Cannot use DM.PulayOnFile=.true.'//&
             'in this version' )
     endif
     write(6,1) 'redata: Write Pulay info on disk',pulfile
  endif
  if (cml_p) then
     call cmlAddParameter(xf=mainXML, name='DM.PulayOnFile',      &
          value=pulfile, dictRef='siesta:pulfile')
  endif

  ! Density Matrix Mixing  (proportion of output DM in new input DM)
  wmix = fdf_get('DM.MixingWeight',0.25_dp)
!!$  if (ionode) then
!!$     write(6,6) 'redata: New DM Mixing Weight',wmix
!!$  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.MixingWeight', &
          value=wmix, dictRef='siesta:wmix', &
          units="cmlUnits:dimensionless" )
  endif

  ! Density Matrix occupancy tolerance
  occtol = fdf_get('DM.OccupancyTolerance',1.e-12_dp)
  if (ionode) then
     write(6,8) 'redata: New DM Occupancy tolerance',occtol
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.OccupancyTolerance', &
          value=occtol, dictRef='siesta:occtol' ,  &
          units="cmlUnits:dimensionless" )
  endif

  ! Perform linear mixing each nkick SCF iterations (to kick system
  ! when it is pinned in a poorly convergent SCF loop)
  nkick = fdf_get('DM.NumberKick',0)
  if (ionode) then
     if (nkick .ge. 1) then
        write(6,5) 'redata: Kick with linear mixing every',nkick,&
             ' iterations'
     else
        write(6,2)'redata: No kicks to SCF'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.NumberKick',     &
          value=nkick, dictRef='siesta:nkick', &
          units="cmlUnits:countable" )
  endif

  ! Density Matrix Mixing each nkick SCF iterations
  wmixkick = fdf_get('DM.KickMixingWeight',0.5_dp)
  if (ionode) then
     write(6,6) 'redata: DM Mixing Weight for Kicks',wmixkick
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.KickMixingWeight',    &
          value=wmixkick, dictRef='siesta:wmixkick',&
          units="cmlUnits:dimensionless" )
  endif


  !-----------------------------------!
  !                                   !
  ! Reading of convergence criteria   !
  !                                   !
  !-----------------------------------!

  ! Harris energy convergence
  ! If this is TRUE, then the following options are defaulted to FALSE.
  converge_Eharr = fdf_get('DM.RequireHarrisConvergence', .false.)
  converge_Eharr = fdf_get('SCF.Harris.Converge', converge_Eharr)
  tolerance_Eharr = fdf_get('DM.HarrisTolerance', 1.e-4_dp*eV, 'Ry' )
  tolerance_Eharr = fdf_get('SCF.Harris.Tolerance', tolerance_Eharr, 'Ry' )
  if ( IONode ) then
     write(6,1) 'redata: Require Harris convergence for SCF', &
          converge_Eharr
     write(6,7) 'redata: Harris energy tolerance for SCF', tolerance_Eharr/eV, ' eV'
  end if

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SCF.Harris.Converge', &
          value=converge_Eharr, &
          dictRef='siesta:ReqHarrisConv' )
     call cmlAddParameter( xf=mainXML, name='SCF.Harris.Tolerance', units='siestaUnits:eV', &
          value=tolerance_Eharr/eV, dictRef='siesta:Harris_tolerance')
  end if
  

  ! Density matrix convergence
  converge_DM = fdf_get('SCF.DM.Converge', .not. converge_Eharr)
  dDtol = fdf_get('DM.Tolerance',1.e-4_dp)
  dDtol = fdf_get('SCF.DM.Tolerance',dDtol)
  if ( IONode ) then
     write(6,1) 'redata: Require DM convergence for SCF', converge_DM
     write(6,11) 'redata: DM tolerance for SCF',dDtol
  end if
  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SCF.DM.Converge', &
          value=converge_DM, &
          dictRef='siesta:ReqDMConv' )
     call cmlAddParameter( xf=mainXML, name='SCF.DM.Tolerance', &
          value=dDtol, dictRef='siesta:dDtol', &
          units='siestaUnits:eAng_3' )
  end if


  ! Energy-density matrix convergence
  converge_EDM = fdf_get('SCF.EDM.Converge', .false.)
  tolerance_EDM = fdf_get('SCF.EDM.Tolerance',1.e-3_dp*eV, 'Ry')
  if ( IONode ) then
     write(6,1) 'redata: Require EDM convergence for SCF', converge_EDM
     write(6,7) 'redata: EDM tolerance for SCF',tolerance_EDM/eV, ' eV'
  end if
  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SCF.EDM.Converge', &
          value=converge_EDM, &
          dictRef='siesta:ReqEDMConv' )
     call cmlAddParameter( xf=mainXML, name='SCF.EDM.Tolerance', &
          value=tolerance_EDM/eV, dictRef='siesta:EDM_tolerance', &
          units='siestaUnits:eVeAng_3' )
  end if


  ! Hamiltonian convergence
  if ( compat_pre_v4_DM_H ) then
     tBool = .false.
  else
     tBool = .not. converge_Eharr
  end if
  converge_H = fdf_get('SCF.H.Converge', tBool)
  dHtol = fdf_get('SCF.H.Tolerance',1.e-3_dp*eV, 'Ry')
  if ( IONode ) then
     write(6,1) 'redata: Require H convergence for SCF', converge_H
     write(6,7) 'redata: Hamiltonian tolerance for SCF', dHtol/eV, ' eV'
  end if
  
  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SCF.H.Converge', &
          value=converge_H, &
          dictRef='siesta:ReqHConv' )
     call cmlAddParameter( xf=mainXML, name='SCF.H.Tolerance',     &
          value=dHtol/eV, dictRef='siesta:dHtol', &
          units='siestaUnits:eV' )
  end if


  ! Free energy convergence
  converge_FreeE = fdf_get('DM.RequireEnergyConvergence', .false.)
  converge_FreeE = fdf_get('SCF.FreeE.Converge',converge_FreeE)
  tolerance_FreeE = fdf_get('DM.EnergyTolerance', 1.e-4_dp*eV, 'Ry' )
  tolerance_FreeE = fdf_get('SCF.FreeE.Tolerance',tolerance_FreeE, 'Ry')
  if ( IONode ) then
     write(6,1) 'redata: Require (free) Energy convergence for SCF', converge_FreeE
     write(6,7) 'redata: (free) Energy tolerance for SCF', tolerance_FreeE/eV, ' eV'
  end if

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SCF.FreeE.Converge', &
          value=converge_FreeE, &
          dictRef='siesta:ReqEnergyConv' )
     call cmlAddParameter( xf=mainXML, name='SCF.FreeE.Tolerance', &
          value=tolerance_FreeE/eV, dictRef='siesta:dEtol', &
          units="siestaUnits:eV" )
  end if
  
  ! Check that there indeed is at least one convergence criteria
  tBool = .false.
  tBool = tBool .or. converge_Eharr
  tBool = tBool .or. converge_FreeE
  tBool = tBool .or. converge_EDM
  tBool = tBool .or. converge_DM
  tBool = tBool .or. converge_H
  if ( .not. tBool ) then
     call die('There is no convergence criteria. Please at least have one.')
  end if
  
  !------------------------------------
  ! Done reading convergence criteria
  !------------------------------------

  ! Monitor forces and stresses during SCF loop
  monitor_forces_in_scf = fdf_get('MonitorForcesInSCF',.false.)
  monitor_forces_in_scf = fdf_get('SCF.MonitorForces',monitor_forces_in_scf)

  !--------------------------------------
  ! Initial spin density: Maximum polarization, Ferro (false), AF (true)
  if ( nspin .eq. 2 ) then
     init_anti_ferro = fdf_get('DM.InitSpin.AF',.false.)
     if ( ionode ) then
        write(6,1) 'redata: Antiferro initial spin density',init_anti_ferro
     endif
     if (cml_p) then
        call cmlAddParameter( xf=mainXML, name='DM.InitSpinAF',   &
             value=init_anti_ferro, dictRef='siesta:inspn')
     end if
  end if

  ! Use Saved Data
  usesaveddata = fdf_get('UseSaveData',.false.)
  if (ionode) then
     write(6,1) 'redata: Using Saved Data (generic)', usesaveddata
  endif

  ! Use continuation DM files
  usesavedm = fdf_get('DM.UseSaveDM',usesaveddata)
  if (ionode) then
     write(6,1) 'redata: Use continuation files for DM',  usesavedm
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.UseSaveDM',            &
          value=usesavedm, dictRef='siesta:usesavedm')
  endif

  ! Neglect Interactions between non-overlapping orbitals ...
  negl = fdf_get('NeglNonOverlapInt',.false.)
  if (ionode) then
     write(6,1) 'redata: Neglect nonoverlap interactions',negl
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='NeglNonOverlapInt', &
          value=negl, dictRef='siesta:negl' )
  endif

  ! Method to Solve LDA Hamiltonian ...
  method = fdf_get('SolutionMethod','diagon')
  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SolutionMethod',        &
          value=method, dictRef='siesta:SCFmethod' )
  endif

  if (leqi(method,'matrix')) then
     isolve = MATRIX_WRITE
     if (ionode)  then
        write(*,3) 'redata: Method of Calculation', 'Matrix write only'
     endif
     
  else if (leqi(method,'diagon')) then
     isolve = SOLVE_DIAGON
     if (ionode)  then
        write(*,3) 'redata: Method of Calculation', 'Diagonalization'
     endif
     
  else if (leqi(method,'ordern')) then
     isolve = SOLVE_ORDERN
     if (ionode) then
        write(*,3) 'redata: Method of Calculation','Order-N'
     endif
     if (nspin .gt. 2) then
        call die( 'redata: You chose the Order-N solution option '// &
             'together with nspin>2.  This is not allowed in '//&
             'this version of siesta' )
     endif
     
  else if (leqi(method,'omm')) then
     isolve = SOLVE_MINIM
     call_diagon_default=fdf_integer('OMM.Diagon',0)
     call_diagon_first_step=fdf_integer('OMM.DiagonFirstStep',call_diagon_default)
     minim_calc_eigenvalues=fdf_boolean('OMM.Eigenvalues',.false.)
     if (ionode) then
        write(*,3) 'redata: Method of Calculation', 'Orbital Minimization Method'
     endif
  else if (leqi(method,"pexsi")) then
#ifdef SIESTA__PEXSI     
     isolve = SOLVE_PEXSI
     if (ionode) then
        call add_citation("10.1088/0953-8984/26/30/305503")
        write(*,3) 'redata: Method of Calculation', 'PEXSI'
     endif
#else
     call die("PEXSI solver is not compiled in. Use -DSIESTA__PEXSI")
#endif
     
  else if (leqi(method,'transi') .or. leqi(method,'transiesta') &
       .or. leqi(method,'negf') ) then
     isolve = SOLVE_TRANSI
     if (ionode) then
        call add_citation("10.1103/PhysRevB.65.165401")
        call add_citation("10.1016/j.cpc.2016.09.022")
        write(*,3) 'redata: Method of Calculation','Transiesta'
     endif
     if ( nspin > 2 ) then
       call die('transiesta does not work for non-collinear or spin-orbit')
     end if
  else
     call die( 'redata: The method of solution must be either '//&
          'Transiesta, '//&
#ifdef SIESTA__PEXSI
          'PEXSI, '//&
#endif
          'OrderN, OMM or Diagon' )
  endif

#ifdef DEBUG
  call write_debug( '    Solution Method: ' // method )
#endif

  ! Electronic temperature for Fermi Smearing ...
  temp = fdf_get('ElectronicTemperature',1.9e-3_dp,'Ry')
  if (ionode .and. isolve.eq.SOLVE_DIAGON) then
     write(6,6) 'redata: Electronic Temperature',temp/Kelvin,' K'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='ElectronicTemperature', &
          value=temp, dictRef='siesta:etemp',       &
          units='siestaUnits:Ry')
  endif

  ! Fix the spin of the system to a given value ; and
  ! value of the Spin of the system (only used if fixspin = TRUE)
  fixspin = fdf_get('FixSpin',.false.)
  fixspin = fdf_get('Spin.Fix',fixspin)
  if (ionode) then
     write(6,1) 'redata: Fix the spin of the system',fixspin 
  endif

  if (fixspin) then
    if (nspin .ne. 2) then
      call die( 'redata: ERROR: You can only fix the spin of '//&
          'the system for collinear spin polarized calculations.' )
    endif
    total_spin = fdf_get('TotalSpin',0.0_dp)
    total_spin = fdf_get('Spin.Total',total_spin)
    if (ionode) then
      write(6,9) 'redata: Total spin of the system (spin value)', total_spin
    endif
  else
    total_spin = 0.0_dp
  endif

  if (cml_p) then 
     call cmlAddParameter( xf=mainXML, name='FixSpin', &
          value=fixspin, dictref='siesta:fixspin' )
     call cmlAddParameter( xf=mainXML, name='TotalSpin', &
          value=total_spin, dictref='siesta:totalspin',&
          units='siestaUnits:eSpin' )
  endif

  ! Order-N solution parameters ...
  !     Maximum number of CG minimization iterations
  ncgmax = fdf_get('ON.MaxNumIter',1000)
  if (ncgmax<1) then
     if (ionode) then
        write(6,2) 'ON.MaxNumIter cannot be less than 1.  Resetting to 1'
     endif
     ncgmax = 1 
  endif

  ! Relative tolerance in total band structure energy
  etol = fdf_get('ON.etol',1.e-8_dp)

  ! Fermi level parameter
  eta(1:2) = 0.0_dp
  eta(1) = fdf_physical('ON.eta',eta(1),'Ry')
  eta(2) = eta(1)
  eta(1) = fdf_physical('ON.eta_alpha',eta(1),'Ry')
  eta(2) = fdf_physical('ON.eta_beta',eta(2),'Ry')

  ! Cutoff radius for Localized Wave Functions
  rcoor = fdf_get('On.RcLWF',9.5_dp,'Bohr')

  ! Use continumation LWF files
  usesavelwf = fdf_get('ON.UseSaveLWF',usesaveddata)

  ! Option on how to build LWF's (disk or functionals)
  lwfopt = fdf_get('ON.functional','kim')
  if (leqi(lwfopt,'files')) then
     ioptlwf = 0
  else if (leqi(lwfopt,'kim')) then
     ioptlwf = 1
  else if (leqi(lwfopt,'ordejon-mauri')) then
     ioptlwf = 2
  else
     call die('redata: wrong ON.funcional option')
  endif

  ! Option to calculate the Chemical potential in O(N)
  ! Option to use the Chemical Potential calculated instead
  ! of the eta variable of the input
  noeta = fdf_get('ON.ChemicalPotentialUse',.false.)
  ! NOTE: This does not yet work in parallel

  if (noeta) then
     ! if so, we must (obviously) calculate the chemical potential
     chebef=.true.
  else
     ! otherwise, we may still want to calculate it but not use it.
     chebef = fdf_get('ON.ChemicalPotential',.false.)
  endif

#ifdef MPI
  if (chebef) then
     call die("ON.ChemicalPotential(Use) options do not work with MPI")
  endif
#endif

  ! Cutoff radius to calculate the Chemical Potential by projection
  rcoorcp = fdf_get( 'ON.ChemicalPotentialRc', 9.5_dp, 'Bohr' )

  ! Temperature of the Fermi distribution to calculate the
  ! Chemical potential by projection
  tcp = fdf_get( 'ON.ChemicalPotentialTemperature', 0.05_dp,'Ry' )
  beta = 1.0_dp/tcp

  ! Order of the Chebishev expansion to calculate the Chemical potential
  pmax = fdf_get('ON.ChemicalPotentialOrder',100)


  if (isolve==SOLVE_ORDERN) then
     if (ionode) then
        write(6,4) 'redata: Maximum number of iterations',ncgmax
        write(6,6) 'redata: Relative tolerance',etol
        if (nspin.eq.2) then
           write(6,6) 'redata: Eta (Fermi level) Alpha spin',eta(1)/eV,' eV'
           write(6,6) 'redata: Eta (Fermi level) Beta spin',eta(2)/eV,' eV'
        else
           write(6,6) 'redata: Eta (Fermi level parameter)',eta(1)/eV,' eV'
        endif
        write(6,6) 'redata: Radius of LWFs',rcoor/Ang,' Ang'
        write(6,1) 'redata: Use continuation files for LWF',usesavelwf
        write(6,3) 'redata: Method to build LWFs',lwfopt
        if (chebef) then
           write(6,1) 'redata: Compute Chemical Potential',chebef
           write(6,2) 'redata: Use the calculated Chemical ..'
           write(6,1) 'redata: ..Potential instead of eta',noeta
           write(6,6) 'redata: Radius to compute the Chem. Pot.',rcoorcp/Ang,' Ang'
           write(6,2) 'redata: Temp. for Fermi distribution ..'
           write(6,6) 'redata: .. to compute the Chem. Pot.',tcp/eV,' eV'
           write(6,4) 'redata: Order of the Chebishev expansion',pmax
        endif
     endif
     if (cml_p) then
        call cmlAddParameter( xf      = mainXML,        &
             name    = 'ON.MaxNumIter',&
             value   = ncgmax,         &
             dictref = 'siesta:ncgmax', &
             units   = "cmlUnits:countable" )

        call cmlAddParameter( xf      = mainXML,       &
             name    = 'ON.etol',     &
             value   = etol,          &
             dictref = 'siesta:etol', &
             units   = "siestaUnits:eV" )
        if (nspin==2) then
           call cmlAddParameter( xf      = mainXML,          &
                name    = 'ON.eta_alpha',   &
                value   = eta(1),           &
                dictref = 'siesta:eta1',    &
                units   = 'siestaUnits:Ry' )

           call cmlAddParameter( xf      = mainXML,          &
                name    = 'ON.eta_beta',    &
                value   = eta(2),           &
                dictref = 'siesta:eta2',    &
                units   = 'siestaUnits:Ry' )
        else
           call cmlAddParameter( xf      = mainXML,         &
                name    = 'ON.eta',        &
                value   = eta(1),          &
                dictref = 'siesta:eta',    &
                units   = 'siestaUnits:Ry')
        endif
        call cmlAddParameter( xf      = mainXML,          &
             name    = 'On.RcLWF',       &
             value   = rcoor,            &
             dictref = 'siesta:rcoor',   &
             units   = 'siestaUnits:Bohr')

        call cmlAddParameter( xf=mainXML,                 &
             name='On.UseSaveLWF',       &
             value=usesavelwf,           &
             dictref='siesta:usesavelwf' )

        call cmlAddParameter( xf      = mainXML,        &
             name    = 'ON.functional',&
             value   = lwfopt,         &
             dictref = 'siesta:lwfopt' )
        if (chebef) then
           call cmlAddParameter( xf      = mainXML,                &
                name    = 'ON.ChemicalPotential', &
                value   = chebef,                 &
                dictref = 'siesta:chebef')

           call cmlAddParameter( xf      = mainXML,                   &
                name    = 'ON.ChemicalPotentialUse', &
                value   = noeta,                     &
                dictref = 'siesta:noeta')

           call cmlAddParameter( xf      = mainXML,                  &
                name    = 'ON.ChemicalPotentialRc', &
                value   = rcoorcp,                  &
                dictref = 'siesta:rcoorcp',         &
                units   = 'siestaUnits:Bohr')

           call cmlAddParameter( xf    = mainXML,                           &
                name  = 'ON.ChemicalPotentialTemperature', &
                value = tcp, dictref='siesta:tcp',         &
                units = 'siestaUnits:Ry' )

           call cmlAddParameter( xf      = mainXML,                     &
                name    = 'ON.ChemicalPotentialOrder', &
                value   = pmax,                        &
                dictref = 'siesta:pmax',               &
                units   = 'cmlUnits:dimensionless')
        endif
     endif
  endif

  ! Dynamics parameters ...
  varcel = fdf_get('MD.VariableCell', .false. )

  ! NB reset below ...
  ! Type of dynamics 

  compat_pre_v4_dynamics = fdf_get('compat-pre-v4-dynamics', .false. )
  if (compat_pre_v4_dynamics) then
     dyntyp = fdf_get('MD.TypeOfRun','verlet')
  else
     dyntyp = fdf_get('MD.TypeOfRun','cg')
  endif

  if (leqi(dyntyp,'cg')) then
     idyn = 0
     usesavecg = fdf_get('MD.UseSaveCG', usesaveddata)
     ! Support the old Broyden switch  for now
     broyden_optim = fdf_get('Optim.Broyden',.false.)
     call deprecated('Optim.Broyden')

  else if (leqi(dyntyp,'broyden')) then
     idyn = 0
     broyden_optim = .true.
  else if (leqi(dyntyp,'fire')) then
     idyn = 0
     fire_optim = .true.
  else if (leqi(dyntyp,'verlet')) then
     idyn = 1
  else if (leqi(dyntyp,'nose')) then
     idyn = 2
  else if (leqi(dyntyp,'parrinellorahman')) then
     idyn = 3
  else if (leqi(dyntyp,'noseparrinellorahman')) then
     idyn = 4
  else if (leqi(dyntyp,'anneal')) then
     idyn = 5
  else if (leqi(dyntyp,'fc')) then
     idyn = 6
  else if (leqi(dyntyp,'phonon')) then
     call die('Dynamics type "PHONON" is no longer supported')
  else if (leqi(dyntyp,'forces').or.leqi(dyntyp,'master')) then
     idyn = 8
#ifdef NCDF_4
    else if (leqi(dyntyp,'explicit')) then
      idyn = 9
#endif
#ifdef SIESTA__FLOOK
   else if (leqi(dyntyp,'lua')) then
      idyn = 10
#endif
  else
     call die('Invalid Option selected - value of MD.TypeOfRun not recognised')
  endif

  ! Maximum number of steps in MD/coordinate optimization
  nmove = fdf_get('MD.NumCGsteps',0)
  nmove = fdf_get('MD.Steps',nmove)

  ! Maximum atomic displacement in one step
  dxmax = fdf_get('MD.MaxCGDispl',0.2_dp,'Bohr')
  dxmax = fdf_get('MD.MaxDispl',dxmax,'Bohr')

  ! Tolerance in the maximum atomic force [0.04 eV/Ang]
  ftol = fdf_get('MD.MaxForceTol', 0.00155574_dp, 'Ry/Bohr')

  ! Tolerance in the maximum residual stress (var cell) [1 GPa]
  strtol = fdf_get('MD.MaxStressTol', 6.79773e-5_dp, 'Ry/Bohr**3')
  strtol = abs(strtol)
  
  GeometryMustConverge = fdf_get('GeometryMustConverge', .false.)

  if (ionode) then
     select case (idyn)
     case(0)
        if (nmove > 0) then
           if (broyden_optim) then
              write(6,3) 'redata: Dynamics option','Broyden coord. optimization'
           elseif (fire_optim) then
              write(6,3) 'redata: Dynamics option', 'FIRE coord. optimization'
           else
              write(6,3) 'redata: Dynamics option','CG coord. optimization'
           endif
           write(6,1) 'redata: Variable cell', varcel
           if (.not. broyden_optim) then
              write(6,1) 'redata: Use continuation files for CG', usesavecg
              write(6,6) 'redata: Max atomic displ per move', dxmax/Ang, ' Ang'
           endif
           write(6,4) 'redata: Maximum number of optimization moves', nmove
           write(6,6) 'redata: Force tolerance', ftol/eV*Ang, ' eV/Ang'
           if (varcel) then
              write(6,6) 'redata: Stress tolerance', &
                   strtol/6.79773e-5_dp, ' GPa'
           endif
           if (cml_p) then
              if (broyden_optim) then
                 call cmlAddParameter( xf   = mainXML,        &
                      name = 'MD.TypeOfRun', &
                      value= 'Broyden' )
              else if (fire_optim) then
                 call cmlAddParameter( xf   = mainXML,        &
                      name = 'MD.TypeOfRun', &
                      value= 'FIRE' )
              else
                 call cmlAddParameter( xf    =mainXML,        &
                      name  ='MD.TypeOfRun', &
                      value ='CG' )

                 call cmlAddParameter( xf    = mainXML,       &
                      name  = 'MD.UseSaveCG',&
                      value = usesavecg )
              endif

              call cmlAddParameter( xf    = mainXML,         &
                   name  = 'MD.NumCGSteps', &
                   value = nmove,           &
                   units = "cmlUnits:countable" )
              call cmlAddParameter( xf    = mainXML,         &
                   name  = 'MD.Steps', &
                   value = nmove,           &
                   units = "cmlUnits:countable" )

              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.MaxCGDispl',   &
                   value = dxmax,             &
                   units = 'siestaUnits:Bohr' )
              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.MaxDispl',   &
                   value = dxmax,             &
                   units = 'siestaUnits:Bohr' )

              call cmlAddParameter( xf=mainXML,                 &
                   name='MD.MaxForceTol',      &
                   value=ftol,                 &
                   units='siestaUnits:Ry_Bohr')
              if (varcel) then
                 call cmlAddParameter( xf=mainXML,                     &
                      name='MD.MaxStressTol',         &
                      value=strtol,                   &
                      units='siestaUnits:Ry_Bohr__3' )
              endif
           endif
        else
           write(6,3) 'redata: Dynamics option','Single-point calculation'
           if (cml_p) then
              call cmlAddParameter( xf   = mainXML,        &
                   name = 'MD.TypeOfRun', &
                   value= 'Single-Point' )
           endif
        endif
     case(1)
        write(6,3) 'redata: Dynamics option', 'Verlet MD run'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,        &
                name  = 'MD.TypeOfRun', &
                value = 'Verlet' )
        endif

     case(2)
        write(6,3) 'redata: Dynamics option', 'Nose thermostat MD run'
        if (cml_p) then
           call cmlAddParameter( xf=mainXML, name='MD.TypeOfRun', value='Nose')
        endif

     case(3)
        write(6,3) 'redata: Dynamics option', 'Parrinello-Rahman MD run'
        if (cml_p) then
           call cmlAddParameter( xf= mainXML,              &
                name= 'MD.TypeOfRun',     &
                value= 'Parrinello-Rahman' )
        endif

     case(4)
        write(6,3) 'redata: Dynamics option', 'Nose-Parrinello-Rahman MD run'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,                &
                name  = 'MD.TypeOfRun',         &
                value = 'Nose-Parrinello-Rahman' )
        endif

     case(5)
        write(6,3) 'redata: Dynamics option', 'Annealing MD run'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,        &
                name  = 'MD.TypeOfRun', &
                value = 'Annealing' )
        endif

     case(6)
        write(6,3) 'redata: Dynamics option', 'Force Constants Matrix Calculation'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,          &
                name  = 'MD.TypeOfRun',   &
                value = 'Force Constants' )
        endif

     case(7)
        ! deprecated

     case(8)
        write(6,3) 'redata: Dynamics option','Force evaluation'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,            &
                name  = 'MD.TypeOfRun',     &
                value = 'Force Evaluation' )
        endif
#ifdef NCDF_4
      case(9)
        write(6,3) 'redata: Dynamics option','Explicit'
        if (cml_p) then
          call cmlAddParameter( xf    = mainXML,            &
                                name  = 'MD.TypeOfRun',     &
                                value = 'Explicit' )
        endif
#endif
#ifdef SIESTA__FLOOK
     case(10)
        write(6,3) 'redata: Dynamics option','LUA'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,        &
                name  = 'MD.TypeOfRun', &
                value = 'LUA' )
        endif
#endif
     end select
  endif

  ! Initial and final time steps for MD
  istart = fdf_get('MD.InitialTimeStep',1)
  if ( fdf_defined('MD.Steps') ) then
    ifinal = fdf_get('MD.FinalTimeStep',max(1,nmove - istart + 1))
  else
    ifinal = fdf_get('MD.FinalTimeStep',1)
  end if

  ! Length of time step for MD
  dt = fdf_get('MD.LengthTimeStep',1._dp,'fs')

  ! Quench Option
  qnch  = fdf_get('MD.Quench',.false.)
  qnch2 = fdf_get('MD.FireQuench',.false.)
  if ((qnch .or. qnch2) .and. (idyn==2 .or. idyn==4)) then 
     call die( 'redata: ERROR: You cannot quench and '//&
          'use a Nose thermostat simultaneously')
  endif

  iquench = 0
  if (qnch) then
     iquench = 1
  endif
  if (qnch2) then
     iquench = 2
  endif

  ! Initial Temperature of MD simulation
  ! (draws random velocities from the Maxwell-Boltzmann distribition
  !  at the given temperature)
  tempinit = fdf_get('MD.InitialTemperature',0.0_dp,'K')

  if (idyn .ge. 1 .and. idyn .le. 5) then
     if (ionode) then
        write(6,4) 'redata: Initial MD time step',istart
        write(6,4) 'redata:   Final MD time step',ifinal
        write(6,6) 'redata: Length of MD time step',dt,' fs'
        write(6,6) 'redata: Initial Temperature of MD run',tempinit,' K'
        if (idyn .ne. 5) then 
           if (qnch2) then
              write(6,1) 'redata: Perform a MD Fire quench',qnch2
           else
              write(6,1) 'redata: Perform a MD quench',qnch
           endif
        endif
     endif

     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,             &
             name  = 'MD.InitialTimeStep',&
             value = istart,              &
             units = 'cmlUnits:countable' )

        call cmlAddParameter( xf    = mainXML,             &
             name  = 'MD.FinalTimeStep',  &
             value = ifinal,              &
             units = 'cmlUnits:countable' )

        call cmlAddParameter( xf=mainXML,              &
             name='MD.LengthTimeStep',&
             value=dt,                &
             units='siestaUnits:fs' )

        call cmlAddParameter( xf=mainXML,                  &
             name='MD.InitialTemperature',&
             value=tempinit,              &
             units='siestaUnits:K' )
        if (idyn/=5) then 
           if(qnch2) then
              call cmlAddParameter( xf    = mainXML,         &
                   name  = 'MD.FireQuench', &
                   value = qnch2 )
           else
              call cmlAddParameter( xf    = mainXML,     &
                   name  = 'MD.Quench', &
                   value = qnch )
           endif
        endif
     endif
  endif

  ! Target Temperature and Pressure
  tt = fdf_get('MD.TargetTemperature',0.0_dp,'K')
  tp = fdf_get('MD.TargetPressure',0.0_dp,'Ry/Bohr**3')
  !
  ! Used for now for the call of the PR md routine if quenching
  if (idyn == 3 .AND. iquench > 0) call set_target_stress()


  ! Mass of Nose variable
  mn = fdf_get('MD.NoseMass',100._dp,'Ry*fs**2')

  ! Mass of Parrinello-Rahman variables
  mpr = fdf_get('MD.ParrinelloRahmanMass',100._dp,'Ry*fs**2')

  if (idyn==2 .or. idyn==4) then
     if (ionode) then
        write(6,6) 'redata: Nose mass',mn/eV,' eV/fs**2'
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,              &
             name  = 'MD.NoseMass',        &
             value = mn,                   &
             units = 'siestaUnits:Ry_fs__2')
     endif
  endif

  if (idyn==3 .or. idyn==4) then
     if (ionode) then
        write(6,6) 'redata: Parrinello-Rahman mass',mpr/eV,' eV/fs**2'
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,                   &
             name  = 'MD.ParrinelloRahmanMass', &
             value = mn,                        &
             units = 'siestaUnits:Ry_fs__2' )
     endif
  endif

  ! Annealing option
  ianneal = 0
  annop = fdf_get( 'MD.AnnealOption','TemperatureAndPressure' )

  if (idyn .eq. 5) then
     if (leqi(annop,'Temperature')) then
        ianneal = 1
     else if (leqi(annop,'Pressure')) then
        ianneal = 2
     else if (leqi(annop,'TemperatureAndPressure')) then
        ianneal = 3
     else
        call die( 'redata: ERROR: With annealing MD, you must '//&
             'choose an appropriate value for MD.AnnealOption' )
     endif

     if (ionode) then
        select case (ianneal)
        case(1)
           write(6,3) 'redata: Annealing Option', 'Temperature'
           if (cml_p) then
              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.AnnealOption', &
                   value = 'Temperature' )
           endif

        case(2)
           write(6,3) 'redata: Annealing Option', 'Pressure'
           if (cml_p) then
              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.AnnealOption', &
                   value = 'Pressure')
           endif

        case(3)
           write(6,3) 'redata: Annealing Option', 'Temperature and Pressure'
           if (cml_p) then
              call cmlAddParameter( xf    = mainXML,                &
                   name  = 'MD.AnnealOption',      &
                   value = 'TemperatureAndPressure')
           endif
        end select
     endif
  endif


  if (idyn==2 .or. idyn==4 .or. (idyn==5 .and. (ianneal ==1 .or. ianneal==3))) then
     if (ionode) then
        write(6,6) 'redata: Target Temperature',tt,' Kelvin'
     endif

     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,                &
             name  = 'MD.TargetTemperature', &
             value = tt,                     &
             units = 'siestaUnits:K' )
     endif
  endif

  if (idyn==3 .or. idyn==4 .or. (idyn==5 .and. (ianneal==2 .or. ianneal==3))) then
     if (ionode) then
        write(6,6) 'redata: Target Pressure', tp/eV*Ang**3, ' eV/Ang**3'
     endif

     if (cml_p) then
        call cmlAddParameter( xf=mainXML,                     &
             name= 'MD.TargetPressure',      &
             value= tp,                      &
             units= 'siestaUnits:Ry_Bohr__3' )
     endif
  endif

  ! Relaxation Time for Annealing
  taurelax = fdf_get( 'MD.TauRelax', 100._dp,'fs' )
  if (idyn==5) then
     if (ionode) then
        write(6,6) 'redata: Annealing Relaxation Time', taurelax,' fs'
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,          &
             name  = 'MD.TauRelax',    &
             value = taurelax,         &
             units = 'siestaUnits:fs' )
     endif
  endif

  ! Estimated Bulk modulus (for Pressure annealing) [100 GPa]
  bulkm = fdf_get( 'MD.BulkModulus',100*6.79773e-5_dp,'Ry/Bohr**3' )
  if (ionode) then
     if (idyn==5 .and. (ianneal==2 .or. ianneal==3)) then
        write(6,6) 'redata: Approx. Bulk Modulus ', bulkm/eV*Ang**3, ' eV/Ang**3'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf    = mainXML,                &
          name  = 'MD.BulkModulus',       &
          value = bulkm,                  &
          units = 'siestaUnits:Ry_Bohr__3')
  endif

  ! Atomic displacement for force constant calculation
  dx = fdf_get('MD.FCDispl',0.04_dp,'Bohr')

  ! First and last atoms to displace for calculation of force constants
  ia1 = fdf_get('MD.FCfirst',1)
  ia2 = fdf_get('MD.FClast',na)

  ! Check that last atom doesn't exceed total number
  if (idyn.eq.6.and.ia2.gt.na) then
     call die( 'redata: ERROR:'//&
          'Last atom index for FC calculation is > number of atoms.')
  endif

  if (idyn==6) then
     if (ionode) then
        write(6,6) 'redata: Atomic displ for force constants',dx/Ang,'  Ang'
        write(6,4) 'redata: First atom to move',ia1
        write(6,4) 'redata: Last atom to move',ia2
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,           &
             name  = 'MD.FCDispl',      &
             value = dx,                &
             units = 'siestaUnits:Bohr' )

        call cmlAddParameter( xf= mainXML,        &
             name= 'MD.FCFirst', &
             value= ia1,         &
             units= 'cmlUnits:countable' )

        call cmlAddParameter( xf= mainXML,       &
             name= 'MD.FCLast', &
             value= ia2,        &
             units= 'cmlUnits:countable' )

     endif
  endif

  ! Variable cell shape? Depending on input and type of dynamics
  varcel = varcel .or. (idyn==3) .or. (idyn==4)             &
       .or. (idyn==5 .and. ianneal==1)           &
       .and. (idyn/=1) .and. (idyn/=2)           &
       .and. (idyn/=6) .and. (idyn/=7)           &
       .and. (idyn/=9) &
       .and. (.not. (idyn==5 .and. ianneal/=1) )

  want_spatial_decomposition = fdf_get('UseSpatialDecomposition', .false.)
  want_domain_decomposition = fdf_get('UseDomainDecomposition', .false.)
#if defined(ON_DOMAIN_DECOMP) || defined(SIESTA__METIS)
#else
#ifdef MPI
  if (want_domain_decomposition) then
     write(*,*)'Need to compile SIESTA__METIS or ON_DOMAIN_DECOMP'
     call die("Need to compile with ON_DOMAIN_DECOMP support")
  endif
#endif
#endif

  ! Read in mixing parameters (SCF)
  call mixers_scf_init( nspin )
  call mixers_scf_print( nspin )

  ! We read in relevant data for ChargeGeometries block
  call read_charge_add( min(2, nspin) , charnet )
  
  ! We read in the relevant data for HartreeGeometries block
  call read_hartree_add( )

  ! Harris Forces?. Then DM.UseSaveDM should be false (use always
  ! Harris density in the first SCF step of each MD step), and
  ! MaxSCFIter should be  2, in the second one the Harris 
  ! forces are computed. Also, should not exit if SCF did 
  ! not converge.

  harrisfun = fdf_get('Harris_functional',.false.)

  ! First read in, then later correct depending on
  ! the other usages
  use_aux_cell = fdf_get('ForceAuxCell', .false.)

  if (harrisfun) then
     mixH = .false.
     mix_charge = .false.
     usesavedm = .false.
     nscf      = 1  ! Note change from tradition, since siesta_forces        
     ! now explicitly separates the "compute_forces"        
     ! phase from the rest of the scf cycle.          
     mix_scf_first = .false.
     SCFMustConverge = .false.
     if ( nspin > 2 ) then
        call die("Harris functional does not work for &
             &non-collinear spin polarization.")
     end if
  endif

  ! Warn the user about deprecated symbols...
  call deprecated('Optim.Broyden.History.Steps')
  call deprecated('Optim.Broyden.Cycle.On.Maxit')
  call deprecated('Optim.Broyden.Variable.Weight')
  call deprecated('Optim.Broyden.Debug')
  call deprecated('Optim.Broyden.Initial.Inverse.Jacobian')


  ! Find some switches 
  writek                = fdf_get( 'WriteKpoints', outlng )
  writef                = fdf_get( 'WriteForces', outlng )

  writedm               = fdf_get( 'WriteDM', .true. )
  write_dm_at_end_of_cycle = fdf_get( 'WriteDM.End.Of.Cycle', writedm )
  writeH                = fdf_get( 'WriteH', .false. )
  write_H_at_end_of_cycle  = fdf_get( 'WriteH.End.Of.Cycle', writeH )

  writedm_cdf           = fdf_get('WriteDM.NetCDF', .false. )
  writedm_cdf_history   = fdf_get('WriteDM.History.NetCDF', .false. )
  writedmhs_cdf         = fdf_get('WriteDMHS.NetCDF', .false. )
  writedmhs_cdf_history = fdf_get('WriteDMHS.History.NetCDF', .false.)
  read_charge_cdf       = fdf_get('SCF.Read.Charge.NetCDF' , .false. )
  read_deformation_charge_cdf = &
       fdf_get('SCF.Read.Deformation.Charge.NetCDF', .false. )

  ! Write the history
  write_tshs_history = fdf_get('Write.TSHS.History', .false.)
  if ( IONode.and.write_tshs_history ) &
       write(*,2) 'redata: Saves TSHS files in MD simulation'
  !write_hs_history = fdf_get('Write.HS.History', .false.)

  if (read_charge_cdf .or. read_deformation_charge_cdf) then
     mix_scf_first = .false.
  endif

  save_initial_charge_density = fdf_get(    &
       'SaveInitialChargeDensity' , .false.)

  analyze_charge_density_only = fdf_get(    &
       'AnalyzeChargeDensityOnly' , .false.)

  tBool = fdf_get( 'UseNewDiagk', .false. )
  if ( tBool ) then
    ctmp = fdf_get('Diag.WFS.Cache', 'CDF')
  else
    ctmp = fdf_get('Diag.WFS.Cache', 'none')
  end if
  if ( leqi(ctmp, 'none') ) then
    diag_wfs_cache = 0
  else if ( leqi(ctmp, 'cdf') ) then
    diag_wfs_cache = 1
  else
    call die('redata: ERROR: Diag.WFS.Cache must be one of none|cdf')
  end if

  writb                  = fdf_get( 'WriteBands', outlng )
  writbk                 = fdf_get( 'WriteKbands', outlng )
  writeig                = fdf_get('WriteEigenvalues', outlng )
  writec                 = fdf_get( 'WriteCoorStep', outlng )
  writmd                 = fdf_get( 'WriteMDhistory', .false. )
  writpx                 = fdf_get( 'WriteMDXmol', .not. writec )
  save_ORB_INDX          = fdf_get( 'WriteOrbitalIndex', .true. )
  ! Do options of graphviz
  ctmp = fdf_get( 'Write.Graphviz', 'none' )
  write_GRAPHVIZ = 0
  if ( leqi(ctmp, 'orb') .or. &
       leqi(ctmp, 'orbital') ) then
     write_GRAPHVIZ = 1
  else if ( leqi(ctmp, 'atom') .or. &
       leqi(ctmp, 'atomic') ) then
     write_GRAPHVIZ = 2
  else if ( leqi(ctmp, 'atom+orb') .or. &
       leqi(ctmp, 'atom+orbital') .or. &
       leqi(ctmp, 'orbital+atom') .or. &
       leqi(ctmp, 'orb+atom') .or. leqi(ctmp, 'both') ) then
     write_GRAPHVIZ = 3
  else if ( leqi(ctmp,'none') ) then
     ! do nothing, correct
  else if ( IONode ) then
     write(*,2) 'Write.Graphviz input could not be understood. &
          &Use [orbital|atom|atom+orbital] in the option.'
  end if
  if ( IONode ) then
     select case ( write_GRAPHVIZ )
     case ( 1 ) 
        write(*,2) 'redata: Save orbital connectivity graph in GRAPHVIZ format'
     case ( 2 )
        write(*,2) 'redata: Save atomic connectivity graph in GRAPHVIZ format'
     case ( 3 )
        write(*,2) 'redata: Save atom+orbital connectivity graphs in GRAPHVIZ format'
     end select
  end if


  writec                 = fdf_get( 'WriteCoorStep', outlng )
  savehs                 = fdf_get( 'SaveHS', .false. )
  initdmaux              = fdf_get( 'ReInitialiseDM', .TRUE. )
  allow_dm_reuse         = fdf_get( 'DM.AllowReuse', .TRUE. )
  allow_dm_extrapolation = fdf_get( 'DM.AllowExtrapolation', .TRUE. )
  DM_history_depth       = fdf_get( 'DM.History.Depth', 1)
  dm_normalization_tol   = fdf_get( 'DM.NormalizationTolerance',1.0d-5)
  normalize_dm_during_scf= fdf_get( 'DM.NormalizeDuringSCF',.true.)
  muldeb                 = fdf_get( 'MullikenInSCF'   , .false.)
  spndeb                 = fdf_get( 'SpinInSCF'   , (nspin>1) )

  ! If no mulliken is requested, set it to false
  if ( mullipop == 0 ) muldeb = .false.
  rijmin                 = fdf_get( 'WarningMinimumAtomicDistance', &
       1.0_dp, 'Bohr' )
  bornz                  = fdf_get( 'BornCharge'   , .false. )
  if (idyn.ne.6) then
     bornz = .false.
  endif
  change_kgrid_in_md           = fdf_get('ChangeKgridInMD', .false.)
  RelaxCellOnly                = fdf_get('MD.RelaxCellOnly', .false.)
  RemoveIntraMolecularPressure = fdf_get( &
       'MD.RemoveIntraMolecularPressure', .false.)
  !
  !   COOP-related flags
  !
  write_coop = fdf_get('COOP.Write', .false.)
  !
  saverho  = fdf_get( 'SaveRho', dumpcharge)
  savedrho = fdf_get( 'SaveDeltaRho',       .false. )
  saverhoxc= fdf_get( 'SaveRhoXC', .false.)
  savevh   = fdf_get( 'SaveElectrostaticPotential', .false. )
  savevna  = fdf_get( 'SaveNeutralAtomPotential', .false. )
  savevt   = fdf_get( 'SaveTotalPotential', .false. )
  savepsch = fdf_get( 'SaveIonicCharge', .false. )
  savebader= fdf_get( 'SaveBaderCharge',  .false.)
  savetoch = fdf_get( 'SaveTotalCharge', savebader )

  !
  !   Siesta2Wannier90 -related flags
  !
  w90_write_mmn = fdf_get( 'Siesta2Wannier90.WriteMmn',   .false. )
  w90_write_unk = fdf_get( 'Siesta2Wannier90.WriteUnk',   .false. )
  w90_write_amn = fdf_get( 'Siesta2Wannier90.WriteAmn',   .false. )
  w90_write_eig = fdf_get( 'Siesta2Wannier90.WriteEig',   .false. )

  w90_processing = ( w90_write_mmn .or. w90_write_unk .or. &
       w90_write_amn .or. w90_write_eig )

  hasnobup   = fdf_defined( 'Siesta2Wannier90.NumberOfBandsUp'   )
  hasnobdown = fdf_defined( 'Siesta2Wannier90.NumberOfBandsDown' )
  hasnob     = fdf_defined( 'Siesta2Wannier90.NumberOfBands'     )

  nobup      = fdf_get( 'Siesta2Wannier90.NumberOfBandsUp',   0)
  nobdown    = fdf_get( 'Siesta2Wannier90.NumberOfBandsDown', 0)
  nob        = fdf_get( 'Siesta2Wannier90.NumberOfBands',     0)

#ifdef NCDF_4
  write_cdf = fdf_get('CDF.Save', .false.)
  ! No compression is by far the fastest
  cdf_comp_lvl = fdf_get('CDF.Compress', 0)
  if ( Nodes > 1 ) then
    cdf_w_parallel = fdf_get('CDF.MPI', .false.)
  else
    cdf_w_parallel = .false.
  end if

  ! Only allow writing the CDF file for FC and non-MD calculations
  ! The MD file should be something different that only contains
  ! these things.
  if ( write_cdf ) then
    if ( idyn == 0 .and. nmove == 0 ) then
      ! non-MD calculation
    else if ( idyn == 6 ) then
      ! FC calculation, the FC calculation is fine
      ! Here we disable saving any real-space grid data
      save_initial_charge_density = .false.
      saverho = .false.
      savedrho = .false.
      saverhoxc = .false.
      savevh = .false.
      savevna = .false.
      savevt = .false.
      savepsch = .false.
      savebader = .false.
      savetoch = .false.
    else
      write_cdf = .false.
    end if
  end if
# ifndef NCDF_PARALLEL
  ! If not compiled with NCDF_PARALLEL, we do not
  ! allow parallel writes.....!!!!
  cdf_w_parallel = .false.
# endif
  if ( cdf_w_parallel ) then
     ! Doing parallel writes does not allow
     ! compression (the offsets cannot be calculated)
     cdf_comp_lvl = 0
  end if
  cdf_r_parallel = fdf_get('CDF.Read.Parallel', cdf_w_parallel )

  if ( IONode ) then
     ! Write out
     write(*,1) 'redata: Save all siesta data in one NC',write_cdf
     if ( write_cdf ) then
        if ( grid_p == dp ) then
           ctmp = fdf_get('CDF.Grid.Precision','double')
           if ( leqi(ctmp,'single') .or. leqi(ctmp,'float') ) then
              write(*,2) 'redata: Grids in NC reduced to single precision'
           end if
        end if
        write(*,4) 'redata: NC compression level',cdf_comp_lvl
        if ( cdf_r_parallel ) then
           write(*,2) 'redata: Reads NC in parallel'
        end if
        if ( cdf_w_parallel ) then
           write(*,2) 'redata: Writes NC in parallel (possibly not working)'
        end if
     end if
  end if
#endif

  if (ionode) then
     write(6,'(2a)') 'redata: ', repeat('*', 71)
  endif
  
  if (cml_p) then
     call cmlEndParameterList(mainXML)
  endif

  ! Print blocks
  call mixers_scf_print_block( )

  RETURN
  !-------------------------------------------------------------------- END
1  format(a,t53,'= ',2x,l1)
2  format(a)
3  format(a,t53,'= ',a)
4  format(a,t53,'= ',i8)
5  format(a,t53,'= ',i8,a)
6  format(a,t53,'= ',f10.4,a)
7  format(a,t53,'= ',f12.6,a)
8  format(a,t53,'= ',f14.12)
9  format(a,t53,'= ',f10.4)
10 format(t55,a)
11 format(a,t53,'= ',f12.6)

CONTAINS
  subroutine deprecated( str )

    implicit none

    character(len=*), intent(in) :: str

    if (ionode) then
       if (fdf_defined(trim(str))) then
          write(6,'(a)') '**Warning: FDF symbol '//trim(str)//&
               ' is deprecated. See the manual.'
       endif
    endif

  end subroutine deprecated

end subroutine read_options


