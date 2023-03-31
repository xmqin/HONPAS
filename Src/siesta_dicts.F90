
module siesta_dicts

  use precision, only : dp
#ifdef SIESTA__FLOOK
  use dictionary
#endif

  implicit none

#ifdef SIESTA__FLOOK

  ! A dictionary for all the options
  type(dictionary_t) :: options

  ! A dictionary for all transiesta options
  type(dictionary_t) :: ts_options

  ! A dictionary for all variables
  type(dictionary_t) :: variables

  private :: dict_variable_add_v_0d
  private :: dict_variable_add_a_1d
  private :: dict_variable_add_b_0d
  private :: dict_variable_add_i_0d
  private :: dict_variable_add_i_1d
  private :: dict_variable_add_d_0d
  private :: dict_variable_add_d_1d
  private :: dict_variable_add_d_2d
  interface dict_variable_add
     module procedure dict_variable_add_v_0d
     module procedure dict_variable_add_a_1d
     module procedure dict_variable_add_b_0d
     module procedure dict_variable_add_i_0d
     module procedure dict_variable_add_i_1d
     module procedure dict_variable_add_d_0d
     module procedure dict_variable_add_d_1d, dict_variable_add_d_2d
  end interface dict_variable_add

contains

  subroutine dict_clean()
    call delete(options,dealloc=.false.)
    call delete(ts_options,dealloc=.false.)
    call delete(variables,dealloc=.false.)
  end subroutine dict_clean

  subroutine dict_populate()
    call dict_populate_options()
    call dict_populate_variables()
  end subroutine dict_populate

  subroutine dict_populate_options()
    use files, only: slabel
    use siesta_options
    use m_mixing_scf, only: scf_mixs
    use m_steps, only: inicoor, fincoor

    integer :: slabel_len

    ! We simply re-create the options, (note the 
    ! de-allocation by "nullification")
    call delete(options, dealloc=.false.)

    slabel_len = len_trim(slabel)
    options = &
         ('Label'.kv.slabel(1:slabel_len))

    ! unluckily the dictionary does not
    ! implement a stringent way of doing characters
    ! by pointers (as we do not know their initial length).

    options = options // &
         ('DM.HistoryDepth'.kvp.DM_history_depth)
    
    ! Output options
    options = options // &
         ('Write.DenChar'.kvp.dumpcharge)
    options = options // &
         ('Write.MullikenPop'.kvp.mullipop)
    options = options // &
         ('Write.HirshfeldPop'.kvp.hirshpop)
    options = options // &
         ('Write.VoronoiPop'.kvp.voropop)

    ! SCF options
    options = options // &
         ('SCF.MinIterations'.kvp.min_nscf)
    options = options // &
         ('SCF.MaxIterations'.kvp.nscf)
    options = options // &
         ('SCF.MixHamiltonian'.kvp.mixH)
    options = options // &
         ('SCF.MixCharge'.kvp.mix_charge)

    ! Old-style options, not really used anymore
    options = options // &
         ('SCF.NumberPulay'.kvp.maxsav)
    options = options // &
         ('SCF.NumberBroyden'.kvp.broyden_maxit)
    options = options // &
         ('SCF.MixingWeight'.kvp.wmix)
    options = options // &
         ('SCF.NumberKick'.kvp.nkick)
    options = options // &
         ('SCF.KickMixingWeight'.kvp.wmixkick)

    ! New mixing options
    ! The first mixer is the mixing 1.
    options = options // &
         ('SCF.Mixer.Weight'.kvp.scf_mixs(1)%w)
    options = options // &
         ('SCF.Mixer.Restart'.kvp.scf_mixs(1)%restart)
    options = options // &
         ('SCF.Mixer.Iterations'.kvp.scf_mixs(1)%n_itt)

    options = options // &
         ('SCF.MonitorForces'.kvp.monitor_forces_in_scf)

    options = options // &
         ('ElectronicTemperature'.kvp.temp)

    ! Convergence criteria
    options = options // &
         ('SCF.Harris.Converge'.kvp.converge_Eharr)
    options = options // &
         ('SCF.Harris.Tolerance'.kvp.tolerance_Eharr)
    options = options // &
         ('SCF.DM.Converge'.kvp.converge_DM)
    options = options // &
         ('SCF.DM.Tolerance'.kvp.dDtol)
    options = options // &
         ('SCF.EDM.Converge'.kvp.converge_EDM)
    options = options // &
         ('SCF.EDM.Tolerance'.kvp.tolerance_EDM)
    options = options // &
         ('SCF.H.Converge'.kvp.converge_H)
    options = options // &
         ('SCF.H.Tolerance'.kvp.dHtol)
    options = options // &
         ('SCF.FreeE.Converge'.kvp.converge_FreeE)
    options = options // &
         ('SCF.FreeE.Tolerance'.kvp.tolerance_FreeE)
    
    options = options // &
         ('MD.MaxDispl'.kvp.dxmax)
    options = options // &
         ('MD.MaxForceTol'.kvp.ftol)
    options = options // &
         ('MD.MaxStressTol'.kvp.strtol)
    options = options // &
         ('MD.FinalTimeStep'.kvp.ifinal)
    options = options // &
         ('MD.FC.Displ'.kvp.dx)
    options = options // &
         ('MD.FC.First'.kvp.ia1)
    options = options // &
         ('MD.FC.Last'.kvp.ia2)
    options = options // &
         ('MD.Temperature.Target'.kvp.tt)
    options = options // &
         ('MD.Relax.CellOnly'.kvp.RelaxCellOnly)
    options = options // &
         ('MD.Relax.Cell'.kvp.varcel)
    options = options // &
         ('MD.Steps.First'.kvp.inicoor)
    options = options // &
         ('MD.Steps.Last'.kvp.fincoor)
    options = options // &
         ('MD.DM.History.Depth'.kvp.DM_history_depth)


    ! All write options    ! fdf-flag
    options = options // & ! SaveHS
         ('Write.HS'.kvp.saveHS)
    options = options // & ! Write.DM
         ('Write.DM'.kvp.writeDM)
    options = options // & ! Write.DM.End.Of.Cycle
         ('Write.EndOfCycle.DM'.kvp.write_DM_at_end_of_cycle)
    options = options // & ! Write.H
         ('Write.H'.kvp.writeH)
    options = options // & ! Write.H.End.Of.Cycle
         ('Write.EndOfCycle.H'.kvp.write_H_at_end_of_cycle)
    options = options // & ! Write.H
         ('Write.Forces'.kvp.writeF)
    options = options // & ! DM.UseSaveDM
         ('Use.DM'.kvp.UseSaveDM)

    options = options // & ! WriteHirshfeldPop
         ('Write.Hirshfeld'.kvp.hirshpop)
    options = options // & ! WriteVoronoiPop
         ('Write.Voronoi'.kvp.voropop)

    ! Options related to the mesh!
    options = options // & ! Required minimum meshcutoff
         ('Mesh.Cutoff.Minimum'.kvp.g2cut)
    options = options // & ! SaveRho
         ('Mesh.Write.Rho'.kvp.saverho)
    options = options // & ! SaveDeltaRho
         ('Mesh.Write.DeltaRho'.kvp.savedrho)
    options = options // & ! SaveRhoXC
         ('Mesh.Write.RhoXC'.kvp.saverhoxc)
    options = options // & ! SaveElectrostaticPotential
         ('Mesh.Write.HartreePotential'.kvp.savevh)
    options = options // & ! SaveNeutralAtomPotential
         ('Mesh.Write.NeutralAtomPotential'.kvp.savevna)
    options = options // & ! SaveTotalPotential
         ('Mesh.Write.TotalPotential'.kvp.savevt)
    options = options // & ! SaveIonicCharge
         ('Mesh.Write.IonicRho'.kvp.savepsch)
    options = options // & ! SaveBaderCharge
         ('Mesh.Write.BaderRho'.kvp.savebader)
    options = options // & ! SaveTotalCharge
         ('Mesh.Write.TotalRho'.kvp.savetoch)

  end subroutine dict_populate_options

  subroutine dict_populate_variables()

    use siesta_geom
    use kpoint_grid, only: kscell, kdispl
    use m_forces
    use m_energies
    use atomlist
    use m_stress

    real(dp), pointer :: r2(:,:)

    ! We simply re-create the options, (note the 
    ! de-allocation by "nullification")
    call delete(variables, dealloc=.false.)

    ! Add geometries (lets do with this for now)
    variables = &
         ('geom.na_u'.kvp.na_u)
    variables = variables // &
         ('geom.cell'.kvp.ucell)
    variables = variables // &
         ('geom.cell_last'.kvp.ucell_last)
    variables = variables // &
         ('geom.vcell'.kvp.vcell)
    variables = variables // &
         ('geom.nsc'.kvp.nsc)
    r2 => xa(:,1:na_u)
    variables = variables // &
         ('geom.xa'.kvp.r2)
    r2 => xa_last(:,1:na_u)
    variables = variables // &
         ('geom.xa_last'.kvp.r2)
    variables = variables // &
         ('geom.va'.kvp.va)

    ! Additional information regarding the
    ! atomic species
    variables = variables // &
         ('geom.species'.kvp.isa(1:na_u))
    variables = variables // &
         ('geom.z'.kvp.iza(1:na_u))
    variables = variables // &
         ('geom.last_orbital'.kvp.lasto(1:na_u))
    variables = variables // &
         ('geom.mass'.kvp.amass)
    variables = variables // &
         ('geom.neutral_charge'.kvp.qa(1:na_u))
    variables = variables // &
         ('geom.orbital_charge'.kvp.Datm(1:no_u))

    ! This is an abstraction made
    ! easy for the user.
    ! The forces that are used internally
    ! are actually the constrained ones.
    ! Hence any constrainst imposed will
    ! be visible to the user via this handle.
    variables = variables // &
         ('geom.fa'.kvp.cfa)
    ! This will let the user interact more
    ! freely by retrieval of all instances of the forces.
    variables = variables // &
         ('geom.fa_pristine'.kvp.fa)
    variables = variables // &
         ('geom.fa_constrained'.kvp.cfa)

    ! Add the stress components to the geometry
    variables = variables // &
         ('geom.stress'.kvp.cstress)
    variables = variables // &
         ('geom.stress_pristine'.kvp.stress)
    variables = variables // &
         ('geom.stress_constrained'.kvp.cstress)


    ! Add energies
    variables = variables // &
         ('E.neutral_atom'.kvp.DEna)
    variables = variables // &
         ('E.electrostatic'.kvp.DUscf)
    variables = variables // &
         ('E.fermi'.kvp.Ef)
    variables = variables // &
         ('E.harris'.kvp.Eharrs)
    variables = variables // &
         ('E.kinetic'.kvp.Ekin)
    variables = variables // &
         ('E.total'.kvp.Etot)
    variables = variables // &
         ('E.exchange_correlation'.kvp.Exc)
    variables = variables // &
         ('E.free'.kvp.FreeE)
    variables = variables // &
         ('E.ions_kinetic'.kvp.Ekinion)
    variables = variables // &
         ('E.ions'.kvp.Eions)
    variables = variables // &
         ('E.band_structure'.kvp.Ebs)
    variables = variables // &
         ('E.spin_orbit'.kvp.Eso)
    variables = variables // &
         ('E.dftu'.kvp.Edftu)
    
    variables = variables // &
         ('E.negf.dN'.kvp.NEGF_DE)
    variables = variables // &
         ('E.negf.harris'.kvp.NEGF_Eharrs)
    variables = variables // &
         ('E.negf.total'.kvp.NEGF_Etot)
    variables = variables // &
         ('E.negf.kinetic'.kvp.NEGF_Ekin)
    variables = variables // &
         ('E.negf.band_structure'.kvp.NEGF_Ebs)

    ! Add the number of charges to the system
    variables = variables // &
         ('charge.electrons'.kvp.qtot)
    variables = variables // &
         ('charge.protons'.kvp.zvaltot)

    ! Add the k-point sampling
    variables = variables // &
         ('BZ.k.Matrix'.kvp.kscell)
    variables = variables // &
         ('BZ.k.Displacement'.kvp.kdispl)

  end subroutine dict_populate_variables

  subroutine dict_repopulate_MD()

    use siesta_geom, only: na_u
    use siesta_geom, only: xa, xa_last, va, isa
    use atomlist, only: no_u, iza, lasto, qa, Datm

    real(dp), pointer :: r1(:), r2(:,:)
    integer, pointer :: i1(:)

    r2 => xa(:,1:na_u)
    call dict_variable_add('geom.xa', r2)
    r2 => xa_last(:,1:na_u)
    call dict_variable_add('geom.xa_last', r2)

    i1 => isa(1:na_u)
    call dict_variable_add('geom.species', i1)
    i1 => iza(1:na_u)
    call dict_variable_add('geom.z', i1)
    i1 => lasto(1:na_u)
    call dict_variable_add('geom.last_orbital', i1)
    r1 => qa(1:na_u)
    call dict_variable_add('geom.neutral_charge', r1)
    r1 => Datm(1:no_u)
    call dict_variable_add('geom.orbital_charge', r1)
    
  end subroutine dict_repopulate_MD

  subroutine dict_repopulate_sparse_matrices()

    use class_dSpData1D, only: val
    use sparse_matrices, only: maxnh
    use sparse_matrices, only: listh, listhptr, numh
    use sparse_matrices, only: Dscf, Escf, H, S, xijo
    use sparse_matrices, only: H_vkb_1D, H_kin_1D

    real(dp), pointer :: r1(:)

    call dict_variable_add('sparse.n_col', numh)
    call dict_variable_add('sparse.list_ptr', listhptr)
    call dict_variable_add('sparse.list_col', listh)
    call dict_variable_add('sparse.nnzs', maxnh)

    call dict_variable_add('sparse.S', S)
    call dict_variable_add('sparse.H', H)
    call dict_variable_add('sparse.DM', Dscf)
    call dict_variable_add('sparse.EDM', Escf)
    call dict_variable_add('sparse.xij', xijo)

    ! Now add the specific matrices
    r1 => val(H_vkb_1D)
    call dict_variable_add('sparse.H_Vkb', r1)
    r1 => val(H_kin_1D)
    call dict_variable_add('sparse.H_kin', r1)
    
  end subroutine dict_repopulate_sparse_matrices

  subroutine dict_variable_add_v_0d(name,val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.true.)
    variables = variables // (name.kv.trim(val))
  end subroutine dict_variable_add_v_0d
  subroutine dict_variable_add_a_1d(name,val)
    character(len=*), intent(in) :: name
    character(len=1), intent(inout), target :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_a_1d
  subroutine dict_variable_add_b_0d(name,val)
    character(len=*), intent(in) :: name
    logical, intent(inout), target :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_b_0d
  subroutine dict_variable_add_i_0d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout), target :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_i_0d
  subroutine dict_variable_add_i_1d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout), target :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_i_1d
  subroutine dict_variable_add_d_0d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_0d
  subroutine dict_variable_add_d_1d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_1d
  subroutine dict_variable_add_d_2d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val(:,:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_2d

#endif

end module siesta_dicts
  
