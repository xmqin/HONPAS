module subs

use precision

public :: ival, manual, manual_dm_creator, orbital, txt2wrd

private

CONTAINS     


      function ival(txt)
        character(len=*), intent(in) :: txt
        integer                      :: ival

        integer  :: tmp, iostat

        ival = -1
!!        read(txt,fmt=*,iostat=iostat) tmp
        read(txt,fmt="(i10)",iostat=iostat) tmp
        ! Sometimes the above read succeeds
        ! for non-numeric strings...
!!!!        if (tmp== 0 .or. iostat /= 0) then
        if (iostat /= 0) then
           ! Not an integer, return -1 flag
           return
        endif

        ival = tmp

      end function ival

      subroutine txt2wrd (txt, wrd, nw, pepe)
      implicit none

      character(len=*), intent(in)   :: txt
      character(len=20), intent(out) :: wrd(:)
      integer, intent(out) :: nw
      integer, intent(in) :: pepe

      character(len=len(txt)+1) taux
      integer  ::  i1, i2, nmaxwords

      nmaxwords = size(wrd)

      taux=txt // " "
      nw=0
      do while (len_trim(taux).ne.0)
        nw=nw+1
        if (nw.gt.nmaxwords) STOP "overflow in txt2wrd"
        i1=verify(taux,' ')
        i2=index(taux(i1+1:),' ')+i1
        wrd(nw)=taux(i1:i2-1)
        taux(i1:i2-1)=repeat(' ',i2-i1)
      enddo
      end subroutine txt2wrd

      subroutine orbital (txt, ia, c, n, l, k)
        implicit none
      character(len=*), intent(in) :: txt
      character(len=*), intent(out) :: c
      integer, intent(out) ::  ia, n, l, k

      character(len=len(txt)) ::  atx, ntx
      integer :: i_, i0, lng, il

      atx=txt
      ntx=txt
      lng=len_trim(txt)
!
!     Error flag set
!
      ia=-1
!
!     Defaults meaning "not specified"
!
      c='  '
      n=-1
      l=-1
      k=-1

      i_=index(txt,'_')
      if (i_.eq.1) return  ! Error

      i0=i_-1
      if (i_.eq.0) i0=lng
      ntx=txt(1:i0)//repeat(' ',20-i0)
      ia=ival(ntx)
      if (ia.eq.-1) then   ! Label specified
!!!AG        if (i0.gt.2) return  ! No more stuff
        ia=0
        c=txt(1:i0)
      endif

      if (i_.eq.0) return

      ! Now search for shell info

      ! Angular momentum
      atx=txt(i_+1:20)//repeat(' ',i_)
      lng=len_trim(atx)
      il=scan(atx,'spdfgh')
      if ((il.ne.0).and.(il.ne.1).and.(scan(atx,'spdfgh',back=.true.).eq.il)) then
        l=index('spdfgh',atx(il:il))-1
      else
        ia=-1   ! error flag
        return
      endif

      ! n quantum number
      ntx=atx(1:il-1)//repeat(' ',20-il+1)
      n=ival(ntx)
      if (n.eq.-1) then  !! error
        ia=-1
        return
      endif

      ! m quantum number (actually index)
      if (il.eq.lng) then
        k=-1
      else
        k=0
        ntx=atx(il+1:20)//repeat(' ',il)
        if (l.eq.1) then
          if (trim(ntx).eq.'y') k=1
          if (trim(ntx).eq.'z') k=2
          if (trim(ntx).eq.'x') k=3
        endif
        if (l.eq.2) then
          if (trim(ntx).eq.'xy') k=1
          if (trim(ntx).eq.'yz') k=2
          if (trim(ntx).eq.'z2') k=3
          if (trim(ntx).eq.'xz') k=4
          if (trim(ntx).eq.'x2-y2') k=5
        endif
        if (k.eq.0) then
          k=ival(ntx)
          if ((k.eq.-1).or.(k.le.0).or.(k.gt.2*l+1)) then
             ! Signal error
            ia=-1
            return
          endif
        endif
      endif

      endsubroutine orbital

      subroutine manual

      write(6,"('* MPROP PROGRAM')")
      write(6,"('  Miquel Llunell, Universitat de Barcelona, 2005')")
      write(6,"('  Alberto Garcia, ICMAB-CSIC, 2007- ')")
      write(6,*)
      write(6,"('    MPROP calculates both DOS projections and COOP curves')")
      write(6,"('    using output files obtained with SIESTA. The atomic orbital (AO)')")
      write(6,"('    sets are defined in an input file (MLabel.mpr).')")
      write(6,"('  ')")
      write(6,*) "Usage: mprop [ options ] MPROP_FILE_BASENAME"
      write(6,*) "Options:"
      write(6,*) "           -h:  print manual                    "
      write(6,*) "           -d:  debug                    "
      write(6,*) "           -l:  print summary of energy information         "
      write(6,*) "   -s SMEAR  :  set value of smearing parameter (default 0.5 eV)"
      write(6,*) "    "
      write(6,*) "   Selection of eigenstates to be used: by eigenvalue range or band index:  "
      write(6,*) "    "
      write(6,*) "   -m Min_e  :  set lower bound of eigenvalue range                    "
      write(6,*) "   -M Max_e  :  set upper bound of eigenvalue range                    "
      write(6,*) "   -b Min_band  :  set minimum band index to be used               "
      write(6,*) "   -B Max_band  :  set maximum band index to be used               "
      write(6,*) "    "
      write(6,*) "   Plotting window and sampling rate. By default, window is based on eigenvalue range"
      write(6,*) "    "
      write(6,*) "   -n NPTS   :  set number of sampling points (default 200)"
      write(6,*) "   -w Ewindow_low  :  set lower energy bound of plotting window                    "
      write(6,*) "   -W Ewindow_high :  set upper energy bound of plotting window                    "
      write(6,*)
      write(6,"('* .mpr FILE STRUCTURE')")
      write(6,"('         SLabel                   # Name of the siesta output files')")
      write(6,"('         DOS/COOP                # Define the curve type to be calculated')")
      write(6,"('    /-[ If DOS selected; as many blocks as projections wanted ]')")
      write(6,"('    |    projection_name         # DOS projection name')")
      write(6,"('    \-   Subset of AO (*)        # Subset of orbitals included')")
      write(6,"('    /-[ If COOP selected; as many blocks as projections wanted ]')")
      write(6,"('    |    curve_name              # COOP curve name')")
      write(6,"('    |    Subset I of AO (*)      # Reference atoms or orbitals')")
      write(6,"('    |    d1 d2                   # Distance range in Angstrom')")
      write(6,"('    \-   Subset II of AO (*)     # Neighbour atoms or orbitals')")
      write(6,"('     (*) See below how to define subsets of AO')")
      write(6,"('     A final line with leading chars  ----  can signal the end of the input')")
      write(6,*)
      write(6,"('* INPUT FILES')")
      write(6,"('    [output files from SIESTA >=  2.4.1]')")
      write(6,"('    SLabel.WFSX and SLabel.HSX (new format)')")
      write(6,*)
      write(6,"('* OUTPUT FORMAT')")
      write(6,*) 
      write(6,*) " SLabel.alldos  :  full-range approximate DOS curve"
      write(6,*) " SLabel.ados    :  specified-range approximate DOS curve"
      write(6,*) " SLabel.intdos  :  full-range integrated-DOS curve"
      write(6,*) " MLabel.CurveName.pdos    :  PDOS curves"
      write(6,*) " MLabel.CurveName.coop    :  COOP curves"
      write(6,*) " MLabel.CurveName.cohp    :  COHP curves"
      write(6,"('    [A control .stt file will always be generated]')")
      write(6,*)
      write(6,"('* PROJECTION AND CURVES NAMES')")
      write(6,"('    Alphanumerical string up to 30 char. with no spaces')")
      write(6,"('* SUBSET OF AO USING ORDER NUMBERS')")
      write(6,"('    List of integer numbers preceeded by a + symbol')")
      write(6,"('    Each number refers to one AO in the final list of AO of SIESTA')")
      write(6,"('    Example: + 23 65 78')")
      write(6,"('* SUBSET OF AO USING ATOM_SHELL NOTATION')")
      write(6,"('    List of atoms and shell groups of AO')")
      write(6,"('    General notation: ATOM_SHELL')")
      write(6,"('     > ATOM:  Atomic symbol refers to all the atoms of that type')")
      write(6,"('              Integer number refers to the N-th atom in unit cell')")
      write(6,"('     > SHELL: Integer1+Letter+Integer2')")
      write(6,"('               > Integer1 refers to the n quantum number')")
      write(6,"('               > Letter   refers to the l quantum number (s,p,d,f,g,h)')")
      write(6,"('               > Integer2 refers to a single AO into the n-l shell')")
      write(6,"('                   Alternatively, alphanumerical strings can be used')")
      write(6,"('                     p-shells   1  y    d-shells   1  xy   4  xz')")
      write(6,"('                                2  z               2  yz   5  x2-y2')")
      write(6,"('                                3  x               3  z2')")
      write(6,"('    Particular cases:')")
      write(6,"('     > Just ATOM is indicated: all the AO of the atom will be included')")
      write(6,"('     > No value for Integer2:  all the AO of the shell will be included')")
      write(6,"('    Example: Ca_3p Al 4_4d3 5 O_2py')")
      stop

      end subroutine manual

      subroutine manual_dm_creator

      write(6,"('* DM_CREATOR PROGRAM')")
      write(6,"('  Alberto Garcia, ICMAB-CSIC, 2009')")
      write(6,*)
      write(6,"('    DM_CREATOR calculates a partial DM from a given energy interval,')")
      write(6,"('    using output files obtained with SIESTA.')")
      write(6,"('  ')")
      write(6,*) "Usage: dm_creator [ options ] SIESTA_SYSTEM_LABEL"
      write(6,*) "Options:"
      write(6,*) "           -h:  print manual                    "
      write(6,*) "           -d:  debug                    "
      write(6,*) "           -l:  print summary of energy information         "
      write(6,*) "   -s SMEAR  :  set value of smearing parameter (default 0.5 eV)"
      write(6,*) "   -m Min_e  :  set lower bound of energy range                    "
      write(6,*) "   -M Max_e  :  set upper bound of energy range                    "
      write(6,*)
      write(6,"('* INPUT FILES')")
      write(6,"('    [output files from SIESTA >=  2.4.1]')")
      write(6,"('    SLabel.WFSX and SLabel.HSX (new format)')")
      write(6,*)
      write(6,"('* OUTPUT FORMAT')")
      write(6,*) 
      write(6,*) " SLabel.alldos  :  full-range approximate DOS curve"
      write(6,*) " SLabel.intdos  :  full-range integrated-DOS curve"
      write(6,*) " DMOUT    :  Partial DM"
      write(6,*) " DM.nc (optional)  :  Partial DM in netcdf form"
      write(6,"('    [A control .stt file will always be generated]')")
      write(6,*)
      stop

      end subroutine manual_dm_creator


end module subs
     
