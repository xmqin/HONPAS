C...............................................................
C
      subroutine test_xv(ii1,nat)
C
C     reads from ii1 (XV file) and returns the number of atoms
C
      implicit none
      integer ii1,ii,jj,nat
      double precision dummy

      rewind ii1
      do ii=1,3
         read (ii1,'(3x,3f18.9)',end=801,err=801) (dummy,jj=1,3)
      enddo
      read (ii1,*,end=802,err=802) nat
      return

  801 write (6,*) ' End/Error reading XV for cell vector ',ii
      stop
  802 write (6,*) ' End/Error reading XV for number of atoms line'
      stop

      end
C
C...............................................................
C
      subroutine read_xv(ii1,nat,ityp,iz,cc_ang,mass,label,coor_ang)
C
C     reads again from ii1 (XV file) to the end,
C     returns cell vectors and coordinates 
C     transformed from Bohr (as in XV) into Angstroem (as for XCrysDen)
C
      implicit none
      integer ii1,nat,nat1,iat,ii,jj,ityp(nat),iz(nat)
      double precision tau(3,3),mass(nat),coor_ang(3,nat),amass(110),
     .                 b2ang,cc_bohr(3,3),cc_ang(3,3),coor_bohr(3)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      character*2 label(nat),alab(110)
C
      data amass(1:40) /  1.00,   4.00,   6.94,   9.01,  10.81,
     .             12.01,  14.01,  16.00,  19.00,  20.18,  ! Ne
     .             23.99,  24.31,  26.98,  28.09,  30.97,
     .             32.07,  35.45,  39.95,  39.10,  40.08,  ! Ca 
     .             44.96,  47.88,  50.94,  52.00,  54.94, 
     .             55.85,  58.93,  58.69,  63.35,  65.39,  ! Zn 
     .             69.72,  72.61,  74.92,  78.96,  79.80, 
     .             83.80,  85.47,  87.62,  88.91,  91.22 / ! Zr 
      data amass(41:110) /
     .             92.91,  95.94,  97.91, 101.07, 102.91,
     .            106.42, 107.87, 112.41, 114.82, 118.71,  ! Sn
     .            121.76, 127.60, 126.90, 131.29, 132.91,
     .            137.33, 138.91, 140.12, 140.91, 144.24,  ! Nd
     .            146.92, 150.36, 151.96, 157.35, 158.93,
     .            162.50, 164.93, 167.26, 168.93, 173.04,  ! Yb
     .            174.97, 178.49, 180.95, 183.84, 186.21,
     .            190.23, 192.22, 195.08, 196.97, 200.59,  ! Hg
     .            204.38, 207.20, 208.98, 208.98, 209.99,
     .            222.02, 223.02, 226.03, 227.03, 232.04,  ! Th
     .            231.04, 238.03, 237.05, 244.06, 243.06,
     .            247.07, 247.07, 251.08, 252.08, 257.10,  ! Fm
     .            258.10, 259.10, 262.11, 261.11, 262.11,
     .            263.12, 262.12, 265.13, 266.13, 271.00 / ! Ds
      data alab / ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', 
     .            'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     .            'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .            'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     .            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .            'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     .            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .            'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', 
     .            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', 
     .            'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     .            'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds' /
C
      rewind (ii1)
C --  read in translation vectors, convert into Ang:
      do ii=1,3
         read  (ii1,101)  (cc_bohr(jj,ii),jj=1,3)
      enddo
      cc_ang = cc_bohr*b2ang

  101 format(3x,3f18.9)
      read (ii1,*) nat1
      if (nat1.ne.nat) then
C       check if the same as returned by test_xv :
        write (6,*) ' Number of atoms from first and second ',
     .              ' reading of XV differ !!'
        stop
      endif
      do iat=1,nat
        read (ii1,103,end=801,err=801) ityp(iat),iz(iat),coor_bohr
  103 format(i3,i6,3f18.9)
        mass(iat)=amass(iz(iat))
        label(iat)=alab(iz(iat))
        do ii=1,3
          coor_ang(ii,iat) = coor_bohr(ii)*b2ang
        enddo
      enddo
      return

  801 write (6,*) ' End/Error reading XV for atom number ',iat
      stop

      end
C
C...............................................................
