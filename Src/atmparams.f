! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module atmparams

      implicit none 
!
!    Hard-wired parameters to be eliminated in the future
!

C INTEGER  NZETMX   : Maximum number of PAOs or polarization orbitals
C                     with the same angular  momentum and 
C                     for the same species.       

         integer, parameter, public  :: nzetmx =  200  

C INTEGER  NKBMX    : Maximum number of Kleinman-Bylander projectors
C                     for each angular momentum

         integer, parameter, public  :: nkbmx  =   3

C INTEGER  NSMX    : Maximum number of semicore shells for each angular
C                    momentum present in the atom ( for normal atom nsmx=0)

         integer, parameter, public  :: nsmx  =    2  
         integer, parameter, public  :: nsemx = 1 + nsmx  

C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                     orbitals,projectors and local neutral-atom 
C                     pseudopotential.

         integer, parameter, public  :: ntbmax =  500  

C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.

         integer, parameter, public  :: lmaxd  =    4  
         integer, parameter, public  :: lmx2   = (lmaxd+1)*(lmaxd+1)  

C INTEGER  NRMAX    : Maximum number of points in the functions read
C                     from file '.vps' or '.psf' (this number is
C                     determined by the parameter nrmax in the
C                     program atm, which generates the files with
C                     the pseudopotential information). The number
C                     of points in the grid can be redefined if the
C                     pseudopotential is reparametrized.
C                     nrmax = 20000 is a typical safe value in this case.
C                     

         integer, parameter, public  :: nrmax  = 20000

         integer, parameter, public  :: maxos=2*nzetmx*lmx2*nsemx

         private

      end module atmparams
