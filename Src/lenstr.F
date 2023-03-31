! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
c----------------------------------------------------------------------
c State: stable          Date:1/1/97
c----------------------------------------------------------------------

        function lenstr( ss )

c-----------------------------------------------------------------------
c Objective: Return length of the string SS until the first white space
c-----------------------------------------------------------------------

        implicit none
c------------------------------------------------------- Input variables
        character*(*) ss

c------------------------------------------------------ Output variables
        integer lenstr

c------------------------------------------------------- Local variables
        integer ii, long

        integer LEN
        intrinsic LEN


c----------------------------------------------------------------- BEGIN
        long = LEN(ss)

        do ii= 1, long
            if (ss(ii:ii) .eq. ' ') then
                lenstr = ii-1
                return
            endif
        enddo

        lenstr = long
c------------------------------------------------------------------- END

        return
        end
