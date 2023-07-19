

      subroutine pseudoXML( ray, npotd, npotu, zion, zratio )

      use FoX_wxml
      use FoX_common

      implicit none

      include 'param.h'
      include 'radial.h'
      include 'ion.h'
      include 'orbital.h'
      include 'charge.h'
      include 'pseudowave.h'

      type(xmlf_t) :: xf

      integer npotd, npotu
      double precision  :: zion, zratio

      character*4      :: polattrib, relattrib, coreattrib
      character*10     :: ray(6)
      character*30 xcfuntype, xcfunparam
      character*30 gridtype, gridunits, gridscale, gridstep, gridnpoint
      character*30 slpotunits, slpotformat, slpotndown, slpotnup
      character*30 slwfnunits, slwfnformat, slwfndown, slwfnup
      integer                        :: ivps, ip
      double precision, allocatable   :: chval(:)

!
!     These determine the format for ASCII files
!
      character(len=*), parameter ::  
     .         fmt_int = "(tr1,i2)" ,                  
     .         fmt_nam = "(tr1,a2,tr1,a2,tr1,a3,tr1,a4)", 
     .         fmt_met = "(tr1,6a10,/,tr1,a70)" ,       
     .         fmt_pot= "(tr1,2i3,i5,3g20.12)" ,      
     .         fmt_rad = "(4(g20.12))"        ,      
     .         fmt_txt = "(tr1,a)"


! Digest and dump the information about the exchange and correlation functional
      select case(icorr) 

        case('ca') 
          xcfuntype    = 'LDA'
          xcfunparam   = 'Ceperley-Alder'

        case('wi') 
          xcfuntype    = 'LDA'
          xcfunparam   = 'Wigner'

        case('hl') 
          xcfuntype    = 'LDA'
          xcfunparam   = 'Hedin-Lundqvist'

        case('gl') 
          xcfuntype    = 'LDA'
          xcfunparam   = 'Gunnarson-Lundqvist'

        case('bh') 
          xcfuntype    = 'LDA'
          xcfunparam   = 'von Barth-Hedin'

        case('pb') 
          xcfuntype    = 'GGA'
          xcfunparam   = 'Perdew-Burke-Ernzerhof'

        case('rp') 
          xcfuntype    = 'GGA'
          xcfunparam   = 'RPBE - Hammer et al'

        case('rv') 
          xcfuntype    = 'GGA'
          xcfunparam   = 'revPBE Zhang+Yang'

        case('bl') 
          xcfuntype    = 'GGA'
          xcfunparam   = 'Becke-Lee-Yang-Parr'

        case('wc') 
          xcfuntype    = 'GGA'
          xcfunparam   = 'Wu-Cohen'

        case('ps') 
          xcfuntype    = 'GGA'
          xcfunparam   = 'Perdew-Burke-Ernzerhof-solid'

      end select

! Digest and dump the information about the pseudopotential flavor
      select case(irel) 

        case('isp') 
          polattrib   = 'yes'
          relattrib   = 'no'

        case('rel') 
          polattrib   = 'yes'
          relattrib   = 'yes'

        case('nrl') 
          polattrib   = 'no'
          relattrib   = 'no'

      end select

! Digest and dump the information about the non-linear core corrections
      select case(nicore) 

        case('pcec') 
          coreattrib  = 'yes'

        case('fcec') 
          coreattrib  = 'yes'

        case('fche') 
          coreattrib  = 'yes'

        case('pche') 
          coreattrib  = 'yes'

        case default
          coreattrib  = 'no'

      end select

! Digest and dump the information about the grid
      gridtype    = 'log'
      gridunits   = 'bohr'
      gridscale   = str(a)
      gridstep    = str(b)
      gridnpoint  = str(nr-1)
      
! Digest and dump the information about the semilocal components
      slpotunits  = 'rydberg'
      slpotformat = 'r*V'
      slpotndown  = str(npotd)
      slpotnup    = str(npotu)

! Digest and dump the information about the pseudowave functions
      slwfnunits  = 'electrons/bohr^(-3/2)'
      slwfnformat = 'u_n,l (r) = 1/r R_n,l (r)'
      slwfndown  = str(npotd)
      slwfnup    = str(npotu)

! Allocate and define the valence charge density
      allocate(chval(1:nr))

      do ip = 2, nr
        chval(ip) = zratio * (cdd(ip)+cdu(ip))
      enddo
                                                                                

! ---------------------------------------------------------------------
                                                                                
      call xml_OpenFile("VPSXML",xf, pretty_print=.true.)

!     call xml_AddXMLDeclaration(xf,"UTF-8")

      call xml_NewElement(xf,"pseudo")

        call xml_NewElement(xf,"header")
          call my_add_attribute(xf,"symbol",nameat)
          call my_add_attribute(xf,"atomic-number",str(znuc))
          call my_add_attribute(xf,"zval",str(zion))
          call my_add_attribute(xf,"creator",ray(1))
          call my_add_attribute(xf,"date",ray(2))
          call my_add_attribute(xf,"flavor",ray(3)//ray(4))
          call my_add_attribute(xf,"relativistic",relattrib)
          call my_add_attribute(xf,"polarized",polattrib)
          call my_add_attribute(xf,"core-corrections",coreattrib)
          call my_add_attribute(xf,"xc-functional-type",xcfuntype)
          call my_add_attribute(xf,"xc-functional-parametrization",
     .                          xcfunparam)
        call xml_EndElement(xf,"header")

        call xml_NewElement(xf,"grid")
          call my_add_attribute(xf,"type",gridtype)
          call my_add_attribute(xf,"units",gridunits)
          call my_add_attribute(xf,"scale",gridscale)
          call my_add_attribute(xf,"step",gridstep)
          call my_add_attribute(xf,"npts",gridnpoint)
        call xml_EndElement(xf,"grid")

        call xml_NewElement(xf,"semilocal")
          call my_add_attribute(xf,"units",slpotunits)
          call my_add_attribute(xf,"format",slpotformat)
          call my_add_attribute(xf,"npots-down",slpotndown)
          call my_add_attribute(xf,"npots-up",slpotnup)
  
! Down pseudopotentials follows

      vpsd: do ivps = 1, lmax
           if (indd(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"vps")
             call my_add_attribute(xf,"principal-n",str(no(indd(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"cutoff",str(rc(ivps)))
             call my_add_attribute(xf,"occupation",str(zo(indd(ivps))))
             call my_add_attribute(xf,"spin","-1")

             call xml_NewElement(xf,"radfunc")
               call xml_NewElement(xf,"grid")
                 call my_add_attribute(xf,"type",gridtype)
                 call my_add_attribute(xf,"units",gridunits)
                 call my_add_attribute(xf,"scale",gridscale)
                 call my_add_attribute(xf,"step",gridstep)
                 call my_add_attribute(xf,"npts",gridnpoint)
               call xml_EndElement(xf,"grid")

               call xml_NewElement(xf,"data")
                 call xml_AddCharacters(xf,viod(ivps,2:nr))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"vps")
         enddo vpsd

! Up pseudopotentials follows

         vpsu: do ivps = 1, lmax
           if (indu(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"vps")
             call my_add_attribute(xf,"principal-n",str(no(indu(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"cutoff",str(rc(ivps)))
             call my_add_attribute(xf,"occupation",str(zo(indu(ivps))))
             call my_add_attribute(xf,"spin","+1")

             call xml_NewElement(xf,"radfunc")
               call xml_NewElement(xf,"grid")
                 call my_add_attribute(xf,"type",gridtype)
                 call my_add_attribute(xf,"units",gridunits)
                 call my_add_attribute(xf,"scale",gridscale)
                 call my_add_attribute(xf,"step",gridstep)
                 call my_add_attribute(xf,"npts",gridnpoint)
               call xml_EndElement(xf,"grid")

               call xml_NewElement(xf,"data")
                 call xml_AddCharacters(xf,viou(ivps,2:nr))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"vps")
        enddo vpsu
        call xml_EndElement(xf,"semilocal")

! Dump of the pseudowave functions
        call xml_NewElement(xf,"pseudowave-functions")
          call my_add_attribute(xf,"units",slwfnunits)
          call my_add_attribute(xf,"format",slwfnformat)
          call my_add_attribute(xf,"n-pseudowave-functions-down",
     .                          slwfndown)
          call my_add_attribute(xf,"n-pseudowave-functions-up",
     .                          slwfnup)
  
! Down pseudowave function follows

        pswfd: do ivps = 1, lmax
           if (indd(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"pswf")
             call my_add_attribute(xf,"principal-n",str(no(indd(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"spin","-1")

             call xml_NewElement(xf,"radfunc")
               call xml_NewElement(xf,"grid")
                 call my_add_attribute(xf,"type",gridtype)
                 call my_add_attribute(xf,"units",gridunits)
                 call my_add_attribute(xf,"scale",gridscale)
                 call my_add_attribute(xf,"step",gridstep)
                 call my_add_attribute(xf,"npts",gridnpoint)
               call xml_EndElement(xf,"grid")

               call xml_NewElement(xf,"data")
                 call xml_AddCharacters(xf,pswfnrd(ivps,2:nr))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"pswf")
        enddo pswfd

! Up pseudowavefunction follows

       pswfu: do ivps = 1, lmax
           if (indu(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"pswf")
             call my_add_attribute(xf,"principal-n",str(no(indu(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"spin","+1")

             call xml_NewElement(xf,"radfunc")
               call xml_NewElement(xf,"grid")
                 call my_add_attribute(xf,"type",gridtype)
                 call my_add_attribute(xf,"units",gridunits)
                 call my_add_attribute(xf,"scale",gridscale)
                 call my_add_attribute(xf,"step",gridstep)
                 call my_add_attribute(xf,"npts",gridnpoint)
               call xml_EndElement(xf,"grid")

               call xml_NewElement(xf,"data")
                 call xml_AddCharacters(xf,pswfnru(ivps,2:nr))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"pswf")
        enddo pswfu
        call xml_EndElement(xf,"pseudowave-functions")

        call xml_NewElement(xf,"valence-charge")
          call xml_NewElement(xf,"radfunc")
            call xml_NewElement(xf,"grid")
              call my_add_attribute(xf,"type",gridtype)
              call my_add_attribute(xf,"units",gridunits)
              call my_add_attribute(xf,"scale",gridscale)
              call my_add_attribute(xf,"step",gridstep)
              call my_add_attribute(xf,"npts",gridnpoint)
            call xml_EndElement(xf,"grid")

            call xml_NewElement(xf,"data")
              call xml_AddCharacters(xf,chval(2:nr))
            call xml_EndElement(xf,"data")
          call xml_EndElement(xf,"radfunc")
        call xml_EndElement(xf,"valence-charge")

        call xml_NewElement(xf,"pseudocore-charge")
          call xml_NewElement(xf,"radfunc")
            call xml_NewElement(xf,"grid")
              call my_add_attribute(xf,"type",gridtype)
              call my_add_attribute(xf,"units",gridunits)
              call my_add_attribute(xf,"scale",gridscale)
              call my_add_attribute(xf,"step",gridstep)
              call my_add_attribute(xf,"npts",gridnpoint)
            call xml_EndElement(xf,"grid")

            call xml_NewElement(xf,"data")
              call xml_AddCharacters(xf,cdc(2:nr))
            call xml_EndElement(xf,"data")
          call xml_EndElement(xf,"radfunc")
        call xml_EndElement(xf,"pseudocore-charge")

        call xml_EndElement(xf,"pseudo")
      call xml_Close(xf)

      deallocate(chval)

      CONTAINS

      subroutine my_add_attribute(xf,name,value)
      type(xmlf_t), intent(inout)   :: xf
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
      end subroutine my_add_attribute

      end subroutine pseudoXML

