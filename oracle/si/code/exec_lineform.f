      program exec_lineform
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LS MPA 2010
! main to call the lineform.f subroutine. Not conceived for "manual"
! use but to be called by line_form.pl perl script.

! units called
! fort.10 all the input quantities (see description below) - produced by ngcl.pl
! fort.14 output

! variables
! abudev(natom), real, varied abundances as passed to 
! abvar, varied abundance read from fort.10
! abvele, integer, atmic number of varied abundance, read from fort.10
! flux_out, real(NWMAX), flux output
! glog, real, stellar gravity
! head_*, character, variuos variables storing header values.
! LF_CTRL, integer, input control bitmask. Is read from fort.10, already in decimal integer form.
! lf_status_out, integer, output lineform bitmask.
! linefile, character*12, parameter set at "linedata.dat" is the local name for the line file. Linking to the
!                         actual file is handled by launching script.
! masw, misw, real, minimum,maximum stepwidth
! ewidth - equivalent width (NOT YET IN MAIN SIU)
! met, real, stellar logarythmic metallicity
! micro, real, stellar microturbulent velocity
! nhead_f, integer parameter, size of head_f_out
! nhead_i, integer parameter, size of head_i_out
! nwmax, integer, max number of wavelength/flux points in output
! nabvar, integer, number of abundance variations to be read from fort.10
! sw_crit, real, stepwidth criterium, read from unit 10
! teff, real, stellar effective temperature
! wlam_cnt_out, integer, number of pixels in output spectrum
! wlam_out, double (NWMAX), wavelengths for output
! wmax,wmin, real: 
!                - start/end wavelength for computation (after call to lineform.f)
!                - central wavelength plus ion ID (e.g. 26.01)(before call to lineform.f)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!20100922 written
!20101001 Maria changed to fixed format
!20120802 many more changes
!
      implicit none
      real(kind=4) :: teff,glog,met,micro,wmin,wmax,abvar, wcenter
      real(kind=4) :: abudev,misw,masw,sw_crit,head_f_out
      real(kind=4) :: ewidth, cos_theta,wstep
      real(kind=4) :: wminnew, wmaxnew, excpot
      real(kind=8) :: wlam_out, flux_out, nunknown
      integer      :: nabvar,i,j,abvele,natom,lf_ctrl,lf_status_out
      integer      :: nhead_i,nhead_f,head_i_out,numwint

      integer(kind=8) :: wlam_cnt_out,NWMAX
 
      character(len=20)  :: head_atmuser,head_linuser
      character(len=30)  :: head_atmdate,head_lindate
      character(len=100) :: head_atmcmt,head_lincmt
      character(len=255) :: head_linref,head_atmref
      character(len=70)  :: linefile,ewout
      parameter(natom=92,linefile="linedata.dat",nhead_i=9)
!      parameter(NWMAX=1000000,numwint=1000)
      parameter(NWMAX=1000000,numwint=1)
!
!     parameter(natom=92,NWMAX=70000,nhead_i=9)     
      parameter(nhead_f=194)
      dimension abudev(natom),head_i_out(nhead_i),head_f_out(nhead_f)
      dimension flux_out(NWMAX),wlam_out(NWMAX)
      dimension wminnew(NUMWINT),wmaxnew(NUMWINT)
      external lineform

! formats
  111 format("# teff, log g, met, micro")
  112 format("# ",f6.1,2x,f5.2,2x,f5.2,2x,f5.2)
  113 format("# wmin, wmax, min_sw, max_sw, sw_crit")
  114 format("# ",f15.4,2x,f10.4,2x,f6.4,2x,f6.4,2x,f5.2)
  115 format("# HEADER:")
  116 format("# ",A20)
  117 format("# ",A30)
  118 format("# ",A100)
  119 format("# ",A255)
  221 format("########################################################")
  222 format("# wave, flux")
  223 format(f15.4,4x,f9.6)                !for absolute fluxes -- e format
!  223 format(f15.4,4x,e12.4)
  224 format(A6,4x,f9.5)
  225 format(f9.2,A7,f6.2,A6,4x,i2,4x,f5.2,4x,f15.5)

!read in
      !linefile = '/afs/mpa/home/mbergema/siu/linedata/master_line.dat'
      !ewout ='/afs/mpa/home/mbergema/siu/ewgrid/results/nlte_ews.log'
      cos_theta = 1.0
      abudev=abudev*0.
      open(10,status='old')
      read(10,*)teff
      read(10,*)glog
      read(10,*)met
      read(10,*)micro
      read(10,*)wmin,wmax
      read(10,*)misw,masw
      read(10,*)lf_ctrl
      read(10,*)sw_crit
!from here on we read the abundance variations
      read(10,*)nabvar
      do i=1,nabvar
        read(10,*)abvele,abudev(abvele)
      enddo
      close(10)

      wcenter = wmin

      write(*,*) 'split wrange into ', numwint, ' sub-intervals'
      wstep = (wmax-wmin)/numwint
      wminnew(1) = wmin
      wmaxnew(1) = wmin+wstep
      do i=2,numwint
       wminnew(i) = wmaxnew(i-1)+0.001
       wmaxnew(i) = wmaxnew(i-1)+wstep
      enddo
      write(*,*) wmin, wmaxnew(numwint)

      open(14,status='unknown')
      write(14,111)
      write(14,112)teff,glog,met,micro
      write(14,113)
      write(14,114)wmin,wmax,misw,masw,sw_crit
      write(14,115)
!     write(14,116)head_atmuser   ! binary stuff
!      write(14,116)head_linuser
!     write(14,117)head_atmdate   ! binary stuff
!      write(14,117)head_lindate
!     write(14,118)head_atmcmt    ! binary
!      write(14,118)head_lincmt
!      write(14,119)head_atmref
!      write(14,119)head_linref
      write(14,221)
      write(14,222)

      call lineform(wmin,wmax,teff,glog,met,micro,cos_theta,abudev,   
     &  lf_ctrl,sw_crit,misw,masw,linefile,lf_status_out,wlam_cnt_out,  
     &  wlam_out,flux_out,ewidth,head_i_out,head_f_out, head_atmdate,
     &  head_atmuser,
     &  head_atmref,head_atmcmt,head_lindate,head_linuser,head_linref,   
     &  head_lincmt, excpot)
       do i=1,wlam_cnt_out
        write(14,223)wlam_out(i),flux_out(i)
       enddo
      close(14)

!!!      do i=1,numwint
!!!        call lineform(wminnew(i),wmaxnew(i),teff,glog,met,micro,
!!!     &    cos_theta,abudev,   
!!!     &    lf_ctrl,sw_crit,misw,masw,linefile,lf_status_out,wlam_cnt_out,  
!!!     &    wlam_out,flux_out,ewidth,head_i_out,head_f_out, head_atmdate,
!!!     &    head_atmuser,head_atmref,head_atmcmt,
!!!     &    head_lindate,head_linuser,head_linref,head_lincmt)
!!!
!!!!        write(14,*) i, wminnew(i),wmaxnew(i)
!!!        write(14,223)wlam_out(1),flux_out(1)
!!!        do j=2,wlam_cnt_out
!!!         if (wlam_out(j).gt.wlam_out(j-1)) then 
!!!          write(14,223)wlam_out(j),flux_out(j)
!!!         endif
!!!        enddo
!!!      enddo

       if (lf_ctrl.EQ.26247215.OR.lf_ctrl.EQ.30441519.OR.
     &     lf_ctrl.EQ.30507055.OR.lf_ctrl.EQ.30441517.OR.
     &     lf_ctrl.EQ.12615695.OR.lf_ctrl.EQ.17858607) then
         open(15,status='unknown')
         write(15,225) wcenter, ' Elow: ',
     &    excpot, ' NLTE: ', abvele,abudev(abvele),ewidth
         close(15)
       else
         open(16,status='unknown')
         write(16,225) wcenter, ' Elow: ', 
     &    excpot, '  LTE: ', abvele,abudev(abvele),ewidth
         close(16)
       endif

      stop
      end
!TBD
!- irradiance is not handled right now, only flux.
