
      subroutine radgen_init(bUseLUT,bGenLUT)

      implicit none

       include "mc_set.inc"
       include "density.inc"
       include "cmpcom.inc"
       include "tailcom.inc"
       include "radgen.inc"
       include "radgenkeys.inc"

      logical bUseLUT,bGenLUT
      CHARACTER*256 tname
      COMMON/PYNUCL/INUMOD,CHANUM,ORDER
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER

! ... force block data modules to be read
      external radata

*     kill_elas_res =0, all rad corrections
*     kill_elas_res =1, ???
*     kill_elas_res =2, no rad corrections
      kill_elas_res=0

*     initialise the rad. correction part !
*     =  2 : fint cernlib extrapolation routine
*     =  1 : 2dim spline
*     =  0 : make lookuptable
*     = -1 : do not lookup table , calculate all events
      if( bUseLUT ) then
        if( bGenLUT ) then
          ixytest = 0
        else
          ixytest = 2
        endif
      else
        ixytest = -1
      endif
     
      if ((NINT(INUMOD).gt.1).and.(NINT(CHANUM).gt.1)) then
         mcSet_TarA = NINT(INUMOD)
         mcSet_TarZ = NINT(CHANUM)
      endif
      
      write(*,*) mcSet_TarA, mcSet_TarZ

*----------------------------
*     ire is target 1,2,3
*     add neutron target
      if(     mcSet_TarA .eq. 1 .and. mcSet_TarZ .eq. 0 ) then
        ire = 0
      elseif( mcSet_TarA .eq. 1 .and. mcSet_TarZ .eq. 1 ) then
        ire = 1
      elseif( mcSet_TarA .eq. 2 .and. mcSet_TarZ .eq. 1 ) then
        ire = 2
      elseif( mcSet_TarA .eq. 3 .and. mcSet_TarZ .eq. 2 ) then
        ire = 3
      elseif( mcSet_TarA .eq. 4 .and. mcSet_TarZ .eq. 2 ) then
        ire = 4
      elseif( mcSet_TarA .eq. 14 .and. mcSet_TarZ .eq. 7 ) then
        ire = 14
      elseif( mcSet_TarA .eq. 20 .and. mcSet_TarZ .eq. 10 ) then
        ire = 20
      elseif( mcSet_TarA .eq. 40 .and. mcSet_TarZ .eq. 20 ) then
        ire = 40
      elseif( mcSet_TarA .eq. 84 .and. mcSet_TarZ .eq. 36 ) then
        ire = 84
      elseif( mcSet_TarA .eq. 131 .and. mcSet_TarZ .eq. 54 ) then
        ire = 131
      elseif( mcSet_TarA .eq. 197 .and. mcSet_TarZ .eq. 79 ) then
        ire = 197
      elseif( mcSet_TarA .eq. 238 .and. mcSet_TarZ .gt. 92 ) then
        ire = 238
      else
        write(*,*)( 'RADGEN_INIT: invalid target selection' )
      endif

C... Using radgen with Pythia the target and beam polarisation can be 
C    hardcoded to zero
      plrun = 0.
      pnrun = 0.
      if(ire.lt.10) then
         tname='radgen/xytab0unp.dat'
      elseif (ire.lt.100) then
         tname='radgen/xytab00unp.dat'
      else
         tname='radgen/xytab000unp.dat'
      endif

      if (ire.lt.10) then
          write(tname(13:13),'(i1)')ire
      elseif (ire.lt.100) then
          write(tname(13:14),'(i2)')ire
      else 
          write(tname(13:15),'(i3)')ire
      endif

*----------------------------
* grid of important regions in theta (7*ntk)
      ntk = 35
*----------------------------
* photonic energy grid
      nrr = 100
*----------------------------
* min. energy in the calo (resolution parameter)
      demin=0.10

*----------------------------
      ap=2.*amp
      amp2=amp**2
      ap2=2.*amp**2
      if(kill_elas_res.eq.1) amc2=4d0

      if(ire.eq.0) then
        amt=.939565d0
        rtara=1d0
        rtarz=0d0
        fermom=0d0
      elseif(ire.eq.1)then
        amt=.938272d0
        rtara=1d0
        rtarz=1d0
        fermom=0d0
      elseif(ire.eq.2)then
        amt=1.87561d0
        rtara=2d0
        rtarz=1d0
        fermom=.07d0
      elseif(ire.eq.3)then
        amt=2.8161d0
        rtara=3d0
        rtarz=2d0
        fermom=.087d0
      elseif(ire.eq.4)then
        amt=3.75567d0
        rtara=4d0
        rtarz=2d0
        fermom=.087d0
      elseif(ire.eq.14)then
        amt=13.1447d0
        rtara=14d0
        rtarz=7d0
        fermom=.100d0
        call fordop
      elseif(ire.eq.20)then
        amt=18.7783d0
        rtara=20d0
        rtarz=10d0
        fermom=.100d0
        call fordop
      elseif(ire.eq.40)then
        amt=37.5567d0
        rtara=40d0
        rtarz=20d0
        fermom=.100d0
        call fordop
      elseif(ire.eq.84)then
        amt=78.8769d0
        rtara=84d0
        rtarz=36d0
        fermom=.105d0
        call fordop
      elseif(ire.eq.131)then
        amt=123.0132d0
        rtara=131d0
        rtarz=54d0
        fermom=.105d0
        call fordop
      elseif(ire.eq.197)then
        amt=184.9922d0
        rtara=197d0
        rtarz=79d0
        fermom=.105d0
        call fordop
      elseif(ire.eq.207)then
        amt=194.3839d0
        rtara=207d0
        rtarz=82d0
        fermom=.105d0
        call fordop
      elseif(ire.eq.238)then
        amt=223.4975d0
        rtara=238d0
        rtarz=92d0
        fermom=.105d0
        call fordop
      endif

*-----------------------------
      if (mcSet_XMin.ge. 1.e-05) then
        ntx=ntdis
        write(*,*) 
     &  ('Radiative corrections for DIS kinematics (xmin=1.e-05')
      elseif (mcSet_XMin.ge. 1.e-09) then
        ntx=ntpho
        write(*,*) 
     &  ('Radiative corrections for photoproduction (xmin=1.e-09)')
        write(*,*) 
     &  ('Elastic and quesielastric contributions set to ZERO !')
        if(ixytest.gt.0 .and. mcSet_XMin.lt.1.e-07) then
          write(*,*) 
     &  ('Xmin in lookup table = 1.e-07 --> extrapolation to 1.e-09')
        endif
        if(ixytest.eq.0 .and. mcSet_XMin.lt.1.e-07) then
          write(*,*) 
     &  ('Xmin in lookup table = 1.e-07 --> change minimum x')
          stop
        endif
      else
        write(*,*) 
     &  ('Minimum x value below minimum value of 10**-9 ! ')
        stop
      endif

* initialize lookup table
      if(ixytest.ge.0) then

* number of x bins depending on x range desired
* finer grid at lower x values
* number of y bins = number of x bins (equidistant in y)
* (needed to avoid double loop - speed optimisation)
      write(*,*) ('*********************************************')
      write(*,*) 
     &  ('Make sure that you did create the correct lookup table')
      write(6,*) 'number of x bins in RC table = ',ntx
      write(6,*) 'size of x bins in RC table depending on x'
      nty=ntx
      write(6,*) 'number of y bins in RC table = ',nty
      write(6,*) 'size of y bins in RC table = ',
     &  (radgen_ymax-radgen_ymin)/(nty-1)
      write(*,*) ('*********************************************')

        if(ixytest.eq.0) then
          write(*,*)('RADGEN_INIT: Creating lookup table '//tname)
        elseif (ixytest.eq.2) then
           write(*,*)
     +       ('RADGEN_INIT: Loading lookup table '//tname)
        endif
        call xytabl(tname,mcSet_EneBeam,plrun,pnrun,ixytest,ire)
      endif

      end
