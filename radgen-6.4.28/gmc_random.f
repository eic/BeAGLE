!-----------------------------------------------------------------
! The point of this set of routines is to replace all potentially
! used random number generators with functions and subroutines
! that utilize a common seed sequence. In this case:
!
!       the CERNLIB RANLUX series
!
! MC programmers should now always use:
!       rndmq to initialize or obtain status
!       rlu to get a single 0:1 random number
!       nra to get a vector of 0:1 random numbers
!       rannor to get 2 Gaussian random numbers
!
! Documentation on RANLUX can be found here:
!     http://wwwinfo.cern.ch/asdoc/shortwrupsdir/v115/top.html 
!-----------------------------------------------------------------
! Initialization and status retrieval routine for random number sequence
!
!       CHOPT = ' '  reset sequence NSEQ to the beginning (seeds 0,0)
!               'S'  set seeds for sequence NSEQ to given values
!               'G'  get the current seeds for the current sequence
!
!       Note1: If ISEQ.le.0, the current (last used) sequence is used.
!-----------------------------------------------------------------

       subroutine rndmq (nseed1, nseed2, nseq, chopt)

       implicit none

       integer LUX_LEVEL
       parameter (LUX_LEVEL=4)

       integer nseed1, nseed2, nseq
       integer iseed1, iseed2, iseq, ilux
       character*(*) chopt
       character*1 c1opt

C MDB 2018-06-17 Remove RNDM. It's redefined in RAPGAP.
! ... force redefined random number generators to be taken from here
C       external rndm, irndm, nran, rannor, ranf, rlu, ranums
       external irndm, nran, rannor, ranf, rlu, ranums

! Parse option string

       c1opt = chopt(1:1)
       if (c1opt.ne.' '.and.c1opt.ne.'S'.and.c1opt.ne.'G') then
         write(*,*)('RNDMQ got unrecognized option')
         stop
       endif

! Take care of the possibilities of resetting the generator

! ... initialize generator to the beginning (seeds 0,0) of the given sequence
        if (c1opt.eq.' ') then
          call rluxgo(LUX_LEVEL,nseq,0,0)

! ... set seeds to given values, after retrieving current sequence number
! ... (and luxury level, why not)
        elseif (c1opt.eq.'S') then
          call rluxat(ilux,iseq,iseed1,iseed2)
          call rluxgo(ilux,iseq,nseed1,nseed2)

! ... retrieve current seeds and hand them back
        elseif (c1opt.eq.'G') then
          call rluxat(ilux,iseq,nseed1,nseed2)
        endif

       return
       end

C MDB 2018-06-17 Remove RNDM. It's redefined in RAPGAP.
C!-----------------------------------------------------------------
C! Replace the obsolete CERNLIB RNDM functions
C
C       real function rndm (dummy)
C
C       implicit none
C
C       real dummy, r
C
C        call ranlux(r,1)
C
C       rndm = r
C
C       return
C       end
C       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       integer function irndm (dummy)

       implicit none

       real dummy, r
       integer i

       equivalence (r,i)

        call ranlux(r,1)
       irndm = i

       return
       end

!-----------------------------------------------------------------
! Replace the obsolete CERNLIB NRAN subroutine

       subroutine nran (r,n)

       implicit none

       integer n
       real r(n)

        call ranlux(r,n)

       return
       end
       
!-----------------------------------------------------------------
! Replace the obsolete CERNLIB RANNOR subroutine

       subroutine rannor (a,b)

       implicit none

       real  a, b, r(2)
       external nran

       call rnormx(r,2,nran)
       a = r(1)
       b = r(2)

       return

       end

!-----------------------------------------------------------------
! Replace the F77 RANF

       real function ranf (dummy)

       implicit none

       real dummy, r

        call ranlux(r,1)    

       ranf = r

       return
       end

!-----------------------------------------------------------------
! Replace the JETSET random number generator

       real function rlu(idummy)

       implicit none

       integer idummy
       real r

        call ranlux(r,1)

       rlu = r

       return
       end

!-----------------------------------------------------------------
! Replace the DIVONNE random number generator

       subroutine ranums (r,n)

       implicit none

       integer n
       real r(n)

        call ranlux(r,n)

       return
       end

