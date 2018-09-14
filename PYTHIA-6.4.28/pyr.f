!-----------------------------------------------------------------
! The point of this set of routines is to replace all potentially
! used random number generators with functions and subroutines
! that utilize a common seed sequence. In this case:
!
!		the CERNLIB RANECU series
!
! MC programmers should now always use:
!	rndmq	to initialize or obtain status
!	rlu	to get a single 0:1 random number
!	nran	to get a vector of 0:1 random numbers
!	rannor	to get 2 Gaussian random numbers
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Replace the PYTHIA random number generator

	double precision function pyr(idummy)

	implicit none

	integer	idummy
	real	r

	call ranlux(r,1)
	pyr = r

	return
	end

