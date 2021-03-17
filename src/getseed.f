      subroutine getseed(seed,verbosity)
C     17-MAR-2021 Initial Version M.D. Baker + B. Schmookler
C
C     Get a random seed that is safely different for multiple 
C     processes being run at nearly the same time.
C
C     If /dev/urandom exists, take tha absolute value of that.
C     Otherwise combine the time (to nearest second) and pid.
C
C     Out: seed (I) - typically 8-10 digits. Positive definite.
C     In: verbosity (Logical) - Print seed & method to screen if .TRUE.
C
      implicit none
      integer seed
      logical verbosity
      integer i, n, un, istat, dt(8), pid, t(2), s
      integer*8 count, tms
      integer*4 today(3), now(3)
      open(newunit=un, file="/dev/urandom", access="stream", 
     & form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
         seed = abs(seed)
         if (verbosity)
     &      write (*,*) "Successfully got seed from /dev/urandom: ",seed
      else
! Fallback to using the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
         call idate(today)      ! today(1)=day, (2)=month, (3)=year
         call itime(now)        ! now(1)=hour, (2)=minute, (3)=second
         pid = getpid()
         seed = today(1)+10*today(2)+today(3)+now(1)+5*now(3)+
     &     10001*pid
         if (verbosity)
     &         write (*,*) "Got seed from time and pid: ",seed
      end if
      return
      end
