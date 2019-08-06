      subroutine ApplyQW(qhat)
      implicit none

      include 'common.f'
      include 'bea_pyqm.inc'

      integer ip,iq,ir,iit ! For do
      integer iloop !current entry number
      integer iEg,iPtF,iqg

      double precision qhat,SupFac,ehat,iet

      double precision th,ph

      double precision inix,iniy,iniz
      double precision ipx,ipy,ipz
      double precision ipg,ipgx,ipgy,ipgz
      double precision tot,ipix,ipiy,ipiz
      double precision ipt,iptx,ipty,iptz
      double precision iE

      double precision iptot,ipl,ptg,plg, cr
      double precision cutoff,sca

      double precision ptot,pt,pt2,mmmm
      double precision N_const, w_n, w_high, w_mean
      double precision w_sum, w_hard, w_soft
      integer n_w, list, ij, i_w
      double precision w_gluon(50)
      integer ijoin(50)
      double precision pgx, pgy, pgz, phe, phi
      double precision phi_final, theta, theta_final, mag

      integer j

      alphas = 1d0/3d0
      iqw = 1
      scor = 1
      ncor = 0
      sfthrd = 1

      cutoff = pyq_iet !energy cut off
C      iEg = 0
C      iPtF = 3
      iEg = PYQ_IEG
      iPtF = PYQ_IPTF
      iet = PYQ_IET
      iqg = 1
      SupFac = 1
      ehat = 0
      cr=0

      iloop = N
      ip = 1
      ij=0

      do while (ip.le.iloop)

        if((K(ip,1).eq.2).or.((K(ip,1).eq.1).and.(ij.ge.1)).and.
     & ((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then


          mmmm = P(ip,5)/P(ip,4)
          call QWComput(qhat,P(ip,1),P(ip,2),P(ip,3),P(ip,4),mmmm,
     & K(ip,2))

          if(abs(K(ip,2)).le.5) then
            cr = 4d0/3d0
          else if(K(ip,2).eq.21) then
            cr=3d0
          endif

          N_const=(2d0*alphas*cr*sqrt(2d0*QW_wc))/(pi)

          if(((P(ip,4)-QW_w).lt.PYQ_IET).and.(QW_w.gt.0)) then
            QW_w=P(ip,4)-PYQ_IET
          endif

        endif




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             no gluons                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          if(iEg.eq.0) then


            w_gluon(1)=QW_w

            if((K(ip,1).eq.2).and.((abs(K(ip,2)).le.5).or.
     &  (K(ip,2).eq.21)).and.(P(ip,4).gt.PYQ_IET).and.(QW_w.gt.0)) then


              ij=ij+1

              call GluonEmission(ip,w_gluon(1),theta_final,phi_final)

              n=n+1
              ij=ij+1

              call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

              P(ip,1)=P(ip,1)-P(N,1)
              P(ip,2)=P(ip,2)-P(N,2)
              P(ip,3)=P(ip,3)-P(N,3)

              n=n-1


            else if((K(ip,1).eq.1).and.(ij.gt.1).and.
     & ((abs(K(ip,2)).le.5).and.(P(ip,4).gt.(PYQ_IET))).and.(QW_w.gt.
     & 0 )) then

              n=n+1
              ij=ij+1

              call GluonEmission(ip,w_gluon(1),theta_final,phi_final)

              call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

              P(ip,1)=P(ip,1)-P(N,1)
              P(ip,2)=P(ip,2)-P(N,2)
              P(ip,3)=P(ip,3)-P(N,3)

              n=n-1
              ij=ij+1
              ij=0

            else

              ij=0

            endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             1 gluon                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


          else if(iEg.eq.1) then


            w_gluon(1)=QW_w

            if(K(ip,1).eq.2) then

              ij = ij+1
              ijoin(ij)=ip

              if(((abs(K(ip,2)).le.5)
     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n=n+1

                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
                  ij=ij+1

                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)

                  ijoin(ij)=n

              endif

            else if(K(ip,1).eq.1) then


              if((abs(K(ip,2)).le.5)
     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then


                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n=n+1
                  ij=ij+1

                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)

                  ijoin(ij)=n

              endif

            ij=ij+1
            ijoin(ij)=ip

            if(ij.ge.3) then
                call PYJOIN(ij,ijoin)
            endif

            ij=0

          endif



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             1 hard gluons + soft                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          else if(iEg.eq.2) then

            if ((QW_w.gt.0).and.((K(ip,1).eq.1).or.(K(ip,1).eq.2))
     & .and.((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then

              w_hard = (sqrt(QW_w)-(1/(4d0*N_const)))**2

              if(w_hard.lt.PYQ_IET) then

                w_hard = QW_w-PYQ_IET

                if(w_hard.lt.PYQ_IET) then
                  QW_w = 0
                  w_hard = 0
                endif

              endif

              w_gluon(1)=w_hard
              w_soft=QW_w - w_hard

            endif

            if(K(ip,1).eq.2) then

              ij = ij+1
              ijoin(ij)=ip

              if(((abs(K(ip,2)).le.5)
     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n=n+1
                  ij=ij+1

                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

                  PYQREC(1)=PYQREC(1)
     & +QW_w*sin(theta_final)*cos(phi_final)-P(N,1)
                  PYQREC(2)=PYQREC(2)
     & +QW_w*sin(theta_final)*sin(theta_final)-P(N,2)
                  PYQREC(3)=PYQREC(3)+QW_w*cos(theta_final)-P(N,3)
                  PYQREC(4)=PYQREC(4)+w_soft

                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)

                  ijoin(ij)=n

              endif

            else if(K(ip,1).eq.1) then

              if((abs(K(ip,2)).le.5)
     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n=n+1
                  ij=ij+1
                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

                  PYQREC(1)=PYQREC(1)
     & +QW_w*sin(theta_final)*cos(phi_final)-P(N,1)
                  PYQREC(2)=PYQREC(2)
     & +QW_w*sin(theta_final)*sin(theta_final)-P(N,2)
                  PYQREC(3)=PYQREC(3)+QW_w*cos(theta_final)-P(N,3)
                  PYQREC(4)=PYQREC(4)+w_soft

                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)

                  ijoin(ij)=n

              endif

            ij=ij+1
            ijoin(ij)=ip

            if(ij.ge.3) then
                call PYJOIN(ij,ijoin)
            endif

            ij=0

          endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            n gluons + soft                                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          else if(iEg.eq.3) then

            if ((QW_w.gt.0).and.((K(ip,1).eq.1).or.(K(ip,1).eq.2))
     & .and.((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then


              n_w=0
              w_mean = (sqrt(QW_w)-(1/(4*N_const)))**2
              w_sum = w_mean
              w_high = QW_w

              do while(w_mean>iet)
                n_w = n_w + 1
c                print*, n_w, w_mean !, w_high 
                w_gluon(n_w)=w_mean
                w_high = w_high-w_mean
                w_mean = (sqrt(w_high)-(1/(4*N_const)))**2
                w_sum = w_sum + w_mean
              enddo

              w_soft = QW_w - w_sum

           endif

            if(K(ip,1).eq.2) then

              ij = ij+1
              ijoin(ij)=ip
              i_w = 1

              if(((abs(K(ip,2)).le.5)
     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

                  mag = QW_w

                  do while(i_w.le.n_w)

                    call GluonEmission(ip,w_gluon(i_w),theta_final,
     & phi_final)

                    n=n+1
                    ij=ij+1

                    call PY1ENT(N,21,w_gluon(i_w),theta_final,phi_final)

                    PYQREC(1)=PYQREC(1)
     & +mag*sin(theta_final)*cos(phi_final)-P(N,1)
                    PYQREC(2)=PYQREC(2)
     & +mag*sin(theta_final)*sin(theta_final)-P(N,2)
                    PYQREC(3)=PYQREC(3)+mag*cos(theta_final)-P(N,3)

                    P(ip,1)=P(ip,1)-P(N,1)
                    P(ip,2)=P(ip,2)-P(N,2)
                    P(ip,3)=P(ip,3)-P(N,3)

                    ijoin(ij)=n

                    mag=mag-w_gluon(i_w)
                    i_w = i_w+1

                  enddo

                  mag = 0
                  PYQREC(4)=PYQREC(4)+w_soft

              endif

            else if(K(ip,1).eq.1) then

              i_w=1

              if((abs(K(ip,2)).le.5)
     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

                  mag = QW_w

                  do while(i_w.le.n_w)

                    call GluonEmission(ip,w_gluon(i_w),theta_final,
     & phi_final)

                    n=n+1
                    ij=ij+1
                    call PY1ENT(N,21,w_gluon(i_w),theta_final,phi_final)

                    PYQREC(1)=PYQREC(1)
     & +mag*sin(theta_final)*cos(phi_final)-P(N,1)
                    PYQREC(2)=PYQREC(2)
     & +mag*sin(theta_final)*sin(theta_final)-P(N,2)
                    PYQREC(3)=PYQREC(3)+mag*cos(theta_final)-P(N,3)

                    P(ip,1)=P(ip,1)-P(N,1)
                    P(ip,2)=P(ip,2)-P(N,2)
                    P(ip,3)=P(ip,3)-P(N,3)

                    ijoin(ij)=n

                    mag=mag-w_gluon(i_w)

                    i_w = i_w+1
                  enddo

                  PYQREC(4)=PYQREC(4)+w_soft
                  mag = 0

              endif

            ij=ij+1
            ijoin(ij)=ip

            if(ij.ge.3) then
                call PYJOIN(ij,ijoin)
            endif

            ij=0

          endif

        endif

        ip=ip+1

      enddo

      end


      subroutine GluonEmission(ip,w_gluon,theta_final,phi_final)

        include 'common.f'
        include 'bea_pyqm.inc'


        double precision w_gluon, theta_final, phi_final
        double precision theta, phi, phe, pgx, pgy, pgz
        integer ip

c        print*, 'ip = ', ip
c        print*, 'w_gluon = ', w_gluon


        P(ip,4)=P(ip,4)-w_gluon

        theta = acos(P(ip,3)/(sqrt(P(ip,1)**2+P(ip,2)**2+P(ip,3)**2)))
        phi = atan2(P(ip,2),P(ip,1))

c        print*, 'thata parton = ', theta
c        print*, 'phi parton = ', phi

        phe = 4*asin(1.)*ranf(0)

c        print*, 'phe = ', phe
c        print*, 'QW_th = ', QW_th


        pgx = (cos(theta)*cos(phi)*sin(QW_th)*cos(phe)-sin(phi)*
     & sin(QW_th)*sin(phe)+sin(theta)*cos(phi)*cos(QW_th))
        pgy = (cos(theta)*sin(phi)*sin(QW_th)*cos(phe)+cos(phi)*
     & sin(QW_th)*sin(phe)+sin(theta)*sin(phi)*cos(QW_th))
        pgz = (-sin(theta)*sin(QW_th)*cos(phe)+cos(theta)*cos(QW_th))


        theta_final = acos((pgz)/(sqrt(pgx**2+pgy**2+pgz**2)))
        phi_final = atan2(pgy, pgx)

c        print*, 'theta_final = ', theta_final
c        print*, 'phi_final = ', phi_final

      end




      subroutine QWComput(qhat,ipx,ipy,ipz,E,mmmm,id)
      implicit none

      include 'common.f'
      include 'bea_pyqm.inc'

      double precision ipx,ipy,ipz,E !input energy momentum of the particle
      double precision partmass !mass of the particle (for conservation purpose)
      double precision radius !distance to the center of the nuclei
      double precision x,y,z !position of the parton
      double precision pp,px,py,pz !integral steps in the space
      double precision integral_step !integral step
      parameter(integral_step = 0.1)
      double precision d !distance already coverd
      integer ipart !0=gluon - otherwise=quark
      double precision cont(1000),disc,step_QW !Variables for energy loss proba
      double precision xx,yy !Variables for energy loss proba
      integer i,nb_step !nb of step for QW calculation
      double precision total !Total of QW for normalization purpose
      double precision randnum !random number to pick the QW
      integer id !id of the parton
      double precision ChiR !Chi sq R
      double precision qhateff
      double precision qhat,ehat

      double precision mmmm !mass/energy of the incoming parton
      integer irej !used for test

      QW_w = 0.
      QW_L = 0.
      QW_wc = 0.
      QW_R  = 0.
      d = 0.
      ehat = 0.
      qhateff = qhat + ehat

      cont=0d+0
      disc=0d+0

ccc Init for qweight
      if (id.eq.21) then
        ipart = 0
      else if(abs(id).lt.7) then
        ipart = 1
      else
        write(*,*) 'Unknown parton with id =',id
      endif
      irw = 0
      nb_step = 200

ccc normalize momentum
      pp = sqrt(ipx**2+ipy**2+ipz**2)
      px = ipx/pp * integral_step
      py = ipy/pp * integral_step
      pz = ipz/pp * integral_step

      x = x_inter
      y = y_inter
      z = z_inter

ccc integration to calculate wc and R
      radius = sqrt(x**2+y**2+z**2)
      do while (radius.lt.20)
        QW_wc = QW_wc + integral_step * d *
     &          density_table(INT(radius/step_size_dens))
        QW_R = QW_R + integral_step *
     &          density_table(INT(radius/step_size_dens))
        d = d + integral_step
        x = x + px
        y = y + py
        z = z + pz
        radius = sqrt(x**2+y**2+z**2)
      enddo

      QW_L = QW_wc / QW_R
      QW_wc = qhateff/density_table(1) * QW_wc
      QW_R = 2 * density_table(1) * QW_wc**2 / QW_R / qhateff

ccccc Convert the units fm -> GeV-1
c      QW_L = QW_L/.1973269
ccccc Convert the units GeV2.fm -> GeV
      QW_wc = QW_wc/.1973269
ccccc Convert the units GeV2.fm2 -> no unit
      QW_R = QW_R / .1973269**2

c      print*,'QW_wc=',QW_wc

ccccc Calculate the energy loss probability
      if(sfthrd.eq.1) step_QW = 2.5/nb_step
      if(sfthrd.eq.2) step_QW = 9.8/nb_step
      yy = E/QW_wc

      total = 0.
      do i=1,nb_step
        xx = step_QW * i
        call qweight(ipart,id,mmmm,dble(QW_R),xx,yy,cont(i),disc)
        total = total + cont(i)*step_QW
      enddo

      total = total + disc
      disc = disc/total

c      print*,'total=',total
c      print*,'disc=',disc
      do i=1,nb_step
        cont(i) = cont(i) / total
      enddo

ccccc Pick randomely a quenching in the table
      if(disc .lt. 1.) then
        randnum = ranf(0)
c        print*, 'total < 1'
        if(randnum.gt.disc) then
          total = disc
          i = 1
          do while (randnum.gt.total)
            total = total + cont(i)*step_QW
            i = i + 1
          enddo
          QW_w = i * step_QW * QW_wc 
        endif
      endif

ccccc Calculate the angle probability
c      if(QW_w .gt. 0) then
c        step_QW = 1./nb_step
c        yy = E/QW_wc
c        xx = QW_w/QW_wc 
c       
c        total = 0.
c        do i=1,nb_step
c          ChiR = (step_QW * i)**2 * QW_R
c          call qweight(ipart,ChiR,xx,yy,cont(i),disc)
cc Do not keep negative probabilities
c          if (cont(i).lt.0) cont(i) = 0
c          total = total + cont(i)*step_QW
c        enddo
c        irej=1
c        if(total.lt.1e-5) irej=0
c        if(irej.eq.0) print*,'total=',total,'QW_w/c=',QW_w,' ',QW_wc
c        do i=1,nb_step
c          cont(i) = cont(i) / total
c        enddo
c        randnum = ranf(0)
c        total = 0.
c        i = 1
c        do while (randnum.gt.total)
c          total = total + cont(i)*step_QW
c          i = i + 1
c        enddo
c        QW_chi = i * step_QW
c        QW_th = asin(QW_chi)
c        if(irej.eq.0) print*,' i=',i
c        if(irej.eq.0) print*,'QW_chi=',QW_chi,' QW_th=',QW_th
c      endif
c      if(isnan(QW_th)) QW_th = pi/2
 
ccccc Scattering angle set to 0
      if(QW_w .gt. 0) then
         QW_th=0.  ! collinear gluon radiation assumption
      endif

      end

************************************************************************
*                          QWEIGHT                                     *
*   Interface to Arleo and Salgado-Wiedemann quenching weights         *
*   by A. Accardi (2004-2007)                                          *
*                                                                      *
*   Includes:                                                          *
*   - qweight    Quenching weight routine                              *
*   - reweight   Normalization for reweighting procedure               *
*                                                                      *
*     * qweight.a   split off modFF.f and adapted to new SW routines   *
*       (08 dec 04)                                                    *
*     * qweight.b   includes reweighting routine                       *
*       (09 mar 07)                                                    *
*                                                                      *
*  Please send me any comments:  aaccardi@nt3.phys.columbia.edu        * 
*                                                                      *
************************************************************************

**************************************************************************
*     Quenching weights interface to Arleo and Salgado-Wiedemann routines
      subroutine qweight(ipart,id,mmmm,rrrr,xx,yy,cont,disc)
*     Programmer: A.Accardi
*     Date: Jul 04  (v1)
*           Dec 04  (v2)
*     Revision: 8 Dec 04
*               9 Mar 07 - cleaned some unused variables
*
*  A. Commentary
*
*     Returns the continuous and discrete part of the quenching weight
*     in various approximations depending on the flags "iqw", "scor" 
*     and "ncor" (i.e, "Quenching weight routine", "size" and "energy" 
*     corrections) in the common block "qw":
*
*      iqw  = 1 --> Salgado-Wiedemann [3]
*             2 --> Arleo [2]
*      scor = 0,1 --> no,yes 
*      ncor = 0,1 --> no,yes 
*
*     Size and energy corrections are available according to the 
*     following table:  
*
*              scor  ncor | decscription      | alphas     | Ref.
*      -------------------+-------------------+------------+-----
*       Arleo   0     0   | asymptotic        |  any       | [2]   
*               0     1   | 1/nu corrections  |            |
*                         |   (no finite size)|            |
*       SW      0     0   | asymptotic size   | 1/3 or 1/2 | [3]
*                         |   (R->infty limit)|            |
*               1     0   | finite size       |            |
*      -------------------+-------------------+------------+-----
*
*     The Salgado-Wiedemann quenching weights come in 2 approximations
*     according to "sfthrd" (unused for Arleo's parametrization):
*   
*      sfthrd = 1 --> multiple soft collisions approximation
*               2 --> single hard scattering (N=1 in opacity expansion)
*      
*     Finally, the value of the size parameter R is stored in "rrrr".
*     If no size correction is desired (scor=0) then R=10d7 is used 
*     instead  of the user supplied value:
*
*       rrrr = value of R  
*
*     NOTES:
*
*     Arleo's quenching weight parametrization was obtained with 
*     alpha_s=1/2. However, for compatibility with Salgado-Wiedemann's 
*     work, who use alpha_s=1/3, we rescale Arleo's weight as follows [1].
*                                                                      
*     Since the mean energy loss (hence the first moment of the 
*     distribution) is directly proportional to alpha_s and these 
*     distributions have to be normalized (zeroth moment), we can rescale 
*     the quenching weights to any alpha_s. We have:  
*
*        D(epsilon, alpha_s_2)                                            
*              = alpha_s_1 / alpha_s_2 * D(alpha_s_1*epsilon/alpha_s_2)
*                                                                      
*     The parametrization in [2] is for alpha=1/2. Then for a generic alpha
*     we have
*                                                                      
*        D(epsilon,alpha) = 0.5/alpha D_Arleo(0.5/alpha*epsilon)
*
*     The same scaling is NOT valid for the SW quenching weights. 
*     Two values of alpha_s = 0.3 and 0.5 are available. For input 
*     values different from these two the routines issues an error and stops. 
*                                                              
*     INPUT VARIABLES
*
*     ipart = (i)  radiating parton (0=gluon - otherwise=quark)
*     rrrr  = (dp) medium finite size parameter R defined in [3]
*                  (inactive when scor=0)
*     xx    = (dp) w/wc, where w  = radiated energy
*                                 wc = 0.5*qhat*L**2 
*                                      coherence gluon energy 
*                                      (see MEDIUM PROPOERTIES)
*     yy    = (dp) nu/wc = parent parton energy normalized to wc
*     cont  = (dp) continous part of the quenching weight
*     disc  = (dp) discrete part of the q.w. (delta function)
*
*     INPUT IN COMMON BLOCK /qw/
*
*     alphas   = (dp) strong coupling at soft scales
*     iqw      = (i)  choice of SW or Arleo's quenching weights
*     scor     = (i)  finite size corrections flag (0=no 1=yes)
*     ncor     = (i)  finite energy corrections flag (0=no 1=yes)
*     sfthard  = (i)  1=soft scatterings  2=hard scatterings
*
*     NOTES ON USAGE:
*
*     1) xx & yy are adimensional
*     2) The quenching weight is normalized as follows:
*        Integral[continous,{x,0,Infinity}] + discrete = 1
*     3) the Salgado-Wiedemann routine works in the range
* 
*          0 < xxxx < 2.59   (soft scatt.)
*          0 < xxxx < 9.87  (hard scatt.)    
*
*          1 < rrrr < 40000  
*
*         For R outside the range an extrapolation is provided.
*         For xxxx outside the range no extrapolation is provided
*         and a segmentation fault ensues.
*
*         - For xxxx > 2.59 (9.87) this interface returns cont=0d0
*           with no explicit warning.
*
*     REFERENCES:

*     [1] F.Arleo, private communication                         
*     [2] F.Arleo, JHEP 0211:044,2002                                
*     [3] C.Salgado and U.Wiedeman, Phys.Rev.Lett.89:092303,2002     
*         C.Salgado and U.Wiedeman, Phys.Rev.D68:014008,2003         
*
*  B. DECLARATIONS 
*
      implicit none
      include 'bea_pyqm.inc'

      integer ipart
      double precision rrrr,xx,yy,cont,disc
      double precision mmmm
      integer id

      double precision rrin,frac,kk

*     FUNCTIONS

      double precision dbarg, dbarq

*     COMMON BLOCKS

*     Quenching weight choice
*     ... iqw             = (i)  1=SW  2=Arleo
*     ... scor            = (i)  finite size corrections flag (0=no 1=yes)
*     ... ncor            = (i)  finite energy corrections flag (0=no 1=yes)
*     ... sfthrd          = (i)  1=soft scatterings  2=hard scatterings
*     ... irw             = (i)  0= no reweighting  1= reweighting
      double precision alphas
            integer iqw,scor,ncor,sfthrd,irw
      common/qw/alphas,iqw,scor,ncor,sfthrd,irw

*     INITIAL PARAMETRS AND DATA
      
*    ... min and max rrrr for S.W. routine
      double precision xxmultmax,xxlinmax
      parameter (xxmultmax=2.59d0,xxlinmax=9.87d0)

*    ... S.W. initialization flag
      logical firstmult,firstlin
      save firstmult,firstlin
      data firstmult,firstlin/.true.,.true./

*
*  C. ACTION
*
c     write(*,*) 'test'
c     write(*,*) ipart,rrrr,xx,yy,cont,disc
c     write(*,*) alphas,iqw,scor,ncor,sfthrd,irw

 
      if ((scor.eq.1).and.(ncor.eq.1)) then
         print*, 'ERROR (qweight): finite size & finite nrg'
     &        //' corrections not yet implemented'
         stop
      end if

*    *** ARLEO's asymptotic medium size
      if (iqw.eq.2) then
*       ... discrete part is identically zero
         disc = 0d0
*       ... with rescaling factor kk for alphas=/=1/2
         if (alphas.le.0d0) then 
            print*, 'ERROR (qweight): alphas < 0'
            stop
         else
            kk = 1d0/(2d0*alphas) 
         end if
*       ... continuous part 
         if (ipart.eq.0) then
            cont =  kk * dbarg(kk*xx,yy)
         else
            cont = kk * dbarq(kk*xx,yy)
         end if

*    *** Salgado-Wiedemann's quenching weights
      else if (iqw.eq.1) then
*       ... initialization
         if ((sfthrd.eq.1).and.(firstmult)) then
            if(pyq_hq.eq.0) then
              call initmult(alphas)
            else if(pyq_hq.eq.1) then
              call  initmassmult
            else if(pyq_hq.eq.2) then
              if(id.le.3) then
                call initmult(alphas)
              else if(id.ge.4) then
                call initmassmult
              end if
            end if
         firstmult = .false.
         end if
         if ((sfthrd.eq.2).and.(firstlin)) then
            call initlin(alphas)
            firstlin = .false.
         end if
*       ... if no size corrections, sets rrin very large
         if (scor.eq.0) then
            rrin = 1d8
         else
            rrin = rrrr
         end if
*       ... calls Salgado-Wiedemann routine
         if (sfthrd.eq.1) then
            if (xx.le.xxmultmax) then
              if(pyq_hq.eq.0) then
                call swqmult(ipart,rrin,xx,cont,disc)
              else if(pyq_hq.eq.1) then
                call qwmassmult(mmmm,rrin,xx,cont,disc)
              else if(pyq_hq.eq.2) then
                if(id.le.3) then
                  call swqmult(ipart,rrin,xx,cont,disc)
                else if(id.ge.4) then
                  call qwmassmult(mmmm,rrin,xx,cont,disc)
                end if
              end if
            else
              if(pyq_hq.eq.0) then
                call swqmult(ipart,rrin,xxmultmax,cont,disc)
              else if(pyq_hq.eq.1) then
                call qwmassmult(mmmm,rrin,xxmultmax,cont,disc)
              else if(pyq_hq.eq.2) then
                if(id.le.3) then
                  call swqmult(ipart,rrin,xxmultmax,cont,disc)
                else if(id.ge.4) then
                  call qwmassmult(mmmm,rrin,xxmultmax,cont,disc)
                end if
              end if
*             ...puts cont=0 if xx exceeds the max value
c               call swqmult(ipart,rrin,xxmultmax,cont,disc)
               cont = 0d0
            end if
         else 
            if (xx.le.xxlinmax) then
               call swqlin(ipart,rrin,xx,cont,disc)
            else 
*             ...puts cont=0 if xx exceeds the max value
               call swqlin(ipart,rrin,xxlinmax,cont,disc)
               cont = 0d0
            end if
         end if

      end if

      return
      end
***************************************************************************
*     ARLEO's quenching weights                                           *
***************************************************************************

c ---------------------------------------------------------------------
      function dbarg(wl,e)

      implicit none 

      double precision dbarg,wl,e,dbarq

      dbarg=4.D0/9.D0*dbarq(4.D0/9.D0*wl,e)

      return
      end

c ---------------------------------------------------------------------
      function dbarq(wl,e)

      implicit none

      double precision dbarq,wl,e,xmu,xsigma

      double precision pi
      parameter(pi=3.1415926D0)

      if (wl.eq.0.D0) then
         dbarq=0.D0
      else
         dbarq=dexp(-(dlog(wl)-xmu(e))**2.D0/(2.D0*xsigma(e)**2.D0))
     +                           /(dsqrt(2.D0*pi)*xsigma(e)*wl)
      endif

      return
      end

c ---------------------------------------------------------------------
      function xmu(e)

      implicit none 

      double precision xmu, e

*     Quenching weight choice
      double precision alphas
            integer iqw,scor,ncor,sfthrd,irw
      common/qw/alphas,iqw,scor,ncor,sfthrd,irw
      
      if (ncor.eq.0) then
*       ... no 1/e corrections  
         xmu=-1.5D0
      else
*       ... with 1/e corrections  
         xmu=-1.5D0+0.81D0*(dexp(-0.2D0/e)-1.D0)
      end if

      return
      end

c ---------------------------------------------------------------------
      function xsigma(e)

      implicit none

      double precision xsigma, e

*     Quenching weight choice
      double precision alphas
            integer iqw,scor,ncor,sfthrd,irw
      common/qw/alphas,iqw,scor,ncor,sfthrd,irw

      if (ncor.eq.0) then
*       ... no 1/e corrections  
         xsigma=0.72D0
      else
*       ... with 1/e corrections  
        xsigma=0.72D0+0.33D0*(dexp(-0.2D0/e)-1.D0)
      end if

      return
      end


***************************************************************************
*     SALGADO-WIEDEMANN quenching weights                                 *
***************************************************************************

C***************************************************************************
C       Quenching Weights for Multiple Soft Scattering
C       	February 10, 2003
C
C       Refs:
C
C  Carlos A. Salgado and Urs A. Wiedemann, hep-ph/0302184.                 
C
C  Carlos A. Salgado and Urs A. Wiedemann Phys.Rev.Lett.89:092303,2002.
C 
*   NOTE: modified by A.Accardi, Dec 2004, to have 
*            gluon --> ipart=0
*            quark --> ipart<>0
*         and to initialize alphas=1/3 or 1/5 at the user's choice
C
C   This package contains quenching weights for gluon radiation in the
C   multiple soft scattering approximation.
C   swqmult returns the quenching weight for a quark (ipart<>0) or 
C   a gluon (ipart=0) traversing a medium with transport coeficient q and
C   length L. The input values are rrrr=0.5*q*L^3 and xxxx=w/wc, where
C   wc=0.5*q*L^2 and w is the energy radiated. The output values are
C   the continuous and discrete (prefactor of the delta function) parts
C   of the quenching weights.
C       
C   In order to use this routine, the files cont.all and disc.all need to be
C   in the working directory. 
C
C   An initialization of the tables is needed by doing call initmult before
C   using swqmult.
C
C   Please, send us any comment:
C
C       urs.wiedemann@cern.ch
C       carlos.salgado@cern.ch
C
C
C-------------------------------------------------------------------

      SUBROUTINE swqmult(ipart,rrrr,xxxx,continuous,discrete)
*
      REAL*8           xx(400), daq(34), caq(34,261), rrr(34)
      COMMON /dataqua/    xx, daq, caq, rrr
*
      REAL*8           xxg(400), dag(34), cag(34,261), rrrg(34)
      COMMON /dataglu/    xxg, dag, cag, rrrg

      REAL*8           rrrr,xxxx, continuous, discrete
      REAL*8           rrin, xxin
      INTEGER          nrlow, nrhigh, nxlow, nxhigh
      REAL*8           rrhigh, rrlow, rfraclow, rfrachigh
      REAL*8           xfraclow, xfrachigh
      REAL*8           clow, chigh
*

      continuous=0.d0
      discrete=0.d0

      rrin = rrrr
      xxin = xxxx
*
      do 666, nr=1,34
         if (rrin.lt.rrr(nr)) then
            rrhigh = rrr(nr)
         else
            rrhigh = rrr(nr-1)
            rrlow = rrr(nr)
            nrlow = nr
            nrhigh = nr-1
            goto 665
         endif
 666     enddo
 665     continue
*
      rfraclow = (rrhigh-rrin)/(rrhigh-rrlow)
      rfrachigh = (rrin-rrlow)/(rrhigh-rrlow)
      if (rrin.gt.10000d0) then
         rfraclow = dlog(rrhigh/rrin)/dlog(rrhigh/rrlow)
         rfrachigh = dlog(rrin/rrlow)/dlog(rrhigh/rrlow)
      endif
*
      if (ipart.ne.0.and.rrin.ge.rrr(1)) then
         nrlow=1
         nrhigh=1
         rfraclow=1
         rfrachigh=0
      endif

      if (ipart.eq.0.and.rrin.ge.rrrg(1)) then
         nrlow=1
         nrhigh=1
         rfraclow=1
         rfrachigh=0
      endif

      if (xxxx.ge.xx(260)) stop
      if (xxxx.ge.xx(260)) go to 245

      nxlow = int(xxin/0.01) + 1
      nxhigh = nxlow + 1
      xfraclow = (xx(nxhigh)-xxin)/0.01
      xfrachigh = (xxin - xx(nxlow))/0.01
*
      if(ipart.ne.0) then
         clow = xfraclow*caq(nrlow,nxlow)+xfrachigh*caq(nrlow,nxhigh)
         chigh = xfraclow*caq(nrhigh,nxlow)+xfrachigh*caq(nrhigh,nxhigh)
      else
         clow = xfraclow*cag(nrlow,nxlow)+xfrachigh*cag(nrlow,nxhigh)
         chigh = xfraclow*cag(nrhigh,nxlow)+xfrachigh*cag(nrhigh,nxhigh)
      endif

      continuous = rfraclow*clow + rfrachigh*chigh


245   continue

      if(ipart.ne.0) then
         discrete = rfraclow*daq(nrlow) + rfrachigh*daq(nrhigh)
      else
         discrete = rfraclow*dag(nrlow) + rfrachigh*dag(nrhigh)
      endif
*
      END

      subroutine initmult(alphas)

      double precision alphas

      REAL*8           xxq(400), daq(34), caq(34,261), rrr(34)
      COMMON /dataqua/    xxq, daq, caq, rrr
*
      REAL*8           xxg(400), dag(34), cag(34,261), rrrg(34)
      COMMON /dataglu/    xxg, dag, cag, rrrg
*
* 14-Dec-2016 MDB Event variable to find datafiles
      CHARACTER*255 ENVDIR
      COMMON /ENVCOM/ ENVDIR
      CHARACTER*255 FILNAM

C      if (nint(alphas*3d0).eq.1) then
C         OPEN(UNIT=20,FILE='PyQM/qweight/cont03.all',
C     &        STATUS='OLD',ERR=90)
C      else if (nint(alphas*2d0).eq.1) then
C         OPEN(UNIT=20,FILE='PyQM/qweight/cont05.all',
C     &        STATUS='OLD',ERR=90)
C      else
C         print*, 'Error (initmult): alphas =/= 1/3 or 1/2'
C         stop
C      end if

      if (nint(alphas*3d0).eq.1) then
         FILNAM=TRIM(ENVDIR)//'/PyQM/qweight/cont03.all'
      else if (nint(alphas*2d0).eq.1) then
         FILNAM=TRIM(ENVDIR)//'/PyQM/qweight/cont05.all'
      else
         print*, 'Error (initmult): alphas =/= 1/3 or 1/2'
         stop
      end if
      OPEN(UNIT=20,FILE=FILNAM,STATUS='OLD',ERR=90)
      do 110 nn=1,261
      read (20,*) xxq(nn), caq(1,nn), caq(2,nn), caq(3,nn),
     +     caq(4,nn), caq(5,nn), caq(6,nn), caq(7,nn), caq(8,nn),
     +     caq(9,nn), caq(10,nn), caq(11,nn), caq(12,nn), 
     +     caq(13,nn),
     +     caq(14,nn), caq(15,nn), caq(16,nn), caq(17,nn), 
     +     caq(18,nn),
     +     caq(19,nn), caq(20,nn), caq(21,nn), caq(22,nn), 
     +     caq(23,nn),
     +     caq(24,nn), caq(25,nn), caq(26,nn), caq(27,nn), 
     +     caq(28,nn),
     +     caq(29,nn), caq(30,nn), caq(31,nn), caq(32,nn), 
     +     caq(33,nn), caq(34,nn)
 110     continue
      do 111 nn=1,261
      read (20,*) xxg(nn), cag(1,nn), cag(2,nn), cag(3,nn),
     +     cag(4,nn), cag(5,nn), cag(6,nn), cag(7,nn), cag(8,nn),
     +     cag(9,nn), cag(10,nn), cag(11,nn), cag(12,nn), 
     +     cag(13,nn),
     +     cag(14,nn), cag(15,nn), cag(16,nn), cag(17,nn), 
     +     cag(18,nn),
     +     cag(19,nn), cag(20,nn), cag(21,nn), cag(22,nn), 
     +     cag(23,nn),
     +     cag(24,nn), cag(25,nn), cag(26,nn), cag(27,nn), 
     +     cag(28,nn),
     +     cag(29,nn), cag(30,nn), cag(31,nn), cag(32,nn), 
     +     cag(33,nn), cag(34,nn)
 111     continue
      close(20)
*
C      if (nint(alphas*3d0).eq.1) then
C         OPEN(UNIT=21,FILE='PyQM/qweight/disc03.all',
C     &        STATUS='OLD',ERR=90)
C      else if (nint(alphas*2d0).eq.1) then
C         OPEN(UNIT=21,FILE='PyQM/qweight/disc05.all',
C     &        STATUS='OLD',ERR=90)
C      else
C         print*, 'Error (initmult): alphas =/= 1/3 or 1/2'
C         stop
C      end if
      if (nint(alphas*3d0).eq.1) then
         FILNAM=TRIM(ENVDIR)//'/PyQM/qweight/disc03.all'
      else if (nint(alphas*2d0).eq.1) then
         FILNAM=TRIM(ENVDIR)//'/PyQM/qweight/disc05.all'
      else
         print*, 'Error (initmult): alphas =/= 1/3 or 1/2'
         stop
      end if
      OPEN(UNIT=21,FILE=FILNAM,STATUS='OLD',ERR=90)
      do 112 nn=1,34
      read (21,*) rrr(nn), daq(nn)
 112     continue
      do 113 nn=1,34
      read (21,*) rrrg(nn), dag(nn)
 113     continue
      close(21)
*
      goto 888
 90   PRINT*, 'input - output error' 
 91   PRINT*, 'input - output error #2' 
 888  continue

      end


C***************************************************************************
C       Quenching Weights for Single Hard Scattering
C               February 20, 2003
C
C       Refs:
C
C  Carlos A. Salgado and Urs A. Wiedemann, hep-ph/0302184.
C
C  Carlos A. Salgado and Urs A. Wiedemann Phys.Rev.Lett.89:092303,2002.
C
C
C   This package contains quenching weights for gluon radiation in the
C   single hard scattering approximation.
C
C   swqlin returns the quenching weight for a quark (ipart=1) or
C   a gluon (ipart=2) traversing a medium with Debye screening mass mu and
C   length L. The input values are rrrr=0.5*mu^2*L^2 and xxxx=w/wc, where
C   wc=0.5*mu^2*L and w is the energy radiated. The output values are
C   the continuous and discrete (prefactor of the delta function) parts
C   of the quenching weights.
C
C   In order to use this routine, the files contlinnew.all and disclinnew.all
C   need to be in the working directory.
C
C   An initialization of the tables is needed by doing call initlin before
C   using swqlin.
C
C   Please, send us any comment:
C
C       urs.wiedemann@cern.ch
C       carlos.salgado@cern.ch
C
C
C-------------------------------------------------------------------

      SUBROUTINE swqlin(ipart,rrrr,xxxx,continuous,discrete)
*
      REAL*8           xx(400), daq(34), caq(34,261), rrr(34)
      COMMON /dataqualin/    xx, daq, caq, rrr
*
      REAL*8           xxg(400), dag(34), cag(34,261), rrrg(34)
      COMMON /dataglulin/    xxg, dag, cag, rrrg

      REAL*8           rrrr,xxxx, continuous, discrete
      REAL*8           rrin, xxin
      INTEGER          nrlow, nrhigh, nxlow, nxhigh
      REAL*8           rrhigh, rrlow, rfraclow, rfrachigh
      REAL*8           xfraclow, xfrachigh
      REAL*8           clow, chigh
*

      continuous=0.d0
      discrete=0.d0

      rrin = rrrr
      xxin = xxxx
*
      do 666, nr=1,34
         if (rrin.lt.rrr(nr)) then
            rrhigh = rrr(nr)
         else
            rrhigh = rrr(nr-1)
            rrlow = rrr(nr)
            nrlow = nr
            nrhigh = nr-1
            goto 665
         endif
 666     enddo
 665     continue
*
      rfraclow = (rrhigh-rrin)/(rrhigh-rrlow)
      rfrachigh = (rrin-rrlow)/(rrhigh-rrlow)
      if (rrin.gt.10000d0) then
         rfraclow = dlog(rrhigh/rrin)/dlog(rrhigh/rrlow)
         rfrachigh = dlog(rrin/rrlow)/dlog(rrhigh/rrlow)
      endif
*
      if (ipart.eq.1.and.rrin.ge.rrr(1)) then
         nrlow=1
         nrhigh=1
         rfraclow=1
         rfrachigh=0
      endif

      if (ipart.ne.1.and.rrin.ge.rrrg(1)) then
         nrlow=1
         nrhigh=1
         rfraclow=1
         rfrachigh=0
      endif

      if (xxxx.ge.xx(260)) go to 245

      nxlow = int(xxin/0.038) + 1
      nxhigh = nxlow + 1
      xfraclow = (xx(nxhigh)-xxin)/0.038
      xfrachigh = (xxin - xx(nxlow))/0.038
*
      if(ipart.eq.1) then
      clow = xfraclow*caq(nrlow,nxlow)+xfrachigh*caq(nrlow,nxhigh)
      chigh = xfraclow*caq(nrhigh,nxlow)+xfrachigh*caq(nrhigh,nxhigh)
      else
      clow = xfraclow*cag(nrlow,nxlow)+xfrachigh*cag(nrlow,nxhigh)
      chigh = xfraclow*cag(nrhigh,nxlow)+xfrachigh*cag(nrhigh,nxhigh)
      endif

      continuous = rfraclow*clow + rfrachigh*chigh

245   continue

      if(ipart.eq.1) then
      discrete = rfraclow*daq(nrlow) + rfrachigh*daq(nrhigh)
      else
      discrete = rfraclow*dag(nrlow) + rfrachigh*dag(nrhigh)
      endif
*
      END

      subroutine initlin(alphas)

      double precision alphas

      REAL*8           xxq(400), daq(34), caq(34,261), rrr(34)
      COMMON /dataqualin/    xxq, daq, caq, rrr
*
      REAL*8           xxg(400), dag(34), cag(34,261), rrrg(34)
      COMMON /dataglulin/    xxg, dag, cag, rrrg

* 14-Dec-2016 MDB Event variable to find datafiles
      CHARACTER*255 ENVDIR
      COMMON /ENVCOM/ ENVDIR
      CHARACTER*255 FILNAM
*
      if (nint(alphas*3d0).eq.1) then
         FILNAM=TRIM(ENVDIR)//'/PyQM/qweight/contlin03.all'
C         OPEN(UNIT=20,FILE='PyQM/qweight/contlin03.all',
C     &        STATUS='OLD',ERR=90)
         OPEN(UNIT=20,FILE=FILNAM,STATUS='OLD',ERR=90)
      else if (nint(alphas*2d0).eq.1) then
*         OPEN(UNIT=20,FILE='contlin05.all',STATUS='OLD',ERR=90)
         print*, 'Error (initlin): alphas=0.5 not yet implemented'
         stop
      else
         print*, 'Error (initlin): alphas =/= 1/3 or 1/2'
         stop
      end if
      do 110 nn=1,261
      read (20,*) xxq(nn), caq(1,nn), caq(2,nn), caq(3,nn),
     +     caq(4,nn), caq(5,nn), caq(6,nn), caq(7,nn), caq(8,nn),
     +     caq(9,nn), caq(10,nn), caq(11,nn), caq(12,nn), 
     +     caq(13,nn),
     +     caq(14,nn), caq(15,nn), caq(16,nn), caq(17,nn), 
     +     caq(18,nn),
     +     caq(19,nn), caq(20,nn), caq(21,nn), caq(22,nn), 
     +     caq(23,nn),
     +     caq(24,nn), caq(25,nn), caq(26,nn), caq(27,nn), 
     +     caq(28,nn),
     +     caq(29,nn), caq(30,nn), caq(31,nn), caq(32,nn), 
     +     caq(33,nn), caq(34,nn)
 110     continue
      do 111 nn=1,261
      read (20,*) xxg(nn), cag(1,nn), cag(2,nn), cag(3,nn),
     +     cag(4,nn), cag(5,nn), cag(6,nn), cag(7,nn), cag(8,nn),
     +     cag(9,nn), cag(10,nn), cag(11,nn), cag(12,nn), 
     +     cag(13,nn),
     +     cag(14,nn), cag(15,nn), cag(16,nn), cag(17,nn), 
     +     cag(18,nn),
     +     cag(19,nn), cag(20,nn), cag(21,nn), cag(22,nn), 
     +     cag(23,nn),
     +     cag(24,nn), cag(25,nn), cag(26,nn), cag(27,nn), 
     +     cag(28,nn),
     +     cag(29,nn), cag(30,nn), cag(31,nn), cag(32,nn), 
     +     cag(33,nn), cag(34,nn)
 111     continue
      close(20)
*
      if (nint(alphas*3d0).eq.1) then
         FILNAM=TRIM(ENVDIR)//'/PyQM/qweight/disclin03.all'
C         OPEN(UNIT=21,FILE='PyQM/qweight/disclin03.all',
C     &        STATUS='OLD',ERR=91)
         OPEN(UNIT=21,FILE=FILNAM,STATUS='OLD',ERR=91)
      else if (nint(alphas*2d0).eq.1) then
*         OPEN(UNIT=21,FILE='disclin05.all',STATUS='OLD',ERR=91)
         print*, 'Error (initlin): alphas=0.5 not yet implemented'
         stop
      else
         print*, 'Error (initlin): alphas =/= 1/3 or 1/2'
         stop
      end if
      do 112 nn=1,34
      read (21,*) rrr(nn), daq(nn)
 112     continue
      do 113 nn=1,34
      read (21,*) rrrg(nn), dag(nn)
 113     continue
      close(21)
*
      goto 888
 90   PRINT*, 'data file conlin03 open error' 
 91   PRINT*, 'data file disclin03 open error #2' 
 888  continue

      end


