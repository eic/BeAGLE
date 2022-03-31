      subroutine ApplyQW(qhat)
      implicit none

      include 'common.f'
      include 'bea_pyqm.inc'

      integer ip,iq,ir,iit ! For do
      integer iloop !current entry number
      integer iEg,iPtF,iqg

      double precision qhat,SupFac,iet

      double precision th,ph

      double precision inix,iniy,iniz,iniE !initial parton 4-momentum
      double precision ipx,ipy,ipz,ipt    !final parton momentum
      double precision iE,iEnew,iEgluon

      double precision iptot,ipl,ptg,plg, cr

      double precision ptot,pt,pt2,mmmm,mass_p
      double precision N_const, w_n, w_high, w_mean, n_w_gluons
      double precision w_sum, w_hard, w_soft
      double precision w_hard_2,w_triplet
      integer n_w, list, ij, i_w
      double precision w_gluon(50)
      integer ijoin(50)
      double precision theta_gluon, phi_gluon
      double precision soft_cut
      double precision E_quark,E_gluon,E_antiq,x_1,x_3
      double precision E_p,E_loss
      double precision MMx,MMy,MMz,MMe
      integer j
      alphas   = 1d0/3d0
      iqw      = 1
      scor     = 1
      ncor     = 0
      sfthrd   = 1       !multiple soft scattering
      pyq_hq   = PYQ_HQ  !New SW calculation for heavy quarks (0(no) or 1(yes))
      soft_cut = 5       !energy cut on soft gluons(GeV)
      iEg      = PYQ_IEG !option gluon model
      iPtF     = PYQ_IPTF!option pt model
      iet      = PYQ_IET !energy cut off(GeV)
      iqg      = 1
      ipt      = 0       !transverse momentum
      cr       = 0
      iloop    = N       !number of partons, given by pythia
      ip       = 1
      ij       = 0       !Parton counter of the string joined after adding a gluon
      ijoin    = 0       !Array with the positions of each parton in the string
      E_p      = 0

      if(qhat.lt.0.00001) return

      do while (ip.le.iloop)
        mass_p =P(ip,5)
cccc   Select partons
       if ((K(ip,1).eq.2).or.(K(ip,1).eq.1))then 
         if(K(ip,1).eq.2) then
cccc     parton counter 
           ij = ij+1
           ijoin(ij)=ip
          endif
        if((abs(K(ip,2)).eq.21).or.(abs(K(ip,2)).le.5))then
cccc   Initialize energies
            w_gluon(1)    = 0.0
            w_gluon(2)    = 0.0
            w_gluon(3)    = 0.0
            w_hard        = 0.0
            w_soft        = 0.0
            w_triplet     = 0.0
            E_p           = 0.0
cccc   initial values        
            inix = P(ip,1)
            iniy = P(ip,2)
            iniz = P(ip,3)
            iniE = P(ip,4)
        
          mmmm = P(ip,5)/P(ip,4)
cccc   calculate QW with a given qhat          
        call QWComput(qhat,P(ip,1),P(ip,2),P(ip,3),P(ip,4),mmmm,
     & K(ip,2))
        if((QW_w.gt.0.00001).and.(iniE.gt.iet))then
cccc  Calculate transverse momentum
           call PTF(ip,iPtf,QW_w,ipt)
           
           call newPartonKinematics(ip,iet,ipt,QW_w,inix,iniy,iniz
     &,iniE,ipx,ipy,ipz,iEnew) 
cccc    Define cr for a quark or gluon
           if(abs(K(ip,2)).le.5) then
             cr = 4d0/3d0
           else if(K(ip,2).eq.21) then
             cr=3d0
           endif
           N_const = 4d0*alphas*alphas*qhat   

cccc     Calculate gluon energy
            E_loss=iniE-P(ip,4)
          if (E_loss.lt.0.0)then
            print*,'ERROR w_gluon < 0'
            exit
          endif

          MMx = inix - P(ip,1)
          MMy = iniy - P(ip,2)
          MMz = iniz - P(ip,3)
          MMe = iniE - P(ip,4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             no gluons radiation                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          if(iEg.eq.0) then

cccc     With this option you lose only the energy of the partons 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             1 hard gluon                                            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          else if(iEg.eq.1) then
            if(E_loss.gt.0.0)then
cccc     calculate gluon kinematics
                call GluonKinematics(ip,E_loss,inix,iniy,iniz,iniE
     &,theta_gluon,phi_gluon)
cccc     add a gluon
                call PY1ENT(N+1,21,E_loss,theta_gluon,phi_gluon)
cccc     parton counter                 
                ij=ij+1
                ijoin(ij)=N
            endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             1 hard gluons + soft                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          else if(iEg.eq.2) then
cccc    Calculate energies
            w_hard = N_const*(QW_L**2)
            if (w_hard.gt.E_loss) w_hard=E_loss
cccc   Constraints for w_hard
            if(w_hard.gt.iet)then 
               if ((iniE-iet).gt.w_hard)then
                  w_gluon(1)= w_hard
               else 
                  w_hard=iniE-iet
                  w_gluon(1)=w_hard
               endif
            else 
                  w_hard = 0
            endif
cccc   Calculate Energy Softs gluons and Triplet              
            w_soft= E_loss - w_hard
            w_gluon(2)=w_soft

            if (w_soft.gt.soft_cut) then
              w_triplet= w_soft-soft_cut
cccc   Calculating energies of q-qbar-g
              call TripletEnergies(w_triplet,E_quark,E_gluon,
     & E_antiq,x_1,x_3)
cccc   Adding a triplet
              call PY3ENT(N+1,2,21,-2,w_triplet,x_1,x_3)
              MMx = MMx - P(N,1) - P(N-1,1) - P(N-2,1)
              MMy = MMy - P(N,2) - P(N-1,2) - P(N-2,2)
              MMz = MMz - P(N,3) - P(N-1,3) - P(N-2,3)
              MMe = MMe - P(N,4) - P(N-1,4) - P(N-2,4)
            endif

            if((w_soft.lt.0).and.(w_triplet.lt.0)) then
               print*,'ERROR unexpected negative values for eloss'
               exit
            endif
cc     calculate gluon kinematics
            if(w_hard.gt.0) then
              call GluonKinematics(ip,w_hard,inix,iniy,iniz,iniE
     &,theta_gluon,phi_gluon)
              call PY1ENT(N+1,21,w_hard,theta_gluon,phi_gluon)
              ij=ij+1
              ijoin(ij)=N

cccc  4-mom going back to the remnant nuclei                
              MMx = MMx -P(N,1)
              MMy = MMy -P(N,2)
              MMz = MMz -P(N,3)
              MMe = MMe -P(N,4)
             endif
               
             if(MMe.lt.0.0) print*,'ERROR unexpected negative energy'

             PYQREC(1) = PYQREC(1)+MMx
             PYQREC(2) = PYQREC(2)+MMy
             PYQREC(3) = PYQREC(3)+MMz
             PYQREC(4) = PYQREC(4)+MMe
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            softs gluons                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           else if(iEg.eq.3) then

cccc   Calculate energies
             w_gluon(1)=E_loss
             w_soft=E_loss

              if (w_soft.gt.soft_cut) then
                w_triplet= w_soft-soft_cut
                w_soft=w_soft-w_triplet
                call TripletEnergies(w_triplet,E_quark,E_gluon,
     & E_antiq,x_1,x_3)
                call PY3ENT(N+1,2,21,-2,w_triplet,x_1,x_3)
              MMx = MMx - P(N,1) - P(N-1,1) - P(N-2,1)
              MMy = MMy - P(N,2) - P(N-1,2) - P(N-2,2)
              MMz = MMz - P(N,3) - P(N-1,3) - P(N-2,3)
              MMe = MMe - P(N,4) - P(N-1,4) - P(N-2,4)
              endif

             if(MMe.lt.0.0) print*,'ERROR unexpected negative energy'

             PYQREC(1) = PYQREC(1)+MMx
             PYQREC(2) = PYQREC(2)+MMy
             PYQREC(3) = PYQREC(3)+MMz
             PYQREC(4) = PYQREC(4)+MMe

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             closing the options 0,1,2 and 3                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccc  closing options iEg=0 or iEg=1 or iEg=2 or iEg=3
           endif
cccc  endif of QW_w>0
        endif 
cccc  endif gluon or quark selection
      endif
cccc     Last Parton of the string
            if(K(ip,1).eq.1) then
cccc     parton counter 
             ij=ij+1
             ijoin(ij)=ip
cccc       Connecting a number of previously defined partons into
cccc       a string configuration
            if(ij.ge.3) then
                call PYJOIN(ij,ijoin)
            endif
cccc      ij and ijoin comeback to zero because it is the last parton KS=1
cccc      of the string
            ijoin = 0
            ij = 0
cccc     closing the particles selection for KS=1 or KS=2
           endif

ccccc  closing initial partons's selection for KS=1 or KS=2
       endif
       ip=ip+1
ccccc closing the particle loop       
      end do
ccccc closing ApplyQW
      end


       subroutine PTF(ip,iPtf,Eloss,ipt)
       include 'common.f'
       include 'bea_pyqm.inc'
       integer ip,iPtf
       double precision ipt
       double precision Eloss        
ccccc    Calculate transverse momentum of final parton IPtf
        if (iPtf.eq.0) then
           ipt=0
        else if(iPtf.eq.1) then
           ipt=qhat*QW_L
        else if(iPtf.eq.2) then
           ipt = ((8d0/3)*Eloss/alphas)/(QW_L**2) !BDMPS mean
        else if(iPtf.eq.3) then
           ipt = (Eloss*sin(QW_th))**2
        endif

       end

       subroutine newPartonKinematics(ip,iet,ipt,Eloss,inix,iniy,iniz
     & ,iniE,ipx,ipy,ipz,iEnew)
       
       include 'common.f'
       include 'bea_pyqm.inc'
       integer ip
       double precision iet
       double precision ipt,Eloss
       double precision th,ph,ipl,iptot2
       double precision iptx,ipty,iptz,tot
       double precision ipgx,ipgy,ipgz,ipg
       double precision inix,iniy,iniz,iniE
       double precision ipix,ipiy,ipiz,itot
       double precision ipx,ipy,ipz,iEnew,ipnew
ccccc  Stock init mom values
       inix = P(ip,1)
       iniy = P(ip,2)
       iniz = P(ip,3)
       iniE = P(ip,4)
ccccc  Normalized init mom values
       itot = sqrt(inix**2+iniy**2+iniz**2)
       ipix = inix/itot
       ipiy = iniy/itot
       ipiz = iniz/itot    
ccccc  mplement ELoss and Pt
       if(P(ip,4)-Eloss.lt.iet) then
C          print*,'Eparton-Eloss < iet', P(ip,4)-Eloss
          th = ranf(0)*2*3.14159265-3.14159265
          ipl = cos(th)*iet !longitudinal parton mom
          ipt = sin(th)*iet !transverse parton mon
       else if (P(ip,4)-Eloss.ge.iet) then
          iptot2 = ((P(ip,4)-Eloss)**2-(P(ip,5))**2)
          if(iptot2.gt.ipt) then
             ipl = iptot2-ipt
             ipt = sqrt(ipt)
             ipl = sqrt(ipl)
          else
             ipl = 0
             ipt = sqrt(iptot2)
          endif
       endif
ccccc     Generate normalized transverse vector
       ph   = 4*asin(1.)*ranf(0)
       iptx = (ipiz-ipiy)*cos(ph)-(ipix*ipiy + ipix*ipiz)*sin(ph)
       ipty = ipix*cos(ph) + (ipix**2+ipiz**2-ipiy*ipiz)*sin(ph)
       iptz = - ipix*cos(ph) + (ipix**2+ipiy**2-ipiy*ipiz)*sin(ph)
       tot  = sqrt(iptx**2+ipty**2+iptz**2)
       iptx = iptx/tot
       ipty = ipty/tot
       iptz = iptz/tot
ccccc  Generate new parton momenta
       ipx = ipt*iptx+ipl*ipix
       ipy = ipt*ipty+ipl*ipiy
       ipz = ipt*iptz+ipl*ipiz
       ipnew = sqrt(ipx**2+ipy**2+ipz**2)
       iEnew = sqrt(P(ip,5)**2+ipx**2+ipy**2+ipz**2)
       
ccccc  Fill Pythia array
ccccc  new parton 4-mom
       if (iEnew.lt.iniE) then
         P(ip,1) = ipx
         P(ip,2) = ipy
         P(ip,3) = ipz
         P(ip,4) = iEnew
       else
          P(ip,1) = inix !initial 4-momentum NO Eloss
          P(ip,2) = iniy
          P(ip,3) = iniz
          P(ip,4) = iniE
          print*, ' NO Energy LOSS there is a problem'
       endif
       end
       
        subroutine GluonKinematics(ip,Egluon,inix,iniy,iniz,iniE
     & ,theta_gluon,phi_gluon)

        include 'common.f'
        include 'bea_pyqm.inc'
        double precision theta_gluon, phi_gluon
        integer ip
        double precision ipgx,ipgy,ipgz,ipg
        double precision inix,iniy,iniz,iniE
        double precision Egluon

ccccc   calculate gluon mom as old mom parton - new mom parton          
        ipgx = inix - P(ip,1)
        ipgy = iniy - P(ip,2)
        ipgz = iniz - P(ip,3)
        ipg = sqrt((ipgx)**2+(ipgy)**2+(ipgz)**2)
       
ccccc   Angles
        theta_gluon = acos(ipgz/ipg)
        phi_gluon = atan2(ipgy,ipgx)
       
        end

cccc   E_total : (Ecm) the total energy of the system.
cccc   X1, X3 : xi = 2Ei/Ecm, i.e. twice the energy fraction taken by the ith parton.
      subroutine TripletEnergies(w_triplet,E_quark,E_gluon,E_antiq,
     & x_1,x_3)
       include 'common.f'
       include 'bea_pyqm.inc'

       double precision w_triplet
       double precision E_quark,E_gluon,E_antiq,E_total
       double precision x_1,x_3

       E_quark=(4d0/10d0)*w_triplet
       E_antiq=(4d0/10d0)*w_triplet
       E_gluon=(2d0/10d0)*w_triplet


       E_total=E_quark+E_antiq+E_gluon
       
       x_1=(2d0/w_triplet)*E_quark
       x_3=(2d0/w_triplet)*E_antiq


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
      double precision qhat

      double precision mmmm !mass/energy of the incoming parton
      integer irej !used for test
      double precision QW_wc_2
      double precision I_QW_wc,I_QW_R
      QW_wc_2 = 0.
      QW_w    = 0.
      QW_L    = 0.
      QW_wc   = 0.
      QW_R    = 0.
      I_QW_wc = 0.
      I_QW_R  = 0.
      d       = 0.
      qhateff = qhat
      ChiR    = 0.

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
ccc  interaction point give by pythia 
      x = x_inter
      y = y_inter
      z = z_inter


ccc integration to calculate wc and R
      radius = sqrt(x**2+y**2+z**2)
      do while (radius.lt.20)
        I_QW_wc = I_QW_wc + integral_step * d *
     &          density_table(INT(radius/step_size_dens))
        I_QW_R = I_QW_R + integral_step *
     &          density_table(INT(radius/step_size_dens))
        d = d + integral_step
        x = x + px
        y = y + py
        z = z + pz
        radius = sqrt(x**2+y**2+z**2)
      enddo
ccc    calculate average of L,R an wc       
       QW_L =(2d0* I_QW_wc) / I_QW_R
       QW_R = 2 * density_table(1) * I_QW_wc**2 / I_QW_R / qhateff
       QW_wc = (qhateff/density_table(1)) * I_QW_wc
 

ccccc Convert the units fm -> GeV-1
c      QW_L = QW_L/.1973269
ccccc Convert the units GeV2.fm -> GeV
      QW_wc = QW_wc/.1973269
ccccc Convert the units GeV2.fm2 -> no unit
      QW_R = QW_R /.1973269**2

ccccc Calculate the energy loss probability
      if(sfthrd.eq.1) step_QW = 2.5/nb_step
      if(sfthrd.eq.2) step_QW = 9.8/nb_step
      yy = E/QW_wc

      total = 0.
      do i=1,nb_step
        xx = step_QW * i
        call qweight(ipart,id,mmmm,QW_R,xx,yy,cont(i),disc)
        total = total + cont(i)*step_QW
      enddo
      

      total = total + disc
      disc = disc/total

      do i=1,nb_step
        cont(i) = cont(i) / total
      enddo
ccccc Pick randomely a quenching in the table
      if(disc .lt. 1.) then
        randnum = ranf(0)
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
      if(QW_w .gt. 0) then
        step_QW = 1./nb_step
        yy = E/QW_wc
        xx = QW_w/QW_wc 
       
        total = 0.
        do i=1,nb_step
          ChiR =(step_QW * i)**2 * QW_R
          call qweight(ipart,id,mmmm,ChiR,xx,yy,cont(i),disc)
cc Do not keep negative probabilities
          if (cont(i).lt.0) cont(i) = 0
          total = total + cont(i)*step_QW
        enddo
        irej=1
        if(total.lt.1e-5) irej=0
        if(irej.eq.0) print*,'total=',total,'QW_w/c=',QW_w,' ',QW_wc
        do i=1,nb_step
          cont(i) = cont(i) / total
        enddo
        randnum = ranf(0)
        total = 0.
        i = 1
        do while (randnum.gt.total)
          total = total + cont(i)*step_QW
          i = i + 1
        enddo
        QW_chi = i * step_QW
        QW_th = asin(QW_chi)
        if(irej.eq.0) print*,' i=',i
        if(irej.eq.0) print*,'QW_chi=',QW_chi,' QW_th=',QW_th
      endif
      if(isnan(QW_th)) QW_th = pi/2 !if QW_Chi is 1
c 56378 continue
c     The continue is temporary (test)
ccccc Scattering angle set to 0
c      if(QW_w .gt. 0) then
c         QW_th=0.  ! collinear gluon radiation assumption
c      endif
c      print*, 'QW_th = ',QW_th
c      print*, 'QW_w =',QW_w

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
C         OPEN(UNIT=22,FILE='PyQM/qweight/disc03.all',
C     &        STATUS='OLD',ERR=90)
C      else if (nint(alphas*2d0).eq.1) then
C         OPEN(UNIT=22,FILE='PyQM/qweight/disc05.all',
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
      OPEN(UNIT=22,FILE=FILNAM,STATUS='OLD',ERR=90)
      do 112 nn=1,34
      read (22,*) rrr(nn), daq(nn)
 112     continue
      do 113 nn=1,34
      read (22,*) rrrg(nn), dag(nn)
 113     continue
      close(22)
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
C         OPEN(UNIT=22,FILE='PyQM/qweight/disclin03.all',
C     &        STATUS='OLD',ERR=91)
         OPEN(UNIT=22,FILE=FILNAM,STATUS='OLD',ERR=91)
      else if (nint(alphas*2d0).eq.1) then
*         OPEN(UNIT=22,FILE='disclin05.all',STATUS='OLD',ERR=91)
         print*, 'Error (initlin): alphas=0.5 not yet implemented'
         stop
      else
         print*, 'Error (initlin): alphas =/= 1/3 or 1/2'
         stop
      end if
      do 112 nn=1,34
      read (22,*) rrr(nn), daq(nn)
 112     continue
      do 113 nn=1,34
      read (22,*) rrrg(nn), dag(nn)
 113     continue
      close(22)
*
      goto 888
 90   PRINT*, 'data file conlin03 open error' 
 91   PRINT*, 'data file disclin03 open error #2' 
 888  continue

      end


