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

      double precision ptot,pt,pt2,mmmm,mass_p
      double precision N_const, w_n, w_high, w_mean, n_w_gluons
      double precision N_const_2, N_const_3
      double precision w_sum, w_hard, w_soft,w_soft_remain
      double precision w_hard_2,w_triplet
      integer n_w, list, ij, i_w
      double precision w_gluon(50)
      integer ijoin(50)
      double precision pgx, pgy, pgz, phe, phi
      double precision phi_final, theta, theta_final, mag
      double precision soft_cut
      double precision E_quark,E_gluon,E_antiq,x_1,x_3
      double precision E_p
      integer j
      integer k_cut_up,k_cut_down

      alphas   = 1d0/3d0
      iqw      = 1
      scor     = 1
      ncor     = 0
      sfthrd   = 1
      soft_cut = 5       !energy cut on soft gluons(GeV)
      cutoff   = pyq_iet !energy cut off (GeV)
      iEg      = PYQ_IEG !option gluon model
      iPtF     = PYQ_IPTF!option pt model
      iet      = PYQ_IET !energy cut(GeV)
      iqg      = 1
      ipt      = 0       !transverse momentum
      
      ehat     = 0
      cr       = 0

      iloop = N
      ip = 1
      ij =0
      k_cut_up=0
      k_cut_down=0
c     print*,'Starting PYQM'
      iet           = 0.25
      w_gluon(1)    = 0.0
      w_gluon(2)    = 0.0
      w_gluon(3)    = 0.0
      w_hard        = 0.0
      w_soft        = 0.0
      w_soft_remain = 0.0
      w_triplet     = 0.0
      E_p           = 0.0
      


      do while (ip.le.iloop)
        mass_p =P(ip,5)
       
c        print*,'ip           --> ',ip
c        print*,'iloop        --> ',iloop
c        print*,'particle ID  --> ',K(ip,2)
c        print*,'just before the IF'
        
c        print*,'MASS=',mass_p
        if((K(ip,1).eq.2).or.((K(ip,1).eq.1)).and.
     & ((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then
c          print*,'after sekection id= ',K(ip,2)
        
          mmmm = P(ip,5)/P(ip,4)
c          print*,'mmmm =',mmmm
          call QWComput(qhat,P(ip,1),P(ip,2),P(ip,3),P(ip,4),mmmm,
     & K(ip,2))
c         print*,'QWComput -> Was called'
c Define cr for a quark or gluon

          if(abs(K(ip,2)).le.5) then
            cr = 4d0/3d0
          else if(K(ip,2).eq.21) then
            cr=3d0
          endif
c Calculate transverse momentum of final parton IPtf
          if (iPtf.eq.0) then
          ipt=0
          else if(iPtf.eq.1) then
          ipt=qhat*QW_L
          else if(iPtf.eq.2) then
          ipt = ((8d0/3)*QW_w/alphas)/(QW_L**2) !BDMPS mean
          else if(iPtf.eq.3) then
          ipt = (QW_w*sin(QW_th))**2
          endif
c          print*,'Transverse momentum=',ipt
         
c Calculate constant for hard energy gluon
          N_const = 4d0*alphas*alphas*qhat
c         print*,'N_const =',N_const

c          print*, 'QW_w Not recalculated = ', QW_w
         if(((P(ip,4)-QW_w).lt.iet).and.(QW_w.gt.0)) then
            QW_w=P(ip,4)-iet
           
c            print*,'Particle E     = ', P(ip,4)
c            print*,'iet            = ', iet
c            print*,'QW_w = E - iet = ', QW_w
          endif
c        print*,'w_gluon(1)         = ',w_gluon(1)

       
       endif
c           print*,'Energy Threshold:',iet
c           print*,'Starting GLUON SELECTION'
           

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             no gluons                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          if(iEg.eq.0) then

            w_gluon(1)=QW_w
c            Select quarks,antiquarks and gluons in an intermidiate state.

            if((K(ip,1).eq.2).and.((abs(K(ip,2)).le.5).or.
     &  (K(ip,2).eq.21)).and.(P(ip,4).gt.PYQ_IET).and.(QW_w.gt.0)) then
              
c             counter partons  
              ij=ij+1
cccccccccc    calculating gluon properties
              call GluonEmission(ip,w_gluon(1),theta_final,phi_final)

              n=n+1
              ij=ij+1

              call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
ccccccccc     P(N,*) gluon 4-mom
ccccccccc     P(ip,*) new parton 4-mom
              P(ip,1)=P(ip,1)-P(N,1)
              P(ip,2)=P(ip,2)-P(N,2)
              P(ip,3)=P(ip,3)-P(N,3)

            else if((K(ip,1).eq.1).and.(ij.gt.1).and.
     & ((abs(K(ip,2)).le.5).and.(P(ip,4).gt.(PYQ_IET))).and.(QW_w.gt.
     & 0 )) then
              n=n+1
              ij=ij+1

              call GluonEmission(ip,w_gluon(1),theta_final,phi_final)

c             Pythia function adding 1 hard-gluon (not apply to this option)

              call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
c             Momentum of the partons less gluon momentum 

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
c             1 hard gluon                                            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


          else if(iEg.eq.1) then

            w_gluon(1)=QW_w
cccccccc  E_p is the new parton energy
            E_p = P(ip,4)-w_gluon(1)

            if(K(ip,1).eq.2) then
ccccccc   partons counter
              ij = ij+1
              ijoin(ij)=ip
cccccccc Select quarks,antiquarks and gluons in an intermidiate state.
              if(((abs(K(ip,2)).le.5)
     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
     & .and.(((P(ip,4).ge.iet).and.(QW_w.gt.0))
     & .and.(E_p.gt.(2d0*iet)))) then

                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n=n+1
cccccccc          CALL PY1ENT(IP,KF,PE,THE,PHI)
cccccccc          Purpose: to add one entry to the event record, i.e. either 
cccccccc          a parton or a particle.
cccccccc          IP=parton/particle flavour code.
cccccccc          KF=parton/particle energy. If PE is smaller than the mass, the parton/particle is taken to be at rest.
cccccccc          THE, PHI : polar and azimuthal angle for the momentum vector of the parton/particle.
                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
                  ij=ij+1

cccccccc          Momentum new parton:
                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)
                  ijoin(ij)=n


              endif

            else if(K(ip,1).eq.1) then


              if((abs(K(ip,2)).le.5)
     & .and.(((P(ip,4).ge.iet).and.(QW_w.gt.0))
     & .and.(E_p.gt.(2d0*iet)))) then


                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n=n+1
                  ij=ij+1
                 

cccccccc          Adding one gluon 
                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
c                 momentum of the parton less gluon momenta
                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)
                  ijoin(ij)=n

              endif
            

            ij=ij+1
            ijoin(ij)=ip
ccccc       ijoin purpuse to be used in PYJOIN
ccccc       Connecting a number of previously defined partons into
ccccc       a string configuration
ccccc       ij:number of entries making up the string formed by PYJOIN
ccccc       ijoin: one dimension array of size at least ij
ccccc       if we find more than 3 partons they will be join with this routine, after
ccccc       that all the partons will have KS=3

            if(ij.ge.3) then
                call PYJOIN(ij,ijoin)
c                
            endif

            ij=0

          endif



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             1 hard gluons + soft                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          else if(iEg.eq.2) then
     
           if ((QW_w.gt.0).and.((K(ip,1).eq.1).or.(K(ip,1).eq.2))
     & .and.((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then

ccccccccccc   Calculate energy HARD gluon
              
               w_hard = N_const*(QW_L**2)
ccccccccccc   Calculate energy Softs gluons and Triplet
ccccccccccc   The energy of the soft gluons cant be greater than soft_cut =5 GeV,
ccccccccccc   if there is energy left between the soft gluons and the cut
ccccccccccc   we create  a triplet.

              if(w_hard.gt.(2d0*iet))then 
                 if (((P(ip,4)-(2d0*iet))).gt.w_hard)then
                    w_gluon(1)= w_hard
c                    print*,'w_hard > 2eit, w_hard =',w_hard
c                    print*, 'w_gluon(1) =',w_gluon(1)
                 else if (((P(ip,4)-iet).lt.w_hard)
     & .and.(P(ip,4).gt.(iet)))then
                    w_hard=P(ip,4)-iet
                    w_gluon(1)=w_hard
c                    print*,'w_hard < Ep-iet, w_hard =',w_hard
c                    print*, 'w_gluon(1) =',w_gluon(1)
                 endif

              else if ((w_hard.lt.(2d0*iet)).and.(w_hard.gt.iet))then
                  w_hard=w_soft_remain
c                  print*, 'w_hard < 2eit & w_hard > iet, w_hard='
c     &, w_hard
c                  print*, 'w_gluon(2) =',w_gluon(2)
c                  print*, 'w_soft     =',w_soft_remain

c                  if((P(ip,4)-iet).gt.w_soft)then
c                     w_gluon(2)=w_soft_remain
c                     print*, 'First loop soft gluons'
c                     print*,'we have a soft g = ',w_soft_remain,'GeV'
c                  endif
               
               
              else if(w_hard.lt.iet)then
c                  print*,'Gluon energy under the threshold'
c                  print*,'Not created'
              endif
              
              if (w_hard.lt.QW_w) then
                 w_soft= (QW_w - w_hard)+w_soft_remain
                 w_gluon(2)=w_soft
                 if (w_soft.gt.soft_cut) then
                    w_triplet= w_soft-soft_cut
                    w_gluon(3)=w_triplet
                    if((w_soft.lt.0).and.(w_triplet.lt.0)) then
                       w_soft = 0
                       w_triplet = 0
c                       print*,'no soft gluons'
                    endif
                  endif
                endif
              
            endif
cccccccc  Energy of the parton
             E_p = P(ip,4)-w_hard
cccccccc  Counting            
             if(w_soft.gt.soft_cut) then
                k_cut_up = k_cut_up+1
             else if(w_soft.lt.soft_cut) then 
                k_cut_down =k_cut_down +1
             endif 
            
cccccccc   Starting new selection of partons for status KS=2 
            if(K(ip,1).eq.2) then

              ij = ij+1
              ijoin(ij)=ip


              if(((abs(K(ip,2)).le.5)
     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
     & .and.(((P(ip,4).ge.iet).and.(QW_w.gt.0)
     & .and.(w_soft.le.soft_cut).and.(w_hard.gt.0.0001)
     & .and.(E_p.gt.(2d0*iet))))) then

                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)
        
                  n=n+1
                  ij=ij+1
                                    
                  
                   call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
                  

cccccccc    Calculating PYQREC: 4-momentum going to the nuclei remnant
cccccccc    Energy of SOFT Gluons and QW left
                  
                  PYQREC(1)=PYQREC(1)
     & +QW_w*sin(theta_final)*cos(phi_final)-P(N,1)
                  PYQREC(2)=PYQREC(2)
     & +QW_w*sin(theta_final)*sin(phi_final)-P(N,2)
                  PYQREC(3)=PYQREC(3)+QW_w*cos(theta_final)-P(N,3)
                  PYQREC(4)=PYQREC(4)+w_soft
cccccccc   Calculating new parton momentum
cccccccc   P(ip,4) is being calculated in Gluon emission

                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)
       

                ijoin(ij)=n
        
              endif
cccccccc Starting new selection of partons for status KS=1
            else if(K(ip,1).eq.1) then

              if((abs(K(ip,2)).le.5)
     & .and.(((P(ip,4).ge.iet).and.(QW_w.gt.0)
     & .and.(w_soft.le.soft_cut).and.(w_hard.gt.0.0001)
     & .and.(E_p.gt.(2d0*iet))))) then

                  call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)
                  
                  n=n+1
                  ij=ij+1

                   call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)
                  
ccccccc    Calculating PYQREC: 4-momentum going back the nuclei remnant
ccccccc    Energy of SOFT Gluons and QW left

                  PYQREC(1)=PYQREC(1)     
     & +QW_w*sin(theta_final)*cos(phi_final)-P(N,1)
                  PYQREC(2)=PYQREC(2)
     & +QW_w*sin(theta_final)*sin(phi_final)-P(N,2)
                  PYQREC(3)=PYQREC(3)+QW_w*cos(theta_final)-P(N,3)
                  PYQREC(4)=PYQREC(4)+w_soft
ccccccc   Calculating new parton momentum                  
                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)
                 
                  ijoin(ij)=n

              endif
ccccccc   Creating the triplet
ccccccc   Calculating energies of q-qbar-g

              call TripletEnergies(N,w_triplet,E_quark,E_gluon,
     & E_antiq,x_1,x_3)
ccccccc   CALL PY3ENT(IP,KF1,KF2,KF3,PECM,X1,X3)
ccccccc   Purpose: to add three entries to the event record, i.e. either 
ccccccc   a 3-parton system or three separate particles.
ccccccc   IP=normally line number for the first parton/particle, with the other two in IP+1 and IP+2
ccccccc   KF1, KF2, KF3: flavour codes for the three partons/particles.
ccccccc   PECM : (Ecm) the total energy of the system.
ccccccc   X1, X3 : xi = 2Ei/Ecm, i.e. twice the energy fraction taken by the ith parton. Thus
ccccccc   x2 = 2 ° x1 ° x3, and need not be given. Note that not all combinations of xi
ccccccc   are inside the physically allowed region.


           if (w_triplet.gt.0) then
              call PY3ENT(N+1,2,21,-2,w_triplet,x_1,x_3)
           endif

              ij=ij+1
              ijoin(ij)=ip
         
ccccccccccc   Joining partons when we add a hard gluon
ccccccccccc   KS = 1 or 2 change to 3
ccccccccccc   In the cases when we dont add a gluon we only change the status

            if(ij.ge.3) then

                call PYJOIN(ij,ijoin)
            endif

            ij=0
          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            softs gluons                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           else if(iEg.eq.3) then

             if ((QW_w.gt.0).and.((K(ip,1).eq.1).or.(K(ip,1).eq.2))
     & .and.((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then

cccccccc   Calculate energies
                w_hard     = N_const*(QW_L**2)
                w_gluon(1) = w_hard


                if (w_hard.lt.QW_w) then
                  w_soft= (QW_w - w_hard)
                  w_gluon(2)=w_soft
                   
                  if (w_soft.gt.soft_cut) then
                      w_triplet= w_soft-soft_cut
                      w_gluon(3)=w_triplet
                      
                      if((w_soft.lt.0).and.(w_triplet.lt.0)) then
                         w_soft = 0
                        w_triplet = 0
                      endif
                   endif
                endif

cccccccc   Starting Selection of partons of KS=2 
              if(K(ip,1).eq.2) then

                 ij       = ij+1
                 ijoin(ij)= ip
                
                 if(((abs(K(ip,2)).le.5)
     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
     & .and.(((P(ip,4).ge.iet).and.(QW_w.gt.0)
     & .and.(w_soft.le.soft_cut).and.(w_hard.gt.0.0001)
     & .and.(E_p.gt.(2d0*iet))))) then
                
                 call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)

                  n  = n+1
                  ij = ij+1
                 
                  call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

cccccccc    Calculating PYQREC: 4-momentum going back to BeAGLE
cccccccc    Energy of SOFT Gluons
                  PYQREC(1)=PYQREC(1)
     & +QW_w*sin(theta_final)*cos(phi_final)-P(N,1)
                  PYQREC(2)=PYQREC(2)
     & +QW_w*sin(theta_final)*sin(phi_final)-P(N,2)
                  PYQREC(3)=PYQREC(3)+QW_w*cos(theta_final)-P(N,3)
                  PYQREC(4)=PYQREC(4)+w_soft
cccccccc   Calculating new parton momentum 
                  P(ip,1)=P(ip,1)-P(N,1)
                  P(ip,2)=P(ip,2)-P(N,2)
                  P(ip,3)=P(ip,3)-P(N,3)
                
                  ijoin(ij)=n
              endif
cccccccc   Starting selection of partons of KS=1
            else if(K(ip,1).eq.1) then
               
                if((abs(K(ip,2)).le.5)
     & .and.(((P(ip,4).ge.iet).and.(QW_w.gt.0)
     & .and.(w_soft.le.soft_cut).and.(w_hard.gt.0.0001)
     & .and.(E_p.gt.(2d0*iet))))) then
                
                   call GluonEmission(ip,w_gluon(1),theta_final,
     & phi_final)
                   n  = n+1
                   ij = ij+1
                   
                   call PY1ENT(N,21,w_gluon(1),theta_final,phi_final)

cccccccc    Calculating PYQREC: 4-momentum going back to BeAGLE
cccccccc    Energy of SOFT Gluons
                   PYQREC(1)=PYQREC(1)
     & +QW_w*sin(theta_final)*cos(phi_final)-P(N,1)
                   PYQREC(2)=PYQREC(2)
     & +QW_w*sin(theta_final)*sin(phi_final)-P(N,2)
                   PYQREC(3)=PYQREC(3)+QW_w*cos(theta_final)-P(N,3)
                   PYQREC(4)=PYQREC(4)+w_soft
ccccccccccc   Calculating new parton momentum
                   P(ip,1)=P(ip,1)-P(N,1)
                   P(ip,2)=P(ip,2)-P(N,2)
                   P(ip,3)=P(ip,3)-P(N,3)
                  
                   ijoin(ij)= n
               endif

                call TripletEnergies(N,w_triplet,E_quark,E_gluon,
     & E_antiq,x_1,x_3)
               if (w_triplet.gt.0) then
                    call PY3ENT(N+1,2,21,-2,w_triplet,x_1,x_3)
                   
               endif

                if(ij.ge.3) then
                  
c                   call PYJOIN(ij,ijoin)
                endif
                ij       = ij+1
                ijoin(ij)= ip
               

                ij = 0
      
            endif

         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            n hard gluons + soft                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c          else if(iEg.eq.3) then

c            if ((QW_w.gt.0).and.((K(ip,1).eq.1).or.(K(ip,1).eq.2))
c     & .and.((abs(K(ip,2)).le.5).or.(K(ip,2).eq.21))) then


c              n_w=0
c              w_mean = (sqrt(QW_w)-(1/(4*N_const)))**2
c              n_w_gluons=2d0*alphas*cr*sqrt((2d0*QW_wc)/(QW_w))
c              print*,n_w_gluons !,n_w_gluons 
cc              w_mean = N_const*QW_wc 
c              w_sum = w_mean
c              w_high = QW_w

c              do while(w_mean>iet)
c                n_w = n_w + 1
c                print*, n_w, w_mean !, w_high 
c                w_gluon(n_w)=w_mean
c                w_high = w_high-w_mean
cc                w_mean = (sqrt(w_high)-(1/(4*N_const)))**2
c                w_mean = N_const*QW_w
c                w_sum = w_sum + w_mean
c              enddo

c              w_soft = QW_w - w_sum

c           endif

c            if(K(ip,1).eq.2) then

c              ij = ij+1
c              ijoin(ij)=ip
c              i_w = 1

c              if(((abs(K(ip,2)).le.5)
c     & .or.((K(ip,2).eq.21).and.(ij.gt.0)))
c     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

c                  mag = QW_w

c                  do while(i_w.le.n_w)

c                    call GluonEmission(ip,w_gluon(i_w),theta_final,
c     & phi_final)
c                    n=n+1
c                    ij=ij+1

c                    call PY1ENT(N,21,w_gluon(i_w),theta_final,phi_final)

c                    PYQREC(1)=PYQREC(1)
c     & +mag*sin(theta_final)*cos(phi_final)-P(N,1)
c                    PYQREC(2)=PYQREC(2)
c     & +mag*sin(theta_final)*sin(theta_final)-P(N,2)
c                    PYQREC(3)=PYQREC(3)+mag*cos(theta_final)-P(N,3)

c                    P(ip,1)=P(ip,1)-P(N,1)
c                    P(ip,2)=P(ip,2)-P(N,2)
c                    P(ip,3)=P(ip,3)-P(N,3)

c                    ijoin(ij)=n

c                   mag=mag-w_gluon(i_w)
c                    i_w = i_w+1

c                  enddo

c                  mag = 0
c                  PYQREC(4)=PYQREC(4)+w_soft

c              endif

c            else if(K(ip,1).eq.1) then

c              i_w=1

c              if((abs(K(ip,2)).le.5)
c     & .and.(((P(ip,4).ge.PYQ_IET).and.(QW_w.gt.0)))) then

c                  mag = QW_w

c                  do while(i_w.le.n_w)

c                    call GluonEmission(ip,w_gluon(i_w),theta_final,
c     & phi_final)

c                    n=n+1
c                    ij=ij+1
c                    call PY1ENT(N,21,w_gluon(i_w),theta_final,phi_final)

c                    PYQREC(1)=PYQREC(1)
c     & +mag*sin(theta_final)*cos(phi_final)-P(N,1)
c                    PYQREC(2)=PYQREC(2)
c     & +mag*sin(theta_final)*sin(phi_final)-P(N,2)
c                    PYQREC(3)=PYQREC(3)+mag*cos(theta_final)-P(N,3)

c                    P(ip,1)=P(ip,1)-P(N,1)
c                    P(ip,2)=P(ip,2)-P(N,2)
c                    P(ip,3)=P(ip,3)-P(N,3)

c                    ijoin(ij)=n

c                    mag=mag-w_gluon(i_w)

c                    i_w = i_w+1
c                  enddo

c                  PYQREC(4)=PYQREC(4)+w_soft
c                  mag = 0

c              endif

c            ij=ij+1
c            ijoin(ij)=ip

c            if(ij.ge.3) then
c                call PYJOIN(ij,ijoin)
c            endif

c            ij=0

          endif
        
        endif 
        print*,'got here - closing selection'

        ip=ip+1

        print*, 'got here - end of loop'

       
      enddo

      end

      subroutine RotationZ(ip,ppx,ppy,ppz)
       include 'common.f'
       include 'bea_pyqm.inc'
       double precision ppx, ppy, ppz, pp
       double precision theta, phi
       integer ip
      
      
cccccc Calculating the momemtum of a particle where z axis is along gamma*
cccccc theta is the angle between gamma* and z axis in the lab frame
cccccc phi is the angle between gamma* and y axis in the lab frame
cccccc QW_th quenching weights angle
cccccc phe is an azimuthal random number characterizing phi for the final states particles
      phi = atan2(P(ip,2),P(ip,1))
      theta =acos(P(ip,3)/(sqrt(P(ip,1)**2+P(ip,2)**2+P(ip,3)**2)))
      pp=sqrt(P(ip,1)**2+P(ip,2)**2+P(ip,3)**2)
      
      ppx=pp*(cos(phe)*sin(QW_th)*cos(phi)*cos(theta)+sin(phe)*
     & sin(theta)*sin(phi)+cos(QW_th)*cos(phi)*sin(theta))
      ppy=pp*(-cos(phe)*sin(QW_th)*sin(phi)*cos(theta)+sin(phe)*
     & sin(QW_th)*cos(phi)-cos(phe)*sin(phi)*cos(theta))
      ppz=pp*(-cos(phe)*sin(QW_th)*sin(theta)+cos(QW_th)*cos(theta))
      
      print*,'Ppx=',ppx
      print*,'Ppy=',ppy
      print*,'Ppz=',ppz
      print*,'Pp=',pp
      end
      

      subroutine GluonEmission(ip,w_gluon,theta_final,phi_final)

        include 'common.f'
        include 'bea_pyqm.inc'


        double precision w_gluon, theta_final, phi_final
        double precision theta, phi, phe, pgx, pgy, pgz
        integer ip
        print*,'inside GluonEmission:'
c        print*, 'ip = ', ip
c        print*, 'w_gluon = ', w_gluon


        P(ip,4)=P(ip,4)-w_gluon
c        print*,'New Ep ->',P(ip,4)

        theta = acos(P(ip,3)/(sqrt(P(ip,1)**2+P(ip,2)**2+P(ip,3)**2)))
        phi = atan2(P(ip,2),P(ip,1))

c        print*, 'theta parton = ', theta
c        print*, 'phi parton = ', phi

        phe = 4*asin(1.)*ranf(0)

c        print*, 'phe = ', phe
c        print*, 'QW_th = ', QW_th

c  Emitted Gluon momentum
        pgx = (cos(theta)*cos(phi)*sin(QW_th)*cos(phe)-sin(phi)*
     & sin(QW_th)*sin(phe)+sin(theta)*cos(phi)*cos(QW_th))
        pgy = (cos(theta)*sin(phi)*sin(QW_th)*cos(phe)+cos(phi)*
     & sin(QW_th)*sin(phe)+sin(theta)*sin(phi)*cos(QW_th))
        pgz = (-sin(theta)*sin(QW_th)*cos(phe)+cos(theta)*cos(QW_th))


        theta_final = acos((pgz)/(sqrt(pgx**2+pgy**2+pgz**2)))
        phi_final = atan2(pgy, pgx)

      end

      subroutine TripletEnergies(ip,w_triplet,E_quark,E_gluon,E_antiq,
     & x_1,x_3)
       include 'common.f'
       include 'bea_pyqm.inc'

       double precision w_triplet
       double precision E_quark,E_gluon,E_antiq,E_total
       double precision x_1,x_3
       integer ip
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
      double precision qhat,ehat

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
      ehat    = 0.
      qhateff = qhat + ehat
      ChiR    = 0.

      cont=0d+0
      disc=0d+0
      print*,'Initialization of QWComput --> inside subroutine'
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
       print*,'id particle =',id
       QW_L =(2d0* I_QW_wc) / I_QW_R
       QW_R = 2 * density_table(1) * I_QW_wc**2 / I_QW_R / qhateff
       QW_wc = (qhateff/density_table(1)) * I_QW_wc
 

ccccc Convert the units fm -> GeV-1
c      QW_L = QW_L/.1973269
ccccc Convert the units GeV2.fm -> GeV
      QW_wc = QW_wc/.1973269
c      QW_wc_2 = QW_wc_2/.1973269 
ccccc Convert the units GeV2.fm2 -> no unit
      QW_R = QW_R /.1973269**2

      
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

      do i=1,nb_step
        cont(i) = cont(i) / total
      enddo
c      print*,'Picking randomely a quenching from the table'
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

c      goto 56378 
ccccc Calculate the angle probability
      print*, 'Angle Probability'
      if(QW_w .gt. 0) then
        step_QW = 1./nb_step
        yy = E/QW_wc
        xx = QW_w/QW_wc 
       
        total = 0.
        do i=1,nb_step
          ChiR = (step_QW * i)**2 * QW_R
          call qweight(ipart,ChiR,xx,yy,cont(i),disc)
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
      print*,'End of QWComput Subroutine'

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


