C...Stuff for PYTHIA 6.4
      double precision PARP(200),PARI(200)
      integer MSTP(200),MSTI(200)
      COMMON/PYPARS/MSTP,PARP,MSTI,PARI
      double precision PARU(200),PARJ(200)
      integer MSTU(200),MSTJ(200)
      COMMON/PYDAT1/MSTU,PARU,MSTJ,PARJ
      double precision P(4000,5),V(4000,5)
      integer N,NPAD,K(4000,5)
      COMMON/PYJETS/N,NPAD,K,P,V
      integer MSEL,MSELPD,MSUB(500),KFIN(2,-40:40)
      double precision CKIN(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB,KFIN,CKIN
      integer NGENPD,NGEN(0:500,3)
      double precision XSEC(0:500,3)
      COMMON/PYINT5/NGENPD,NGEN,XSEC

ccccc Table for density distribution
      double precision density_table(2000),step_size_dens
      double precision quantity_table(2000),init_dens
      common/Density/density_table,step_size_dens,quantity_table,
     &               init_dens

ccccc Quenching Weight variables
      integer QW_nb
      double precision QW_wc,QW_R,QW_L,QW_w,QW_qhat,QW_chi,QW_th
      common/QuenWei/QW_nb,QW_wc,QW_R,QW_L,QW_w,QW_qhat,QW_chi,QW_th

ccccc Interaction position
      double precision x_inter,y_inter,z_inter
      double precision pos_radius,pos_theta,pos_phi
      common/InteracPos/x_inter,y_inter,z_inter,
     &      pos_radius,pos_theta,pos_phi

ccccc AA routine variables
      double precision alphas
            integer iqw,scor,ncor,sfthrd,irw
      common/qw/alphas,iqw,scor,ncor,sfthrd,irw

ccccc Miscellanous
c random number generator from CERNLIB
      real ranf

c Flags and config values
c     integer iTg,iFM,iDens,iQuenching,iSim,nevent 
c     integer nkin,nucleon,specId,iColl,iqg,iEg,iPtF
c     double precision rFM,E0,qhat,ehat,EColl,FMlimit,SupFac
c     common/flags/iTg,iFM,iDens,iQuenching,iSim,nevent,
c    &             nkin,nucleon,specId,rFM,E0,qhat,ehat,iColl,EColl,
c    &             FMlimit,iqg,iEg,SupFac,iPtF

      double precision pi
      data pi/3.1415926535/

