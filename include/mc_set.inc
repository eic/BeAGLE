       common /common_mc_set/
     +  mcSet_EneBeam,
     +  mcSet_TarA,
     +  mcSet_TarZ,
     +  mcSet_Q2Min,
     +  mcSet_Q2Max,
     +  mcSet_YMin,
     +  mcSet_YMax,
     +  qedrad,
     +  iModel,
     +  genSet_FStruct,
     +  genSet_R,
     +  mcSet_PTarget,
     +  mcSet_XMin,
     +  mcSet_XMax,
     +  massp,
     +  masse,
     +  ebeam,
     +  pbeam
       save /common_mc_set/

       integer
     +  qedrad,
     +  iModel,
     +  mcSet_TarA,
     +  mcSet_TarZ

       real
     +  mcSet_EneBeam,
     +  mcSet_Q2Min,
     +  mcSet_Q2Max,
     +  mcSet_YMin,
     +  mcSet_YMax,
     +  mcSet_XMin,
     +  mcSet_XMax,
     +  massp,
     +  masse,
     +  ebeam,
     +  pbeam

       character*4
     +  genSet_FStruct,
     +  genSet_R,
     +  mcSet_PTarget

        common /mcevnt/
     +          weight,
     +          genq2, gennu, genx, geny, genw2,
     +          genthe, genphi, geneprim, ebeamEnucl
       save /mcevnt/

        double precision
     +          weight,
     +          genq2, gennu, genx, geny, genw2,        ! vertex kinematics
     +          genthe, genphi, geneprim, ebeamEnucl    ! scattered lepton
