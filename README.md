This is the readme file for BeAGLE

This directory contains the following files:
***********************************************************************
in src/	
=======================================================================
	dpmjet3.0-5F-new.f	DPMJET code with interface to FLUKA and PYTHIA
	phojet1.12-35c4.f	PHOJET CODE version 1.12
	dpm_pythia.f		interface from dpmjet to PYTHIA
	azmass.f		lookup table for M(A,Z)	in amu
	kickit.f		Elastic intrinsic parton kt kick/nucleon recoil
	ktkick.f		Rolls an intrinsic kt similar to PYREMN
	shdmap.f		Function dipole sigma vs. shadowing R_A
	ulangl.f		Azimuthal angle given x,y (from sPHENIX)
	any other codes necessary for pythia or (radgen in the future)

in include/
=======================================================================
	header file used in pythia

in PYTHIA-6.4.28/
=======================================================================
	the main source code of pythia

in PyQM/
=======================================================================
	the main source code of quenching effect for pythia from 
	Alberto Accardi, Raphael Dupre

in top dir
=======================================================================
	user-eA3.0-5.f    main program to call these subroutines

	Makefile				makefile to install DPMJET

	nuclear.bin			Needed for fluka evaporation

	ep.inp/eAu.inp		input file for dpmjet
	
	read_ep/read_eA		auxiliary input parameter to initialize pthia




Usage:
***********************************************************************
Preconfiguration:
===================
	In order to use this version of hybrid model, FLUKA and LHAPDF must
	be installed to take care of the nuclear evaporation and nuclear PDF.
	For FLUKA:
	~~~~~~~~~~~~~~~~~~~~~
	- download and install FLUKA (www.fluka.org), 
	- define a variable FLUPRO pointing to the FLUKA-directory 
	- link libflukahp.a when creating the DPMJET-executable (taken care of by the Makefile)
	- copy nuclear.bin either into the directory where you run DPMJET or create a link to it 

	For LHAPDF:
	~~~~~~~~~~~~~~~~~~~~~
	- download and install LHAPDF, this code is tested with LHAPDF-5.8.6
	- link to the libLHAPDF.a when creating the executable file
	- define an environment variable $LHAPATH pointing to the LHAPDF sets path

	The linking step of the FLUKA and LHAPDF can be done in the Makefile.
	You need to specify the path to your liblukaph.a and libLHAPDF.a there.


Compile:
===================
	To compile this code, you need to use the Makefile script in the current directory.
	After you have set up your environment, use the make utility 'make all'.

	Tips: It is suggested that all the programs you will use are in the 32-bit version,
	since the link to fluka requires the compiling to be done compatiable with a 32-bit
	program.

	To clean all the .o file, exe file and core dump, use clean utility 'make clean'


Running:
===================
	To run this code, for ep use the command './dpmjetHybrid < ep.inp > ep.log'
	for eA, just replace the ep.inp by eA.inp, './dpmjetHybrid < eA.inp > eA.log'
	Remember that corresponding 'read_ep' and 'read_eA' parameter list should be put
	in the same directory, and the beam particle information needs to be set
	by ep.inp/eA.inp while output file name (when using pythia model) is given in
	read_ep/read_eA for ep or eA running respectively.

	More details about options you can change in the input file and the output file
	structure can be found in the "Doumentation" in this directory.



Created by liang zheng <liangzhphy@gmail.com> on Feb,29 2012
