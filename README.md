This is the readme file for BeAGLE

The primary documentation for BeAGLE resides at:
	https://wiki.bnl.gov/eic/index.php/BeAGLE

This directory contains the following files:
***********************************************************************
in src/	
=======================================================================
Multi-routine files:
	dpmjet3.0-5F-new.f	DPMJET code with interface to FLUKA and PYTHIA
	phojet1.12-35c4.f	PHOJET CODE version 1.12
	dpm_pythia.f		interface from dpmjet to PYTHIA
	dpm_rapgap.f		interface from dpmjet to RAPGAP
Single-routine files:
	anear.f			lookup table for nearest A with valid nPDF
	azmass.f		lookup table for M(A,Z)	in amu
	dcalc.f			Calculates "d", effective distance in nucleus
	deutfix.f		Ad hoc handling of e+D 4-momentum conservation
	ishadr.f		Logical function to detect hadrons from pid 
	kickit.f		Elastic intrinsic parton kt kick/nucleon recoil
	nbary.f	 		Returns 3* baryon # for pid (i.e. net quark #)
	nrbinom.f		Roll a # according to a binomial distribution
	pfshift.f		Apply missing pF to Pythia/Rapgap subevent
	reinit.f		Switch DPMJET/Pythia/Rapgap between e+n & e+p
	rmkick.f 		Rolls an intrinsic kt similar to PYREMN 
	shdmap.f		Function dipole sigma vs. shadowing R_A
	any other codes necessary for pythia or (radgen in the future)
Missing/obsolete/temporary routines:
	ktkick.f ??		Old name for rmkick.f??
	rgsteer.f		Allows us to just run vanilla RAPGAP in BeAGLE
	ulangl.f		Duplicate of pyangl.f. Just use pyangl.f.

in include/
=======================================================================
	header files 

in PYTHIA-6.4.28/
=======================================================================
	the main source code of pythia

in PyQM/
=======================================================================
	the main source code of quenching effect for pythia from 
	Alberto Accardi, Raphael Dupre, with Mathieu Ehrhart

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
	In order to use this version of hybrid model, FLUKA, LHAPDF and RAPGAP
        must be installed to take care of the nuclear breakup, nuclear PDF, and
        optional improved description of diffraction respectively.

	For FLUKA:
	~~~~~~~~~~~~~~~~~~~~~
	- download and install FLUKA (www.fluka.org), 
	- define a variable FLUPRO pointing to the FLUKA-directory 
	- link libflukahp.a when creating the DPMJET-executable (taken care 
	  of by the Makefile)
	- copy nuclear.bin either into the directory where you run DPMJET or 
          create a link to it 

	For LHAPDF:
	~~~~~~~~~~~~~~~~~~~~~
	- download and install LHAPDF, this code is tested with LHAPDF-5.8.6
	- link to the libLHAPDF.a when creating the executable file
	- define an environment variable $LHAPATH pointing to the LHAPDF sets path

	The linking step of the FLUKA and LHAPDF can be done in the Makefile.
	You need to specify the path to your liblukaph.a and libLHAPDF.a there.

	For RAPGAP, we have an awkward setup at the moment. We build RAPGAP 
        elsewhere, and import a modified version of some of the libraries
        to ...BeAGLE/RAPGAP_3.302/lib. 	The include files are already contained
        in BeAGLE/RAPGAP_3.302/include.

	For RAPGAP:
	~~~~~~~~~~~~~~~~~~~~~
	- Download RAPGAP 
          (https://rapgap.hepforge.org/downloads/?f=rapgap-3.302.tar.gz)
        - For pythia, you can download pythia6428.tgz from RAPGAP or just point to ours.
        - Either way, copy the Pythia 6.4.28 pydata.f into .../rapgap-3.302/src/rapgap/
          The one in there is from 6.4.13 for some reason
        - Install/build RAPGAP (ARIADNE & HZTOOL are optional)
	- Copy libar4.a, libbases.a, librapgap33.a to .../BeAGLE/RAPGAP-3.302/lib/
	- Issue the command: ar d librapgap33.a sfecfe.o
	  NOTE: This is to avoid a linking clash. SFECFE is copied from Fluka which we
          are already linking.
        - The Makefile should link correctly to all of these libraries.

Compile:
===================
	To compile this code, you need to use the Makefile script in the current directory.
	After you have set up your environment, use the make utility 'make all'.

	Note: The code is now 64-bit. 

	To clean all the .o file, exe file and core dump, use clean utility 'make clean'


Running:
===================
	The code requires an environment variable $BEAGLESYS to point to the
        top level BeAGLE directory where the executable and a collection of
        input data files reside.

	To run this code, use the command './BeAGLE < myinputfile.inp > myinputfile.log'
	Your input (DPMJET control) file should also reference a Pythia parameter file
	which needs to be in the same directory.

	More details about options you can change in the input file can be found at:
	https://wiki.bnl.gov/eic/index.php/BeAGLE



