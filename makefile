# DETELECTRO - Determinant Quantum Monte Carlo for Hubbard-like models
# Created by Francisco Brito October 23, 2018

#	The program is prepared to handle a variety of geometries (listed below).

#		Input the number corresponding to the desired geometry
#		(1)		1D Periodic Chain
#		(2) 	1D Open Chain
#		(3) 	2D Periodic Square Lattice
#		(4) 	2D Open Square Lattice
#		(5) 	2D Periodic Rectangular Lattice
#		(6) 	2D Open Rectangular Lattice
#		(7) 	2D Periodic Triangular Lattice
#		(8) 	2D Nanoribbon Triangular Lattice
#		(9) 	2D Periodic Honeycomb Lattice
#		(10)	2D Honeycomb Nanoribbon
#		(11)	2D Honeycomb Strained Nanoribbon
#		(12)	2D Minimal model of a periodic MoS2 sample (Liu2013)
#		(13)	2D Minimal model of a MoS2 nanoribbon (Liu2013)
#		(14)	2D Minimal model of a MoS2 nanoribbon with Strain (Liu2013)
#		(15)	2D Minimal model of a periodic WS2 sample (Liu2013)
#		(16)	2D Minimal model of a WS2 nanoribbon (Liu2013)
#		(17)	2D Minimal model of a WS2 nanoribbon with Strain (Liu2013)
#		(18)	2D Minimal model of a periodic MoSe2 sample (Liu2013)
#		(19)	2D Minimal model of a MoSe2 nanoribbon (Liu2013)
#		(20)	2D Minimal model of a MoSe2 nanoribbon with Strain (Liu2013)
#		(21)	2D Minimal model of a periodic WSe2 sample (Liu2013)
#		(22)	2D Minimal model of a WSe2 nanoribbon (Liu2013)
#		(23)	2D Minimal model of a WSe2 nanoribbon with Strain (Liu2013)
#		(24)	2D Minimal model of a periodic MoTe2 sample (Liu2013)
#		(25)	2D Minimal model of a MoTe2 nanoribbon (Liu2013)
#		(26)	2D Minimal model of a MoTe2 nanoribbon with Strain (Liu2013)
#		(27)	2D Minimal model of a periodic WTe2 sample (Liu2013)
#		(28)	2D Minimal model of a WTe2 nanoribbon (Liu2013)
#		(29)	2D Minimal model of a WTe2 nanoribbon with Strain (Liu2013)

#	DEFAULT PARAMETERS

# Size of the system (orbital + real space)
nsites=2
#width of the ribbon
ny=4
# Trotter error
dt_inv=16
# Inverse temperature
beta=2
# Frequency of recomputing G
green_afresh_freq=4
# Toggle prints
verbose=0
# Choose source file
source=main_equal_time

# Set parameters of the simulation here.
CXX = g++ -DNSITES=$(nsites) -DDT_INV=$(dt_inv) -DBETA=$(beta)\
 -DGREEN_AFRESH_FREQ=$(green_afresh_freq) -DVERBOSE=$(verbose) -DNY=$(ny)

include_dir=./includes

CXXFLAGS = -Wall -g -O3 -std=c++11 -I$(include_dir)

simulation: src/$(source).o
ifeq ($(verbose),1)
	@echo ""
	@echo "		DETELECTRO - Determinant Quantum Monte Carlo for Hubbard-like models"
	@echo ""
	@echo "			Created by Francisco Brito (2018)"
	@echo ""
	@echo "The code has compiled successfully. To change the number of sites, \
	inverse Trotter error, inverse temperature, or the frequency of recomputing \
	the Green's functions, type:"
	@echo ""
	@echo "make clean"
	@echo ""
	@echo "make nsites=<Number of sites> dt_inv=<Inverse Trotter Error> \
	beta=<Inverse Temperature> green_afresh_freq=<Frequency of Recomputing G>"
	@echo ""
	@echo "To run a simulation, simply type ./simulation followed by its arguments:"
	@echo ""
	@echo "./simulation <t> <U> <mu> <geom>\
	 <Total Number of Sweeps (Space-Time)>\
	 <Number of Warm-up Sweeps (Space-Time)> \
	 <Number of Auto-correlation Sweeps (Space-Time)>"
	@echo ""
	@echo "where t is the (tight-binding) normalization, U is the on-site interaction, \
mu is the chemical potential (in particle-hole symmetric form, mu -> mu + U /2 )"
	@echo "The program is prepared to handle a variety of geometries (listed below).\
	 Input the number corresponding to the desired geometry:"
	@echo ""
	@echo "(1)		1D Periodic Chain"
	@echo ""
	@echo "(2) 		1D Open Chain"
	@echo ""
	@echo "(3) 		2D Periodic Square Lattice"
	@echo ""
	@echo "(4) 		2D Open Square Lattice"
	@echo ""
	@echo "(5) 		2D Periodic Rectangular Lattice"
	@echo ""
	@echo "(6) 		2D Open Rectangular Lattice"
	@echo ""
	@echo "(7) 		2D Periodic Triangular Lattice"
	@echo ""
	@echo "(8) 		2D Nanoribbon Triangular Lattice"
	@echo ""
	@echo "(9) 		2D Periodic Honeycomb Lattice"
	@echo ""
	@echo "(10)		2D Honeycomb Nanoribbon"
	@echo ""
	@echo "(11)		2D Honeycomb Nanoribbon w/ Strain"
	@echo ""
	@echo "(12)		2D Minimal model of a periodic MoS2 sample\
	 (Liu et al., Phys Rev B 88, 085433, 2013 )\
	 nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites."
	@echo ""
	@echo "(13)		2D Minimal model of a MoS2 nanoribbon\
	 (Liu et al., Phys Rev B 88, 085433, 2013 )\
	 nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites."
	@echo ""
	@echo "(14)		2D Minimal model of a MoS2 nanoribbon with strain\
	 (Liu et al., Phys Rev B 88, 085433, 2013 )\
	 nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites."
	@echo ""
	@echo "(15-17)		Same with WS2."
	@echo ""
	@echo "(18-20)		Same with MoSe2."
	@echo ""
	@echo "(21-23)		Same with WSe2."
	@echo ""
	@echo "(24-26)		Same with MoTe2."
	@echo ""
	@echo "(27-29)		Same with WTe2."
	@echo ""
	@echo "The Ny parameter is only meaningful\
	 for geometry options 5, 6, 8, 10-29."
	@echo ""
endif

	$(CXX) $(CXXFLAGS) -o simulation src/$(source).o


src/$(source).o: src/$(source).cpp $(include_dir)/matrixgen.h\
	 $(include_dir)/green.h $(include_dir)/QDT.h $(include_dir)/TDQ.h

clean:
	rm -f simulation src/*.o
	rm -f ./temp-data/*
