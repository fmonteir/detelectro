//
//  main_equal_time.cpp
//
//
//  Created by Francisco Brito on 18/10/2018.
//
//  This program simulates the Hubbard model for an arbitrary geometry lattice
//  using auxiliary field (or determinant) Quantum Monte Carlo: in particular,
//  the BSS algorithm. This is version is optimized to compute observables that
//  require measuring equal-time Green's functions
//  The used notation is based on the lecture notes "Numerical Methods for
//  Quantum Monte Carlo Simulations of the Hubbard Model by Zhaojun Bai,
//  Wenbin Chen, Richard Scalettar, and Ichitaro Yamazaki (2009)
//

//  Total number of "sites" (actual spatial sites plus orbitals)
#ifndef NSITES
#define NSITES 2
#endif

//  Inverse Trotter error
#ifndef DT_INV
#define DT_INV 16
#endif

//  Inverse temperature
#ifndef BETA
#define BETA 2
#endif

//  How often to calculate Green's functions afresh
//  (measured in number of imaginary-time slices)
#ifndef GREEN_AFRESH_FREQ
#define GREEN_AFRESH_FREQ 4
#endif

//  Output information about the progress of the run
#ifndef VERBOSE
#define VERBOSE 0
#endif

//  Includes

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <array>
#include <algorithm>
#include <functional>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "unsupported/Eigen/KroneckerProduct"
#include "matrixgen.h"
#include "green.h"

int main(int argc, char **argv)
{
    if ( argc != 9) //  t, U, mu, geom, Ny,
                    //  # sweeps, # warm-up sweeps, # auto-correlation sweeps
    {
        std::cout << "Not enough arguments given to simulation."
        << std::endl << std::endl;

        return -1;
    }

    double t = atof(argv[1]);  //  tight binding parameter
    double U = atof(argv[2]);  //  on-site interaction
    double mu = atof(argv[3]);  //  chemical potential
    int geom = atof(argv[4]);  //   geometry (see makefile)
    int Ny = atof(argv[5]);    //   parameter defining the width of the sample
    int totalMCSweeps = atof(argv[6]);  //  number of sweeps
    int W = atof(argv[7]);  //  number of warm-up sweeps
    int A = atof(argv[8]);  //  number of auto-correlation sweeps

    double dt = 1. / DT_INV;  //  Trotter error. The error scales as dt^2
    const int L = (int)(BETA * DT_INV);  //  # slices
    const int Lbda = (int)(L / GREEN_AFRESH_FREQ);
    //  Lbda = # intervals in which the product of B's is divided to stabilize.
    double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;
    //  HS transformation parameter (to order dtau^2)

    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    std::mt19937 gen;  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    //  Ensure that we generate uncorrelated random numbers
    std::random_device r;
    std::array<int, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    gen.seed(seq);
    double decisionMaker; int totalMCSteps = totalMCSweeps * NSITES * L;


    // -- INITIALIZATION ---


    //  HOPPING MATRIX
    Geometry< NSITES > K;
    //  1D CHAIN PBC
    if (geom == 1)
    {
        K.oneDimensionalChainPBC(t, dt, mu);
    }
    //  1D CHAIN OBC
    if (geom == 2)
    {
        K.oneDimensionalChainOBC(t, dt, mu);
    }
    //  SQUARE LATTICE PBC
    if (geom == 3)
    {
        K.twoDimensionalRectanglePBC(sqrt(NSITES), t, dt, mu);
    }
    //  SQUARE LATTICE OBC
    if (geom == 4)
    {
        K.twoDimensionalRectangleOBC(sqrt(NSITES), t, dt, mu);
    }
    //  RECTANGULAR LATTICE PBC
    if (geom == 5)
    {
        K.twoDimensionalRectanglePBC(Ny, t, dt, mu);
    }
    //  RECTANGULAR LATTICE OBC
    if (geom == 6)
    {
        K.twoDimensionalRectangleOBC(Ny, t, dt, mu);
    }
    // //  TRIANGULAR LATTICE PBC
    // if (geom == 7)
    // {
    //     K.twoDimensionalTrianglePBC(t, dt);
    // }
    // //  TRIANGULAR LATTICE OBC
    // if (geom == 8)
    // {
    //     K.triangleNanoribbon(t, dt);
    // }
    //  HONEYCOMB LATTICE PBC
    if (geom == 9)
    {
        K.hcPBC();
        K.computeExponential(t, dt);
    }
    //  NANORIBBON
    if (geom == 10)
    {
        K.hcNanoribbon(Ny);
        K.computeExponential(t, dt);
    }
    //  STRAINED NANORIBBON
    if (geom == 11)
    {
        double Delta = 0.3;
        K.hcStrainedNanoribbon(Ny, Delta);
        K.computeExponential(t, dt);
    }
    //  3 ORBITAL TIGHT BIDING MODEL ON THE M-ATOM TRIANGULAR LATTICE
    //  MODEL OF A TMD NANORIBBON (SEE Liu2013)
    //  MoS2
    if (geom == 12)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)."
            << std::endl << std::endl;

            return -1;
        }
        double params[] = {1.046, 2.104, -0.184, 0.401, 0.507, 0.218, 0.338, 0.057};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();
        K.computeExponential(t, dt);
    }
    if (geom == 13)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.046, 2.104, -0.184, 0.401, 0.507, 0.218, 0.338, 0.057};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);
        K.computeExponential(t, dt);
    }

    if (geom == 14)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.046, 2.104, -0.184, 0.401, 0.507, 0.218, 0.338, 0.057};
        K.setParamsThreeOrbitalTB(params);
        double Delta = 0.4;
        K.tmdNanoribbonStrained(Ny, Delta);
        K.computeExponential(t, dt);
    }
    //  WS2
    if (geom == 15)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)."
            << std::endl << std::endl;

            return -1;
        }
        double params[] = {1.130, 2.275, -0.206, 0.567, 0.536, 0.286, 0.384, -0.061};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();
        K.computeExponential(t, dt);
    }
    if (geom == 16)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.130, 2.275, -0.206, 0.567, 0.536, 0.286, 0.384, -0.061};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);
        K.computeExponential(t, dt);
    }

    if (geom == 17)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.130, 2.275, -0.206, 0.567, 0.536, 0.286, 0.384, -0.061};
        K.setParamsThreeOrbitalTB(params);
        double Delta = 0.4;
        K.tmdNanoribbonStrained(Ny, Delta);
        K.computeExponential(t, dt);
    }
    //  MoSe2
    if (geom == 18)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)."
            << std::endl << std::endl;

            return -1;
        }
        double params[] = {0.919, 2.065, -0.188, 0.317, 0.456, 0.211, 0.290, 0.130};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();
        K.computeExponential(t, dt);
    }
    if (geom == 19)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.919, 2.065, -0.188, 0.317, 0.456, 0.211, 0.290, 0.130};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);
        K.computeExponential(t, dt);
    }

    if (geom == 20)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.919, 2.065, -0.188, 0.317, 0.456, 0.211, 0.290, 0.130};
        K.setParamsThreeOrbitalTB(params);
        double Delta = 0.4;
        K.tmdNanoribbonStrained(Ny, Delta);
        K.computeExponential(t, dt);
    }
    //  WSe2
    if (geom == 21)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)."
            << std::endl << std::endl;

            return -1;
        }
        double params[] = {0.943, 2.179, -0.207, 0.457, 0.486, 0.263, 0.329, 0.034};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();
        K.computeExponential(t, dt);
    }
    if (geom == 22)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.943, 2.179, -0.207, 0.457, 0.486, 0.263, 0.329, 0.034};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);
        K.computeExponential(t, dt);
    }

    if (geom == 23)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.943, 2.179, -0.207, 0.457, 0.486, 0.263, 0.329, 0.034};
        K.setParamsThreeOrbitalTB(params);
        double Delta = 0.4;
        K.tmdNanoribbonStrained(Ny, Delta);
        K.computeExponential(t, dt);
    }
    //  MoTe2
    if (geom == 24)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)."
            << std::endl << std::endl;

            return -1;
        }
        double params[] = {0.605, 1.972, -0.169, 0.228, 0.390, 0.207, 0.239, 0.252};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();
        K.computeExponential(t, dt);
    }
    if (geom == 25)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.605, 1.972, -0.169, 0.228, 0.390, 0.207, 0.239, 0.252};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);
        K.computeExponential(t, dt);
    }

    if (geom == 26)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.605, 1.972, -0.169, 0.228, 0.390, 0.207, 0.239, 0.252};
        K.setParamsThreeOrbitalTB(params);
        double Delta = 0.4;
        K.tmdNanoribbonStrained(Ny, Delta);
        K.computeExponential(t, dt);
    }
    //  WTe2
    if (geom == 27)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)."
            << std::endl << std::endl;

            return -1;
        }
        double params[] = {0.606, 2.102, -0.175, 0.342, 0.410, 0.233, 0.270, 0.190};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();
        K.computeExponential(t, dt);
    }
    if (geom == 28)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.606, 2.102, -0.175, 0.342, 0.410, 0.233, 0.270, 0.190};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);
        K.computeExponential(t, dt);
    }

    if (geom == 29)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {0.606, 2.102, -0.175, 0.342, 0.410, 0.233, 0.270, 0.190};
        K.setParamsThreeOrbitalTB(params);
        double Delta = 0.4;
        K.tmdNanoribbonStrained(Ny, Delta);
        K.computeExponential(t, dt);
    }

    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Configuration< L , NSITES > * h = new Configuration< L , NSITES >;
    h->genHsMatrix();

    //  GENERATE THE B-MATRICES.
    OneParticlePropagators< NSITES, L > * Bup =
      new OneParticlePropagators< NSITES, L >;
    OneParticlePropagators< NSITES, L > * Bdown=
      new OneParticlePropagators< NSITES, L >;
    Bup->fillMatrices( true, nu, h->matrix(), K.BpreFactor(dt, mu) );
    Bdown->fillMatrices( false, nu, h->matrix(), K.BpreFactor(dt, mu) );

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green< NSITES, L, Lbda> * Gup = new Green< NSITES, L, Lbda>;
    Green< NSITES, L, Lbda> * Gdown = new Green< NSITES, L, Lbda>;
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
    //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    // Gup->computeGreenNaive(Bup->list(), L - 1);
    // Gdown->computeGreenNaive(Bdown->list(), L - 1);
    Gup->storeVDU( Bup->list() ); Gdown->storeVDU( Bdown->list() );
    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[W * L];
    double LOGweight = 0.;

    double electronDensities = 0;
    double doubleOcs = 0;
    double energies = 0;
    Eigen::MatrixXd magCorrZZs =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();

    // double sign = std::copysign(1, Gup->matrix().determinant()
    //   * Gdown->matrix().determinant() );
    double sign = 1;
    double meanSign = 0;

    double electronDensity;
    double doubleOc;
    double energy;
    Eigen::MatrixXd magCorrZZ = Eigen::Matrix<double, NSITES, NSITES>::Zero();

    double nEl = 0;
    double nUp_nDw = 0;
    double Hkin = 0;
    Eigen::MatrixXd SiSjZ =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE,
    //  THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0; int i = 0; int latticeSweepUntilAfresh = 0; int sweep = 0; int step;


    // --- MC LOOP ---

    if (VERBOSE == 1)
    {
      std::cout << "\nMC loop started. Progress:\n";
    }

    for (step = 0; step < totalMCSteps; step++)
    {
        if (VERBOSE == 1)
        {
             // DISPLAY PROGRESS OF THE RUN.
            if ( (step + 1)  % (totalMCSteps/8) == 0 )
            {
                std::cout << (step + 1) * 1. / totalMCSteps * 100 << " %"
                  << std::endl << std::endl;
                std::cout << "Average Sign: " << meanSign
                  << std::endl << std::endl;
                std::cout << "Log Weight: " << LOGweight
                  << std::endl << std::endl;
            }
        }
        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h->get(l, i) * nu ) - 1 );
        alphaDown = ( exp( 2 * h->get(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup->get(i, i) ) );
        dDown = ( 1 + alphaDown  * ( 1 - Gdown->get(i, i) ) );
        //  SAMPLING: METROPOLIS OR HEAT BATH
        accRatio = fabs( dUp * dDown );
        // accRatio = fabs( dUp * dDown / ( 1 + dUp * dDown ) );

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= accRatio )
        {
            //  KEEP TRACK OF WEIGHT
            LOGweight += log( fabs( dUp ) ) + log( fabs ( dDown ) );
            sign *= std::copysign(1, dUp * dDown );
            //  FLIP A SPIN
            h->flip(l, i);
            //  UPDATE Bs
            Bup->update(l, i, alphaUp); Bdown->update(l, i, alphaDown);
            //  RANK-ONE UPDATE -> O(N^2)
            Gup->update(alphaUp, dUp, i); Gdown->update(alphaDown, dDown, i);
        }


        // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. ---


        if (i < NSITES - 1)
        {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE.
            i += 1;
        }
        else
        {
            //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH.
            latticeSweepUntilAfresh += 1;
            //  --- MEASUREMENTS ---
            if ( sweep < W )
            {
              //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
              weights[sweep * L + l] = LOGweight;
            }
            //  STORE ELECTRON DENSITY, DOUBLE OCCUPANCY, AND SPIN-SPIN CORRELATIONS.
            electronDensity = 0.; doubleOc = 0.; energy = 0.;
            for (int x = 0; x < NSITES; x++)
            {
                electronDensity -= ( Gup->get(x, x) + Gdown->get(x, x) );
                doubleOc += - Gup->get(x, x) - Gdown->get(x, x) + Gup->get(x, x)
                * Gdown->get(x, x);
                magCorrZZ(x, x) = ( Gup->get(x, x) + Gdown->get(x, x) )
                  - 2 * Gup->get(x, x) * Gdown->get(x, x);
                energy += 2 * ( Gup->get(x, x) + Gdown->get(x, x) ) * (t * K.get(x, x) + mu);
                for (int y = 0; y < x; y++)
                {
                    magCorrZZ(x, y) =
                      - ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                      - ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->get(y, y) )
                      + ( 1 - Gup->get(x, x) ) * ( 1 - Gup->get(y, y) )
                      + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                      - Gup->get(y, x) * Gup->get(x, y)
                      - Gdown->get(y, x) * Gdown->get(x, y);
		                magCorrZZ(y, x) = magCorrZZ(x, y);
                    energy += ( Gup->get(x, y) + Gdown->get(x, y)
                      + Gup->get(y, x) + Gdown->get(y, x) ) * t * K.get(x, y);
                }
            }
            electronDensity /= NSITES; electronDensity += 2;
            doubleOc /= NSITES; doubleOc += 1;
            energy /= NSITES;

            electronDensities +=
              ( electronDensity * sign - electronDensities ) / ( l + 1 ) ;
            doubleOcs +=
              ( doubleOc * sign - doubleOcs ) / ( l + 1 ) ;
            magCorrZZs +=
              (magCorrZZ * sign - magCorrZZs ) / ( l + 1 );
            energies +=
              ( energy * sign - energies ) / ( l + 1 ) ;

            //  DEAL WITH THE GREEN'S FUNCTIONS.


            //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == GREEN_AFRESH_FREQ)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
                // Gup->computeGreenNaive(Bup->list(), l);
                // Gdown->computeGreenNaive(Bdown->list(), l);
                //  Uncomment to compute the product in the stabilized,
                //  but slightly inefficient way
                // Gup->computeStableGreenNaiveR(Bup->list(), l);
                // Gdown->computeStableGreenNaiveR(Bdown->list(), l);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    Gup->storeUDV(Bup->list(), l, GREEN_AFRESH_FREQ);
                    Gdown->storeUDV(Bdown->list(), l, GREEN_AFRESH_FREQ);
                    //  This is the standard way described in
                    //  "Stable simulations of models of interacting electrons"
                    Gup->computeStableGreen(l, GREEN_AFRESH_FREQ);
                    Gdown->computeStableGreen(l, GREEN_AFRESH_FREQ);
                }
                else
                {
                    Gup->storeVDU( Bup->list() ); Gdown->storeVDU( Bdown->list() );
                    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();
                }
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                Gup->wrap( Bup->matrix(l) ); Gdown->wrap( Bdown->matrix(l) );
            }
            if (l < L - 1)
            {
                l += 1; i = 0;
            }
            else
            {
                if ( (sweep >= W) )
                {
                    if ( sweep % A == 0 )
                    {
                      meanSign += ( sign - meanSign ) / ( ( sweep - W ) / A + 1 );
                      nEl += ( electronDensities - nEl )
                       / ( (sweep - W)/A + 1 ) ;
                      nUp_nDw += ( doubleOcs - nUp_nDw )
                       / ( (sweep - W)/A + 1 ) ;
                      SiSjZ += ( magCorrZZs - SiSjZ )
                       / ( (sweep - W)/A + 1 ) ;
                      Hkin += ( energies - Hkin )
                       / ( (sweep - W)/A + 1 ) ;

                    }
                    electronDensities = 0.; doubleOcs = 0.;
                    magCorrZZs = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                    energies = 0.;

                }
                //  MOVE SWEEP COUNTER
                sweep += 1;
                l = 0; i = 0;

            }
        }
    }   //  END OF MC LOOP.

    write(meanSign, sweep, W, A, nEl, nUp_nDw,
      Hkin, U, NSITES, dt, BETA, L, t, mu, GREEN_AFRESH_FREQ, Lbda,
      geom, Ny, weights, SiSjZ);

    delete[] weights;
    delete Gup; delete Gdown; delete h; delete Bup; delete Bdown;

    return 0;
}
