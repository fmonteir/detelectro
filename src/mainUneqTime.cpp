//
//  mainUneqTime.cpp
//
//
//  Created by Francisco Brito on 17/09/2018.
//
//  This program simulates the Hubbard model for an arbitrary geometry lattice
//  using auxiliary field (or determinant) Quantum Monte Carlo: in particular, the BSS algorithm.
//  The used notation is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki (2009)
//

//  DEFAULT SIMULATION PARAMETERS FOR MINIMAL EXAMPLE.
//  Compare with J. E. Hirsch - Phys Rev B 28 7, 1983
//  For U = 4, we get <n_up n_dw > -> 0.1384 (exact)

#ifndef NSITES
#define NSITES 2 //  # sites
#endif

#ifndef DT_INV
#define DT_INV 16 //  Inverse Trotter error
#endif

#ifndef BETA
#define BETA 2  //  inverse temperature
#endif

#ifndef GREEN_AFRESH_FREQ
#define GREEN_AFRESH_FREQ 4  //   how often to calculate Green's functions afresh (in # im-time slices)
#endif

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "unsupported/Eigen/KroneckerProduct"
#include "matrixgen.h"
#include "green.h"


int main(int argc, char **argv)
{
    if ( argc != 9) //  U, mu, # sweeps, # warm-up sweeps, geometry
    {
        return -1;
    }

    Eigen::IOFormat CleanFmt(10, 0, ", ", "\n", "", "");

    double t = atof(argv[1]);  //  tight binding parameter
    double U = atof(argv[2]);  //  on-site interaction
    double mu = atof(argv[3]);  //  chemical potential
    int geom = atof(argv[4]);  //   geometry (for example, 1 stands for a 1D chain with PBCs - see makefile)
    int Ny = atof(argv[5]);    //   geometry parameter to define width of the simulated sample
    int totalMCSweeps = atof(argv[6]);  //  number of sweeps
    int W = atof(argv[7]);  //  number of warm-up sweeps
    int A = atof(argv[8]);  //  number of auto-correlation sweeps

    double dt = 1. / DT_INV;  //  Trotter error, or time subinterval width. error scales as dt^2
    const int L = (int)(BETA * DT_INV);  //  # slices
    //  Lbda = # intervals in which the product of B's is divided to stabilize.
    const int Lbda = L / GREEN_AFRESH_FREQ;
    //  HS transformation parameter (to order dtau^2)
    double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;

    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 1;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
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
    //  TRIANGULAR LATTICE PBC
    if (geom == 6)
    {

    }
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
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.046, 2.104, 0.401, 0.507, 0.218, 0.338, 0.057};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();

        // FOR DEBUGGING: TEST MATRIX CREATED BY THE PROGRAM IN OTHER CODES
        // std::ofstream TMDhopping("temp-data/tmd-hopping.csv");
        // if (TMDhopping.is_open())
        // {
        //     TMDhopping << K.matrix() << '\n';
        // }
        // TMDhopping.close();
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

        // FOR DEBUGGING: TEST MATRIX CREATED BY THE PROGRAM IN OTHER CODES
        // std::ofstream TMDhoppingNano("temp-data/tmd-hopping-nanoribbon.csv");
        // if (TMDhoppingNano.is_open())
        // {
        //     TMDhoppingNano << K.matrix() << '\n';
        // }
        // TMDhoppingNano.close();
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
    Gup->initializeUneqs(); Gdown->initializeUneqs();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[W * L];
    double LOGweight = 0.;

    double electronDensities = 0;
    double doubleOcs = 0;
    double energies = 0;
    double zzMags = 0;
    Eigen::MatrixXd magCorrZZs =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenFunctionUps =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenFunctionDowns =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd magCorrXXs =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd uneqMagCorrZZs =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd uneqMagCorrXXs =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();

    // double sign = std::copysign(1, Gup->matrix().determinant()
    //   * Gdown->matrix().determinant() );
    double sign = 1;
    double meanSign = 0;

    double electronDensity;
    double doubleOc;
    double energy;
    double zzMag;
    Eigen::MatrixXd magCorrZZ = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd magCorrXX = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd uneqMagCorrZZ = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd uneqMagCorrXX = Eigen::Matrix<double, NSITES, NSITES>::Zero();

    double nEl = 0;
    double nUp_nDw = 0;
    double Hkin = 0;
    double zzAFstFactor = 0;
    Eigen::MatrixXd SiSjZ =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd SiSjX =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd intSiTSjZ =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd intSiTSjX =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenUp = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenDown = Eigen::Matrix<double, NSITES, NSITES>::Zero();

    double nElSq = 0;
    double nUp_nDwSq = 0;
    double HkinSq = 0;
    double zzAFstFactorSq = 0;
    Eigen::MatrixXd SiSjZSq =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd SiSjXSq =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd intSiTSjZSq =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd intSiTSjXSq =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenUpSq = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenDownSq = Eigen::Matrix<double, NSITES, NSITES>::Zero();

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE,
    //  THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0; int i = 0; int latticeSweepUntilAfresh = 0; int sweep = 0; int step;


    // --- MC LOOP ---


    std::cout << "\nMC loop started. Progress:\n";
    for (step = 0; step < totalMCSteps; step++)
    {
        //  DISPLAY PROGRESS OF THE RUN.
        if ( (step + 1)  % (totalMCSteps/8) == 0 )
        {
            std::cout << (step + 1) * 1. / totalMCSteps * 100 << " %"
              << std::endl << std::endl;
            std::cout << "Average Sign: " << meanSign
              << std::endl << std::endl;
            std::cout << "Log Weight: " << LOGweight
              << std::endl << std::endl;
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
            electronDensity = 0.; doubleOc = 0.; zzMag = 0.; energy = 0.;
            for (int x = 0; x < NSITES; x++)
            {
                electronDensity -= ( Gup->get(x, x) + Gdown->get(x, x) );
                doubleOc += - Gup->get(x, x) - Gdown->get(x, x) + Gup->get(x, x)
                * Gdown->get(x, x);
                magCorrZZ(x, x) = ( Gup->get(x, x) + Gdown->get(x, x) )
                  - 2 * Gup->get(x, x) * Gdown->get(x, x);
                // magCorrXX(x, x) = - 2 * Gup->get(x, x) * Gdown->get(x, x);
                zzMag += magCorrZZ(x, x);
                energy += 2 * ( Gup->get(x, x) + Gdown->get(x, x) ) * t * K.get(x, x);
                if (l == 0)
                {
                    uneqMagCorrZZ(x, x) = ( Gup->zero(x, x) + Gdown->zero(x, x) )
                      - 2 * Gup->zero(x, x) * Gdown->zero(x, x);
                    // uneqMagCorrXX(x, x) = - 2 * Gup->zero(x, x) * Gdown->zero(x, x);
                }
                else
                {
                    uneqMagCorrZZ(x, x) = ( 1 - Gup->get(x, x) ) * ( 1 - Gup->zero(x, x) )
                     - ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->zero(x, x) )
                     - ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->zero(x, x) )
                     + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->zero(x, x) )
                     - Gup->uneqBackward(x, x) * Gup->uneqForward(x, x)
                     - Gdown->uneqBackward(x, x)  * Gdown->uneqForward(x, x);
                    // uneqMagCorrXX(x, x) =
                    // - Gup->uneqBackward(x, x) * Gdown->uneqForward(x, x)
                    // - Gdown->uneqBackward(x, x)  * Gup->uneqForward(x, x);
                }
                for (int y = 0; y < x; y++)
                {
                    magCorrZZ(x, y) =
                      - ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                      - ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->get(y, y) )
                      + ( 1 - Gup->get(x, x) ) * ( 1 - Gup->get(y, y) )
                      + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                      - Gup->get(y, x) * Gup->get(x, y)
                      - Gdown->get(y, x) * Gdown->get(x, y);
                    // magCorrXX(x, y) = - Gdown->get(x, y) * Gup->get(y, x) -
                    //   Gdown->get(y, x) * Gup->get(x, y);
                    magCorrZZ(y, x) = magCorrZZ(x, y);
                    energy += ( Gup->get(x, y) + Gdown->get(x, y)
                      + Gup->get(y, x) + Gdown->get(y, x) ) * t * K.get(x, y);
                    if ( geom == 1 or geom == 2 )
                    {
                        zzMag += 2 * pow(-1, x - y ) * magCorrZZ(x, y);
                    }

                    if ( geom == 3 or geom == 4 or geom == 5 or geom == 6 )
                    {
                        if ( ( x + ( ( x - x % int (sqrt(NSITES)) )
                          / int (sqrt(NSITES)) ) % 2 ) % 2
                          == ( y + ( ( y - y % int (sqrt(NSITES)) )
                          / int (sqrt(NSITES)) ) % 2 ) % 2 )
                        {
                            zzMag += 2 * magCorrZZ(x, y);
                        }
                        else
                        {
                            zzMag -= 2 * magCorrZZ(x, y);
                        }
                    }
                    if (l == 0)
                    {
                        uneqMagCorrZZ(x, y) =
                          - ( 1 - Gup->zero(x, x) ) * ( 1 - Gdown->zero(y, y) )
                          - ( 1 - Gdown->zero(x, x) ) * ( 1 - Gup->zero(y, y) )
                          + ( 1 - Gup->zero(x, x) ) * ( 1 - Gup->zero(y, y) )
                          + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gdown->zero(y, y) )
                          + ( 1 - Gup->zero(y, x) ) * Gup->zero(x, y)
                          + ( 1 - Gdown->zero(y, x) ) * Gdown->zero(x, y);
                        uneqMagCorrZZ(y, x) = uneqMagCorrZZs(x, y);
                        // uneqMagCorrXX(x, x) = - 2 * Gup->zero(x, y) * Gdown->zero(y, x);
                    }
                    else
                    {
                        uneqMagCorrZZ(x, y) =
                          ( 1 - Gup->get(x, x) ) * ( 1 - Gup->zero(y, y) )
                        - ( 1 - Gup->get(x, x) )  * ( 1 - Gdown->zero(y, y) )
                        - ( 1 - Gdown->get(x, x) )  * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gdown->get(x, x) )  * ( 1 - Gdown->zero(y, y) )
                        - Gup->uneqBackward(y, x) * Gup->uneqForward(x, y)
                        - Gdown->uneqBackward(y, x) * Gdown->uneqForward(x, y);
                        uneqMagCorrZZ(y, x) =
                          ( 1 - Gup->get(y, y) )  * ( 1 - Gup->zero(x, x) )
                        - ( 1 - Gup->get(y, y) )  * ( 1 - Gdown->zero(x, x) )
                        - ( 1 - Gdown->get(y, y) )  * ( 1 - Gup->zero(x, x) )
                        + ( 1 - Gdown->get(y, y) )  * ( 1 - Gdown->zero(x, x) )
                        - Gup->uneqBackward(x, y) * Gup->uneqForward(y, x)
                        - Gdown->uneqBackward(x, y) * Gdown->uneqForward(y, x);
                        // uneqMagCorrXX(x, y) = - Gup->uneqBackward(x, y) * Gdown->uneqForward(x, y)
                        // - Gdown->uneqBackward(x, y) * Gup->uneqForward(x, y);
                    }
                }
            }

            GreenFunctionUps +=
              ( Gup->matrix() * sign - GreenFunctionUps ) / ( l + 1 ) ;
            GreenFunctionDowns +=
              ( Gdown->matrix() * sign - GreenFunctionDowns ) / ( l + 1 ) ;
            electronDensity /= NSITES; electronDensity += 2;
            doubleOc /= NSITES; doubleOc += 1;
            zzMag /= NSITES; energy /= NSITES;

            electronDensities +=
              ( electronDensity * sign - electronDensities ) / ( l + 1 ) ;
            doubleOcs +=
              ( doubleOc * sign - doubleOcs ) / ( l + 1 ) ;
            magCorrZZs +=
              (magCorrZZ * sign - magCorrZZs ) / ( l + 1 );
            // magCorrXXs +=
            //   (magCorrXX * sign - magCorrXXs ) / ( l + 1 );
            uneqMagCorrZZs +=
              (uneqMagCorrZZ * sign - uneqMagCorrZZs ) / ( l + 1 );
            // uneqMagCorrXXs +=
            //   (uneqMagCorrXX * sign - uneqMagCorrXXs ) / ( l + 1 );
            zzMags +=
              (zzMag * sign - zzMags ) / ( l + 1 );
            energies +=
              ( energy * sign - energies ) / ( l + 1 ) ;


            //  DEAL WITH THE GREEN'S FUNCTIONS.


                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == GREEN_AFRESH_FREQ)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
//                Gup->computeGreenNaive(Bup->list(), l);
//                Gdown->computeGreenNaive(Bdown->list(), l);
                //  Uncomment to compute the product in the stabilized, but slightly inefficient way
//                Gup->computeStableGreenNaiveR(Bup->list(), l);
//                Gdown->computeStableGreenNaiveR(Bdown->list(), l);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    Gup->storeUDV(Bup->list(), l, GREEN_AFRESH_FREQ);
                    Gdown->storeUDV(Bdown->list(), l, GREEN_AFRESH_FREQ);
                    //  This is the standard way described in "Stable simulations
                    //  of models of interacting electrons"
//                    Gup->computeStableGreen(l, GREEN_AFRESH_FREQ);
//                    Gdown->computeStableGreen(l, GREEN_AFRESH_FREQ);
                    //  Using the BlockOfGreens Method, we can obtain
                    //  time-displaced Green's as well
                    Gup->computeBlockOfGreens(l, GREEN_AFRESH_FREQ);
                    Gdown->computeBlockOfGreens(l, GREEN_AFRESH_FREQ);
                }
                else
                {
                    Gup->storeVDU( Bup->list() ); Gdown->storeVDU( Bdown->list() );
                    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();
                    Gup->initializeUneqs(); Gdown->initializeUneqs();
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
                      GreenUp += ( GreenFunctionUps - GreenUp )
                       / ( (sweep - W)/A + 1 ) ;
                      GreenDown += ( GreenFunctionDowns - GreenDown )
                       / ( (sweep - W)/A + 1 ) ;
                      meanSign += ( sign - meanSign ) / ( ( sweep - W ) / A + 1 );
                      nEl += ( electronDensities - nEl )
                       / ( (sweep - W)/A + 1 ) ;
                      nUp_nDw += ( doubleOcs - nUp_nDw )
                       / ( (sweep - W)/A + 1 ) ;
                      SiSjZ += ( magCorrZZs - SiSjZ )
                       / ( (sweep - W)/A + 1 ) ;
                      // SiSjX += ( magCorrXXs - SiSjX )
                      //  / ( (sweep - W)/A + 1 ) ;
                      zzAFstFactor += ( zzMags - zzAFstFactor )
                       / ( (sweep - W)/A + 1 ) ;
                      Hkin += ( energies - Hkin )
                       / ( (sweep - W)/A + 1 ) ;
                      intSiTSjZ += ( uneqMagCorrZZs * BETA - intSiTSjZ )
                       / ( (sweep - W)/A + 1 ) ;
                      // intSiTSjX += ( uneqMagCorrXXs * BETA - intSiTSjX )
                      //  / ( (sweep - W)/A + 1 ) ;

                      GreenUpSq += ( GreenFunctionUps.unaryExpr(&matSq) - GreenUpSq )
                       / ( (sweep - W)/A + 1 ) ;
                      GreenDownSq += ( GreenFunctionDowns.unaryExpr(&matSq) - GreenDownSq )
                       / ( (sweep - W)/A + 1 ) ;
                      nElSq += ( pow(electronDensities, 2) - nElSq )
                       / ( (sweep - W)/A + 1 ) ;
                      nUp_nDwSq += ( pow(doubleOcs, 2) - nUp_nDwSq )
                       / ( (sweep - W)/A + 1 ) ;
                      SiSjZSq += ( magCorrZZs.unaryExpr(&matSq) - SiSjZSq )
                       / ( (sweep - W)/A + 1 ) ;
                      // SiSjXSq += ( magCorrXXs.unaryExpr(&matSq) - SiSjXSq )
                      //  / ( (sweep - W)/A + 1 ) ;
                      intSiTSjZSq += ( ( BETA * uneqMagCorrZZs).unaryExpr(&matSq)  - intSiTSjZSq )
                       / ( (sweep - W)/A + 1 ) ;
                      // intSiTSjXSq += ( ( BETA * uneqMagCorrXXs).unaryExpr(&matSq) - intSiTSjXSq )
                      //  / ( (sweep - W)/A + 1 ) ;
                      zzAFstFactorSq += ( pow(zzMags, 2) - zzAFstFactorSq )
                       / ( (sweep - W)/A + 1 ) ;
                      HkinSq += ( pow(energies, 2) - HkinSq )
                       / ( (sweep - W)/A + 1 ) ;
                     }
                     electronDensities = 0.; doubleOcs = 0.;
                     magCorrZZs = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                     GreenFunctionUps = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                     GreenFunctionDowns = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                     // magCorrXXs = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                     uneqMagCorrZZs = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                     // uneqMagCorrXXs = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                     zzMags = 0.;
                     energies = 0.;
                }
                //  MOVE SWEEP COUNTER
                sweep += 1;
                l = 0; i = 0;
            }
        }
    }   //  END OF MC LOOP.

    //  Normalize to mean sign
    nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign; zzAFstFactor /= meanSign;
    Hkin /= meanSign; GreenUp /= meanSign; GreenDown /= meanSign;
    nElSq /= meanSign; nUp_nDwSq /= meanSign; SiSjZSq /= meanSign; zzAFstFactorSq /= meanSign;
    HkinSq /= meanSign; GreenUpSq /= meanSign; GreenDownSq /= meanSign;

    std::cout << "Simulation ended" << std::endl << std::endl;
    std::cout << "ds / <s>: " << sqrt( 1 - pow(meanSign, 2) )
     / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) / meanSign << std::endl << std::endl;
    std::cout << "nEl: " << nEl << " +- " <<
     sqrt( nElSq - pow(nEl, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << std::endl << std::endl;
    std::cout << "nUp_nDw: " << nUp_nDw << " +- " <<
     sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << std::endl << std::endl;
    std::cout << "< m^2 >: " << nEl - 2 * nUp_nDw << " +- " <<
     (nEl - 2 * nUp_nDw) * ( sqrt( nElSq - pow(nEl, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) / nEl
     + 2 * sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) / nUp_nDw ) << std::endl << std::endl;
    std::cout << "Hkin: " << Hkin << " +- " <<
     sqrt( HkinSq - pow(Hkin, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << std::endl << std::endl;
    std::cout << "Hint: " << U * nUp_nDw << " +- " <<
     U * sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << std::endl << std::endl;
    std::cout << "E: " << Hkin + U * nUp_nDw + U / 2 * nEl << " +- " <<
     sqrt( HkinSq - pow(Hkin, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) +
     U * sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) )<< std::endl << std::endl;

    //  SAVE OUTPUT.
    std::ofstream file0("temp-data/simulationParameters.csv");
    if (file0.is_open())
    {
      file0 << std::left << std::setw(50) << "Number of sites," << NSITES << '\n';
      file0 << std::left << std::setw(50) << "dt," << dt << '\n';
      file0 << std::left << std::setw(50) << "beta," << BETA << '\n';
      file0 << std::left << std::setw(50) << "L," << L << '\n';
      file0 << std::left << std::setw(50) << "t," << t << '\n';
      file0 << std::left << std::setw(50) << "U," << U << '\n';
      file0 << std::left << std::setw(50) << "mu," << mu << '\n';
      file0 << std::left << std::setw(50) << "totalMCSweeps," << totalMCSweeps << '\n';
      file0 << std::left << std::setw(50) << "Frequency of recomputing G,"
        << GREEN_AFRESH_FREQ << '\n';
      file0 << std::left << std::setw(50)
        << "Number of multiplied Bs after stabilization," << Lbda << '\n';
      file0 << std::left << std::setw(50) << "Geometry," << geom << '\n';
      file0 << std::left << std::setw(50) << "Ny," << Ny << '\n';
    } file0.close();
    //  STORE MEASUREMENTS
    std::ofstream file1("temp-data/Log-weights.csv");
    std::ofstream file2("temp-data/MeasurementsScalars.csv");
    std::ofstream file3("temp-data/EqTimeSzCorrelations.csv");
    std::ofstream file4("temp-data/EqTimeSzCorrelationsError.csv");
    std::ofstream file5("temp-data/UneqTimeSzCorrelations.csv");
    std::ofstream file6("temp-data/UneqTimeSzCorrelationsError.csv");
    // std::ofstream file7("temp-data/EqTimeSxCorrelations.csv");
    // std::ofstream file8("temp-data/EqTimeSxCorrelationsError.csv");
    // std::ofstream file9("temp-data/UneqTimeSxCorrelations.csv");
    // std::ofstream file10("temp-data/UneqTimeSxCorrelationsError.csv");
    std::ofstream file11("temp-data/GreenUp.csv");
    std::ofstream file12("temp-data/GreenUpError.csv");
    std::ofstream file13("temp-data/GreenDown.csv");
    std::ofstream file14("temp-data/GreenDownError.csv");
    if ( file1.is_open() and file2.is_open() and
      file3.is_open() and file4.is_open() and
      file5.is_open() and file6.is_open() and
     //  and
     // file7.is_open() and file8.is_open() and
     // file9.is_open() and file10.is_open()
      file11.is_open() and file12.is_open() and
      file13.is_open() and file14.is_open()
    )
    {
        file1 << std::left << std::setw(50) << "Configuration log weight" << '\n';
        for (int s = 0; s < W; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file1 << std::left << std::setw(50) << weights[s * L + slice] << '\n';
            }
        }
        file2 << std::left << std::setw(50) << "Electron density <n>,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nEl << '\n';
        file2 << std::left << std::setw(50) << "d<n>,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << sqrt( nElSq - pow(nEl, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << '\n';
        file2 << std::left << std::setw(50) << "Double occupancy <n+ n->,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "d<n+ n->,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << '\n';
        file2 << std::left << std::setw(50) << "ZZ AF Structure Factor,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << zzAFstFactor << '\n';
        file2 << std::left << std::setw(50) << "< m^2 >,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nEl - 2 * nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "d< m^2 >,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << (nEl - 2 * nUp_nDw) * ( sqrt( nElSq - pow(nEl, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) / nEl
         + 2 * sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) / nUp_nDw ) << '\n';
        file2 << std::left << std::setw(50) << "Hkin,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << Hkin << '\n';
        file2 << std::left << std::setw(50) << "dHkin,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << sqrt( HkinSq - pow(Hkin, 2) ) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) << '\n';
        file2 << std::left << std::setw(50) << "Hint,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << U * nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "E,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << Hkin + U * nUp_nDw + U / 2 * nEl << '\n';
        file2 << std::left << std::setw(50) << "Average sign,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << meanSign << '\n';
        file3 << std::left << std::setw(50) << "<Sz_i Sz_j >" << '\n';
        file3 << std::setprecision(10) << SiSjZ.format(CleanFmt) << '\n';
        file4 << std::left << std::setw(50) << "d<Sz_i Sz_j >" << '\n';
        file4 <<
         std::setprecision(10) << ( ( SiSjZSq - SiSjZ.unaryExpr(&matSq) )
         .unaryExpr(&matSqrt) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) )
         .format(CleanFmt) << '\n';
        file5 << std::left << std::setw(50) << "int_0^beta dt <Sz_i (t) Sz_j (0) >" << '\n';
        file5 << std::setprecision(10) << intSiTSjZ.format(CleanFmt) << '\n';
        file6 << std::left << std::setw(50) << "d int_0^beta dt <Sz_i (t) Sz_j (0) >" << '\n';
        file6 <<
         std::setprecision(10) << ( ( intSiTSjZSq - intSiTSjZ.unaryExpr(&matSq) )
         .unaryExpr(&matSqrt) / sqrt( (totalMCSweeps - W) / A - 1 ) )
         .format(CleanFmt) << '\n';
        // file7 << std::left << std::setw(50) << "<Sx_i Sx_j >" << '\n';
        // file7 << std::setprecision(10) << SiSjX.format(CleanFmt) << '\n';
        // file8 << std::left << std::setw(50) << "d<Sx_i Sx_j >" << '\n';
        // file8 <<
        //  std::setprecision(10) << ( ( SiSjXSq - SiSjX.unaryExpr(&matSq) )
        //  .unaryExpr(&matSqrt) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) )
        //  .format(CleanFmt) << '\n';
        // file9 << std::left << std::setw(50) << "int_0^beta dt <Sx_i (t) Sx_j (0) >" << '\n';
        // file9 << std::setprecision(10) << intSiTSjX.format(CleanFmt) << '\n';
        // file10 << std::left << std::setw(50) << "d int_0^beta dt <Sx_i (t) Sx_j (0) >" << '\n';
        // file10 <<
        //  std::setprecision(10) << ( ( intSiTSjXSq - intSiTSjX.unaryExpr(&matSq) )
        //  .unaryExpr(&matSqrt) / sqrt( (totalMCSweeps - W) / A - 1 ) )
        //  .format(CleanFmt) << '\n';
        file11 << std::left << std::setw(50) << "Gup" << '\n';
        file11 << std::setprecision(10) << GreenUp.format(CleanFmt) << '\n';
        file12 << std::left << std::setw(50) << "dGup" << '\n';
        file12 <<
         std::setprecision(10) << ( ( GreenUpSq - GreenUp.unaryExpr(&matSq) )
         .unaryExpr(&matSqrt) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) )
         .format(CleanFmt) << '\n';
        file13 << std::left << std::setw(50) << "Gdown" << '\n';
        file13 << std::setprecision(10) << GreenDown.format(CleanFmt) << '\n';
        file14 << std::left << std::setw(50) << "dGdown" << '\n';
        file14 <<
         std::setprecision(10) << ( ( GreenDownSq - GreenDown.unaryExpr(&matSq) )
         .unaryExpr(&matSqrt) / sqrt( ( (totalMCSweeps - W) / A - 1 ) ) )
         .format(CleanFmt) << '\n';
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();
    file6.close();
    // file7.close();
    // file8.close();
    // file9.close();
    // file10.close();
    file11.close();
    file12.close();
    file13.close();
    file14.close();

    delete[] weights;
    delete Gup; delete Gdown; delete h; delete Bup; delete Bdown;

    return 0;
    }
