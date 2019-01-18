//
//  matrixgen.h
//
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef matrixgen_h
#define matrixgen_h

#ifndef NORB
#define NORB 3  //   number of orbitals
#endif

double matSq(double x) // element-wise square of a matrix
{
    return std::pow(x, 2);
}

double matSqrt(double x) // element-wise square root of a matrix
{
    return std::sqrt(x);
}

void write(double meanSign, int sweep, int W, int A, double nEl, double nUp_nDw,
  double Hkin, double U, int nSites, double dt, double beta, int L, double t,
  double mu, int green_afresh_freq, int Lbda, int geom, int Ny, double * weights, Eigen::MatrixXd SiSjZ,
  double * corrs, int totalMCSweeps)
{
    //  Normalize to mean sign
    nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign;
    Hkin /= meanSign;

    int precision = 10;
    Eigen::IOFormat CleanFmt(precision, 0, ", ", "\n", "", "");

    if (VERBOSE == 1)
    {
        std::cout << "Writing results" << std::endl << std::endl;
        std::cout << "<s>: " << meanSign << std::endl << std::endl;
        std::cout << "ds / <s>: " << sqrt( 1 - pow(meanSign, 2) )
         / sqrt( ( (sweep - W) / A - 1 ) ) / meanSign
          << std::endl << std::endl;
        std::cout << "nEl: " << nEl << std::endl << std::endl;
        std::cout << "Hkin: " << Hkin << std::endl << std::endl;
        std::cout << "U nUp_nDw: " << U * nUp_nDw << std::endl << std::endl;
    }


    //  SAVE OUTPUT.
    std::ofstream file0("temp-data/simulationParameters.csv");
    if (file0.is_open())
    {
      file0 << "Number of sites NSITES," << NSITES << '\n';
      file0 << "Trotter Error dt," << dt << '\n';
      file0 << "Inverse Temperature BETA," << BETA << '\n';
      file0 << "Number of Imaginary-time Slices L," << L << '\n';
      file0 << "Hopping Normalization t," << t << '\n';
      file0 << "On-site interaction U," << U << '\n';
      file0 << "mu," << mu << '\n';
      file0 << "totalMCSweeps," << sweep << '\n';
      file0 << "Frequency of recomputing G,"
        << GREEN_AFRESH_FREQ << '\n';
      file0 << "Number of multiplied Bs after stabilization," << Lbda << '\n';
      file0 << "Geometry," << geom << '\n';
      file0 << "Ny," << Ny << '\n';
    } file0.close();
    //  STORE MEASUREMENTS
    std::ofstream file1("temp-data/Log-weights.csv");
    std::ofstream file2("temp-data/MeasurementsScalars.csv");
    std::ofstream file3("temp-data/EqTimeSzCorrelations.csv");
    std::ofstream file4("temp-data/corrs.csv");
    if ( file1.is_open() and file2.is_open()
     and file3.is_open() and file4.is_open() )
    {
        file1 << std::left << std::setw(50) << "Configuration log weight" << '\n';
        file4 << std::left << std::setw(50) << "Correlation" << '\n';
        for (int s = 0; s < W; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file1 << std::left << std::setw(50) << weights[s * L + slice] << '\n';
                file4 << std::left << std::setw(50) << corrs[s * L + slice] << '\n';
            }
        }
        for (int s = W; s < totalMCSweeps; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file4 << std::left << std::setw(50) << corrs[s * L + slice] << '\n';
            }
        }
        file2 << std::left << std::setw(50) << "Electron density <n>,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nEl << '\n';
        file2 << std::left << std::setw(50) << "Double occupancy <n+ n->,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "Hkin,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << Hkin << '\n';
        file2 << std::left << std::setw(50) << "Average sign,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << meanSign << '\n';
        file3 << std::left << std::setw(50) << "<SiSj>" << '\n';
        file3 << std::setprecision(10) << SiSjZ.format(CleanFmt) << '\n';
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
}

void writeTMDNR(double meanSign, int sweep, int W, int A, double nEl, double nUp_nDw,
  double Hkin, double U, int nSites, double dt, double beta, int L, double t,
  double mu, int green_afresh_freq, int Lbda, int geom, int Ny, double * weights, Eigen::MatrixXd SiSjZ,
double * corrs, int totalMCSweeps)
{
    //  Normalize to mean sign
    nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign;
    Hkin /= meanSign;

    int Nx = nSites / Ny / NORB;
    Eigen::MatrixXd spin_corr;
    spin_corr = Eigen::MatrixXd::Zero(Ny * NORB, nSites);

    for (int a = 0; a < NORB; a++)
        for (int x1 = 0; x1 < Nx; x1++)
            for (int y1 = 0; y1 < Ny; y1++)
                for (int b = 0; b < NORB; b++)
                    for (int x2 = 0; x2 < Nx; x2++)
                        for (int y2 = 0; y2 < Ny; y2++)
                        {
                            spin_corr( NORB * y1 + a ,
                            NORB * ( Nx * y2 + abs(x2 - x1) ) + b ) +=
                            SiSjZ( NORB * ( Nx * y1 + x1 ) + a ,
                            NORB * ( Nx * y2 + x2 ) + b ) / 2 / Nx;

                            spin_corr( NORB * y1 + a ,
                            NORB * ( Nx * y2 + ( Nx - abs(x2 - x1) ) % Nx ) + b ) +=
                            SiSjZ( NORB * ( Nx * y1 + x1 ) + a ,
                            NORB * ( Nx * y2 + x2 ) + b ) / 2 / Nx;
                        }

    int precision = 10;
    Eigen::IOFormat CleanFmt(precision, 0, ", ", "\n", "", "");

    if (VERBOSE == 1)
    {
        std::cout << "Writing results" << std::endl << std::endl;
        std::cout << "<s>: " << meanSign << std::endl << std::endl;
        std::cout << "ds / <s>: " << sqrt( 1 - pow(meanSign, 2) )
         / sqrt( ( (sweep - W) / A - 1 ) ) / meanSign
          << std::endl << std::endl;
        std::cout << "nEl: " << nEl << std::endl << std::endl;
        std::cout << "Hkin: " << Hkin << std::endl << std::endl;
        std::cout << "U nUp_nDw: " << U * nUp_nDw << std::endl << std::endl;
    }


    //  SAVE OUTPUT.
    std::ofstream file0("temp-data/simulationParameters.csv");
    if (file0.is_open())
    {
      file0 << "Number of sites NSITES," << NSITES << '\n';
      file0 << "Trotter Error dt," << dt << '\n';
      file0 << "Inverse Temperature BETA," << BETA << '\n';
      file0 << "Number of Imaginary-time Slices L," << L << '\n';
      file0 << "Hopping Normalization t," << t << '\n';
      file0 << "On-site interaction U," << U << '\n';
      file0 << "mu," << mu << '\n';
      file0 << "totalMCSweeps," << sweep << '\n';
      file0 << "Frequency of recomputing G,"
        << GREEN_AFRESH_FREQ << '\n';
      file0 << "Number of multiplied Bs after stabilization," << Lbda << '\n';
      file0 << "Geometry," << geom << '\n';
      file0 << "Ny," << Ny << '\n';
    } file0.close();
    //  STORE MEASUREMENTS
    std::ofstream file1("temp-data/Log-weights.csv");
    std::ofstream file2("temp-data/MeasurementsScalars.csv");
    std::ofstream file3("temp-data/EqTimeSzCorrelations.csv");
    std::ofstream file4("temp-data/corrs.csv");
    if ( file1.is_open() and file2.is_open()
     and file3.is_open() and file4.is_open() )
    {
        file1 << std::left << std::setw(50) << "Configuration log weight" << '\n';
        file4 << std::left << std::setw(50) << "Correlation" << '\n';
        for (int s = 0; s < W; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file1 << std::left << std::setw(50) << weights[s * L + slice] << '\n';
                file4 << std::left << std::setw(50) << corrs[s * L + slice] << '\n';
            }
        }
        for (int s = W; s < totalMCSweeps; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file4 << std::left << std::setw(50) << corrs[s * L + slice] << '\n';
            }
        }
        file2 << std::left << std::setw(50) << "Electron density <n>,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nEl << '\n';
        file2 << std::left << std::setw(50) << "Double occupancy <n+ n->,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "Hkin,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << Hkin << '\n';
        file2 << std::left << std::setw(50) << "Average sign,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << meanSign << '\n';
        file3 << std::left << std::setw(50) << "<SiSj>" << '\n';
        file3 << std::setprecision(precision)
        << spin_corr.format(CleanFmt) << '\n';
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
}

void writeALL(double meanSign, int sweep, int W, int A, double nElSq, double nEl,
 double nUp_nDw, double nUp_nDwSq, double Hkin,
 double HkinSq, double U, int nSites, double dt, double beta, int L, double t,
 double mu, int green_afresh_freq, int Lbda, int geom, int Ny, double * weights,
 Eigen::MatrixXd SiSjZ, Eigen::MatrixXd SiSjZSq)
{
    //  Normalize to mean sign
    nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign;
    Hkin /= meanSign;
    nElSq /= meanSign; nUp_nDwSq /= meanSign; SiSjZSq /= meanSign;
    HkinSq /= meanSign;

    Eigen::IOFormat CleanFmt(10, 0, ", ", "\n", "", "");

    if (VERBOSE == 1)
    {
        std::cout << "Writing results" << std::endl << std::endl;
        std::cout << "<s>: " << meanSign << std::endl << std::endl;
        std::cout << "ds / <s>: " << sqrt( 1 - pow(meanSign, 2) )
         / sqrt( ( (sweep - W) / A - 1 ) ) / meanSign
          << std::endl << std::endl;
        std::cout << "nEl: " << nEl << " +- " <<
         sqrt( nElSq - pow(nEl, 2) ) / sqrt( ( (sweep - W) / A - 1 ) )
          << std::endl << std::endl;
        std::cout << "nUp_nDw: " << nUp_nDw << " +- " <<
         sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (sweep - W) / A - 1 ) )
          << std::endl << std::endl;
        std::cout << "Hkin: " << Hkin << " +- " <<
         sqrt( HkinSq - pow(Hkin, 2) ) / sqrt( ( (sweep - W) / A - 1 ) )
          << std::endl << std::endl;
        std::cout << "Hint: " << U * nUp_nDw << " +- " <<
         U * sqrt( nUp_nDwSq - pow(nUp_nDw, 2) )
          / sqrt( ( (sweep - W) / A - 1 ) ) << std::endl << std::endl;
    }


    //  SAVE OUTPUT.
    std::ofstream file0("temp-data/simulationParameters.csv");
    if (file0.is_open())
    {
      file0 << "Number of sites NSITES," << NSITES << '\n';
      file0 << "Trotter Error dt," << dt << '\n';
      file0 << "Inverse Temperature BETA," << BETA << '\n';
      file0 << "Number of Imaginary-time Slices L," << L << '\n';
      file0 << "Hopping Normalization t," << t << '\n';
      file0 << "On-site interaction U," << U << '\n';
      file0 << "mu," << mu << '\n';
      file0 << "totalMCSweeps," << sweep << '\n';
      file0 << "Frequency of recomputing G,"
        << GREEN_AFRESH_FREQ << '\n';
      file0 << "Number of multiplied Bs after stabilization," << Lbda << '\n';
      file0 << "Geometry," << geom << '\n';
      file0 << "Ny," << Ny << '\n';
    } file0.close();
    //  STORE MEASUREMENTS
    std::ofstream file1("temp-data/Log-weights.csv");
    std::ofstream file2("temp-data/MeasurementsScalars.csv");
    std::ofstream file3("temp-data/EqTimeSzCorrelations.csv");
    std::ofstream file4("temp-data/EqTimeSzCorrelationsError.csv");
    if ( file1.is_open() and file2.is_open()
     and file3.is_open() and file4.is_open() )
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
        << sqrt( nElSq - pow(nEl, 2) ) / sqrt( ( (sweep - W) / A - 1 ) ) << '\n';
        file2 << std::left << std::setw(50) << "Double occupancy <n+ n->,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "d<n+ n->,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << sqrt( nUp_nDwSq - pow(nUp_nDw, 2) ) / sqrt( ( (sweep - W) / A - 1 ) ) << '\n';
        file2 << std::left << std::setw(50) << "Hkin,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << Hkin << '\n';
        file2 << std::left << std::setw(50) << "dHkin,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << sqrt( HkinSq - pow(Hkin, 2) )
         / sqrt( ( (sweep - W) / A - 1 ) ) << '\n';
        file2 << std::left << std::setw(50) << "Hint,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << U * nUp_nDw << '\n';
        file2 << std::left << std::setw(50) << "Average sign,";
        file2 << std::left << std::setw(50) << std::setprecision(10)
        << meanSign << '\n';
        file3 << std::left << std::setw(50) << "<Sz_i Sz_j >" << '\n';
        file3 << std::setprecision(10) << SiSjZ.format(CleanFmt) << '\n';
        file4 << std::left << std::setw(50) << "d<Sz_i Sz_j >" << '\n';
        file4 <<
         std::setprecision(10) << ( ( SiSjZSq - SiSjZ.unaryExpr(&matSq) )
         .unaryExpr(&matSqrt) / sqrt( ( (sweep - W) / A - 1 ) ) )
         .format(CleanFmt) << '\n';
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
}

template<int N>
class Geometry
{
    Eigen::MatrixXd B;
    Eigen::MatrixXd Hoppings;
    //  the minimal model of Liu2013 has 9 parameters
    //  the hopping t0 is taken as the argument t, which in the
    //  uniform hopping case is multiplied by all elements
    //  of the hopping matrix
    double e1; double e2;
    double t0; double t1; double t2;
    double t11; double t12; double t22;
    double lbd;
    Eigen::Matrix<double, 3, 3> Ezero;
    Eigen::Matrix<double, 3, 3> Eone;
    Eigen::Matrix<double, 3, 3> Etwo;
    Eigen::Matrix<double, 3, 3> Ethree;
    Eigen::Matrix<double, 3, 3> Efour;
    Eigen::Matrix<double, 3, 3> Efive;
    Eigen::Matrix<double, 3, 3> Esix;
public:
    //  For the 1D and square lattices, we use the checkerboard breakup.
    //  Otherwise, we compute the matrix exponential using Eigen.
    void oneDimensionalChainPBC(double t, double dt, double mu);
    void oneDimensionalChainOBC(double t, double dt, double mu);
    void twoDimensionalRectanglePBC(int Ny, double t, double dt, double mu);
    void twoDimensionalRectangleOBC(int Ny, double t, double dt, double mu);
    void computeExponential(double t, double dt);
    Eigen::MatrixXd matrix();
    double get(int x, int y);
    void trianglePBC();
    void triangleNanoribbon(int Ny);
    void hcPBC();
    void hcNanoribbon(int Ny); //  Ny = width of the ribbon
    void hcStrainedNanoribbon(int Ny, double Delta); //  Ny = width of the ribbon
    void setParamsThreeOrbitalTB(double threeOrbitalTBparameters[8]);
    void tmdPBC();
    void tmdNanoribbon(int Ny); //  Ny = width of the ribbon
    void tmdNanoribbonStrained(int Ny, double Delta); //  Ny = width of the ribbon
    Eigen::MatrixXd BpreFactor(double dt, double mu);
    Geometry() : B(N, N), Hoppings(N, N) {
    };
};

template<int N>
void Geometry<N>::oneDimensionalChainPBC(double t, double dt, double mu)
{
   Hoppings = Eigen::Matrix<double, N, N>::Zero();
   //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
   Hoppings(0, 1) += 1.;
   Hoppings(0, N - 1) += 1.;
   Hoppings(N - 1, 0) += 1.;
   Hoppings(N - 1, N - 2) += 1.;
   //  Set the remaining ones
   for (int i = 1; i < N - 1; i++)
   {
       Hoppings(i, i - 1) += 1; Hoppings(i, i + 1) += 1;
   }
    Eigen::MatrixXd exp_k1 = Eigen::Matrix<double, N, N>::Zero();
    Eigen::MatrixXd exp_k2 = Eigen::Matrix<double, N, N>::Zero();
    for (int i = 1; i < N / 2; i++)
    {
        exp_k1(2 * i, 2 * i) = cosh( t * dt );
        exp_k1(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1(0, 0) = cosh( t * dt );
    exp_k1(1, 1) = cosh( t * dt );
    exp_k1(0, 1) = sinh( t * dt );
    exp_k1(1, 0) = sinh( t * dt );
    exp_k2(0, 0) = cosh( t * dt / 2 );
    exp_k2(N - 1, N - 1) = cosh( t * dt / 2 );
    exp_k2(0, N - 1) = sinh( t * dt / 2 );
    exp_k2(N - 1, 0) = sinh( t * dt / 2 );
    B = exp_k2 * exp_k1 * exp_k2;
}

template<int N>
void Geometry<N>::oneDimensionalChainOBC(double t, double dt, double mu)
{
    Hoppings = Eigen::Matrix<double, N, N>::Zero();
    //  Set the elements of the hopping matrix that define OBC
    //  corresponding to the ends of the 1D chain
    Hoppings(0, 1) += 1.;
    //  B(0, N - 1) += 1.; //   This hopping only occurs for PBC
    //  B(N - 1, 0) += 1.; //   So does this one
    Hoppings(N - 1, N - 2) += 1.;
    //  Set the remaining ones
    for (int i = 1; i < N - 1; i++)
    {
       Hoppings(i, i - 1) += 1; Hoppings(i, i + 1) += 1;
    }
    Eigen::MatrixXd exp_k1 = Eigen::Matrix<double, N, N>::Zero();
    Eigen::MatrixXd exp_k2 = Eigen::Matrix<double, N, N>::Zero();
    for (int i = 1; i < N / 2; i++)
    {
        exp_k1(2 * i, 2 * i) = cosh( t * dt );
        exp_k1(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1(0, 0) = cosh( t * dt );
    exp_k1(1, 1) = cosh( t * dt );
    exp_k1(0, 1) = sinh( t * dt );
    exp_k1(1, 0) = sinh( t * dt );
    exp_k2(0, 0) = 1;
    exp_k2(N - 1, N - 1) = 1;
    B = exp_k2 * exp_k1 * exp_k2;
}

template<int N>
void Geometry<N>::twoDimensionalRectanglePBC(int Ny, double t, double dt, double mu)
{
    int Nx = N / Ny;

    Eigen::MatrixXd HopX = Eigen::MatrixXd::Zero(Nx, Nx);
    Eigen::MatrixXd HopY =  Eigen::MatrixXd::Zero(Ny, Ny);

    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    HopX(0, 1) += 1.;
    HopX(0, Nx - 1) += 1.;
    HopX(Nx - 1, 0) += 1.;
    HopX(Nx - 1, Nx - 2) += 1.;
    //  Set the remaining ones
    for (int i = 1; i < Nx - 1; i++)
    {
        HopX(i, i - 1) += 1; HopX(i, i + 1) += 1;
    }

    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    HopY(0, 1) += 1.;
    HopY(0, Ny - 1) += 1.;
    HopY(Ny - 1, 0) += 1.;
    HopY(Ny - 1, Ny - 2) += 1.;
    //  Set the remaining ones
    for (int i = 1; i < Ny - 1; i++)
    {
        HopY(i, i - 1) += 1; HopY(i, i + 1) += 1;
    }

    Hoppings = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(Nx, Nx) , HopY ) +
     Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(Ny, Ny) , HopX );

    Eigen::MatrixXd exp_k1x = Eigen::MatrixXd::Zero(Nx, Nx);
    Eigen::MatrixXd exp_k2x = Eigen::MatrixXd::Zero(Nx, Nx);
    for (int i = 1; i < Nx / 2; i++)
    {
        exp_k1x(2 * i, 2 * i) = cosh( t * dt );
        exp_k1x(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1x(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1x(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2x(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2x(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1x(0, 0) = cosh( t * dt );
    exp_k1x(1, 1) = cosh( t * dt );
    exp_k1x(0, 1) = sinh( t * dt );
    exp_k1x(1, 0) = sinh( t * dt );
    exp_k2x(0, 0) = cosh( t * dt / 2 );
    exp_k2x(Nx - 1, Nx - 1) = cosh( t * dt / 2 );
    exp_k2x(0, Nx - 1) = sinh( t * dt / 2 );
    exp_k2x(Nx - 1, 0) = sinh( t * dt / 2 );

    Eigen::MatrixXd exp_k1y = Eigen::MatrixXd::Zero(Ny, Ny);
    Eigen::MatrixXd exp_k2y = Eigen::MatrixXd::Zero(Ny, Ny);
    for (int i = 1; i < Ny / 2; i++)
    {
        exp_k1y(2 * i, 2 * i) = cosh( t * dt );
        exp_k1y(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1y(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1y(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2y(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2y(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1y(0, 0) = cosh( t * dt );
    exp_k1y(1, 1) = cosh( t * dt );
    exp_k1y(0, 1) = sinh( t * dt );
    exp_k1y(1, 0) = sinh( t * dt );
    exp_k2y(0, 0) = cosh( t * dt / 2 );
    exp_k2y(Ny - 1, Ny - 1) = cosh( t * dt / 2 );
    exp_k2y(0, Ny - 1) = sinh( t * dt / 2 );
    exp_k2y(Ny - 1, 0) = sinh( t * dt / 2 );

    exp_k1x = exp_k2x * exp_k1x * exp_k2x;
    exp_k1y = exp_k2y * exp_k1y * exp_k2y;

    B = Eigen::kroneckerProduct(exp_k1y , exp_k1x);
}

template<int N>
void Geometry<N>::twoDimensionalRectangleOBC(int Ny, double t, double dt, double mu)
{
    int Nx = N / Ny;

    Eigen::MatrixXd HopX = Eigen::MatrixXd::Zero(Nx, Nx);
    Eigen::MatrixXd HopY =  Eigen::MatrixXd::Zero(Ny, Ny);

    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    HopX(0, 1) += 1.;
    HopX(Nx - 1, Nx - 2) += 1.;
    //  Set the remaining ones
    for (int i = 1; i < Nx - 1; i++)
    {
        HopX(i, i - 1) += 1; HopX(i, i + 1) += 1;
    }

    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    HopY(0, 1) += 1.;
    HopY(Ny - 1, Ny - 2) += 1.;
    //  Set the remaining ones
    for (int i = 1; i < Ny - 1; i++)
    {
        HopY(i, i - 1) += 1; HopY(i, i + 1) += 1;
    }

    Hoppings = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(Nx, Nx) , HopY ) +
     Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(Ny, Ny) , HopX );

    Eigen::MatrixXd exp_k1x = Eigen::MatrixXd::Zero(Nx, Nx);
    Eigen::MatrixXd exp_k2x = Eigen::MatrixXd::Zero(Nx, Nx);
    for (int i = 1; i < Nx / 2; i++)
    {
        exp_k1x(2 * i, 2 * i) = cosh( t * dt );
        exp_k1x(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1x(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1x(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2x(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2x(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1x(0, 0) = cosh( t * dt );
    exp_k1x(1, 1) = cosh( t * dt );
    exp_k1x(0, 1) = sinh( t * dt );
    exp_k1x(1, 0) = sinh( t * dt );
    exp_k2x(0, 0) = 1;
    exp_k2x(Nx - 1, Nx - 1) = 1;

    Eigen::MatrixXd exp_k1y = Eigen::MatrixXd::Zero(Ny, Ny);
    Eigen::MatrixXd exp_k2y = Eigen::MatrixXd::Zero(Ny, Ny);
    for (int i = 1; i < Ny / 2; i++)
    {
        exp_k1y(2 * i, 2 * i) = cosh( t * dt );
        exp_k1y(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1y(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1y(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2y(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2y(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1y(0, 0) = cosh( t * dt );
    exp_k1y(1, 1) = cosh( t * dt );
    exp_k1y(0, 1) = sinh( t * dt );
    exp_k1y(1, 0) = sinh( t * dt );
    exp_k2y(0, 0) = 1;
    exp_k2y(Ny - 1, Ny - 1) = 1;

    exp_k1x = exp_k2x * exp_k1x * exp_k2x;
    exp_k1y = exp_k2y * exp_k1y * exp_k2y;

    B = Eigen::kroneckerProduct(exp_k1y , exp_k1x);
}

template<int N>
void Geometry<N>::computeExponential(double t, double dt)
{
    B = (t * dt * B).exp();
}

template<int N>
void Geometry<N>::trianglePBC()
{

}

template<int N>
void Geometry<N>::triangleNanoribbon(int Ny)
{

}


template<int N>
void Geometry<N>::hcPBC()
{
    B = Eigen::Matrix<double, N, N>::Zero();
    int Ny = sqrt(N / 2);
    int Nx = N / Ny / 2;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  SUBLATTICE A
            if (y == Ny - 1)
            {
                if (x == 0)
                {
                    B(Nx * y + x, Nx * Ny + Nx * y) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1;
                    B(x, Nx * Ny + Nx - 1) = 1; // additional
                    B(Nx * Ny + Nx * y, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1;
                    B(Nx * Ny + Nx - 1, x) = 1; //  additional
                }
                else
                {
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1;
                    B(x, Nx * Ny + x - 1) = 1; // additional
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1;
                    B(Nx * Ny + x - 1, x) = 1; //  additional
                }
            }
            else
            {
                if (x == 0)
                {
                    B(Nx * y + x, Nx * Ny + Nx * y) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * (y + 1) + Nx - 1) = 1;
                    B(Nx * Ny + Nx * y, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * (y + 1) + Nx - 1, Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * ( y + 1 ) + x - 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * ( y + 1 ) + x - 1, Nx * y + x) = 1;
                }
            }
            //  SUBLATTICE B
            if (y == 0)
            {
                if (x == Nx - 1)
                {
                    B(Nx * Ny + Nx * y + x, Nx * y) = 1;
                    B(Nx * Ny + Nx * ( Ny - 1 ) + x, Nx * (Ny - 1) ) = 1; // additional
                    B(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1;
                    B(Nx * y, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * (Ny - 1), Nx * Ny + Nx * ( Ny - 1 ) + x) = 1; // additional
                    B(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1;
                    B(Nx * Ny + Nx * ( Ny - 1 ) + x, Nx * ( Ny - 1 ) + x + 1) = 1; // additional
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * ( Ny - 1 ) + x + 1, Nx * Ny + Nx * ( Ny - 1 ) + x) = 1; // additional
                }
            }
            else
            {
                if (x == Nx - 1)
                {
                    B(Nx * Ny + Nx * y + x, Nx * y) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) ) = 1;
                    B(Nx * y, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * ( y - 1 ) , Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) + x + 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * ( y - 1 ) + x + 1, Nx * Ny + Nx * y + x) = 1;
                }
            }
        }
    }
    Hoppings = B;
}

template<int N>
void Geometry<N>::hcNanoribbon(int Ny)
{   //  Ny = width of the ribbon
    B = Eigen::Matrix<double, N, N>::Zero();
    int Nx = N / Ny / 2;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  SUBLATTICE A
            if (y == Ny - 1)
            {
                if (x == 0)
                {
                    B(Nx * y + x, Nx * Ny + Nx * y) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1;
                    B(Nx * Ny + Nx * y, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1;
                }
            }
            else
            {
                if (x == 0)
                {
                    B(Nx * y + x, Nx * Ny + Nx * y) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * (y + 1) + Nx - 1) = 1;
                    B(Nx * Ny + Nx * y, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * (y + 1) + Nx - 1, Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * ( y + 1 ) + x - 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * ( y + 1 ) + x - 1, Nx * y + x) = 1;
                }
            }
            //  SUBLATTICE B
            if (y == 0)
            {
                if (x == Nx - 1)
                {
                    B(Nx * Ny + Nx * y + x, Nx * y) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1;
                    B(Nx * y, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1;
                }
            }
            else
            {
                if (x == Nx - 1)
                {
                    B(Nx * Ny + Nx * y + x, Nx * y) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) ) = 1;
                    B(Nx * y, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * ( y - 1 ) , Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) + x + 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1;
                    B(Nx * ( y - 1 ) + x + 1, Nx * Ny + Nx * y + x) = 1;
                }
            }
        }
    }
    Hoppings = B;
}

template<int N>
void Geometry<N>::hcStrainedNanoribbon(int Ny, double Delta)
{   //  Ny = width of the ribbon
    B = Eigen::Matrix<double, N, N>::Zero();
    int Nx = N / Ny / 2;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  SUBLATTICE A
            if (y == Ny - 1)
            {
                if (x == 0)
                {
                    B(Nx * y + x, Nx * Ny + Nx * y) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1 - Delta;
                    B(Nx * Ny + Nx * y, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1 - Delta;
                }
                else
                {
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1 - Delta;
                }
            }
            else
            {
                if (x == 0)
                {
                    B(Nx * y + x, Nx * Ny + Nx * y) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * (y + 1) + Nx - 1) = 1;
                    B(Nx * Ny + Nx * y, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * (y + 1) + Nx - 1, Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * ( y + 1 ) + x - 1) = 1;
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * ( y + 1 ) + x - 1, Nx * y + x) = 1;
                }
            }
            //  SUBLATTICE B
            if (y == 0)
            {
                if (x == Nx - 1)
                {
                    B(Nx * Ny + Nx * y + x, Nx * y) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1 - Delta;
                    B(Nx * y, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1 - Delta;
                }
                else
                {
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1 - Delta;
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1 - Delta;
                }
            }
            else
            {
                if (x == Nx - 1)
                {
                    B(Nx * Ny + Nx * y + x, Nx * y) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) ) = 1;
                    B(Nx * y, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * ( y - 1 ) , Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    B(Nx * Ny + Nx * y + x, Nx * y + x) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1 - Delta;
                    B(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) + x + 1) = 1;
                    B(Nx * y + x, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1 - Delta;
                    B(Nx * ( y - 1 ) + x + 1, Nx * Ny + Nx * y + x) = 1;
                }
            }
        }
    }
    Hoppings = B;
}

template<int N>
void Geometry<N>::setParamsThreeOrbitalTB(double threeOrbitalTBparameters[8])
{
    t0 = fabs(threeOrbitalTBparameters[2]);
    e1 = threeOrbitalTBparameters[0] / t0;
    e2 = threeOrbitalTBparameters[1] / t0;
    t1 = threeOrbitalTBparameters[3] / t0;
    t2 = threeOrbitalTBparameters[4] / t0;
    t11 = threeOrbitalTBparameters[5] / t0;
    t12 = threeOrbitalTBparameters[6] / t0;
    t22 = threeOrbitalTBparameters[7] / t0;

    t0 = -1;

    Ezero << e1,   0.,   0.,
             0.,   e2,   0.,
             0.,   0.,   e2;

    Eone << t0,      t1,     t2,
            -t1,    t11,    t12,
            t2,     -t12,   t22;

    Efour << t0,    -t1,     t2,
             t1,    t11,   -t12,
             t2,    t12,   t22;

    Etwo << t0,                 t1 / 2 - sqrt(3) / 2 * t2,          - sqrt(3) / 2 * t1 - t2 / 2,
-t1 / 2 - sqrt(3) / 2 * t2,   ( t11 + 3 * t22 ) / 4,              sqrt(3) / 4 * ( t22 - t11 ) - t12,
sqrt(3) / 2 * t1 - t2 / 2,      sqrt(3) / 4 * ( t22 - t11 ) + t12,  ( 3 * t11 + t22 ) / 4;

    Efive << t0,                -t1 / 2 - sqrt(3) / 2 * t2,             sqrt(3) / 2 * t1 - t2 / 2,
t1 / 2 - sqrt(3) / 2 * t2,   ( t11 + 3 * t22 ) / 4,              sqrt(3) / 4 * ( t22 - t11 ) + t12,
-sqrt(3) / 2 * t1 - t2 / 2,      sqrt(3) / 4 * ( t22 - t11 ) - t12,  ( 3 * t11 + t22 ) / 4;

    Ethree << t0,                -t1 / 2 + sqrt(3) / 2 * t2,             -sqrt(3) / 2 * t1 - t2 / 2,
t1 / 2 + sqrt(3) / 2 * t2,   ( t11 + 3 * t22 ) / 4,              -sqrt(3) / 4 * ( t22 - t11 ) + t12,
sqrt(3) / 2 * t1 - t2 / 2,      -sqrt(3) / 4 * ( t22 - t11 ) - t12,  ( 3 * t11 + t22 ) / 4;

    Esix << t0,                t1 / 2 + sqrt(3) / 2 * t2,             sqrt(3) / 2 * t1 - t2 / 2,
-t1 / 2 + sqrt(3) / 2 * t2,   ( t11 + 3 * t22 ) / 4,              -sqrt(3) / 4 * ( t22 - t11 ) - t12,
-sqrt(3) / 2 * t1 - t2 / 2,      -sqrt(3) / 4 * ( t22 - t11 ) + t12,  ( 3 * t11 + t22 ) / 4;

}

template<int N>
void Geometry<N>::tmdPBC()
{
    B = Eigen::Matrix<double, N, N>::Zero();
    int Nx = sqrt(N / NORB);
    int Ny = Nx;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  Diagonal term

            B.block(
            ( Nx * y + x ) * NORB , ( Nx * y + x ) * NORB , NORB , NORB )
            = Ezero;

            //  E1

            B.block(
            ( Nx * y + x ) * NORB , ( Nx * y + ( (x + 1) % Nx ) ) * NORB , NORB , NORB )
            = Eone;

            //  E4

            B.block(
            ( Nx * y + ( (x + 1) % Nx ) ) * NORB , ( Nx * y + x ) * NORB , NORB , NORB )
            = Efour;

            if ( y == 0 )
            {
                B.block(
                x * NORB , ( Nx + x ) * NORB , NORB , NORB )
                = Esix;

                B.block(
                ( Nx + x ) * NORB , x * NORB , NORB , NORB )
                = Ethree;

                B.block(
                x * NORB, (Nx * (Ny - 1) + x)*NORB, NORB, NORB )
                = Ethree;

                B.block(
                (Nx * (Ny - 1) + x)*NORB, x * NORB, NORB, NORB )
                = Esix;

                B.block(
                x * NORB, ( Nx * (Ny - 1) + (x + 1) % Nx ) * NORB, NORB, NORB)
                = Etwo;

                B.block(
                ( Nx * (Ny - 1) + (x + 1) % Nx ) * NORB, x * NORB, NORB, NORB)
                = Efive;

                if ( x == 0 )
                {
                    B.block(
                    x * NORB , ( Nx + (Nx - 1 ) ) * NORB , NORB , NORB )
                    = Efive;

                    B.block(
                    ( Nx + (Nx - 1 ) ) * NORB, x * NORB , NORB , NORB )
                    = Etwo;
                }
                else
                {
                    B.block(
                    x * NORB, ( Nx + ( x - 1 ) ) * NORB, NORB, NORB )
                    = Efive;

                    B.block(
                    ( Nx + ( x - 1 ) ) * NORB, x * NORB, NORB, NORB )
                    = Etwo;
                }
            }
            else
            {
                if ( y == Ny - 1 )
                {
                    B.block(
                    ( Nx * ( Ny - 1 ) + x ) * NORB,
                    ( Nx * ( Ny - 2 ) + ( x + 1 ) % Nx ) * NORB, NORB, NORB )
                    = Etwo;

                    B.block(
                    ( Nx * ( Ny - 2 ) + ( x + 1 ) % Nx ) * NORB,
                    ( Nx * ( Ny - 1 ) + x ) * NORB, NORB, NORB )
                    = Efive;

                    B.block(
                    ( Nx * ( Ny - 1 ) + x ) * NORB,
                    ( Nx * ( Ny - 2 ) + x ) * NORB, NORB, NORB )
                    = Ethree;

                    B.block(
                    ( Nx * ( Ny - 2 ) + x ) * NORB,
                    ( Nx * ( Ny - 1 ) + x ) * NORB, NORB, NORB )
                    = Esix;

                    B.block(
                    (Nx * (Ny - 1) + x) * NORB,
                    x * NORB, NORB, NORB )
                    = Esix;

                    B.block(
                    x * NORB,
                    (Nx * (Ny - 1) + x) * NORB,
                    NORB, NORB )
                    = Ethree;

                    if (x == 0)
                    {
                        B.block(
                        x * NORB, (Nx - 1) * NORB, NORB, NORB)
                        = Efive;

                        B.block(
                        (Nx - 1) * NORB, x * NORB, NORB, NORB)
                        = Etwo;
                    }
                    else
                    {
                        B.block(
                        ( Nx * y + x )* NORB, (x-1)*NORB, NORB, NORB)
                        = Efive;

                        B.block(
                        (x-1)*NORB, ( Nx * y + x )* NORB, NORB, NORB)
                        = Etwo;
                    }

                }
                else
                {
                    B.block(
                    ( Nx * y + x ) * NORB,
                    ( Nx * ( y - 1 ) + ( x + 1 ) % Nx ) * NORB, NORB, NORB )
                    = Etwo;

                    B.block(
                    ( Nx * ( y - 1 ) + ( x + 1 ) % Nx ) * NORB,
                    ( Nx * y + x ) * NORB, NORB, NORB )
                    = Efive;

                    B.block(
                    ( Nx * y + x ) * NORB,
                    ( Nx * ( y - 1 ) + x ) * NORB, NORB, NORB )
                    = Ethree;

                    B.block(
                    ( Nx * ( y - 1 ) + x ) * NORB,
                    ( Nx * y + x ) * NORB, NORB, NORB )
                    = Esix;

                    if ( x == 0 )
                    {
                        B.block(
                        ( Nx * y + x ) * NORB,
                        ( Nx * ( y + 1 ) + Nx - 1 ) * NORB, NORB, NORB )
                        = Efive;

                        B.block(
                        ( Nx * ( y + 1 ) + Nx - 1 ) * NORB,
                        ( Nx * y + x ) * NORB, NORB, NORB )
                        = Etwo;
                    }
                    else
                    {
                        B.block(
                        ( Nx * y + x ) * NORB,
                        ( Nx * ( y + 1 ) + x - 1 ) * NORB, NORB, NORB )
                        = Efive;

                        B.block(
                        ( Nx * ( y + 1 ) + x - 1 ) * NORB,
                        ( Nx * y + x ) * NORB, NORB, NORB )
                        = Etwo;
                    }
                }
            }

        }
    }
    B = - B;
    Hoppings = B;
}


template<int N>
void Geometry<N>::tmdNanoribbon(int Ny)
{
    B = Eigen::Matrix<double, N, N>::Zero();
    int Nx = N / Ny / NORB;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  Diagonal term

            B.block(
            ( Nx * y + x ) * NORB , ( Nx * y + x ) * NORB , NORB , NORB )
            = Ezero;

            //  E1

            B.block(
            ( Nx * y + x ) * NORB , ( Nx * y + ( (x + 1) % Nx ) ) * NORB , NORB , NORB )
            = Eone;

            //  E4

            B.block(
            ( Nx * y + ( (x + 1) % Nx ) ) * NORB , ( Nx * y + x ) * NORB , NORB , NORB )
            = Efour;

            if ( y == 0 )
            {
                B.block(
                x * NORB , ( Nx + x ) * NORB , NORB , NORB )
                = Esix;

                B.block(
                ( Nx + x ) * NORB , x * NORB , NORB , NORB )
                = Ethree;

                if ( x == 0 )
                {
                    B.block(
                    x * NORB , ( Nx + (Nx - 1 ) ) * NORB , NORB , NORB )
                    = Efive;

                    B.block(
                    ( Nx + (Nx - 1 ) ) * NORB, x * NORB , NORB , NORB )
                    = Etwo;
                }
                else
                {
                    B.block(
                    x * NORB, ( Nx + ( x - 1 ) ) * NORB, NORB, NORB )
                    = Efive;

                    B.block(
                    ( Nx + ( x - 1 ) ) * NORB, x * NORB, NORB, NORB )
                    = Etwo;
                }
            }
            else
            {
                if ( y == Ny - 1 )
                {
                    B.block(
                    ( Nx * ( Ny - 1 ) + x ) * NORB,
                    ( Nx * ( Ny - 2 ) + ( x + 1 ) % Nx ) * NORB, NORB, NORB )
                    = Etwo;

                    B.block(
                    ( Nx * ( Ny - 2 ) + ( x + 1 ) % Nx ) * NORB,
                    ( Nx * ( Ny - 1 ) + x ) * NORB, NORB, NORB )
                    = Efive;

                    B.block(
                    ( Nx * ( Ny - 1 ) + x ) * NORB,
                    ( Nx * ( Ny - 2 ) + x ) * NORB, NORB, NORB )
                    = Ethree;

                    B.block(
                    ( Nx * ( Ny - 2 ) + x ) * NORB,
                    ( Nx * ( Ny - 1 ) + x ) * NORB, NORB, NORB )
                    = Esix;
                }
                else
                {
                    B.block(
                    ( Nx * y + x ) * NORB,
                    ( Nx * ( y - 1 ) + ( x + 1 ) % Nx ) * NORB, NORB, NORB )
                    = Etwo;

                    B.block(
                    ( Nx * ( y - 1 ) + ( x + 1 ) % Nx ) * NORB,
                    ( Nx * y + x ) * NORB, NORB, NORB )
                    = Efive;

                    B.block(
                    ( Nx * y + x ) * NORB,
                    ( Nx * ( y - 1 ) + x ) * NORB, NORB, NORB )
                    = Ethree;

                    B.block(
                    ( Nx * ( y - 1 ) + x ) * NORB,
                    ( Nx * y + x ) * NORB, NORB, NORB )
                    = Esix;

                    if ( x == 0 )
                    {
                        B.block(
                        ( Nx * y + x ) * NORB,
                        ( Nx * ( y + 1 ) + Nx - 1 ) * NORB, NORB, NORB )
                        = Efive;

                        B.block(
                        ( Nx * ( y + 1 ) + Nx - 1 ) * NORB,
                        ( Nx * y + x ) * NORB, NORB, NORB )
                        = Etwo;
                    }
                    else
                    {
                        B.block(
                        ( Nx * y + x ) * NORB,
                        ( Nx * ( y + 1 ) + x - 1 ) * NORB, NORB, NORB )
                        = Efive;

                        B.block(
                        ( Nx * ( y + 1 ) + x - 1 ) * NORB,
                        ( Nx * y + x ) * NORB, NORB, NORB )
                        = Etwo;
                    }
                }
            }
        }
    }
    B = - B;
    Hoppings = B;
}

template<int N>
void Geometry<N>::tmdNanoribbonStrained(int Ny, double Delta)
{
    B = Eigen::Matrix<double, N, N>::Zero();
    int Nx = N / Ny / NORB;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  Diagonal term

            B.block(
            ( Nx * y + x ) * NORB , ( Nx * y + x ) * NORB , NORB , NORB )
            = Ezero;

            //  E1

            B.block(
            ( Nx * y + x ) * NORB , ( Nx * y + ( (x + 1) % Nx ) ) * NORB , NORB , NORB )
            = Eone - Delta * Eigen::Matrix<double, NORB, NORB>::Identity();

            //  E4

            B.block(
            ( Nx * y + ( (x + 1) % Nx ) ) * NORB , ( Nx * y + x ) * NORB , NORB , NORB )
            = Efour - Delta * Eigen::Matrix<double, NORB, NORB>::Identity();

            if ( y == 0 )
            {
                B.block(
                x * NORB , ( Nx + x ) * NORB , NORB , NORB )
                = Esix - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                B.block(
                ( Nx + x ) * NORB , x * NORB , NORB , NORB )
                = Ethree - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                if ( x == 0 )
                {
                    B.block(
                    x * NORB , ( Nx + (Nx - 1 ) ) * NORB , NORB , NORB )
                    = Efive - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx + (Nx - 1 ) ) * NORB, x * NORB , NORB , NORB )
                    = Etwo - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();
                }
                else
                {
                    B.block(
                    x * NORB, ( Nx + ( x - 1 ) ) * NORB, NORB, NORB )
                    = Efive - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx + ( x - 1 ) ) * NORB, x * NORB, NORB, NORB )
                    = Etwo - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();
                }
            }
            else
            {
                if ( y == Ny - 1 )
                {
                    B.block(
                    ( Nx * ( Ny - 1 ) + x ) * NORB,
                    ( Nx * ( Ny - 2 ) + ( x + 1 ) % Nx ) * NORB, NORB, NORB )
                    = Etwo - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx * ( Ny - 2 ) + ( x + 1 ) % Nx ) * NORB,
                    ( Nx * ( Ny - 1 ) + x ) * NORB, NORB, NORB )
                    = Efive - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx * ( Ny - 1 ) + x ) * NORB,
                    ( Nx * ( Ny - 2 ) + x ) * NORB, NORB, NORB )
                    = Ethree - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx * ( Ny - 2 ) + x ) * NORB,
                    ( Nx * ( Ny - 1 ) + x ) * NORB, NORB, NORB )
                    = Esix - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();
                }
                else
                {
                    B.block(
                    ( Nx * y + x ) * NORB,
                    ( Nx * ( y - 1 ) + ( x + 1 ) % Nx ) * NORB, NORB, NORB )
                    = Etwo - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx * ( y - 1 ) + ( x + 1 ) % Nx ) * NORB,
                    ( Nx * y + x ) * NORB, NORB, NORB )
                    = Efive - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx * y + x ) * NORB,
                    ( Nx * ( y - 1 ) + x ) * NORB, NORB, NORB )
                    = Ethree - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    B.block(
                    ( Nx * ( y - 1 ) + x ) * NORB,
                    ( Nx * y + x ) * NORB, NORB, NORB )
                    = Esix - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                    if ( x == 0 )
                    {
                        B.block(
                        ( Nx * y + x ) * NORB,
                        ( Nx * ( y + 1 ) + Nx - 1 ) * NORB, NORB, NORB )
                        = Efive - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                        B.block(
                        ( Nx * ( y + 1 ) + Nx - 1 ) * NORB,
                        ( Nx * y + x ) * NORB, NORB, NORB )
                        = Etwo - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();
                    }
                    else
                    {
                        B.block(
                        ( Nx * y + x ) * NORB,
                        ( Nx * ( y + 1 ) + x - 1 ) * NORB, NORB, NORB )
                        = Efive - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();

                        B.block(
                        ( Nx * ( y + 1 ) + x - 1 ) * NORB,
                        ( Nx * y + x ) * NORB, NORB, NORB )
                        = Etwo - Delta / 4 * Eigen::Matrix<double, NORB, NORB>::Identity();
                    }
                }
            }
        }
    }
    B = - B;
    Hoppings = B;
}

template<int N>
Eigen::MatrixXd Geometry<N>::matrix()
{
    return Hoppings;
}

template<int N>
double Geometry<N>::get(int x, int y)
{
    return Hoppings(x, y);
}

template<int N>
Eigen::MatrixXd Geometry<N>::BpreFactor(double dt, double mu)
{
    return exp(dt * mu) * B;
}

template<int L, int N>
class Configuration
{
    Eigen::MatrixXd HSfield;
public:
    void genHsMatrix();
    Eigen::MatrixXd matrix();
    double get(int x, int y);
    void flip(int l, int i);
    Configuration() : HSfield(L, N) {
    };
};

template<int L, int N>
void Configuration<L, N>::genHsMatrix()
{
    //  Generate the HS field matrix
    int l;
    int i;
    HSfield = Eigen::Matrix<double, L, N>::Random(L,N);
    for (l = 0; l < L; l++)
    {
        for (i = 0; i < N; i++)
        {
            if ( HSfield(l, i) < 0 )
            {
                HSfield(l, i) = -1;
            }
            else
            {
                HSfield(l, i) = 1;
            }
        }
    }
}

template<int L, int N>
Eigen::MatrixXd Configuration<L, N>::matrix()
{
    return HSfield;
}

template<int L, int N>
double Configuration<L, N>::get(int x, int y)
{
    return HSfield(x, y);
}

template<int L, int N>
void Configuration<L, N>::flip(int l, int i)
{
    HSfield(l, i) *= -1;
}

template<int N, int L>
class OneParticlePropagators
{
    Eigen::MatrixXd B[L];
public:
    void fillMatrices(bool spin, double nu, Eigen::MatrixXd h, Eigen::MatrixXd BpreFactor);
    Eigen::MatrixXd matrix(int l);
    Eigen::MatrixXd * list();
    void update(int l, int i, double alpha);
};

template<int N, int L>
void OneParticlePropagators<N, L>::fillMatrices(bool spin, double nu,
   Eigen::MatrixXd h, Eigen::MatrixXd BpreFactor)
{
    int i;
    int l;
    for (l = 0; l < L; l++)
    {
        B[l] = Eigen::Matrix<double, N, N>::Zero();

        if (spin == true)
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = exp ( nu * h(l, i) ) ;
            }
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = exp ( - 1. * nu * h(l, i) );
            }
        }
        B[l] = BpreFactor * B[l] ;
    }
}

template<int N, int L>
void OneParticlePropagators<N, L>::update(int l, int i, double alpha)
{
    B[l].col(i) *= ( alpha + 1 );
}

template<int N, int L>
Eigen::MatrixXd OneParticlePropagators<N, L>::matrix(int l)
{
    return B[l];
}

template<int N, int L>
Eigen::MatrixXd * OneParticlePropagators<N, L>::list()
{
    return B;
}

#endif /* matrixgen_h */
