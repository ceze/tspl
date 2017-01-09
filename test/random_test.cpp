/*****************************************************************************
 *                               random_test.cpp
 *
 * Random number generator testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <iomanip>
#include <random.h>
#include <statistics.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 100;


int main()
{
    int seed = 37;
    int low = 0, high = 100;
    Type mu = 0.0, sigma = 1.0;
    Type beta = 2.0;
    Type lambda = 4.0;
    Type p = 0.5;

    Vector<Type> rs(N);
    Vector<int> irs(N);
    Vector<double> tmp(N);

//    Random rand(seed);
//    cout << "Random Number Generator:" << endl;
//    for( int i=0; i<N; ++i )
//    {
//        cout << rand.random() << "\t";
//        if( !((i+1)%4) )
//            cout << endl;
//    }
//    cout << endl << endl;

    cout << "Uniform distribution: U ~ ( " << low << ", " << high << " )" << endl;
    for( int i=0; i<N; ++i )
        cout << randu( seed, low, high ) << "\t";
    cout << endl;
    irs = randu( seed, low, high, N );
    cout << "mean     (theoretical   generated):   "
         << mean(irs) << "\t" << (high-low)/2 << endl;
    cout << "variance (theoretical   generated):   "
         << var(irs) << "\t" << (high-low)*(high-low)/12
         << endl << endl << endl;

    cout << "Normal distribution: N ~ ( " << mu << ", " << sigma << " )" << endl;
    cout << setiosflags(ios::fixed) << setprecision(4);
    for( int i=0; i<N; ++i )
        cout << randn( seed, mu, sigma ) << "\t\t";
    cout << endl;
    rs = randn( seed, mu, sigma, N );
    cout << "mean     (theoretical   generated):   "
         << mean(rs) << "\t" << mu << endl;
    cout << "variance (theoretical   generated):   "
         << var(rs) << "\t" << sigma*sigma
         << endl << endl << endl;

    cout << "Exponential distribution: E ~ ( " << beta << " )" << endl;
    cout << setiosflags(ios::fixed) << setprecision(4);
    for( int i=0; i<N; ++i )
        cout << rande( seed, beta ) << "\t\t";
    cout << endl;
    rs = rande( seed, beta, N );
    cout << "mean     (theoretical   generated):   "
         << mean(rs) << "\t" << beta << endl;
    cout << "variance (theoretical   generated):   "
         << var(rs) << "\t" << beta*beta
         << endl << endl << endl;

    cout << "Rayleigh distribution: R ~ ( " << sigma << " )" << endl;
    cout << setiosflags(ios::fixed) << setprecision(4);
    for( int i=0; i<N; ++i )
        cout << randr( seed, sigma ) << "\t\t";
    cout << endl;
    rs = randr( seed, sigma, N );
    cout << "mean     (theoretical   generated):   "
         << mean(rs) << "\t" << sigma*sqrt(PI/2.0) << endl;
    cout << "variance (theoretical   generated):   "
         << var(rs) << "\t" << (2-PI/2)*sigma*sigma
         << endl << endl << endl;

    cout << "Poisson distribution: B ~ ( " << p << " )" << endl;
    for( int i=0; i<N; ++i )
        cout << randp( seed, lambda ) << "\t\t";
    cout << endl;
    irs = randp( seed, lambda, N );
    for( int i=0; i<N; ++i )
        tmp[i] = irs[i];
    cout << "mean     (theoretical   generated):   "
         << mean(tmp) << "\t" << lambda << endl;
    cout << "variance (theoretical   generated):   "
         << var(tmp) << "\t" << lambda
         << endl << endl << endl;

    cout << "Bernoulli distribution: B ~ ( " << p << " )" << endl;
    for( int i=0; i<N; ++i )
        cout << randb( seed, p ) << "\t\t";
    cout << endl;
    irs = randb( seed, p, N );
    for( int i=0; i<N; ++i )
        tmp[i] = irs[i];
    cout << "mean     (theoretical   generated):   "
         << mean(tmp) << "\t" << p << endl;
    cout << "variance (theoretical   generated):   "
         << var(tmp) << "\t" << p*(1-p)
         << endl << endl;

    return 0;
}
