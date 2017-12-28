/*
// Copyright (C) 2016, HydroComplexity Group
// All rights reserved.
//
// Distributed Hydrologicc and Regional Analysis (DHARA) Model
// DHARA model is made available as a restricted, non-exclusive, 
// non-transferable license for education and research purpose only, 
// and not for commercial use. See the LICENSE.txt for more details.
//
// Author: levuvietphong@gmail.com (Phong Le)
*/

#include "main.h"
#include "_bool.h"
#include "_defs.h"


/*
 * FUNC: MaxFind()
 * ---------------
 * 
 * Returns the maximum value between 2 numbers in which they must be in 
 * the same data type (double, float, int, etc).
 *
 *      a       : first number
 *      b       : first number
 *      returns : maximum of a and b
 */
double MaxFind ( const double a, const double b )
{
    return (a < b) ? b : a;   // Ternary conditional comparison
}



/*
 * FUNC: QuadraticSolver()
 * -----------------------
 * Solve equation x^2 + b x + c = 0 and returns the number of roots
 * The coefficient 'a' is normalized to 1, only two arguments are required
 * 
 *     ce     : a pointer to an array of equation coefficients.
 *     roots  : a pointer to roots to be stored
 *     return : number of solution
 */
unsigned int QuadraticSolver (double *ce,double *roots )
{
    double discriminant;

    // Calculate the discriminant
    discriminant = 0.25 * ce[0] * ce[0] - ce[1];

    // discriminant is negative, no solution found
    if (discriminant < 0.)
        return 0;

    // discriminant is non-negative, two solutions found (can be identical)
    roots[0]= -0.5*ce[0] + sqrt(discriminant);
    roots[1]= -ce[0] - roots[0];
    return 2;
}



/*
 * FUNC: CubicSolver()
 * -----------------------
 * Solve equation x^3 + b x^2 + c x + d = 0 and returns the number of roots
 * The coefficient 'a' is normalized to 1, only three arguments are required
 * 
 *     ce     : a pointer to an array of equation coefficients.
 *     roots  : a pointer to roots to be stored
 *     return : number of solution
 */
unsigned int CubicSolver (double *ce, double *roots )
{

    // Transform the equation so that b = 0.
    // We will get the depressed cubic x^3 + p t + q = 0
    // This is a Tschirnhaus transformation and x is shifted x = t - b/(3a)
    unsigned int ret = 0;
    double shift = (1./3)*ce[0];
    double p = ce[1] - shift*ce[0];
    double q = ce[0] * ( (2./27)*ce[0]*ce[0]-(1./3)*ce[1] ) + ce[2];
    double discriminant = (1./27)*p*p*p+(1./4)*q*q;

    // If value of p is closed to 0, one solution found
    if ( fabs(p)< 1.0e-75)
    {
        ret = 1;
        *roots = (q>0) ? -pow(q,(1./3)) : pow(-q,(1./3));
        *roots -= shift;
        return ret;
    }

    // discriminant is positive, one solution found
    if(discriminant > 0)
    {
        double ce2[2] = {q, -1./27*p*p*p},u3[2];
        ret = QuadraticSolver(ce2,u3);
        if (! ret )
        {   //should not happen
            cerr<<"CubicSolver::Error CubicSolver("<<ce[0]<<' '<<ce[1]<<' '<<ce[2]<<")\n";
        }

        ret = 1;
        double u, v;
        u = (q<=0) ? pow(u3[0], 1./3): -pow(-u3[1],1./3);
        v = (-1./3)*p/u;
        *roots = u + v - shift;
        return ret;
    }

    // discriminant is non-positive, three solutions found
    ret=3;
    complex<double> u(q,0),rt[3];
    u=pow(-0.5*u-sqrt(0.25*u*u+p*p*p/27),1.0/3.0);
    rt[0]=u-p/(3.*u)-shift;
    complex<double> w(-0.5,sqrt(3.)/2);
    rt[1]=u*w-p/(3.*u*w)-shift;
    rt[2]=u/w-p*w/(3.*u)-shift;

    roots[0]=rt[0].real();
    roots[1]=rt[1].real();
    roots[2]=rt[2].real();
    return ret;
}


/*
 * FUNC: quarticSolver()
 * -----------------------
 * Solve equation x^4 + ce[0] x^3 + ce[1] x^2 + ce[2] x + ce[3] = 0 and returns the 
 * number of roots. The coefficient 'a' is normalized to 1, only three arguments are required
 * 
 *     ce     : a pointer to an array of equation coefficients.
 *     roots  : a pointer to roots to be stored
 *     return : number of solution
 */
unsigned int quarticSolver(double * ce, double *roots)
{
    unsigned int ret=0;
    double shift=0.25*ce[0];
    double shift2=shift*shift;
    double a2=ce[0]*ce[0];
    double p= ce[1] - (3./8)*a2;
    double q= ce[2] + ce[0]*((1./8)*a2 - 0.5*ce[1]);
    double r= ce[3] - shift*ce[2] + (ce[1] - 3.*shift2)*shift2;

    if (fabs(q) <= 1.0e-75)
    {   // Biquadratic equations
        double discriminant= 0.25*p*p -r;
        if (discriminant < 0.) {
            return 0;
        }

        double t2[2];
        t2[0]=-0.5*p-sqrt(discriminant);
        t2[1]= -p - t2[0];

        if ( t2[0] >= 0. )
        {   // four real roots
            roots[0]=sqrt(t2[0])-shift;
            roots[1]= -sqrt(t2[0])-shift;
            roots[2]=sqrt(t2[1])-shift;
            roots[3]= -sqrt(t2[1])-shift;
            return 4;
        }

        if ( t2[1] >= 0.) { // two real roots
            roots[0]=sqrt(t2[1])-shift;
            roots[1]= -roots[0]-shift;
            return 2;
        }

        return 0;
    }

    if ( fabs(r)< 1.0e-75 )
    {
        double cubic[3]= {0.,p,q};
        roots[0]=0.;
        ret=1+CubicSolver(cubic,roots+1);
        for(unsigned int i=0; i<ret; i++) roots[i] -= shift;
        return ret;
    }

    double cubic[3]= {2.*p,p*p-4.*r,-q*q},croots[3];
    ret = CubicSolver(cubic,croots);
    if (ret==1)
    {   //one real root from cubic
        if (croots[0]< 0.)
        {//this should not happen
            cerr<<"Quartic Error:: Found one real root for cubic, but negative\n";
            return 0;
        }

        double sqrtz0=sqrt(croots[0]);
        double ce2[2];
        ce2[0]= -sqrtz0;
        ce2[1]=0.5*(p+croots[0])+0.5*q/sqrtz0;
        ret=QuadraticSolver(ce2,roots);

        if (! ret )
        {
            ce2[0] = sqrtz0;
            ce2[1] = 0.5*(p+croots[0])-0.5*q/sqrtz0;
            ret    = QuadraticSolver(ce2,roots);
        }

        ret=2;
        for(unsigned int i=0; i<ret; i++)
            roots[i] -= shift;

        return ret;
    }

    if ( croots[0]> 0. && croots[1] > 0. )
    {
        double sqrtz0=sqrt(croots[0]);
        double ce2[2];
        ce2[0] = -sqrtz0;
        ce2[1] = 0.5*(p+croots[0])+0.5*q/sqrtz0;
        ret    = QuadraticSolver(ce2,roots);
        ce2[0] = sqrtz0;
        ce2[1] = 0.5*(p+croots[0])-0.5*q/sqrtz0;
        ret    = QuadraticSolver(ce2,roots+2);
        ret    = 4;

        for(unsigned int i=0; i<ret; i++)
            roots[i] -= shift;

        return ret;
    }

    return 0;
}

// Using bisection, find the root of a function func known to lie between x1 and x2.
// The root, returned as rtbis, will be refined until its accuracy is =/xacc.
double rtbis(double (*func)(double, double, double, double, double, double, double, double, double,
             double, double, double, double, double),
             double Rabs, double Ta1, double Ts1, double ea1, double pa1, double U1, double z1,
             double dzs1, double psis1_MPa, double vonk, double z0, double TC1, double epss,
             double x1, double x2, double xacc)
{
    int j;
    double dx,f,fmid,xmid,rtb;

    f=(*func)(x1, Rabs, Ta1, Ts1, ea1, pa1, U1, z1, dzs1, psis1_MPa, vonk, z0, TC1, epss);
    fmid=(*func)(x2, Rabs, Ta1, Ts1, ea1, pa1, U1, z1, dzs1, psis1_MPa, vonk, z0, TC1, epss);

    if (f*fmid >= 0.0)
        printf("Root must be bracketed for bisection in rtbis\n");

    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    for (j=1; j<=Jmax_Bisec; j++)
    {
        fmid=(*func)(xmid=rtb+(dx *= 0.5), Rabs, Ta1, Ts1, ea1, pa1, U1, z1, dzs1, psis1_MPa,
                vonk, z0, TC1, epss);

        if (fmid <= 0.0)
            rtb=xmid;

        if (fabs(dx) < xacc || fmid == 0.0)
            return rtb;
    }

    printf("Too many bisections in rtbis\n");
    return 0.0;
}


// --------------------------------------------------------------------
// MaxValue()
//    Returns the maximum value of an array
// --------------------------------------------------------------------
double MaxValue(VectorXd& Array_in, int SIZE)
{
    int ind;                  // Index for loop
    double MaxValue;          // Maximum value that will be returned.

    // MaxValue is initially set to the first element
    MaxValue = abs(Array_in[0]);

    // Loop run through the all elements in the array
    for (ind = 1; ind < SIZE; ++ind)
    {
        // If an element is larger than MaxValue, update MaxValue to this element
        if ( abs(Array_in[ind]) > MaxValue )
            MaxValue = abs(Array_in[ind]);
    }

    return MaxValue;          // Return the maximum value
}

// Sub diagonal of the matrix
// Main diagonal of the matrix
// Supper diagonal of the matrix
// Right hand side of the linear system
// Vector solution of the linear system
// Size of the linear system
void ThomasAlgorithm (VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& d, VectorXd& f, int N)
{
    double m; // temporaty variable used for substitution
    vector<double> c_star(N, 0.0);
    vector<double> d_star(N, 0.0);

    // Modify the first-row coefficients
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    // Forward sweep and modification
    for (int i = 1; i < N; i++) 
    {
        m = 1.0 / (b[i] - a[i] * c_star[i-1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
    }

    // Backward substitution
    f[N-1] = d_star[N-1];
    for (int i = N-2; i >= 0; i--) 
        f[i] = d_star[i] - c_star[i] * f[i+1];
}

/*
void CanopyGrid(CanopyClass *canopies, VectorXd& LADnorm, VerticalCanopyClass *vertcanopies)
{
    int nl_can = canopies->nl_can;
    double hcan = canopies->hcan;

    for (int i=0; i< nl_can; i++)
    {
        vertcanopies->zhc[i] = (i+1) * (hcan / nl_can);
    }

    vertcanopies->dzc = 0.5 * vertcanopies->zhc[0];

    for (int i=0; i< nl_can; i++)
    {
        vertcanopies->znc[i] = vertcanopies->zhc[i] - vertcanopies->dzc;
    }

}

void SoilGrid( ProjectClass *project, VerticalSoilClass *vertsoils, SoilClass *Soil)
{
    int nl_soil = soils->nl_soil;
    double rhod;
    double scalek = Soil->scalek;
    double smpmin = Soil->smpmin;
    double sand = Soil->sand;
    double clay = Soil->clay;

    vertsoils->zhs[0] = project->dz_meter;
    for (int i=1; i<nl_soil; i++)
    {
        vertsoils->zhs[i] = vertsoils->zhs[i] + project->dz_meter;
    }

    for (int i=0; i<nl_soil; i++)
    {
        vertsoils->zns[i] = vertsoils->zhs[i] - 0.5 * project->dz_meter;
        vertsoils->dzs[i] = project->dz_meter;
        vertsoils->sand[i] = sand;
        vertsoils->clay[i] = clay;
        vertsoils->HC_sol[i] = (2.128*sand + 2.385*clay) / (sand+clay) * 1e6;
        vertsoils->porosity[i] = Soil->theta_S;
        vertsoils->psi0[i] = -10 * pow( 10, 1.88-0.0131*sand );
        vertsoils->bsw[i] = 2.91 + 0.159*clay;
        vertsoils->TK_sol[i] = (8.80*sand + 2.92*clay) / (sand+clay);
        rhod = 2700*(1 - vertsoils->porosity[i]);
        vertsoils->TK_dry[i]  = (0.135*rhod + 64.7) / (2700 - 0.947*rhod);
        vertsoils->HKsat[i] = Ksat[i];
        vertsoils->theta_dry[i] = Soil->theta_R;
    }
}
*/

void PH_Dist(CanopyClass *canopies, VerticalCanopyClass *vertcanopies,
    PhotosynthesisClass *photosynthesis, Ref<VectorXd> Vz)
{
    int nl_can = canopies->nl_can;
    double kn = photosynthesis->kn_canopy;
    double *LAIcum = new double[nl_can+1];

    // Copy array structs to eigen vectors
    VectorXd LAIz = Map<VectorXd>(vertcanopies->LAIz, nl_can);

    LAIcum[0] = 0;
    for (int i=1; i<=nl_can; i++)
        LAIcum[i] = LAIcum[i-1] + LAIz[nl_can - i];

    for (int i=0; i<nl_can; i++)
    {
        if (i < nl_can-1)
            Vz[i] = exp(-kn * LAIcum[nl_can-i-1]);
        else
            Vz[i] = 1.0;
    }
}

void First_Order_Closure_U(ForcingClass *forcings, CanopyClass *canopies,
    VerticalCanopyClass *vertcanopies, MicroEnvironmentClass *microenviron,
    Ref<VectorXd> Um, Ref<VectorXd> Km)
{
    double hcan = canopies->hcan;
    double Cd = microenviron->Cd;
    double alph = microenviron->alph;
    double l_mix, Utop, Ulow, dx, eps1, dzc;
    int nl_can, N;

    Utop   = forcings->u;
    dzc    = vertcanopies->dzc;
    nl_can = canopies->nl_can;
    N      = nl_can + 1;
    eps1   = 0.5;

    VectorXd z(N), y(N), dU(N), tau(N), co(N);
    VectorXd a1(N), a2(N), a3(N), aa(N), bb(N), cc(N), dd(N);
    VectorXd upd(N), dia(N), lod(N), U(N), Un(N), LAD(N);

    // Copy array structs to eigen vectors
    VectorXd znc  = Map<VectorXd>(vertcanopies->znc, nl_can);
    VectorXd LADz = Map<VectorXd>(vertcanopies->LADz, nl_can);

    co = VectorXd::Zero(N);

    // Mixing Length
    l_mix = alph * hcan;
    z << 0, znc;
    Ulow = 0;
    U = VectorXd::LinSpaced(N, Ulow, Utop);
    LAD << 0, LADz;

    // Start iterative solution
    double err     = 1e9;
    double epsilon = 1e-2;

    while (err > epsilon) 
    {
        for (int i=1; i<N; i++)
            y[i] = (U[i] - U[i-1]) / dzc;

        y[0] = y[1];
        dU   = y;
        Km   = l_mix * l_mix * y.cwiseAbs();
        tau  = -Km.cwiseProduct(y);

        // Set up coefficients for ODE in Eqn 28
        a1 = -Km;
        for (int i=1; i<N; i++)
            a2[i] = -(Km[i] - Km[i-1]) / dzc;

        a2[0] = a2[1];
        a3    = 0.5 * Cd * LAD.cwiseProduct(U);
        dx    = dzc;

        // Set the elements of the Tri-diagonal Matrix
        upd     = ( a1/(dx*dx) + a2/(2*dx) );
        dia     = ( -a1*2/(dx*dx) + a3 );
        lod     = ( a1/(dx*dx) - a2/(2*dx) );
        co[0]   = Ulow;
        co[N-1] = Utop;
        aa      = lod;
        bb      = dia;
        cc      = upd;
        dd      = co;
        aa[0]   = 0;
        aa[N-1] = 0;
        cc[0]   = 0;
        cc[N-1] = 0;
        bb[0]   = 1.0;
        bb[N-1] = 1.0;

        // Use the Thomas Algorithm to solve the tridiagonal matrix
        ThomasAlgorithm (aa, bb, cc, dd, Un, N);
        err = (Un - U).cwiseAbs().maxCoeff();
        U = ( eps1*Un + (1.0-eps1) * U );
    }
  
    for (int i=1; i<N; i++) 
        Um[i-1] = U[i];
}

void Precip_Interception(ForcingClass *forcings, CanopyClass *canopies,
    VerticalCanopyClass *vertcanopies, VerticalSoilClass *vertsoils,
    RadiationClass *radiation, Ref<VectorXd> Sh2o_prof, Ref<VectorXd> Smaxz,
    Ref<VectorXd> wetfrac, Ref<VectorXd> dryfrac, double *ppt_ground)
{
    int nl_can = canopies->nl_can;
    double ppt = forcings->ppt;
    double Smax = canopies->Smax;
    double Ffact = canopies->Ffact;
    double pptintfact = canopies->pptintfact;
    double clump = radiation->clump;
    double pptint, pptint_pot;

    Map<VectorXd> LAIz(vertcanopies->LAIz, nl_can);
    Map<VectorXd> znc(vertcanopies->znc, nl_can);
    Smaxz = Smax * LAIz;

    // Attenuate precipitation through canopy
    if ( ppt > 0 )
    {
        // Loop from top to bottom of canopy
        for (int zz=nl_can-1; zz >=0; zz--)
        {
            if ( Sh2o_prof(zz) < Smaxz(zz) )
            {
                pptint_pot = ( 1 - exp(-pptintfact*clump*LAIz(zz)) ) * ppt;
                
                if ( (Sh2o_prof(zz) + pptint_pot) <= Smaxz(zz) ) 
                {
                    Sh2o_prof(zz) += pptint_pot;
                    ppt -= pptint_pot;
                } else {
                    pptint = Smaxz(zz) - Sh2o_prof(zz);
                    Sh2o_prof(zz) = Smaxz(zz);
                    ppt -= pptint;
                }
            }
        }
    }

    *ppt_ground = ppt;
    wetfrac     = Ffact * Sh2o_prof.cwiseQuotient(Smaxz);
    dryfrac     = VectorXd::Ones(nl_can) - wetfrac;
}


double Diffuse_Fraction(double zendeg, int doy, double SWdn)
{
    double fdiff, Scs, solelev, So, Sg_o_So, sinbeta, R, K;
    Scs = 1370;
    solelev = (90 - zendeg) * pinum/180;
    So = Scs * ( 1+0.033*cos((360*doy/365.)*pinum/180.) ) * sin(solelev);

    if (zendeg >= 90)
        So = numeric_limits<double>::quiet_NaN();
    
    Sg_o_So = SWdn/So;

    // Fraction Diffuse
    sinbeta = sin(solelev);
    R = 0.847 - 1.61 * sinbeta + 1.04 * sinbeta * sinbeta;
    K = (1.47 - R) / 1.66;

    if ( Sg_o_So <= 0.22 )
    {
        fdiff = 1;
    } else if ( Sg_o_So <= 0.35 ) {
        fdiff =  1 - 6.4*(Sg_o_So - 0.22)*(Sg_o_So - 0.22);
    } else if ( Sg_o_So <=K ) {
        fdiff = 1.47 - 1.66*Sg_o_So;
    } else {
        fdiff = R;
    }

    return fdiff;    
}


void SW_Attenuation(
    Ref<VectorXd> sun_abs, Ref<VectorXd> shade_abs, Ref<VectorXd> fsun,
    Ref<VectorXd> fshade, Ref<VectorXd> diffdn, Ref<VectorXd> diffup,
    Ref<VectorXd> LAIz, double *soil_abs_ptr, double *radabs_tot_ptr,
    double *radlost_ptr, double *radremain_ptr, double beam_top, double diff_top,
    double trans, double refl, double refl_soil, double clump, double Kbm,
    double Kdf, int count, int nl_can)
{
    double soil_inc = 0;
    double radlost, radabs_tot, radremain, soil_abs;
    soil_abs = *soil_abs_ptr;
    radlost = *radlost_ptr;

    //-------------------------------------------
    // DOWNWARD rADIATION . . . . . . . . . . . .
    //-------------------------------------------
    int tind = nl_can-1;
    int bind = 0;

    double taub, taud, beam_int;
    double reflbup    = 0.;
    double refldup    = 0.;
    double diff_int   = 0.;
    double diffuplost = 0.;
    double beam_inc   = beam_top;
    double LAIc       = 0.0;
    diffdn[tind] += diff_top;

    for (int i=tind; i>=0; i--) 
    {
        LAIc = LAIc + LAIz[i];
        fsun[i] = exp(-Kbm * clump * LAIc);
        fshade[i] = 1.0 - fsun[i];

        // . . . Beam Radiaion . . . . . . . . . . .
        if (Kbm>0 && beam_top>0 && count==0)
        {
            taub = exp( -Kbm * clump * LAIz[i] );
            beam_int = beam_inc - taub*beam_inc;
            beam_inc = taub * beam_inc;
            sun_abs[i] = ( 1-refl-trans ) * beam_int;

            // intercepted beam that is transmitted
            if (i == bind) 
            {
                soil_inc = trans * beam_int + beam_inc;
            } else {
                diffdn[i-1] = trans * beam_int;
            }

            // Intercepted beam that is reflected
            if (i < tind) 
            {
                diffup[i+1] = refl * beam_int;
            } else {
                reflbup = refl * beam_int;
            }
        } else {
            reflbup = 0;
        }

        // . . . Diffuse Radiaion . . . . . . . . . . .
        taud = exp( -Kdf * clump * LAIz[i] );
        diff_int = diffdn[i] - taud * diffdn[i];

        // Downward transmission
        if (i==bind)
        {
            soil_inc += trans * diff_int + taud*diffdn[i];
        } else {
            diffdn[i-1] += trans * diff_int + taud*diffdn[i];
        }

        // Upward reflection
        if (i < tind) 
        {
            diffup[i+1] += refl * diff_int;
        } else {
            refldup = refl * diff_int;
        }

        // Absorbed fraction
        sun_abs[i] += ( 1-refl-trans ) * diff_int * fsun[i];
        shade_abs[i] += ( 1-refl-trans ) * diff_int * fshade[i];
    }

    soil_abs += ( 1-refl_soil ) * soil_inc;
    diffup[0] += refl_soil * soil_inc;
    soil_inc = 0;
    diffdn = diffdn * 0.0;

    //-------------------------------------------
    // UPWARD DIFFUSE rADIATION . . . . . . . . .
    //-------------------------------------------
    for (int i=bind; i<=tind; i++)
    {
        taud = exp( -Kdf * clump * LAIz[i] );
        diff_int = diffup[i] * (1 - taud);

        // Upward tranmission
        if (i < tind) 
        {
            diffup[i+1] += trans * diff_int + taud * diffup[i];
        } else {
            diffuplost = trans * diff_int + taud * diffup[i];
        }

        // Downward tranmission
        if (i==bind) 
        {
            soil_inc += refl * diff_int;
        } else {
            diffdn[i-1] += refl * diff_int;
        }

        // Absorbed fraction
        sun_abs[i] += ( 1-refl-trans ) * diff_int * fsun[i];
        shade_abs[i] += ( 1-refl-trans ) * diff_int * fshade[i];
    }


    diffup = diffup * 0;
    soil_abs += ( 1-refl_soil ) * soil_inc;
    diffup[0] = refl_soil * soil_inc;

    // Absorbed radiation
    radabs_tot = sun_abs.sum() + shade_abs.sum() + soil_abs;
    // Lost radiation
    radlost += reflbup + refldup + diffuplost;
    // radiation Remaining in System
    radremain = diffdn.sum() + diffup.sum();

    // Return via pointers
    *radabs_tot_ptr = radabs_tot;
    *radremain_ptr  = radremain;
    *soil_abs_ptr   = soil_abs;
    *radlost_ptr    = radlost;
}


void ShortWaveradiation(
    ForcingClass *forcings, CanopyClass *canopies,
    VerticalCanopyClass *vertcanopies, VerticalSoilClass *vertsoils,
    ConstantClass *constants, RadiationClass *radiation, Ref<VectorXd> fsun,
    Ref<VectorXd> fshade, Ref<VectorXd> LAIsun, Ref<VectorXd> LAIshade,
    Ref<VectorXd> SWabs_sun, Ref<VectorXd> SWabs_shade, Ref<VectorXd> PARabs_sun,
    Ref<VectorXd> PARabs_shade, Ref<VectorXd> NIRabs_sun, Ref<VectorXd> NIRabs_shade,
    Ref<VectorXd> taud, double *SWabs_soil, double *SWout)
{
    int nl_can = canopies->nl_can;
    int doy, count;
    double SWin, zendeg, zenrad, fdiff, Pa, Kdf, Kbm, Wm2toumol, umoltoWm2;
    double PARin, PARtop, PARtop_beam, PARtop_diff;
    double NIRtop, NIRtop_beam, NIRtop_diff;
    double transmiss, xx, clump, trans_PAR, trans_NIR, refl_PAR, refl_NIR, refl_soil;
    double percdiff, beam_top, diff_top;
    double radlost, radabs_tot, radremain, radtot, radin;
    double PARabs_soil, NIRabs_soil, PARout, NIRout;

    // Mapping pointer array to Eigen vectors
    Map<VectorXd> LAIz(vertcanopies->LAIz, nl_can);
    Map<VectorXd> diffdn(vertcanopies->diffdn, nl_can);
    Map<VectorXd> diffup(vertcanopies->diffup, nl_can);

    SWin      = forcings->rg;
    zendeg    = forcings->zen;
    doy       = forcings->doy;
    Pa        = forcings->pa;

    transmiss = radiation->transmiss;
    xx        = radiation->xx;
    clump     = radiation->clump;
    trans_PAR = radiation->trans_PAR;
    refl_PAR  = radiation->refl_PAR;
    trans_NIR = radiation->trans_NIR;
    refl_NIR  = radiation->refl_NIR;
    refl_soil = radiation->refl_soil;
    Kdf       = radiation->Kdf;
    Wm2toumol = constants->Wm2toumol;
    umoltoWm2 = constants->umoltoWm2;

    // Convert Zenith in degree to radian
    zenrad = zendeg * pinum/180.0;

    // No reference supported
    //if (zendeg > 89) 
    //{
    //    SWin = 0;
    //    PARin = 0;
    //}

    PARtop = SWin * 0.45 * Wm2toumol;
    NIRtop = SWin * 0.55;

    // Diffuse fraction
    fdiff = Diffuse_Fraction(zendeg, doy, SWin);
    // No reference supported
    //if (zendeg > 80)
    //    fdiff = 1;

    // the maximum fdiff should be 1
    if (fdiff > 1)
        fdiff = 1;

    PARtop_beam = PARtop * (1 - fdiff);
    PARtop_diff = PARtop - PARtop_beam;
    NIRtop_beam = NIRtop * (1 - fdiff);
    NIRtop_diff = NIRtop - NIRtop_beam;

    // BEAM EXTINCTION COEFF FOR ELLIPSOIDAL LEAF DISTRIBUTION (EQN. 15.4 C&N)
    // or Equation (24) in Drewry et al, 2009 - Part B: Online Supplement.
    Kbm = sqrt( xx*xx + pow(tan(zenrad),2.0) ) / ( xx + 1.774 * pow(xx+1.182, -0.733) );

    PARabs_sun.setZero();
    PARabs_shade.setZero();
    NIRabs_sun.setZero();
    NIRabs_shade.setZero();
    fsun.setZero();
    fshade.setZero();
    LAIsun.setZero();
    LAIshade.setZero();
    diffdn.setZero();
    diffup.setZero();

    PARabs_soil = 0.0;
    NIRabs_soil = 0.0;
    radlost     = 0.0;
    PARout      = 0.0;
    NIRout      = 0.0;

    if (PARtop < 1) 
    {
        LAIshade = LAIz;
        fshade = fshade + VectorXd::Ones(nl_can);
    } else {
        // Iterate to Solve PAR Absorption Profile
        count = 0;
        percdiff = 1;
        radin = PARtop;
        beam_top = PARtop_beam;
        diff_top = PARtop_diff;

        while (percdiff > 0.01)
        {
            SW_Attenuation(PARabs_sun, PARabs_shade, fsun, fshade, diffdn, diffup, LAIz,
                    &PARabs_soil, &radabs_tot, &radlost, &radremain, beam_top, diff_top,
                    trans_PAR, refl_PAR, refl_soil, clump, Kbm, Kdf, count, nl_can);

            beam_top = 0.0;
            diff_top = 0.0;
            radtot   = radabs_tot + radlost;
            percdiff = (radin - radtot) / radin;

            count += 1;
            if (count > 5){
                printf("Count > 5 in PAR loop!!!\n");
                break;
            }
        }

        LAIsun = fsun.cwiseProduct(LAIz);
        LAIshade = fshade.cwiseProduct(LAIz);
        PARout = radlost;


        // Iterate to Solve NIR Absorption Profile
        diffdn.setZero();
        diffup.setZero();
        radlost = 0;
        count = 0;
        percdiff = 1;
        radin = NIRtop;
        beam_top = NIRtop_beam;
        diff_top = NIRtop_diff;

        while (percdiff > 0.01)
        {
            SW_Attenuation(NIRabs_sun, NIRabs_shade, fsun, fshade, diffdn, diffup, LAIz,
                    &NIRabs_soil, &radabs_tot, &radlost, &radremain, beam_top, diff_top,
                    trans_NIR, refl_NIR, refl_soil, clump, Kbm, Kdf, count, nl_can);

            beam_top = 0.0;
            diff_top = 0.0;
            radtot   = radabs_tot + radlost;
            percdiff = ( radin - radtot ) / radin;

            count += 1;
            if (count > 5)
            {
                printf("Count > 5 in NIR loop!!!\n");
                break;
            }
        }
        NIRout = radlost;
    }

    for (int i=0; i<nl_can; i++)
        taud[i] = exp ( -Kdf * clump *LAIz[i] );

    // Total SW outgoing [W/m^2]
    *SWout       = PARout * umoltoWm2 + NIRout;

    // Shortwave Absorption Profiles
    SWabs_sun    = PARabs_sun * umoltoWm2 + NIRabs_sun;
    SWabs_shade  = PARabs_shade*umoltoWm2 + NIRabs_shade;

    // Shortwave Absorbed by Soil
    *SWabs_soil  = NIRabs_soil + PARabs_soil * umoltoWm2;  // [W/m^2]

    PARabs_sun   = PARabs_sun.cwiseQuotient(LAIsun);
    PARabs_shade = PARabs_shade.cwiseQuotient(LAIshade);

    // To prevent numerical error.
    for (int i = 0; i<nl_can; i++) {
        if (fsun[i] < 0.0001)
            fsun[i] = 0;
      
        if (fshade[i] < 0.0001)
            fshade[i] = 0;
      
        if (LAIsun[i] <= 0.0001) {
            LAIsun[i] = 0;
            PARabs_sun[i] = 0;
        }

        if (LAIshade[i] < 0.0001) {
            LAIshade[i] = 0;
            PARabs_shade[i] = 0;
        }

        if (SWabs_sun[i] < 0.0001)
            SWabs_sun[i] = 0;

        if (SWabs_shade[i] < 0.0001)
            SWabs_shade[i] = 0;

        if (PARabs_sun[i] < 0.0001)
            PARabs_sun[i] = 0;

        if (PARabs_shade[i] < 0.0001)
            PARabs_shade[i] = 0;

        if (NIRabs_sun[i] < 0.0001)
            NIRabs_sun[i] = 0;
      
        if (NIRabs_shade[i] < 0.0001)
            NIRabs_shade[i] = 0;
    }
}


void LW_Attenuation( Ref<VectorXd> LWabs_can, Ref<VectorXd> diffdn,
    Ref<VectorXd> diffup, Ref<VectorXd> LWcan_out, Ref<VectorXd> LWsun_out,
    Ref<VectorXd> LWshade_out, Ref<VectorXd> Tl_sun, Ref<VectorXd> Tl_shade,
    Ref<VectorXd> LAI_sun, Ref<VectorXd> LAI_shade, Ref<VectorXd> LAIz,
    Ref<VectorXd> fsun, Ref<VectorXd> fshade, double *radlost_ptr,
    double *LWabs_soil_ptr, double *LW_soil_out_ptr, double LW_sky,
    double Tatop, double Tsoil, double epsv, double epss, double clump,
    double Kdf, int count, int nl_can)
{
    double taud, radlost;
    double soil_inc      = 0.;
    double boltz         = 5.6697 * 1e-8;
    double LW_int        = 0.;
    double LW_flux       = 0.;
    double LW_flux_sun   = 0.;
    double LW_flux_shade = 0.;
    double LW_soil_out   = 0.;
    double LWabs_soil    = 0.;
    

    LWabs_soil = *LWabs_soil_ptr;
    radlost = *radlost_ptr;

    int tind = nl_can - 1;
    int bind = 0;

    
    ///////////////////////////////////
    // DOWNWARD rADIATION            //
    ///////////////////////////////////
    diffdn[tind] += LW_sky;
    for (int i=tind; i>=0; i--)
    {
        // Absorbed Downward
        taud = exp (-Kdf * clump * LAIz[i]);
        LW_int = (1 - taud) * diffdn[i];
        LWabs_can[i] += epsv * LW_int;

        // Downward Flux
        if (i > bind) 
        {
            diffdn[i-1] += taud * diffdn[i] + (1-epsv)*LW_int;
        } else {
            soil_inc += taud * diffdn[i] + (1-epsv) * LW_int;
        }

        ////////////////////////////////////////////////
        // THERMAL CONTRIBUTION FROM FOLIAGE AND SOIL //
        ////////////////////////////////////////////////
        if (count ==0)
        {
            LW_flux_sun = fsun[i]*(1-taud)*epsv*boltz*pow(Tl_sun[i]+273.15, 4);
            LW_flux_shade = fshade[i]*(1-taud)*epsv*boltz*pow(Tl_shade[i]+273.15, 4);
            LW_flux = LW_flux_sun + LW_flux_shade;

            LWcan_out[i] = 2.0 * LW_flux;
            LWsun_out[i] = 2.0 * LW_flux_sun;
            LWshade_out[i] = 2.0 * LW_flux_shade;

            // Upward Flux
            if (i < tind)
            {
                diffup[i+1] += LW_flux;
            } else {
                radlost += LW_flux;
            }

            // Downward Flux
            if (i > bind)
            {
                diffdn[i-1] += LW_flux;
            } else {
                soil_inc += LW_flux;
            }

            // Soil Emission
            diffup[0]   = epss * boltz * pow(Tsoil+273.15, 4);
            LW_soil_out = epss * boltz * pow(Tsoil+273.15, 4);

        } else {
            LWcan_out.setZero();
            LW_soil_out = 0;
        }
        
        diffdn[i] = 0;    // Downward flux has been absorbed or tranmitted
    }

    LWabs_soil += epss * soil_inc;
    diffup[0] += (1-epss) * soil_inc;

    
    ///////////////////////////////////
    // UPWARD RADIATION              //
    ///////////////////////////////////
    for (int i=bind; i<=tind; i++)
    {
        // Absorbed upward
        taud = exp( -Kdf * clump * LAIz[i]);
        LW_int = diffup[i] - taud * diffup[i];
        LWabs_can[i] += epsv * LW_int;

        // Upward Flux
        if (i < tind){
            diffup[i+1] += taud * diffup[i] + (1-epsv) * LW_int;
        } else {
            radlost += taud * diffup[i] + (1-epsv) * LW_int;
        }
        diffup[i] = 0;
    }

    // Return via pointers
    *LW_soil_out_ptr = LW_soil_out;
    *LWabs_soil_ptr  = LWabs_soil;
    *radlost_ptr     = radlost;
}


void LongWaveradiation(ForcingClass *forcings, CanopyClass *canopies,
    VerticalCanopyClass *vertcanopies, VerticalSoilClass *vertsoils,
    ConstantClass *constants, RadiationClass *radiation, Ref<VectorXd> LWabs_can,
    Ref<VectorXd> LWabs_sun, Ref<VectorXd> LWabs_shade, double *LWabs_soil,
    double *LWup, Ref<VectorXd> LWcan_emit, Ref<VectorXd> LWsun_emit,
    Ref<VectorXd> LWshade_emit, double *LWemit_soil, int LW_eq)
{
    double LWin, LW_sky, LW_top, LWtot = 0.;
    double zendeg, zenrad, Tatop, eatop;
    double percdiff, Kdf, clump, epsv, epss, epsa, boltz, Ts;
    double radlost, radremain, LW_soil_out;
    int count, maxiters;
    int nl_can = canopies->nl_can;

    VectorXd fsun(nl_can), fshade(nl_can), diffdn(nl_can), diffup(nl_can);
    VectorXd LAIsun(nl_can), LAIshade(nl_can), LAIz(nl_can);
    VectorXd Tl_sun(nl_can), Tl_shade(nl_can);
    VectorXd LWcan_out(nl_can), LWsun_out(nl_can), LWshade_out(nl_can);

    zendeg = forcings->zen;
    LWin   = forcings->lwdn;
    Tatop  = forcings->ta;
    eatop  = forcings->ea;

    Kdf    = radiation->Kdf;
    clump  = radiation->clump;
    epsv   = radiation->epsv;
    epss   = radiation->epss;
    boltz  = constants->boltz;
    Ts     = vertsoils->Ts[0];

    // Copy array structs to eigen vectors
    fsun     = Map<VectorXd>(vertcanopies->fsun, nl_can);
    fshade   = Map<VectorXd>(vertcanopies->fshade, nl_can);
    LAIz     = Map<VectorXd>(vertcanopies->LAIz, nl_can);
    LAIsun   = Map<VectorXd>(vertcanopies->LAIsun, nl_can);
    LAIshade = Map<VectorXd>(vertcanopies->LAIshade, nl_can);
    Tl_sun   = Map<VectorXd>(vertcanopies->Tl_sun, nl_can);
    Tl_shade = Map<VectorXd>(vertcanopies->Tl_shade, nl_can);
    diffdn   = Map<VectorXd>(vertcanopies->diffdn, nl_can);
    diffup   = Map<VectorXd>(vertcanopies->diffup, nl_can);

    // Convert Zenith in degree to radian
    zenrad = zendeg * pinum/180.0;

    // Downward Longwave from sky
    if (isnan(LWin))
    {
        if (LW_eq == 0) {
            epsa = 1.72 * pow(eatop / (Tatop + 273.15), 1.0/7.0);
            LW_sky = epsa * boltz * pow(Tatop + 273.15, 4);
        }
        else {
            // See http://www.met.wau.nl/Courses/Micrometcourse/Modules/Longwave/ModulePage3.html
            double c1 = 0.53;
            double c2 = 0.067;
            double Lc = 60;
            double e_atm;
            e_atm = c1 + c2 * pow(eatop, 1.0/2.0);
            LW_sky = e_atm * boltz * pow(Tatop + 273.15, 4) + Lc;
        }
    } else {
        LW_sky = LWin;
    }

    LWabs_can.setZero();
    diffdn.setZero();
    diffup.setZero();
    *LWabs_soil = 0.0;
    radlost = 0.0;

    // Iterate to Solve LW Absorption Profile
    count    = 0;
    percdiff = 1;
    maxiters = 20;
    LW_top   = LW_sky;

    while (percdiff > 0.01) 
    {
        LW_Attenuation(LWabs_can, diffdn, diffup, LWcan_out, LWsun_out, LWshade_out, Tl_sun,
                Tl_shade, LAIsun, LAIshade, LAIz, fsun, fshade, &radlost, LWabs_soil,
                &LW_soil_out, LW_top, Tatop, Ts, epsv, epss, clump, Kdf, count, nl_can);

        if (count == 0) 
        {
            LWtot        = LW_sky + LWcan_out.sum() + LW_soil_out;
            LWcan_emit   = LWcan_out;
            LWsun_emit   = LWsun_out;
            LWshade_emit = LWshade_out;
            *LWemit_soil = LW_soil_out;
        }

        LW_top    = 0;
        radremain = diffdn.sum() + diffup.sum();
        percdiff  = radremain / LWtot;

        count += 1;
        if (count > maxiters) 
        {
              printf(" TOO MANY ITERATIONS IN LW LOOP!!!\n");
              break;
        }
    }

    *LWup       = radlost;
    LWabs_sun   = LWabs_can.cwiseProduct(fsun);
    LWabs_shade = LWabs_can.cwiseProduct(fshade);

    // Return single variable via pointer
    //*vertsoils->LWabs_soil = LWabs_soil;
    //*vertsoils->LWemit_soil = LWemit_soil;
    //*vertcanopies->LWout = LWup;
}


double Polynomial_Root(double a, double b, double c)
{
    double x1 = 0.0;
    double x2 = 0.0;
    double delta = b*b - 4*a*c;

    if (delta < 0) 
    {
        cout << "No real root found!!!" << endl;
    }
    else
    {
        x1 = (-b + sqrt(delta))/(2*a);
        x2 = (-b - sqrt(delta))/(2*a);
    }

    if (x1 > x2)
        x1 = x2;
    
    return x1;
}

void photosynthesis_C3(PhotosynthesisClass *photosynthesis, VerticalCanopyClass *vertcanopies, 
                       ConstantClass *constants, Ref<VectorXd> Ph, Ref<VectorXd> An, int nl_can,
                       int sunlit)
{
    double Vcmax25, Jmax25, Rd25, beta, O, R_J, R;
    double Ko, Kc, phiPSIImax, thetaPSII, Q2, J, Jc, Jj, Js, Jp, tt, bb;

    // Mapping array to eigen vector
    Map<VectorXd> Vz(vertcanopies->Vz, nl_can);

    // Initializing eigen vector;
    VectorXd Qabs(nl_can), Tl(nl_can), Ci(nl_can);
    VectorXd Vcmax25_vec(nl_can), Jmax25_vec(nl_can), TlK(nl_can), Rd(nl_can);
    VectorXd Phtype(nl_can), gamstar(nl_can), Wc(nl_can), Wj(nl_can), Vcmax(nl_can), Jmax(nl_can);

    // Copy array structs to eigen vectors
    if (sunlit == 1)
    {      
        Qabs = Map<VectorXd>(vertcanopies->PARabs_sun, nl_can);
        Tl = Map<VectorXd>(vertcanopies->Tl_sun, nl_can);
        Ci = Map<VectorXd>(vertcanopies->Ci_sun, nl_can);
    } else {      
        Qabs = Map<VectorXd>(vertcanopies->PARabs_shade, nl_can);
        Tl = Map<VectorXd>(vertcanopies->Tl_shade, nl_can);
        Ci = Map<VectorXd>(vertcanopies->Ci_shade, nl_can);
    }

    R_J = constants->R;
    Vcmax25 = photosynthesis->Vcmax25_C3 * photosynthesis->Vcmax25_fact;
    Jmax25 = photosynthesis->Jmax25_C3 * photosynthesis->Vcmax25_fact;
    Rd25 =photosynthesis->Rd25 * photosynthesis->Vcmax25_fact;
    beta = photosynthesis->beta_ph_C3;  
    O = photosynthesis->Oi;
    R = R_J / 1000.;

    for (int i=0; i<nl_can; i++)
    {
        TlK[i] = Tl[i] + 273.15;
        Vcmax25_vec[i] = Vz[i] * Vcmax25;
        Jmax25_vec[i] = Vz[i] * Jmax25;
        gamstar[i] = exp(19.02 - 37.83/(R*TlK[i]));
        Ko = exp(20.30 - 36.38/(R*TlK[i]));
        Kc = exp(38.05 - 79.43/(R*TlK[i]));    
        Rd[i] = Rd25 * exp(18.72 - 46.39/(R*TlK[i]));
        Vcmax[i] = Vcmax25_vec[i] * exp(26.35 - 65.33/(R*TlK[i]));

        phiPSIImax = 0.352 + 0.022*Tl[i] - 0.00034*Tl[i]*Tl[i];
        Q2 = Qabs[i] * phiPSIImax * beta;
        thetaPSII = 0.76 + 0.018*Tl[i] - 0.00037*Tl[i]*Tl[i];
        Jmax[i] = Jmax25_vec[i] * exp(17.57 - 43.54/(R*TlK[i]));
        J = (Q2 + Jmax[i] - sqrt((Q2+Jmax[i])*(Q2+Jmax[i]) - 4*thetaPSII*Q2*Jmax[i])) / (2*thetaPSII);

        // Constraint for the negative phiPSIImax
        if (phiPSIImax < 0){
            phiPSIImax = 0;
            J = 0;
        }
        if (thetaPSII < 0){
            thetaPSII = 0;
            J = 0;
        }

        Wc[i] = (Vcmax[i] * Ci[i]) / (Ci[i] + Kc*(1+O/Ko));
        Wj[i] = (J * Ci[i]) / (4.5*Ci[i] + 10.5*gamstar[i]);

        // Limiting Rates
        Jc = (1. - gamstar[i] / Ci[i]) * Wc[i]; // Rubisco-Limited Rate [umol/m^2 leaf area/s] 
        Jj = (1. - gamstar[i] / Ci[i]) * Wj[i]; // Light-Limited Rate [umol/m^2 leaf area/s] 
        Js = Vcmax[i] / 2;                      // Sucrose-Limited Rate [umol/m^2 leaf area/s]

        // Solve quadratics from [Collatz et al, ] to account for co-limitation between rates
        tt = 0.98;
        bb = 0.96;
        Jp = ( (Jc+Jj) - sqrt( (-(Jc+Jj))*(-(Jc+Jj)) - 4*tt*Jc*Jj ) ) / (2*tt);
        Ph[i] = ( (Jp+Js) - sqrt( (-(Jp+Js))*(-(Jp+Js)) - 4*bb*Jp*Js ) ) / (2*bb);
        
        // Constrain for the interal CO2
        if ((Ci[i] <= gamstar[i]) & (Qabs[i] >0)){
            Ph[i] = Rd[i];
        }

        An[i] =  Ph[i] - Rd[i];                 // Photosynthetic minus Leaf Respiration flux from ecosystem [umol CO2/ m^2 leaf / s] 

        if (Jc < Jj && Jc < Js) 
            Phtype[i] = 1;
        else if (Jj <= Jc && Jj <= Js)
              Phtype[i] = 2;
        else
            Phtype[i] = 3;
    }
}

void photosynthesis_C4(PhotosynthesisClass *photosynthesis, VerticalCanopyClass *vertcanopies,
                       Ref<VectorXd> Ph, Ref<VectorXd> An, int nl_can, int sunlit)
{
    double Vmaxs, Q10s, VT, RT, KT, aa, bb, cc, M;
    double Vmax, kk, Q10, theta, beta, Rd, al;

    // Mapping array to eigen vector
    Map<VectorXd> Vz(vertcanopies->Vz, nl_can);

    // Initializing eigen vector;
    VectorXd Qabs(nl_can), Tl(nl_can), Ci(nl_can);

    // Copy array structs to eigen vectors
    if (sunlit == 1)
    {
        Qabs = Map<VectorXd>(vertcanopies->PARabs_sun, nl_can);
        Tl   = Map<VectorXd>(vertcanopies->Tl_sun, nl_can);
        Ci   = Map<VectorXd>(vertcanopies->Ci_sun, nl_can);
    } else {
        Qabs = Map<VectorXd>(vertcanopies->PARabs_shade, nl_can);
        Tl   = Map<VectorXd>(vertcanopies->Tl_shade, nl_can);
        Ci   = Map<VectorXd>(vertcanopies->Ci_shade, nl_can);
    }

    Vmax  = photosynthesis->Vmax_C4;
    kk    = photosynthesis->kk_C4;
    Q10   = photosynthesis->Q10_C4;
    theta = photosynthesis->theta_C4;
    beta  = photosynthesis->beta_C4;
    Rd    = photosynthesis->Rd_C4;
    al    = photosynthesis->al_C4;

    for (int i=0; i<nl_can; i++)
    {
        // Equation 5B in Appendix B, pp.537
        Vmaxs = Vz[i] * Vmax;
        Q10s  = pow( Q10, (Tl[i]-25) * 0.1 );
        VT    = (Vmaxs * Q10s) / ( (1 + exp(0.3*(13-Tl[i]))) * (1+exp(0.3*(Tl[i]-36))) );
        RT    = Rd * Q10s / ( 1 + exp( 1.3*(Tl[i]-55) ) );
        KT    = kk * Q10s;

        // Equation 2B in Appendix B, pp.537
        aa    = theta;
        bb    = -( VT + al * Qabs[i] );
        cc    = VT * al * Qabs[i];

        M     = Polynomial_Root(aa, bb, cc);

        // Equation 3B in Appendix B, pp.537
        aa    = beta;
        bb    = -(M + KT * Ci[i]);
        cc    = M * KT * Ci[i];
        Ph[i] = Polynomial_Root(aa, bb, cc);    // leaf photosynthesis [umol/m^2 leaf area/s]
        An[i] = Ph[i] - RT;                     // net leaf CO2 uptake rate [umol/m^2 leaf area/s]
    }
}

void BLC_Nikolov(PhotosynthesisClass *photosynthesis, CanopyClass *canopies, 
                 VerticalCanopyClass *vertcanopies, Ref<VectorXd> gbv, Ref<VectorXd> gbh,
                 int nl_can, int sunlit)
{
    double gbv_forced, gbv_free, ce, cf, Tv_diff;
    double ld, lw, gsv, Tak, Tlk, Pa, ea, esTl, eb;
    int leaftype;

    // Initialize eigen vectors;
    VectorXd Tl_in(nl_can), gsv_in(nl_can);

    // Mapping array to eigen vectors
    Map<VectorXd> Ta_in(vertcanopies->TAz, nl_can);
    Map<VectorXd> Pa_in(vertcanopies->PAz, nl_can);
    Map<VectorXd> ea_in(vertcanopies->EAz, nl_can);
    Map<VectorXd> U(vertcanopies->Uz, nl_can);

    ld       = canopies->ld;
    lw       = canopies->lw;
    leaftype = canopies->leaftype;

    // Copy array structs to eigen vectors
    if (sunlit == 1)
    {
        Tl_in = Map<VectorXd>(vertcanopies->Tl_sun, nl_can);
        gsv_in = Map<VectorXd>(vertcanopies->gsv_sun, nl_can);
    } else {
        Tl_in = Map<VectorXd>(vertcanopies->Tl_shade, nl_can);
        gsv_in = Map<VectorXd>(vertcanopies->gsv_shade, nl_can);
    }

    for (int i=0; i<nl_can; i++)
    {
        // Unit conversion
        gsv = gsv_in[i] / 41.4;            // [m/s]
        Tak = Ta_in[i] + 273.15;           // [K]
        Tlk = Tl_in[i] + 273.15;           // [K]
        Pa  = Pa_in[i] * 1000;             // [Pa]
        ea  = ea_in[i] * 1000;             // [Pa]

        esTl = 1000 * 0.611 * exp( (17.502*Tl_in[i]) / (Tl_in[i] + 240.97) );

        if (leaftype == 1)
        {
            cf = 1.6361 * 1e-3;
        } else {
            cf = 0.8669 * 1e-3;
        }
        gbv_forced = cf * pow(Tak, 0.56) * pow((Tak+120) * (U[i]/ld/Pa), 0.5);


        if (leaftype == 1)
        {
            ce = 1.6361 * 1e-3;
        } else {
            ce = 0.8669 * 1e-3;
        }
        gbv_free = gbv_forced;
        eb = (gsv * esTl + gbv_free * ea) / (gsv + gbv_free);


        Tv_diff = (Tlk / (1 - 0.378*eb/Pa)) - (Tak / (1 - 0.378*ea/Pa));
        gbv_free = ce * pow(Tlk, 0.56) * sqrt( (Tlk+120)/Pa ) * pow(abs(Tv_diff)/lw, 0.25);

        gbv_forced *= 41.4;
        gbv_free *= 41.4;

        if (gbv_forced > gbv_free){
            gbv[i] = gbv_forced;
        } else {
            gbv[i] = gbv_free;
        }

        gbh[i] = 0.924 * gbv[i];        
    }
}


void Leaf_Water_Potential(VerticalSoilClass *vertsoils, VerticalCanopyClass *vertcanopies,
                          ConstantClass *constants, StomaConductClass *stomaconduct,
                          Ref<VectorXd> psil_MPA, int nl_can)
{
    double rpp_wgt, rpp_wgt_MPa, Rp, grav, dtime, mmH2OtoMPa, rho_kg;

    // Initialize eigen vectors;
    VectorXd TR(nl_can), TR_m(nl_can), TR_sun(nl_can), TR_shade(nl_can);
    VectorXd znc(nl_can);

    // Copy array structs to eigen vectors
    TR_sun = Map<VectorXd>(vertcanopies->TR_sun, nl_can);                             // transpiration PER UNIT LEAF AREA sunlit [mm/s/unit LAI] 
    TR_shade = Map<VectorXd>(vertcanopies->TR_shade, nl_can);                         // transpiration PER UNIT LEAF AREA shade [mm/s/unit LAI] 
    znc = Map<VectorXd>(vertcanopies->znc, nl_can);                                   // Height of canopy levels [m]

    rpp_wgt = vertsoils->rpp_wgt[0];                                                  // root pressure potential weighted by root distribution [mm]
    Rp = stomaconduct->Rp;                                                            // plant resistance to water flow [MPa / m / s]
    grav = constants->grav;                                                           // Gravity Acceleration [m / s^2]
    dtime = constants->dtime;                                                         // Time Step [s]
    mmH2OtoMPa = constants->mmH2OtoMPa;                                               // Conversion Factor from mmH2O to MPa

    rho_kg = 1;                                                                       // [kg / m^3]
    TR = TR_sun + TR_shade;                                                           // [W/LAI/s]
    TR_m = TR / 1000.0;                                                               // [m/s/unit LAI]
    rpp_wgt_MPa = rpp_wgt * mmH2OtoMPa;                                               // [MPa] 
    psil_MPA = rpp_wgt_MPa*VectorXd::Ones(nl_can) - TR*Rp - (rho_kg*grav*znc)/1e6;    // Leaf Water Potential [MPa] 
}


void Tuzet_Function(
    StomaConductClass *stomaconduct, Ref<VectorXd> fsv, Ref<VectorXd> psil_MPa,
    int nl_can)
{
  double sf, psif;
  sf = stomaconduct->sf;
  psif = stomaconduct->psif;

  for (int i=0; i<nl_can; i++){
    fsv[i] = (1 + exp(sf * psif)) / (1 + exp(sf * (psif - psil_MPa[i])));
  }
}


void Ball_Berry(
    StomaConductClass *stomaconduct, VerticalCanopyClass *vertcanopies,
    Ref<VectorXd> gsv, Ref<VectorXd> Ci, Ref<VectorXd> An, Ref<VectorXd> fsv,
    Ref<VectorXd> gbv, int sunlit, int nl_can)
{
  double mslope, bint, Cs, estarTl, ei, es;
  Map<VectorXd> Ca(vertcanopies->CAz, nl_can);
  Map<VectorXd> ea(vertcanopies->EAz, nl_can);
  VectorXd Tl_C(nl_can), Hs(nl_can);

  // Copy array structs to eigen vectors
  if (sunlit == 1){
    Tl_C = Map<VectorXd>(vertcanopies->Tl_sun, nl_can);
    gsv = Map<VectorXd>(vertcanopies->gsv_sun, nl_can);
    Ci = Map<VectorXd>(vertcanopies->Ci_sun, nl_can);
  } else {
    Tl_C = Map<VectorXd>(vertcanopies->Tl_shade, nl_can);
    gsv = Map<VectorXd>(vertcanopies->gsv_shade, nl_can);
    Ci = Map<VectorXd>(vertcanopies->Ci_shade, nl_can);
  }

  mslope = stomaconduct->mslope;
  bint = stomaconduct->bint;

  for (int i=0; i<nl_can; i++){
    Cs = Ca[i] - (1.37*An[i])/gbv[i];
    estarTl = 0.611 * exp( 17.502 * Tl_C[i]/(Tl_C[i]+240.97) );
    ei = estarTl;
    es = (gsv[i]*ei + gbv[i]*ea[i]) / (gsv[i] + gbv[i]);
    Hs[i] = es / estarTl;
    gsv[i] = fsv[i] * mslope * An[i] * Hs[i] / Cs + bint;
    Ci[i] = Ca[i] - (1.37 * An[i]) / gbv[i] - (1.6 * An[i]) / gsv[i];
    if (An[i] <= 0){
      Ci[i] = Ca[i];
      gsv[i] = bint;
    }
  }
}


void LEB_Quartic(CanopyClass *canopies, VerticalCanopyClass *vertcanopies,
                 ConstantClass *constants, RadiationClass *radiation, Ref<VectorXd> Tl,
                 Ref<VectorXd> gv, Ref<VectorXd> H, Ref<VectorXd> LE, int sunlit,
                 int nl_can, int dry)
{
    double epsv, LWfact, Hfact, LEfact, Lv, boltz, cp;
    double c1, c2, c3, c4, c5, aa, bb, dd, ee, estarTl, Me, maxval = 0.0;
    double a1, ce[4], roots[4];
    VectorXd taud(nl_can), eaz(nl_can), Taz(nl_can), Paz(nl_can);
    VectorXd An(nl_can), Rabs(nl_can), LAIfrac(nl_can), leaffrac(nl_can);
    VectorXd gbh(nl_can), gsv(nl_can), gbv(nl_can);
    int num_real, ind;
    char lebwet[5] = "WET";
    char lebdry[5] = "DRY";
    char leb[5];

    // Copy array structs to eigen vectors
    if (sunlit == 1)
    {
        Rabs = Map<VectorXd>(vertcanopies->Rabs_sun, nl_can);
        An = Map<VectorXd>(vertcanopies->An_sun, nl_can);
        gsv = Map<VectorXd>(vertcanopies->gsv_sun, nl_can);
        gbv = Map<VectorXd>(vertcanopies->gbv_sun, nl_can);
        gbh = Map<VectorXd>(vertcanopies->gbh_sun, nl_can);
        leaffrac = Map<VectorXd>(vertcanopies->fsun, nl_can);
        LAIfrac = Map<VectorXd>(vertcanopies->LAIsun, nl_can);
    }
    else
    {
        Rabs = Map<VectorXd>(vertcanopies->Rabs_shade, nl_can);
        An = Map<VectorXd>(vertcanopies->An_shade, nl_can);
        gsv = Map<VectorXd>(vertcanopies->gsv_shade, nl_can);
        gbv = Map<VectorXd>(vertcanopies->gbv_shade, nl_can);
        gbh = Map<VectorXd>(vertcanopies->gbh_shade, nl_can);
        leaffrac = Map<VectorXd>(vertcanopies->fshade, nl_can);
        LAIfrac = Map<VectorXd>(vertcanopies->LAIshade, nl_can);
    }


    // Mapping array to eigen vectors
    taud = Map<VectorXd>(vertcanopies->taud, nl_can);
    eaz = Map<VectorXd>(vertcanopies->EAz, nl_can);
    Taz = Map<VectorXd>(vertcanopies->TAz, nl_can);
    Paz = Map<VectorXd>(vertcanopies->PAz, nl_can);

    // Get parameters
    epsv = radiation->epsv;
    LWfact = canopies->LWfact;
    Hfact = canopies->Hfact;
    LEfact = canopies->LEfact;

    Lv = constants->Lv;
    boltz = constants->boltz;
    cp = constants->cp_mol;

    if (dry == 0)
    {
        gbh = gbh * 0.0;
    }


    for (int i=0; i<nl_can; i++)
    {
        // Energy stored in biochemical reactions (Me)
        Me = 0.506 * An[i];

        // Total vapor conductance (gv)
        if (dry == 1)
        {
            gv[i] = (gbv[i] * gsv[i]) / (gbv[i] + gsv[i]);
            strcpy( leb, lebdry);
        }
        else
        {
            gv[i] = gbv[i];
            strcpy( leb, lebwet);
        }

        aa = 273.15;
        bb = LWfact * epsv * boltz * (1-taud[i]) * leaffrac[i] / LAIfrac[i];
        dd = Hfact * cp * gbh[i];
        ee = LEfact * gv[i] * Lv / Paz[i];

        //printf("gbh = %10.8f \n", gbh[i]);

        c1 = (5.82436 * 1e-4) / 1000;
        c2 = (1.5842*1e-2) / 1000;
        c3 = 1.55186 / 1000;
        c4 = 44.513596 / 1000;
        c5 = 607.919 / 1000;

        // Terms of Quartic
        a1 = ee * c1 + bb;                                                            // * Tl^4
        ce[0] = (ee * c2 + 4 * aa * bb) / a1;                                         // * Tl^3
        ce[1] = (ee * c3 + 6 * bb * pow(aa, 2)) / a1;                                 // * Tl^2
        ce[2] = (ee * c4 + dd + 4 * bb * pow(aa,3)) / a1;                             // * Tl
        ce[3] = (ee * (c5-eaz[i])-dd * Taz[i]+bb * pow(aa, 4) - Rabs[i] + Me) / a1;   // constant term

        num_real = quarticSolver(ce, roots);
        if (num_real == 0)
        {   // No real root found
            // printf("WARNING: NO REAL ROOTS FOUND IN LEAF ENERGY BALANCE -- %s \n", leb);
            Tl[i] = Taz[i];
        }
        else
        {
            // Find positive solution first
            ind = 0;
            for (int j=0; j<num_real; j++)
            {
                if (roots[j] > 0)
                {
                    ind += 1;
                    maxval = roots[j];
                }
            }

            if (ind == 0)
            {
                //printf("WARNING: NO REAL POSITIVE ROOTS FOUND IN LEAF ENERGY BALANCE -- %s\n", leb);
                Tl[i] = Taz[i];
            }
            else
            {
                // Find largest positive
                for (int j=0; j<num_real; j++)
                {
                    if (maxval < roots[j] && roots[j] > 0)
                    {
                        maxval = roots[j];
                        ind = j;
                    }
                }
                Tl[i] = maxval;
            }
        }

        estarTl = 0.611 * exp( 17.502 * Tl[i]/(Tl[i]+240.97) );
        LE[i] = LEfact * ( Lv*gv[i]/Paz[i] ) * ( estarTl - eaz[i] );
        H[i] = Hfact * cp * gbh[i] * (Tl[i] - Taz[i]);
    }
}


void LeafSolution(ForcingClass *forcings, CanopyClass *canopies, VerticalCanopyClass *vertcanopies,
                  VerticalSoilClass *vertsoils, PhotosynthesisClass *photosynthesis,
                  ConstantClass *constants, RadiationClass *radiation,
                  StomaConductClass *stomaconduct, Ref<VectorXd> Ph, Ref<VectorXd> An,
                  Ref<VectorXd> Ci, Ref<VectorXd> gsv, Ref<VectorXd> Tl, Ref<VectorXd> LE,
                  Ref<VectorXd> TR, Ref<VectorXd> Evap_mm, Ref<VectorXd> H, Ref<VectorXd> psil,
                  Ref<VectorXd> fsv, Ref<VectorXd> Ch2o_mm, Ref<VectorXd> gbv, Ref<VectorXd> gbh,
                  int sunlit, int t)
{
    double Lv_g, relax, relaxval, maxchange, converged;
    double md1, md2;
    int ph_type, maxiters, cnt, dry;
    int nl_can = canopies->nl_can;

    ph_type = photosynthesis->ph_type;
    Lv_g = constants->Lv_g;

    // Mapping array to eigen vectors
    Map<VectorXd> CAz(vertcanopies->CAz, nl_can);
    Map<VectorXd> TAz(vertcanopies->TAz, nl_can);
    Map<VectorXd> wetfrac(vertcanopies->wetfrac, nl_can);
    Map<VectorXd> dryfrac(vertcanopies->dryfrac, nl_can);

    // Initializing eigen vector for local usages
    VectorXd Ph_prev(nl_can), An_prev(nl_can), Ci_prev(nl_can), Tl_prev(nl_can);
    VectorXd gsv_prev(nl_can), gsvdiff(nl_can), gsvdiffprev(nl_can);
    VectorXd Ci_quot(nl_can), gsv_quot(nl_can), Tl_quot(nl_can);
    VectorXd Tl_dry(nl_can), LE_dry(nl_can), H_dry(nl_can), gv_dry(nl_can);
    VectorXd Tl_wet(nl_can), LE_wet(nl_can), H_wet(nl_can), gv_wet(nl_can);
    VectorXd Ch2o_mm_dry(nl_can), Ch2o_mm_wet(nl_can), Evap_wm2(nl_can);

    relax = 0;
    relaxval = 0.25;
    maxchange = 0.25;
    maxiters = 50;
    converged = 0;
    cnt = 0;

    while (converged == 0)
    {
        // PHOTOSYNTHESIS
        if (ph_type == 1) 
        {
            photosynthesis_C3(photosynthesis, vertcanopies, constants, Ph, An, nl_can, sunlit);
        } else {
            photosynthesis_C4(photosynthesis, vertcanopies, Ph, An, nl_can, sunlit);
        }

        if (cnt > 0)
        {
            Ph = Ph - relax * (Ph-Ph_prev);
            An = An - relax * (An-An_prev);
        }

        // LEAF BOUNDARY LAYER CONDUCTANCES
        BLC_Nikolov(photosynthesis, canopies, vertcanopies, gbv, gbh, nl_can, sunlit);

        Leaf_Water_Potential(vertsoils, vertcanopies, constants, stomaconduct, psil, nl_can);

        Tuzet_Function(stomaconduct, fsv, psil, nl_can);

        Ball_Berry(stomaconduct, vertcanopies, gsv, Ci, An, fsv, gbv, sunlit, nl_can);
        if (cnt > 0)
        {
            gsv = gsv - relax*(gsv-gsv_prev);
            Ci = Ci - relax*(Ci-Ci_prev);
        }

        if (sunlit == 1)
        {
            Map<VectorXd>(vertcanopies->gsv_sun, nl_can) = gsv;            // Stomatal Conductance [mol/m^2/s]
            Map<VectorXd>(vertcanopies->Ci_sun, nl_can) = Ci;
        }
        else
        {
            Map<VectorXd>(vertcanopies->gsv_shade, nl_can) = gsv;          // Stomatal Conductance [mol/m^2/s]
            Map<VectorXd>(vertcanopies->Ci_shade, nl_can) = Ci;
        }

        // LEAF ENERGY BALANCE - DRY LEAF FRACTION
        dry = 1;
        LEB_Quartic(canopies, vertcanopies, constants, radiation, Tl_dry, gv_dry, H_dry, LE_dry,
                sunlit, nl_can, dry);

        // Copy eigen vectors back to array structs
        if (sunlit == 1)
            Map<VectorXd>(vertcanopies->TR_sun, nl_can) = LE_dry/Lv_g;    // [mm/s/LAI]
        else
            Map<VectorXd>(vertcanopies->TR_shade, nl_can) = LE_dry/Lv_g;  // [mm/s/LAI]

        TR = LE_dry / Lv_g;                                               // [mm/s/LAI]

        // LEAF ENERGY BALANCE - WET LEAF FRACTION
        dry = 0;
        LEB_Quartic(canopies, vertcanopies, constants, radiation, Tl_wet, gv_wet, H_wet, LE_wet,
                sunlit, nl_can, dry);

        // Mean Leaf Temperature
        Tl = Tl_dry.cwiseProduct(dryfrac) + Tl_wet.cwiseProduct(wetfrac);

        // Copy eigen vectors back to array structs
        if (sunlit == 1)
            Map<VectorXd>(vertcanopies->Tl_sun, nl_can) = Tl;
        else
            Map<VectorXd>(vertcanopies->Tl_shade, nl_can) = Tl;

        Map<VectorXd>(vertcanopies->Tl, nl_can) = Tl;

        // Test for solution divergence
        if (cnt > 0) {
            gsvdiff = gsv - gsv_prev;

            // Check for solution divergence
            //Ci_quot = (Ci-Ci_prev).cwiseQuotient(Ci_prev);
            //gsv_quot = (gsv-gsv_prev).cwiseQuotient(gsv_prev);
            //Tl_quot = (Tl-Tl).cwiseQuotient(Tl_prev);
            // They should be the absolute value.
            Ci_quot = (Ci - Ci_prev).cwiseQuotient(Ci_prev).cwiseAbs();
            gsv_quot = (gsv - gsv_prev).cwiseQuotient(gsv_prev).cwiseAbs();
            Tl_quot = (Tl - Tl).cwiseQuotient(Tl_prev).cwiseAbs();

            if (Ci_quot.maxCoeff() > maxchange || gsv_quot.maxCoeff() > maxchange ||
                Tl_quot.maxCoeff() > maxchange)
            {
                // Rewind calculation and set relaxation on
                gsv = gsv_prev;
                Ci = Ci_prev;
                Tl = Tl_prev;
                relax = relaxval;
            }
            else if (relax > 0)
            {
                relax = 0;
            }

            if (cnt > 2)
            {
                md1 = (gsvdiffprev.cwiseAbs()).maxCoeff();
                md2 = (gsvdiff.cwiseAbs()).maxCoeff();

                if ( md1 > (md2-0.01*md2) && md1 < (md2+0.01*md2) )
                    relax = relaxval;
            }
            gsvdiffprev = gsvdiff;
        }

        // Test convergence
        if (cnt > 0) {
            //Ci_quot = ((Ci-Ci_prev).cwiseAbs()).cwiseQuotient(Ci_prev);
            //gsv_quot = ((gsv-gsv_prev).cwiseAbs()).cwiseQuotient(gsv_prev);
            //Tl_quot = ((Tl-Tl_prev).cwiseAbs()).cwiseQuotient(Tl_prev);

            // They should be the absolute value.
            Ci_quot = ((Ci - Ci_prev).cwiseQuotient(Ci_prev)).cwiseAbs();
            gsv_quot = ((gsv - gsv_prev).cwiseQuotient(gsv_prev)).cwiseAbs();
            Tl_quot = ((Tl - Tl_prev).cwiseQuotient(Tl_prev)).cwiseAbs();

            if (Ci_quot.maxCoeff() < 0.01 && gsv_quot.maxCoeff() < 0.01 &&
                Tl_quot.maxCoeff() < 0.01)
            {
                converged = 1;
            }
        }

        // Update convergence check variables
        Ph_prev = Ph;
        An_prev = An;
        Ci_prev = Ci;
        gsv_prev = gsv;
        Tl_prev = Tl;

        if ( cnt>maxiters && converged==0 )
        {
            //cout<<"*** TOO MANY INTERATIONS IN LEAF MODEL, "<< t << endl;
            break;
        }
        cnt += 1;
    }

    H = H_dry.cwiseProduct(dryfrac);    // only dry fraction can produce H
    LE = LE_dry.cwiseProduct(dryfrac) + LE_wet.cwiseProduct(wetfrac);

    // Compute Evaporation
    Evap_wm2 = LE_wet.cwiseProduct(wetfrac);
    for (int i=0; i<nl_can; i++)
    {
        if (Evap_wm2[i] < 0)
          Evap_wm2[i] = 0;
    }
    Evap_mm = Evap_wm2 / Lv_g;

    // Compute Condensation
    Ch2o_mm_dry = LE_dry / Lv_g;
    Ch2o_mm_wet = LE_wet / Lv_g;

    for (int i=0; i<nl_can; i++)
    {
        if (Ch2o_mm_dry[i] > 0)
          Ch2o_mm_dry[i] = 0;

        if (Ch2o_mm_wet[i] > 0)
          Ch2o_mm_wet[i] = 0;
    }

    Ch2o_mm = Ch2o_mm_dry + Ch2o_mm_wet;
}


void SEB_Remainder_return(double Ts, double Rabs, double Ta1, double Ts1, double ea1, double pa,
                          double U1, double z1, double dzs1, double psis1_MPa, double vonk,
                          double z0, double TC1, double epss, double *Hs, double *RH, double *LEs,
                          double *Gs, double *LWups, double *remain)
{
    double density_dry_air = 1.2923;
    double rhoa = density_dry_air  * 1000 / 28.97;
    double cp = 29.3;
    double Lv = 44000;
    double boltz = 5.6697 * 1e-8;
    double Vw = 18;
    double R = 8.3143;

    double D = U1 * vonk*vonk / (log(z1/z0) * log(z1/z0));
    *Hs = cp * rhoa * D * (Ts - Ta1);
    double esatTs = 0.611 * exp( 17.502 * Ts / (Ts + 240.97) );
    *RH = exp( psis1_MPa * Vw / R / (Ts + 273.15) );
    // Constraint
    if (*RH > 1){
        *RH = 1;
    }
    if (*RH < 0){
        *RH = 0;
    }
    *LEs = Lv * rhoa * D * (0.622/pa) * (esatTs * (*RH) - ea1);

    // Typo in Drewry et al, 2009
    //*Gs = TC1 * (Ts-Ts1) / dzs1;
    *Gs = TC1 * (Ts - Ts1) / dzs1 / 2;
    *LWups = epss * boltz * pow(Ts + 273.15, 4);
    *remain = Rabs - *Hs - *LEs - *Gs - *LWups;
}


double SEB_Remainder(double Ts, double Rabs, double Ta1, double Ts1, double ea1, double pa,
                     double U1, double z1, double dzs1, double psis1_MPa, double vonk, double z0,
                     double TC1, double epss)
{
    double density_dry_air = 1.2923;
    double rhoa = density_dry_air  * 1000 / 28.97;
    double cp = 29.3;
    double Lv = 44000;
    double boltz = 5.6697 * 1e-8;
    double Vw = 18;
    double R = 8.3143;

    double D = U1 * vonk*vonk / (log(z1/z0) * log(z1/z0));
    double Hs = cp * rhoa * D * (Ts - Ta1);
    double esatTs = 0.611 * exp( 17.502 * Ts / (Ts + 240.97) );
    double RH = exp( psis1_MPa * Vw / R / (Ts + 273.15) );
    // Constraint
    if (RH > 1){
        RH = 1;
    }
    if (RH < 0){
        RH = 0;
    }
    double LEs = Lv * rhoa * D * (0.622/pa) * (esatTs * RH - ea1);

    // Typo in Drewry et al, 2009
    //double Gs = TC1 * (Ts-Ts1) / dzs1;
    double Gs = TC1 * (Ts-Ts1) / dzs1 / 2;
    double LWups = epss * boltz * pow(Ts + 273.15, 4);
    double remain = Rabs - Hs - LEs - Gs - LWups;
    return remain;
}


void Soil_Surface_Fluxes(VerticalCanopyClass *vertcanopies, CanopyClass *canopies,
                         VerticalSoilClass *vertsoils, SoilClass *soils, RadiationClass *radiation,
                         ConstantClass *constants, double *Tsurf, double *Hs, double *LEs,
                         double *Gs, double *RH, double *LWups)
{
    double Ta1, Ts1, ea1, pa1, U1, Rabs, volliq1, psis1, z1, dzs1, TC1, porsl1;
    double epss, z0, vonk, Lv, boltz, cp_mol, R, rho_dry_air, mmH2OtoMPa, psis1_MPa, remain;

    Ta1         = vertcanopies->TAz[0];
    ea1         = vertcanopies->EAz[0];
    pa1         = vertcanopies->PAz[0];
    U1          = vertcanopies->Uz[0];

    Rabs        = *vertsoils->Totabs_soil;
    volliq1     = vertsoils->volliq[0];
    psis1       = vertsoils->smp[0];

    z1          = vertcanopies->znc[0];
    dzs1        = vertsoils->dzs[0];
    TC1         = vertsoils->TK_sol[0];
    porsl1      = vertsoils->porosity[0];
    Ts1         = vertsoils->Ts[0];

    epss        = radiation->epss;
    z0          = soils->z0;
    vonk        = constants->vonk;

    Lv          = constants->Lv;
    boltz       = constants->boltz;
    cp_mol      = constants->cp_mol;
    R           = constants->R;
    rho_dry_air = constants->rho_dry_air;
    mmH2OtoMPa  = constants->mmH2OtoMPa;
    psis1_MPa   = psis1 * mmH2OtoMPa;

    // Setup SEB_Remainder equation
    double (*SEB_Rem) (double, double Rabs, double Ta1, double Ts1, double ea1, double pa1,
            double U1, double z1, double dzs1, double psis1_MPa, double vonk, double z0,
            double TC1, double epss);
    SEB_Rem = &SEB_Remainder;

    // Using Bisection method to find the root of Soil Energy Balance Eqn;
    // This function is equivalent to fzero used in Matlab
    *Tsurf = rtbis(SEB_Rem, Rabs, Ta1, Ts1, ea1, pa1, U1, z1, dzs1, psis1_MPa, vonk, z0, TC1,
            epss, Ta1-50, Ta1+50, 1e-12);

    SEB_Remainder_return(*Tsurf, Rabs, Ta1, Ts1, ea1, pa1, U1, z1, dzs1, psis1_MPa, vonk, z0, TC1,
            epss, Hs, RH, LEs, Gs, LWups, &remain);
}



void Order_1_Closure_All(Ref<VectorXd> Ca_in, Ref<VectorXd> z_in, Ref<VectorXd> K_in,
                         Ref<VectorXd> SS_in, double dz, double soilf, double ctz, int nl_can)
{
    int ct_ind, N, cnt, max_iters;
    double err, epsilon, eps1;

    ct_ind = 0;
    for (int i=0; i<nl_can; i++)
    {
        if (ctz <= z_in[i]){
          ct_ind = i;
          break;
        }
    }

    if (ct_ind == 0)
    ct_ind = nl_can;

    VectorXd Ca(ct_ind), K(ct_ind), SS(ct_ind), z(ct_ind);
    VectorXd Conc(ct_ind), Ca_prev(ct_ind);
    VectorXd a1(ct_ind), a2(ct_ind), a3(ct_ind), a4(ct_ind);
    VectorXd aa(ct_ind), bb(ct_ind), cc(ct_ind), dd(ct_ind);
    VectorXd upd(ct_ind), dia(ct_ind), lod(ct_ind), co(ct_ind), Cn(ct_ind);

    // Prevent numerical error in MLCan
    VectorXd Conc_fixed(nl_can);

    for (int i=0; i<ct_ind; i++){
    Ca[i] = Ca_in[i];
    K[i] = K_in[i];
    SS[i] = SS_in[i];
    z[i] = z_in[i];
    }
    N = ct_ind;

    cnt = 1;
    err = 1e6;
    epsilon = 1e-2;
    max_iters = 20;

    while (err > epsilon && cnt < max_iters)
    {
        if (cnt != 1)
            Ca_prev = Ca;

        Conc = Ca;
        // Set up coefficients for ODE - Eqn (31) Drewry et al., 2009 - part B
        a1 = K;
        for (int i=1; i<N; i++)
        {
            a2[i] = (K[i] - K[i-1]) / dz;
        }
        a2[0] = a2[1];
        a3 = 0*z;
        a4 = -SS;

        // Set the elements of the Tri-diagonal Matrix
        upd = ( a1/(dz*dz) + a2/(2*dz) );
        dia = ( -a1*2/(dz*dz) + a3 );
        lod = ( a1/(dz*dz) - a2/(2*dz) );
        co = a4;

        aa = lod;
        bb = dia;
        cc = upd;
        dd = co;

        aa[0] = 0;
        bb[0] = 1.0;
        cc[0] = -1.0;
        dd[0] = soilf * dz / (K[0] + 0.00001);

        aa[N-1] = 0;
        bb[N-1] = 1.0;
        cc[N-1] = 0;
        dd[N-1] = Ca[N-1];

        // Use the Thomas Algorithm to solve the tridiagonal matrix
        ThomasAlgorithm (aa, bb, cc, dd, Cn, N);
        eps1 = 0.5;
        Conc = (eps1 * Cn + (1-eps1)*Conc);
        Ca = Conc;

        if (cnt != 1)
            err = (Ca - Ca_prev).cwiseAbs().maxCoeff();
        cnt += 1;
    }

    for (int i = 0; i<nl_can; i++){
        if (i < ct_ind){
            Conc_fixed[i] = Conc[i];
        }
        else {
            Conc_fixed[i] = Conc[ct_ind - 1];
        }
    }
    Ca_in = Conc_fixed;

    if (cnt >= max_iters)
        printf("*** Closure Max Iters!!!\n");
}

void Micro_Environment(VerticalCanopyClass *vertcanopies, CanopyClass *canopies,
                       VerticalSoilClass *vertsoils, ConstantClass *constants,
                       Ref<VectorXd> CAz, Ref<VectorXd> EAz, Ref<VectorXd> TAz, int nl_can)
{
    double Fc_soil, LE_soil, H_soil, dzc, hcan, Lv, cp_mol;
    double psy, Sv_soil, Sh_soil;

    VectorXd An_sun(nl_can), An_shade(nl_can), LE_sun(nl_can), LE_shade(nl_can);
    VectorXd H_sun(nl_can), H_shade(nl_can), LAIsun(nl_can), LAIshade(nl_can);
    VectorXd znc(nl_can), Km(nl_can+1), molar_density(nl_can), PAz(nl_can);
    VectorXd q(nl_can), Sc(nl_can), Sv(nl_can), Sh(nl_can);

    dzc      = vertcanopies->dzc;
    hcan     = canopies->hcan;
    Fc_soil  = *vertsoils->Fc_soil;
    LE_soil  = *vertsoils->LE_soil;
    H_soil   = *vertsoils->H_soil;
    Lv       = constants->Lv;
    cp_mol   = constants->cp_mol;

    Km       = Map<VectorXd>(vertcanopies->Km, nl_can+1);
    PAz      = Map<VectorXd>(vertcanopies->PAz, nl_can);
    An_sun   = Map<VectorXd>(vertcanopies->An_sun, nl_can);
    An_shade = Map<VectorXd>(vertcanopies->An_shade, nl_can);
    LE_sun   = Map<VectorXd>(vertcanopies->LE_sun, nl_can);
    LE_shade = Map<VectorXd>(vertcanopies->LE_shade, nl_can);
    H_sun    = Map<VectorXd>(vertcanopies->H_sun, nl_can);
    H_shade  = Map<VectorXd>(vertcanopies->H_shade, nl_can);
    LAIsun   = Map<VectorXd>(vertcanopies->LAIsun, nl_can);
    LAIshade = Map<VectorXd>(vertcanopies->LAIshade, nl_can);
    znc      = Map<VectorXd>(vertcanopies->znc, nl_can);

    psy           = 6.66 * 1e-4;
    molar_density = 44.6 * 273.15 * PAz.cwiseQuotient(101.3 * (TAz+273.15*VectorXd::Ones(nl_can)));
    Sc            = (An_sun.cwiseProduct(LAIsun) + An_shade.cwiseProduct(LAIshade))/dzc;
    Sc            = (-Sc).cwiseQuotient(molar_density);
    Order_1_Closure_All(CAz, znc, Km, Sc, dzc, Fc_soil, hcan, nl_can);

    q       = (EAz.cwiseQuotient(PAz)).cwiseProduct(molar_density);
    Sv      = (LE_sun.cwiseProduct(LAIsun) + LE_shade.cwiseProduct(LAIshade))/dzc;
    Sv      = Sv / Lv;
    Sv_soil = LE_soil / Lv;
    Order_1_Closure_All(q, znc, Km, Sv, dzc, Sv_soil, hcan, nl_can);
    EAz = (q.cwiseQuotient(molar_density)).cwiseProduct(PAz);

    Sh      = (H_sun.cwiseProduct(LAIsun) + H_shade.cwiseProduct(LAIshade))/dzc;
    Sh      = (Sh/cp_mol).cwiseQuotient(molar_density);
    Sh_soil = H_soil / cp_mol / molar_density[0];
    Order_1_Closure_All(TAz, znc, Km, Sh, dzc, Sh_soil, hcan, nl_can);
}

// Fix water mass balance in canopy
void Evap_Condensation_Adjust(VerticalCanopyClass *vertcanopies, CanopyClass *canopies,
                              VerticalSoilClass *vertsoils, Ref<VectorXd> Sh2o_prof,
                              double *Sh2o_can, Ref<VectorXd> wetfrac, Ref<VectorXd> dryfrac, int nl_can)							  
{
	double Ffact, ppt_ground;
    VectorXd znc(nl_can), H2oinc(nl_can), Smaxz(nl_can);
    VectorXd Ch2o_prof(nl_can), Evap_prof(nl_can);
    
    double dripout;
    VectorXd dripv(nl_can);

    znc        = Map<VectorXd>(vertcanopies->znc, nl_can);
    Ch2o_prof  = Map<VectorXd>(vertcanopies->Ch2o_prof, nl_can);
    Evap_prof  = Map<VectorXd>(vertcanopies->Evap_prof, nl_can);
    Smaxz      = Map<VectorXd>(vertcanopies->Smaxz, nl_can);
    Ffact      = canopies->Ffact;
    ppt_ground = *vertsoils->ppt_ground;
    H2oinc     = Ch2o_prof - Evap_prof;
    
    dripout = 0;
    dripv = VectorXd::Zero(nl_can);

    for (int i=nl_can-1; i>=0; i--)
    {
        Sh2o_prof[i] += H2oinc[i] + dripv[i];
        if (Sh2o_prof[i] < 0) 
        {
            Evap_prof[i] += Sh2o_prof[i];
            Sh2o_prof[i] = 0;
        }
        else if (Sh2o_prof[i] >= Smaxz[i])
        {
            if (i == 0)
                dripout = Sh2o_prof[i] - Smaxz[i];
            else
                dripv[i - 1] = Sh2o_prof[i] - Smaxz[i];
            Sh2o_prof[i] = Smaxz[i];
        }
    }

    Map<VectorXd>(vertcanopies->Evap_prof, nl_can) = Evap_prof;
    *Sh2o_can = Sh2o_prof.sum();

    ppt_ground = ppt_ground + dripout;
    *vertsoils->ppt_ground = ppt_ground;                  // Rate of rainfall hit the ground [mm]
    wetfrac = Ffact * (Sh2o_prof.cwiseQuotient(Smaxz));
    dryfrac = VectorXd::Ones(nl_can) - wetfrac;
}


void CanopyModel(ProjectClass *project, SwitchClass *Switches, ConstantClass *constants,
                 CanopyClass *canopies, VerticalCanopyClass *vertcanopies, SoilClass *soils,
                 VerticalSoilClass *vertsoils, TimeForcingClass *timeforcings, 
                 ForcingClass *forcings, RadiationClass *radiation, 
                 PhotosynthesisClass *photosynthesis, StomaConductClass *stomaconduct,
                 RespirationClass *respiration, MicroEnvironmentClass *microenviron, 
                 OutputClass *outmlcan, int t, int rank, int size)
{
    int converged_LW, cnt_LW, maxiters;
    int turb_on = Switches->Turbulence;
    int LW_eq = Switches->LWequation;
    int nl_can  = canopies->nl_can;
    double ppt_ground, percdiff, Fc_soil, H_soil, LE_soil, G, RH_soil, T_surf;
    double LWout, LWabs_soil, LWemit_soil, Totabs_soil, Rnrad_soil, LWups;
    double SWout, SWabs_soil, Rnrad_eco, Evap_can, Ch2o_can, Sh2o_can;
    double An_can, LE_can, H_can, TR_can, Rnrad_can;

    double Tsoil   = vertsoils->Ts[0];
    double Ro      = respiration->Ro;
    double Q10     = respiration->Q10;
    double dtime   = constants->dtime;

    forcings->doy  = timeforcings->doy[t];
    forcings->rg   = timeforcings->rg[t];
    forcings->pa   = timeforcings->pa[t];
    forcings->lwdn = timeforcings->lwdn[t];
    forcings->zen  = timeforcings->zen[t];
    forcings->u    = timeforcings->u[t];
    forcings->ppt  = timeforcings->ppt[t];
    forcings->ta   = timeforcings->ta[t];
    forcings->ea   = timeforcings->ea[t];
    forcings->ca   = project->co2concentration;

    for (int i=0; i<nl_can; i++)
    {
        vertcanopies->LAIz[i]     = timeforcings->lai[t] * vertcanopies->LADnorm[i];
        vertcanopies->LADz[i]     = vertcanopies->LAIz[i] / vertcanopies->dzc;
        vertcanopies->TAz[i]      = forcings->ta;
        vertcanopies->CAz[i]      = forcings->ca;
        vertcanopies->EAz[i]      = forcings->ea;
        vertcanopies->PAz[i]      = forcings->pa;
        vertcanopies->Uz[i]       = forcings->u;
        vertcanopies->TR_sun[i]   = 0.0;
        vertcanopies->TR_shade[i] = 0.0;
    }
  
    // Create local Eigen Vectors for operations
    VectorXd Totabs_sun(nl_can), Totabs_shade(nl_can);
    VectorXd diffprof(nl_can), percdiffprof(nl_can), LWabs_prev(nl_can);
    VectorXd Ci_sun(nl_can), gsv_sun(nl_can), Tl_sun(nl_can);
    VectorXd Ci_shade(nl_can), gsv_shade(nl_can), Tl_shade(nl_can);

    // Copy array structs to eigen vectors
    gsv_sun   = Map<VectorXd>(vertcanopies->gsv_sun, nl_can);
    Ci_sun    = Map<VectorXd>(vertcanopies->Ci_sun, nl_can);
    Tl_sun    = Map<VectorXd>(vertcanopies->Tl_sun, nl_can);
    gsv_shade = Map<VectorXd>(vertcanopies->gsv_shade, nl_can);
    Ci_shade  = Map<VectorXd>(vertcanopies->Ci_shade, nl_can);
    Tl_shade  = Map<VectorXd>(vertcanopies->Tl_shade, nl_can);

    
    // Mapping Eigen vectors to Structs
    Map<VectorXd> SWabs_sun(vertcanopies->SWabs_sun, nl_can);
    Map<VectorXd> SWabs_shade(vertcanopies->SWabs_shade, nl_can);
    Map<VectorXd> Rabs_sun_lai(vertcanopies->Rabs_sun, nl_can);
    Map<VectorXd> Rabs_shade_lai(vertcanopies->Rabs_shade, nl_can);
    Map<VectorXd> fsun(vertcanopies->fsun, nl_can);
    Map<VectorXd> fshade(vertcanopies->fshade, nl_can);
    Map<VectorXd> LAIsun(vertcanopies->LAIsun, nl_can);
    Map<VectorXd> LAIshade(vertcanopies->LAIshade, nl_can);
    Map<VectorXd> PARabs_sun(vertcanopies->PARabs_sun, nl_can);
    Map<VectorXd> PARabs_shade(vertcanopies->PARabs_shade, nl_can);
    Map<VectorXd> NIRabs_sun(vertcanopies->NIRabs_sun, nl_can);
    Map<VectorXd> NIRabs_shade(vertcanopies->NIRabs_shade, nl_can);
    Map<VectorXd> taud(vertcanopies->taud, nl_can);

    Map<VectorXd> LWabs_can(vertcanopies->LWabs_can, nl_can);
    Map<VectorXd> LWabs_sun(vertcanopies->LWabs_sun, nl_can);
    Map<VectorXd> LWabs_shade(vertcanopies->LWabs_shade, nl_can);
    Map<VectorXd> LWcan_out(vertcanopies->LWcan_out, nl_can);
    Map<VectorXd> LWshade_out(vertcanopies->LWshade_out, nl_can);
    Map<VectorXd> LWsun_out(vertcanopies->LWsun_out, nl_can);
    Map<VectorXd> LWemit_can(vertcanopies->LWemit_can, nl_can);
    Map<VectorXd> LWemit_sun(vertcanopies->LWemit_sun, nl_can);
    Map<VectorXd> LWemit_shade(vertcanopies->LWemit_shade, nl_can);

    Map<VectorXd> Ph_sun(vertcanopies->Ph_sun, nl_can);
    Map<VectorXd> Ph_shade(vertcanopies->Ph_shade, nl_can);
    Map<VectorXd> Vz(vertcanopies->Vz, nl_can);
    Map<VectorXd> Uz(vertcanopies->Uz, nl_can);
    Map<VectorXd> Km(vertcanopies->Km, nl_can+1);

    Map<VectorXd> CAz(vertcanopies->CAz, nl_can);
    Map<VectorXd> EAz(vertcanopies->EAz, nl_can);
    Map<VectorXd> TAz(vertcanopies->TAz, nl_can);

    Map<VectorXd> Sh2o_prof(vertcanopies->Sh2o_prof, nl_can);
    Map<VectorXd> Smaxz(vertcanopies->Smaxz, nl_can);
    Map<VectorXd> wetfrac(vertcanopies->wetfrac, nl_can);
    Map<VectorXd> dryfrac(vertcanopies->dryfrac, nl_can);

    Map<VectorXd> An_sun(vertcanopies->An_sun, nl_can);
    Map<VectorXd> An_shade(vertcanopies->An_shade, nl_can);
    Map<VectorXd> fsv_sun(vertcanopies->fsv_sun, nl_can);
    Map<VectorXd> fsv_shade(vertcanopies->fsv_shade, nl_can);
    Map<VectorXd> LE_sun(vertcanopies->LE_sun, nl_can);
    Map<VectorXd> LE_shade(vertcanopies->LE_shade, nl_can);
    Map<VectorXd> H_sun(vertcanopies->H_sun, nl_can);
    Map<VectorXd> H_shade(vertcanopies->H_shade, nl_can);
    Map<VectorXd> TR_sun(vertcanopies->TR_sun, nl_can);
    Map<VectorXd> TR_shade(vertcanopies->TR_shade, nl_can);
    Map<VectorXd> Evap_sun(vertcanopies->Evap_sun, nl_can);
    Map<VectorXd> Evap_shade(vertcanopies->Evap_shade, nl_can);
    Map<VectorXd> psil_sun(vertcanopies->psil_sun, nl_can);
    Map<VectorXd> psil_shade(vertcanopies->psil_shade, nl_can);
    Map<VectorXd> gbv_sun(vertcanopies->gbv_sun, nl_can);
    Map<VectorXd> gbv_shade(vertcanopies->gbv_shade, nl_can);
    Map<VectorXd> gbh_sun(vertcanopies->gbh_sun, nl_can);
    Map<VectorXd> gbh_shade(vertcanopies->gbh_shade, nl_can);
    Map<VectorXd> Ch2o_sun(vertcanopies->Ch2o_sun, nl_can);
    Map<VectorXd> Ch2o_shade(vertcanopies->Ch2o_shade, nl_can);

    Map<VectorXd> Rnrad_sun(vertcanopies->Rnrad_sun, nl_can);
    Map<VectorXd> Rnrad_shade(vertcanopies->Rnrad_shade, nl_can);
    Map<VectorXd> Evap_prof(vertcanopies->Evap_prof, nl_can);
    Map<VectorXd> Ch2o_prof(vertcanopies->Ch2o_prof, nl_can);


    // Vertical distribution of photosynthetic capacity
    PH_Dist(canopies, vertcanopies, photosynthesis, Vz);

    // Wind Profile
    First_Order_Closure_U(forcings, canopies, vertcanopies, microenviron, Uz, Km);

    // Canoy Precipitation Interception
    Precip_Interception(forcings, canopies, vertcanopies, vertsoils, radiation, Sh2o_prof, Smaxz,
            wetfrac, dryfrac, &ppt_ground);
    *vertsoils->ppt_ground = ppt_ground;

    // Shortwave radiation Profile
    ShortWaveradiation(forcings, canopies, vertcanopies, vertsoils, constants, radiation, fsun,
            fshade, LAIsun, LAIshade, SWabs_sun, SWabs_shade, PARabs_sun, PARabs_shade,
            NIRabs_sun, NIRabs_shade, taud, &SWabs_soil, &SWout);

    *vertcanopies->SWout = SWout;
    *vertsoils->SWabs_soil = SWabs_soil;

    // Longwave Convergence Loop
    converged_LW = 0;
    cnt_LW = 0;
    
    // The more iteration, the more unstable solution due to nonlinearity. 
    // Or direct solution by inposing ground heat flux might help.
    maxiters = 5;   //20;
    percdiff = 0.5; //0.01;

    while (converged_LW == 0)
    {
        // LONGWAVE rADIATION ABSORPTION
        LongWaveradiation(forcings, canopies, vertcanopies, vertsoils,
            constants, radiation, LWabs_can, LWabs_sun, LWabs_shade, &LWabs_soil,
            &LWout, LWemit_can, LWemit_sun, LWemit_shade, &LWemit_soil, LW_eq);

        // TOTAL ABSORBED rADIATION [W/m^2 ground]
        Totabs_sun = LWabs_sun + SWabs_sun;
        Totabs_shade = LWabs_shade + SWabs_shade;

        // TOTAL ABSORBED rADIATION PER UNIT LEAF AREA [W/m^2 leaf area]
        Rabs_sun_lai = Totabs_sun.cwiseQuotient(LAIsun);
        Rabs_shade_lai = Totabs_shade.cwiseQuotient(LAIshade);

        // SOIL ABSORBED ENERGY
        Totabs_soil = SWabs_soil + LWabs_soil;
        Rnrad_soil = Totabs_soil - LWemit_soil;
        *vertsoils->Totabs_soil = Totabs_soil;

        //==================================================================
        //                   SHADED CANOPY SOLUTION
        //   Calculations performed per [m^2 ground area], and canopy fluxes
        //   are calculated by integrating over the shaded leaf area
        //==================================================================
        double sunlit = 0;
        LeafSolution(forcings, canopies, vertcanopies, vertsoils, photosynthesis, constants,
                radiation, stomaconduct, Ph_shade, An_shade, Ci_shade, gsv_shade, Tl_shade,
                LE_shade, TR_shade, Evap_shade, H_shade, psil_shade, fsv_shade, Ch2o_shade,
                gbv_shade, gbh_shade, sunlit, t);

        //==================================================================
        //                       SUNLIT CANOPY SOLUTION
        //   Calculations performed per [m^2 leaf area], and canopy fluxes
        //   are calculated by integrating over the sunlit leaf area
        //==================================================================
        if (fsun.sum()==0)
        {   // under nocturnal conditions all leaf area is
            // considered to be shaded
            An_sun.setZero();
            LE_sun.setZero();
            H_sun.setZero();
            gsv_sun = gsv_shade;
            gbv_sun  = gbv_shade;
            gbh_sun = gbh_shade;
            Ci_sun = Ci_shade;
            Tl_sun = Tl_shade;
            psil_sun = psil_shade;
            fsv_sun = fsv_shade;
            Evap_sun.setZero();
            TR_sun.setZero();
            Ch2o_sun.setZero();
        }
        else
        {
            sunlit = 1;
            LeafSolution(forcings, canopies, vertcanopies, vertsoils, photosynthesis, constants,
                    radiation, stomaconduct, Ph_sun, An_sun, Ci_sun, gsv_sun, Tl_sun, LE_sun,
                    TR_sun, Evap_sun, H_sun, psil_sun, fsv_sun, Ch2o_sun, gbv_sun, gbh_sun,
                    sunlit, t);
        }
        Fc_soil = Ro * pow(Q10, (Tsoil - 10)/10);

        Soil_Surface_Fluxes(vertcanopies, canopies, vertsoils, soils, radiation, constants,
                &T_surf, &H_soil, &LE_soil, &G, &RH_soil, &LWups);

        *vertsoils->Fc_soil = Fc_soil;
        *vertsoils->T_surf = T_surf;
        *vertsoils->LE_soil = LE_soil;
        *vertsoils->H_soil = H_soil;
        *vertsoils->G = G;

        if (turb_on == 1)
            Micro_Environment(vertcanopies, canopies, vertsoils, constants, CAz, EAz, TAz, nl_can);

        // Test Longwave Convergence
        cnt_LW += 1;
        if (cnt_LW > 1)
        {
            diffprof = LWabs_can - LWabs_prev;
            // This should be a absolute value
            //percdiffprof = diffprof.cwiseQuotient(LWabs_prev);
            percdiffprof = diffprof.cwiseQuotient(LWabs_prev).cwiseAbs();     
            if (percdiffprof.maxCoeff() < percdiff)
                converged_LW = 1;
        }
        LWabs_prev = LWabs_can;

        if (cnt_LW > maxiters && converged_LW == 0)
            break;
    }

    // Net radiation
    Rnrad_eco   = (forcings->rg - SWout) + (forcings->lwdn - LWout);
    Rnrad_sun   = SWabs_sun + LWabs_sun  - LWemit_sun;
    Rnrad_shade = SWabs_shade + LWabs_shade  - LWemit_shade;

    // H2O storage on foliage: Precipitation and Condensation
    Evap_prof   = (Evap_sun.cwiseProduct(LAIsun) + Evap_shade.cwiseProduct(LAIshade)) * dtime;
    Evap_can    = Evap_prof.sum();
    Ch2o_prof   = -(Ch2o_sun.cwiseProduct(LAIsun) + Ch2o_shade.cwiseProduct(LAIshade)) * dtime;
    Ch2o_can    = Ch2o_prof.sum();

    *vertcanopies->Evap_can = Evap_can;
    *vertcanopies->Ch2o_can = Ch2o_can;

    Evap_Condensation_Adjust(vertcanopies, canopies, vertsoils, Sh2o_prof, &Sh2o_can, wetfrac, dryfrac, nl_can);

    *vertcanopies->Sh2o_can = Sh2o_can;

    Ch2o_prof = Map<VectorXd>(vertcanopies->Ch2o_prof, nl_can);
    Ch2o_can = Ch2o_prof.sum();                                  // Condensation water in the canopy [mm] 
    *vertcanopies->Ch2o_can = Ch2o_can;

    Evap_prof = Map<VectorXd>(vertcanopies->Evap_prof, nl_can);
    Evap_can = Evap_prof.sum();                                  // Total canopy evaporation [mm] 
    *vertcanopies->Evap_can = Evap_can;
	//-------------------------------

    An_can    = (An_sun.cwiseProduct(LAIsun)+An_shade.cwiseProduct(LAIshade)).sum();
    LE_can    = (LE_sun.cwiseProduct(LAIsun)+LE_shade.cwiseProduct(LAIshade)).sum();
    H_can     = (H_sun.cwiseProduct(LAIsun)+H_shade.cwiseProduct(LAIshade)).sum();
    TR_can    = (TR_sun.cwiseProduct(LAIsun)+TR_shade.cwiseProduct(LAIshade)).sum();
    Rnrad_can = Rnrad_sun.sum() + Rnrad_shade.sum();

    *vertcanopies->An_can    = An_can;
    *vertcanopies->LE_can    = LE_can;
    *vertcanopies->H_can     = H_can;
    *vertcanopies->TR_can    = TR_can;
    *vertcanopies->Rnrad_can = Rnrad_can;
    //*vertsoils->E_soil       = max(*vertcanopies->LE_can/constants->Lv_g, 0.0); // Soil evaporation [mm/s]=[g/m^2/s]
    *vertsoils->E_soil       = *vertcanopies->LE_can/constants->Lv_g; // Soil evaporation (positive) or condensation (negative) [mm/s]=[g/m^2/s]

    outmlcan->An_can[t]      = An_can;
    outmlcan->LE_can[t]      = LE_can;
    outmlcan->H_can[t]       = H_can;
    outmlcan->TR_can[t]      = TR_can;                                          // Transpiration 
    outmlcan->Rnrad_can[t]   = Rnrad_can;

    for (int i=0; i<nl_can; i++)
    {
        outmlcan->LAI_sun[t*nl_can+i]   = LAIsun[i];
        outmlcan->LAI_shade[t*nl_can+i] = LAIshade[i];
        outmlcan->An_sun[t*nl_can+i]    = An_sun[i]*LAIsun[i];
        outmlcan->An_shade[t*nl_can+i]  = An_shade[i]*LAIshade[i];
        outmlcan->LE_sun[t*nl_can+i]    = LE_sun[i]*LAIsun[i];
        outmlcan->LE_shade[t*nl_can+i]  = LE_shade[i]*LAIshade[i];
        outmlcan->H_sun[t*nl_can+i]     = H_sun[i]*LAIsun[i];
        outmlcan->H_shade[t*nl_can+i]   = H_shade[i]*LAIshade[i];
        outmlcan->TR_sun[t*nl_can+i]    = TR_sun[i]*LAIsun[i];
        outmlcan->TR_shade[t*nl_can+i]  = TR_shade[i]*LAIshade[i];
        outmlcan->Evap_prof[t*nl_can+i] = Evap_prof[i];
    }

    // Canopy water balance
    if (t == 0) 
    {
        *vertcanopies->mbw_can = 0;
    }
    else {
        *vertcanopies->mbw_can = (-Evap_can / dtime + Ch2o_can / dtime + forcings->ppt / dtime - *vertsoils->ppt_ground / dtime - (Sh2o_can - *vertcanopies->Sh2o_can_prev) / dtime) * dtime; // [mm]  
    }
    outmlcan->mbw_can[t] = *vertcanopies->mbw_can;
    *vertcanopies->Sh2o_can_prev = Sh2o_can;
}


/*
void RootConductivities(SwitchClass *Switches, VerticalCanopyClass *vertcanopies,
                        SoilClass *soils, VerticalSoilClass *vertsoils, Ref<VectorXd> krad,
                        Ref<VectorXd> kax, double *etr, int nl_soil)
{
    int rhc      = Switches->RootHydrauConduct;
    int hr       = Switches->HydraulicRedistribution;
    int rtcond   = Switches->RootConductivity;
    //double TR    = *vertcanopies->TR_can;
    double K_rad = soils->K_rad;
    double K_axs = soils->K_axs;
    VectorXd kaxs(nl_soil);

    // Copy array structs to eigen vectors
    VectorXd smp          = Map<VectorXd>(vertsoils->smp, nl_soil);
    VectorXd rpp          = Map<VectorXd>(vertsoils->rpp, nl_soil);
    VectorXd volliq       = Map<VectorXd>(vertsoils->volliq, nl_soil);
    VectorXd z            = Map<VectorXd>(vertsoils->zns, nl_soil) * 1000;
    VectorXd dz           = Map<VectorXd>(vertsoils->dzs, nl_soil) * 1000;
    VectorXd rootfr       = Map<VectorXd>(vertsoils->rootfr, nl_soil);
    VectorXd thetadry     = Map<VectorXd>(vertsoils->theta_dry, nl_soil);
    VectorXd eff_porosity = Map<VectorXd>(vertsoils->porosity, nl_soil);

    // Compute the radial and axial conductivities of the roots
    krad = (K_rad*rootfr.cwiseProduct(volliq)).cwiseQuotient(eff_porosity);
    if (rtcond == 1)
    {
        kaxs = rootfr / rootfr[0] * K_axs;
        kaxs = kaxs.cwiseQuotient(dz) * dz[0];
        kax = kaxs.cwiseQuotient(dz);
    }
    else
    {
        kax = (K_axs*rootfr.cwiseQuotient(dz)).cwiseProduct(volliq.cwiseQuotient(eff_porosity));
    }

    // For the case where the root hydraulic conductivity is allowed to increase
    // with depth, a linear increasing effect is considered
    if (rhc == 1)
    {
        krad += z.cwiseProduct(krad)/1000;
        kax += z.cwiseProduct(kax)/1000;
    }

    if (hr == 0)
    {
        for (int i=0; i<nl_soil; i++)
        {
            if (smp[i] < rpp[i])
                krad[i] = 0.0;
        }
        if (krad.maxCoeff() == 0.0)
            *etr = 0;
    }
}


void RootModel( VerticalSoilClass *vertsoils, Ref<VectorXd> rpp, double etr, int nl_soil)
{
    int i;
    double den, den1, den2;
    VectorXd aa(nl_soil), bb(nl_soil), cc(nl_soil), dd(nl_soil), ff(nl_soil);

    // Copy array structs to eigen vectors
    VectorXd krad = Map<VectorXd>(vertsoils->krad, nl_soil);
    VectorXd kax  = Map<VectorXd>(vertsoils->kax, nl_soil);
    VectorXd z    = Map<VectorXd>(vertsoils->zns, nl_soil) * 1000;  // unit [mm]
    VectorXd smp  = Map<VectorXd>(vertsoils->smp, nl_soil);

    // For the top soil layer
    i     = 0;
    den   = z[i+1] - z[i];
    aa[i] = 0;
    bb[i] = kax[i]/den + krad[i];
    cc[i] = -kax[i]/den;
    dd[i] = krad[i]* smp[i] - etr - kax[i];

    // For the middile soil layers
    for (i=1; i<nl_soil-1; i++)
    {
        den1  = z[i] - z[i-1];
        den2  = z[i+1] - z[i];
        aa[i] = -kax[i-1]/den1;
        bb[i] = kax[i-1]/den1 + kax[i]/den2 + krad[i];
        cc[i] = -kax[i]/den2;
        dd[i] = krad[i]* smp[i] + kax[i-1] - kax[i];
    }

    // For the bottom soil layer
    i     = nl_soil-1;
    den   = z[i] - z[i-1];
    aa[i] = -kax[i-1]/den;
    bb[i] = kax[i-1]/den + krad[i];
    cc[i] = 0;
    dd[i] = krad[i]* smp[i] + kax[i-1];

    ThomasAlgorithm(aa, bb, cc, dd, ff, nl_soil);

    rpp = ff;
    if (krad.maxCoeff() == 0.0)
        rpp = VectorXd::Zero(nl_soil);
}

void vanGenuchten(SoilClass *soils, VectorXd& C, VectorXd& K, VectorXd& Ksat, VectorXd& theta,
                  VectorXd& h, int SIZE)
{
    double theta_S = soils->theta_S;
    double theta_R = soils->theta_R;
    double alpha = soils->alpha;
    double n = soils->nv;
    double m = (n - 1.0)/n;

    double Se;
    double h_, theta_;

    for (int i=0; i<SIZE; i++)
    {
        h_ = h[i] * 0.1;            // [cm] Convert to centimeter

        // . . .Compute the volumetric moisture content [eqn 21] . . .
        if (h_ < 0)
        {
            theta[i] = (theta_S - theta_R) / pow(1.0 + pow((alpha*(-h_)),n), m) + theta_R;
        } else{
            theta[i] = theta_S;
        }

        // . . .Compute the effective saturation [eqn 2] . . .
        Se = (theta[i] - theta_R)/(theta_S - theta_R);  // [-]

        // . . .Compute the hydraulic conductivity [eqn 8] . . .[ Convert to unit: mm/s ]
        K[i] = Ksat[i]*sqrt(Se)*(1.0 - pow(1.0-pow(Se,1.0/m),m))*(1.0 - pow(1.0-pow(Se,1.0/m),m));

        // . . .Compute the specific moisture storage (derivative of eqn 21: C = d(theta)/dh . . .
        if (h_ < 0)
        {
            C[i] = -alpha * n * -1 * (1.0/n-1.0)*pow(alpha*abs(h_),n-1) * (theta_R-theta_S) 
                    * pow(pow(alpha*abs(h_),n)+1,1.0/n-2.0) * 0.10; //[ Unit: 1/mm ]
        } else {
            C[i] = 0.0;
        }
    }
}

void vanGenuchten_inverse(SoilClass *soils, VectorXd& theta, VectorXd& h, int SIZE)
{
    double theta_S = soils->theta_S;
    double theta_R = soils->theta_R;
    double alpha = soils->alpha;
    double n = soils->nv;
    double m = (n - 1.0)/n;

    for (int i=0; i<SIZE; i++)
    {
        if (theta[i] < theta_S){
            // in [mm]
            h[i] = -(1/alpha) * pow(pow((theta_S-theta_R)/(theta[i]-theta_R),1/m)-1.0,1/n) * 10;  
        } else {
            h[i] = 0;
        }
    }
}

void Set_Matrix(SoilClass *soils, VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& d, 
                VectorXd& Cnp1m, VectorXd& Knp1m, VectorXd& hnp1m, VectorXd& hn, 
                VectorXd& thetanp1m, VectorXd& thetan, VectorXd& TR, VectorXd& krad,
                VectorXd& rpp, VectorXd& dz, int P, int TopBound, double htop, double qin,
                int BotBound, double hbottom, double qout, double dt)
{
    double Ss = soils->Ss;
    double poros = soils->theta_S;
    double Knp1m_up;
    double Knp1m_down;

    for (int k=0; k<P; k++)
    {
        if (k==0 || k==P-1)
        {
            Knp1m_up = Knp1m[k];
            Knp1m_down = Knp1m[k];
        } else {
            Knp1m_up = (Knp1m[k+1]-Knp1m[k])/(dz[k]+dz[k+1])*dz[k] + Knp1m[k];
            Knp1m_down = (Knp1m[k]-Knp1m[k-1])/(dz[k]+dz[k-1])*dz[k-1] + Knp1m[k-1];
        }

        if (k==0 || k==P-1)
        {
            b[k] = (Cnp1m[k]+Ss*thetanp1m[k]/poros)/(dt) + krad[k]/dz[k];
            a[k] = 0.0;
            c[k] = 0.0;
        } else {
            b[k] = (Cnp1m[k]+Ss*thetanp1m[k]/poros)/(dt) + krad[k]/dz[k]
                  + (1/dz[k])*(2*Knp1m_up/(dz[k+1]+dz[k]) + 2*Knp1m_down/(dz[k]+dz[k-1]));
            a[k] = -2.0*Knp1m_down/(dz[k]+dz[k-1])/dz[k];
            c[k] = -2.0*Knp1m_up/(dz[k+1]+dz[k])/dz[k];
        }
        d[k] = Cnp1m[k]/(dt)*hnp1m[k] + (Ss*thetanp1m[k]/poros)*hn[k]/(dt)
             - (1.0/dz[k])*(Knp1m_up - Knp1m_down) - (thetanp1m[k] - thetan[k])/(dt)
             + krad[k]/dz[k] * rpp[k];

        if (TopBound == 0) 
        {
            if (k==0)
            {
                b[k] = 1;
                d[k] = htop;
            }
            if (k==1)
            {
                a[k] = 0;
                d[k] += 2.0*Knp1m_down/(dz[k]+dz[k-1])/dz[k] * htop;
            }
        } else {
            if (k==0)
            {
                b[k] += 2.0*Knp1m_up/(dz[k+1]+dz[k]) / dz[k];
                c[k] = -2.0*Knp1m_up/(dz[k+1]+dz[k]) / dz[k];
                d[k] += ( -Knp1m[k] + qin ) / dz[k];
            }
        }

        if (BotBound == 0)
        {
            if (k==P-1)
            {
                b[k] = 1;
                d[k] = hbottom;
            }
            
            if (k==P-2)
            {
                c[k] = 0;
                d[k] += 2.0*Knp1m_up/(dz[k+1]+dz[k])/dz[k] * hbottom;
            }
        } else {
            if (k==P-1)
            {
                b[k] += 2.0*Knp1m_down/(dz[k]+dz[k-1])/dz[k];
                a[k] = -2.0*Knp1m_down/(dz[k]+dz[k-1])/dz[k];
                d[k] += ( Knp1m[k] + qout ) / dz[k];
            }
        }
    }
}


void SoilMoistureModel(SoilClass *soils, VerticalSoilClass *vertsoils,
    VectorXd& hn, VectorXd& thetan, VectorXd& dsmp, VectorXd& dwat, VectorXd& Ksat,
    VectorXd& TR, double *qss, VectorXd& krad, VectorXd& rpp, double PPT, double ET,
    double *PH, double *htb, double *diff, int *tbc, int *case_out, int BotBound, double hbottom,
    double qout, double dt, int t, int P)
{
    int maxiter = soils->maxiter;
    double psimin = soils->psimin;
    double air_dry = soils->air_dry;
    double theta_S = soils->theta_S;
    double stop_tol = soils->stop_tol;

    int stop_flag, niter, TopBound, cases;
    double deficit, qpotential, qin, htop;
    VectorXd hnp1m(P), hnp1mp1(P), thetanp1m(P), Knp1m(P), Cnp1m(P);
    VectorXd a_z(P), b_z(P), c_z(P), d_z(P), f_z(P);
    VectorXd dz = Map<VectorXd>(vertsoils->dzs, P) * 1000;

    hnp1m     = hn;
    stop_flag = 0;
    niter     = 0;
    qin       = 0;

    while (stop_flag == 0 && niter < maxiter)
    {
        htop = *htb;
        if (*PH > psimin)
        {
            // Ponding conditions,  Dirichlet BC type
            TopBound = 0;
            htop     = *PH + PPT + ET;
            cases    = 1;
        } else {
            // Non-ponding condition, including 3 sub-cases
            if (*PH > 0)
            {
                // 1. Saturated condition,  Dirichlet BC type
                TopBound = 0;
                htop     = *PH + PPT + ET;
                cases    = 2;
            } else {
                // 2. Unsaturated or Air-dry condition
                if (htop > air_dry)
                {
                    // Unsaturated condition, nagative pressure but greater than air-dry threshold
                    deficit    = (theta_S - thetan[0]) * dz[0];
                    qpotential = *PH + PPT + ET;
                    if (qpotential > deficit) 
                    {
                        // More rainfall than soil water deficit, Dirichlet BC type
                        TopBound = 0;
                        htop     = hn[0] + qpotential - deficit;
                        cases    = 3;
                    } else {
                        // Less rainfall than soil water deficit
                        if (htop + PPT + ET < Ksat[0]*dt) 
                        {
                            // Rain rate < Ksat, Neumann BC type
                            TopBound = 1;
                            qin      = (*PH + PPT + ET) / dt;    // infiltration rate [mm/s]
                            htop     = 0.0;
                            cases    = 4;
                        } else {
                            // Rain rate > Ksat, Dirichlet BC type
                            TopBound = 0;
                            htop     = *PH + PPT + ET - Ksat[0]*dt;
                            cases    = 5;
                        }
                    }
                } else {
                    // 3. Air-dry condition, Dirichlet BC type
                    TopBound = 0;
                    htop = air_dry + 2000;
                    cases = 6;
                }
            }
        }

        vanGenuchten(soils, Cnp1m, Knp1m, Ksat, thetanp1m, hnp1m, P);
        
        //qout = 0.0;                 // no-flow
        qout = -Knp1m[P-1];         // free-flow
        Set_Matrix(soils, a_z, b_z, c_z, d_z, Cnp1m, Knp1m, hnp1m, hn, thetanp1m, thetan, TR,
                krad, rpp, dz, P, TopBound, htop, qin, BotBound, hbottom, qout, dt);
        
        ThomasAlgorithm(a_z, b_z, c_z, d_z, hnp1mp1, P);
        niter += 1;

        f_z = hnp1mp1 - hnp1m;
        if (MaxValue(f_z, P) < stop_tol)
        {
            stop_flag = 1;
            vanGenuchten(soils, Cnp1m, Knp1m, Ksat, thetanp1m, hnp1mp1, P);
        } else {
            hnp1m = hnp1mp1;
        }

        *diff = qin;

        if (TopBound == 1){
            *qss = qin;
            htop = hn[0] - dz[0] + *qss * dz[0] / Knp1m[0];
        } else {
            *qss = -Knp1m[0] * (hn[0] - htop - dz[0]) / dz[0];
            htop -= *qss*dt;
        }
    }

    if (htop > *PH + PPT)
        htop = *PH + PPT;

    *PH = htop;
    *htb = htop;
    if (*PH < 0)
        *PH = 0;

    *tbc = TopBound;
    *case_out = cases;
    if (*PH < 0)   
        *PH = 0;

    // *PH = 0;
    dwat = thetanp1m - thetan;
    dsmp = hnp1m - hn;
}


void RootSoilModel(SwitchClass *Switches, ConstantClass *constants, SoilClass *soils,
                   VerticalSoilClass *vertsoils, VerticalCanopyClass *vertcanopies, int t)
{
    double converge, difprev, diff;
    int cnt, tbc, case_out;

    int nl_soil = soils->nl_soil;
    int maxiter = soils->maxiter;
    int hr      = Switches->HydraulicRedistribution;
    int rhc     = Switches->RootHydrauConduct;

    double dt         = constants->dtime;
    double etr        = vertcanopies->TR_can[0];
    double K_rad      = soils->K_rad;
    double K_axs      = soils->K_axs;
    double kpar_ax    = soils->kpar_ax;
    double Ss         = soils->Ss;
    double qss        = vertsoils->qss[0];
    double PH         = vertsoils->PH[0];
    double htb        = vertsoils->htb[0];
    double E_soil     = -vertsoils->E_soil[0]*dt;
    double ppt_ground = vertsoils->ppt_ground[0];

    // Create local Eigen Vectors for operations
    VectorXd z(nl_soil), dz(nl_soil), zhs(nl_soil), zi(nl_soil);
    VectorXd rootfr(nl_soil), rpp(nl_soil), smp(nl_soil), dsmp(nl_soil);
    VectorXd dwat(nl_soil), theta(nl_soil);
    VectorXd krad(nl_soil), kax(nl_soil), Ksat(nl_soil);
    VectorXd TR_soil(nl_soil), TRodz(nl_soil);

    // Copy array structs to eigen vectors
    z      = Map<VectorXd>(vertsoils->zns, nl_soil) * 1000;
    zhs    = Map<VectorXd>(vertsoils->zhs, nl_soil) * 1000;
    dz     = Map<VectorXd>(vertsoils->dzs, nl_soil) * 1000;
    zi     = Map<VectorXd>(vertsoils->zhs, nl_soil);
    rootfr = Map<VectorXd>(vertsoils->rootfr, nl_soil);
    rpp    = Map<VectorXd>(vertsoils->rpp, nl_soil);
    smp    = Map<VectorXd>(vertsoils->smp, nl_soil);
    theta  = Map<VectorXd>(vertsoils->volliq, nl_soil);
    Ksat   = Map<VectorXd>(vertsoils->HKsat, nl_soil);
    RootConductivities(Switches, vertcanopies, soils, vertsoils, krad, kax, &etr, nl_soil);

    // Copy eigen vectors back to array structs
    Map<VectorXd>(vertsoils->krad, nl_soil) = krad;     // [1/s]
    Map<VectorXd>(vertsoils->kax, nl_soil)  = kax;      // [mm/s]

    TR_soil = etr * rootfr;
    TRodz = -TR_soil.cwiseQuotient(dz);

    converge = 10e10;
    cnt = 0;
    while (converge > 0.1) 
    {
        cnt++;
        if (cnt > 100)
            break;

        RootModel(vertsoils, rpp, etr, nl_soil);         // rpp is in [mm]
        Map<VectorXd>(vertsoils->rpp, nl_soil) = rpp;

        RootConductivities(Switches, vertcanopies, soils, vertsoils, krad, kax, &etr, nl_soil);

        SoilMoistureModel(soils, vertsoils, smp, theta, dsmp, dwat, Ksat, TRodz, &qss, krad, rpp,
                ppt_ground, E_soil, &PH, &htb, &diff, &tbc, &case_out, 1, -40000.0, 0.0, dt, t,
                nl_soil);

        RootConductivities(Switches, vertcanopies, soils, vertsoils, krad, kax, &etr, nl_soil);

        if (cnt == 1)
        {
            converge = 1000;
        } else {
            converge = abs(abs(etr - ((smp-rpp).cwiseProduct(krad)).sum()) - difprev) / difprev;
        }
        difprev = abs(etr - ((smp-rpp).cwiseProduct(krad)).sum());
        converge = 0.0;
    }

    smp   += dsmp;
    theta += dwat;
    // Copy eigen vectors back to array structs
    Map<VectorXd>(vertsoils->smp, nl_soil)    = smp;
    Map<VectorXd>(vertsoils->volliq, nl_soil) = theta;
    Map<VectorXd>(vertsoils->dsmp, nl_soil)   = dsmp;
    Map<VectorXd>(vertsoils->dwat, nl_soil)   = dwat;
    Map<VectorXd>(vertsoils->rpp, nl_soil)    = rpp;
    Map<VectorXd>(vertsoils->krad, nl_soil)   = krad;    // [1/s]
    Map<VectorXd>(vertsoils->kax, nl_soil)    = kax;      // [mm/s]

    *vertsoils->qss   = qss;
    *vertsoils->PH    = PH;
    *vertsoils->htb   = htb;
    *vertsoils->tbc   = tbc;
    *vertsoils->diff  = diff;
    *vertsoils->cases = case_out;
}


void Set_Heat_Matrix(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& d, VectorXd& Ts, 
                     VectorXd& lambda, VectorXd& cv, VectorXd& dz, int P, int TopBound, 
                     double Ttop, double g_in, int BotBound, double Tbottom, double g_out, 
                     double dt)
{
    double Kd_up;
    double Kd_down;
    VectorXd Kd = lambda.cwiseQuotient(cv);

    for (int k=0; k<P; k++)
    {
        if (k==0 || k==P-1)
        {
            Kd_up = Kd[k];
            Kd_down = Kd[k];
        } else {
            Kd_up = (Kd[k+1]-Kd[k])/(dz[k]+dz[k+1])*dz[k] + Kd[k];
            Kd_down = (Kd[k]-Kd[k-1])/(dz[k]+dz[k-1])*dz[k-1] + Kd[k-1];
        }

        if (k==0 || k==P-1)
        {
            b[k] = 1.0/dt;
            a[k] = 0.0;
            c[k] = 0.0;
        } else {
            b[k] = 1/(dt) + (1/dz[k])*(2*Kd_up/(dz[k+1]+dz[k]) + 2*Kd_down/(dz[k]+dz[k-1]));
            a[k] = -2.0*Kd_down/(dz[k]+dz[k-1])/dz[k];
            c[k] = -2.0*Kd_up/(dz[k+1]+dz[k])/dz[k];
        }
        d[k] = 1/(dt) * Ts[k];


        if (TopBound == 0)
        {
            if (k==0)
            {
                b[k] = 1;
                d[k] = Ttop;
            }
            
            if (k==1)
            {
                a[k] = 0;
                d[k] += 2.0*Kd_down/(dz[k]+dz[k-1])/dz[k] * Ttop;
            }
        } else {
            if (k==0)
            {
                b[k] += 2.0*Kd_up/(dz[k+1]+dz[k]) / dz[k];
                c[k] = -2.0*Kd_up/(dz[k+1]+dz[k]) / dz[k];
                d[k] += 2.0*Kd_up/(dz[k+1]+dz[k]) * g_in / lambda[k];
            }
        }

        if (BotBound == 0) 
        {
            if (k==P-1)
            {
                b[k] = 1;
                d[k] = Tbottom;
            }

            if (k==P-2)
            {
                c[k] = 0;
                d[k] += 2.0*Kd_up/(dz[k+1]+dz[k])/dz[k] * Tbottom;
            }
        } else {
            if (k==P-1)
            {
                b[k] += 2.0*Kd_down/(dz[k]+dz[k-1])/dz[k];
                a[k] = -2.0*Kd_down/(dz[k]+dz[k-1])/dz[k];
                d[k] += -2.0*Kd_down/(dz[k]+dz[k-1]) * g_out / lambda[k];
            }
        }
    }
}

// ----------------------------------------------------------------------------
// SoilHeatModel()
//    This function implements the numerical solution for soil heat transport
//    described in Oleson et al (2004) for the Community Land Model (CLM).
//    Top BC: Hg is heat flux into top of soil column.
//    Bottom BC: zero heat flux at bottom of soil column.
//
//    The following equation is solved:
//        dT    d        dT
//      C -- = -- ( TK * -- )
//        dt   dz        dz
// ----------------------------------------------------------------------------
void SoilHeatModel(ConstantClass *constants, SoilClass *soils, VerticalSoilClass *vertsoils, int t)
{
    int P          = soils->nl_soil;            // Number of soil layers [-]
    double Tf      = soils->Tf;                 // Freezing temperature [K]
    double HC_liq  = soils->HC_liq;             // Heat capacity of liquid water [J/kg/K]
    double HC_ice  = soils->HC_ice;             // Heat capacity of ice [J/kg/K]
    double rho_liq = soils->rho_liq;            // Density of water [kg/m^3]
    double rho_ice = soils->rho_ice;            // Density of ice [kg/m^3]
    double TK_liq  = soils->TK_liq;             // Thermal conductivity of liquid water [W/m/K]
    double TK_ice  = soils->TK_ice;             // Thermal conductivity of ice [W/m/K]
    double Hg      = *vertsoils->G;             // Heat flux goes into soil [W/m^2]
    double dt      = constants->dtime;          // Time step interval [s]

    VectorXd TKsoil(P), TKsoil_h(P), TK_sat(P), Ke(P), Sr(P), fact(P), One(P);
    VectorXd Ts_new(P), theta_ice(P), cpv(P), aa(P), bb(P), cc(P), rr(P);

    One = VectorXd::Ones(P);
    // Copy array structs to eigen vectors
    VectorXd znode     = Map<VectorXd>(vertsoils->zns, P);
    VectorXd zlayer    = Map<VectorXd>(vertsoils->zhs, P);
    VectorXd dz        = Map<VectorXd>(vertsoils->dzs, P);
    VectorXd TK_dry    = Map<VectorXd>(vertsoils->TK_dry, P);
    VectorXd TK_sol    = Map<VectorXd>(vertsoils->TK_sol, P);
    VectorXd HC_sol    = Map<VectorXd>(vertsoils->HC_sol, P);
    VectorXd poros     = Map<VectorXd>(vertsoils->porosity, P);
    VectorXd theta_liq = Map<VectorXd>(vertsoils->volliq, P);
    VectorXd Ts        = Map<VectorXd>(vertsoils->Ts, P);

    theta_ice = VectorXd::Zero(P);

    // Volumetric Heat Capacity [J / m^3 / K]  (eqn 6.67)
    cpv = HC_sol.cwiseProduct(One - poros) + (theta_liq * rho_liq * HC_liq) 
          + (theta_ice * rho_ice * HC_ice);

    // Thermal Conductivities at nodes (eqn 6.60)
    // Saturated Thermal Conductivity [W / m / K]
    for (int k=0; k<P; k++) 
    {
        if (Ts[k] >= Tf) 
        {
            // non-frozen layers
            TK_sat[k] = pow(TK_sol[k], 1.0 - poros[k]) * pow(TK_liq, poros[k]);
        } else {
            // frozen layers
            TK_sat[k] = pow(TK_sol[k], 1.0 - poros[k]) * pow(TK_liq, poros[k])
                        * pow(TK_ice, poros[k]-theta_liq[k]);
        }
    }

    Sr = (theta_liq + theta_ice).cwiseQuotient(poros);
    if (MaxValue(Sr, P) > 1) 
        cout << "ERROR in Soil Heat Transport: Sr > 1 at time: " << t << endl;

    // Kersten Number (6.63)
    for (int k=0; k<P; k++) 
    {
        if (Ts[k] >= Tf) 
        {   // non-frozen layers
            Ke[k] = log10(Sr[k]) + 1.0;
        } else {            // frozen layers
            Ke[k] = Sr[k];
        }
        
        if (Ke[k] < 0)
            Ke[k] = 0;
    }

    // Soil Thermal Conductivity [W / m / K] (6.58)
    TKsoil = Ke.cwiseProduct(TK_sat) + (One - Ke).cwiseProduct(TK_dry);
    for (int k=0; k<P; k++) 
    {
        if (Sr[k] <= 1e-7)
            TKsoil[k] = TK_dry[k];
    }

    int TopBound     = 1;
    int BotBound     = 1;
    double Ttop      = 0.0;
    double Tbottom   = 0.0;
    double  qout_bot = 0.0;

    Set_Heat_Matrix(aa, bb, cc, rr, Ts, TKsoil, cpv, dz, P, TopBound, Ttop, Hg, BotBound, Tbottom,
            qout_bot, dt);

    // Call Tri-diagonal solver using Thomas algorithm
    ThomasAlgorithm(aa, bb, cc, rr, Ts_new, P);

    // Copy eigen vectors back to array structs
    Map<VectorXd>(vertsoils->Ts, P) = Ts_new;
}
*/