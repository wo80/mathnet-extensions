// Based on Fortran code RADAU5
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;

    using S_fp = System.Action<int, double, double[], double[]>;
    using J_fp = System.Action<int, double, double[], double[]>;

    /// <summary>
    /// Numerical solution of a stiff (or differential algebraic) system of first order
    /// ordinary differential equations
    /// 
    ///     M* Y'=F(X,Y).
    ///     
    /// The system can be (linearly) implicit (mass-matrix M != I) or explicit (M = I).
    /// 
    /// The method used is an implicit Runge-Kutta method (RADAU IIA) of order 5
    /// with step size control and continuous output.
    ///
    /// Authors: E. Hairer and G. Wanner
    ///
    /// This code is part of the book:
    ///         E. Hairer and G. Wanner
    ///         Solving Ordinary Differential Equations II.
    ///         Stiff and differential-algebraic problems. (2nd edition)
    ///         Springer-Verlag (1996)
    /// </summary>
    public class Radau5
    {
        const double t11 = 0.091232394870892942792;
        const double t12 = -0.14125529502095420843;
        const double t13 = -0.030029194105147424492;
        const double t21 = 0.24171793270710701896;
        const double t22 = 0.20412935229379993199;
        const double t23 = 0.38294211275726193779;
        const double t31 = 0.96604818261509293619;

        const double ti11 = 4.325579890063155351;
        const double ti12 = 0.33919925181580986954;
        const double ti13 = 0.54177053993587487119;
        const double ti21 = -4.1787185915519047273;
        const double ti22 = -0.32768282076106238708;
        const double ti23 = 0.47662355450055045196;
        const double ti31 = -0.50287263494578687595;
        const double ti32 = 2.5719269498556054292;
        const double ti33 = -0.59603920482822492497;

        const double SQRT6 = 2.449489742783178098197284;

        const double c1 = (4.0 - SQRT6) / 10.0;
        const double c2 = (SQRT6 + 4.0) / 10.0;

        const double dd1 = -(SQRT6 * 7.0 + 13.0) / 3.0;
        const double dd2 = (SQRT6 * 7.0 - 13.0) / 3.0;
        const double dd3 = -1.0 / 3.0;

        const double c1m1 = c1 - 1.0;
        const double c2m1 = c2 - 1.0;
        const double c1mc2 = c1 - c2;
        
        // a = (12 - cbrt(81) + cbrt(9)) / 60 = 0.16255558520216131612608569820029
        // b = (cbrt(81) + cbrt(9)) * sqrt(3) / 60 = 0.1849493244071407842750912237438
        // c = a * a + b * b = 0.0606305708790749288701403085596
        
        const double alph = 2.68108287362775213389579074321111; // a / c
        const double beta = 3.0504301992474105694263776247876; // b / c
        const double u1 = 3.6378342527444957322084185135778; // 30 / (cbrt(81) + 6 - cbrt(9))

        bool caljac;

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        double _xsol, _hsol;

        S_fp fcn;
        J_fp jac;

        double rtol, atol;

        int n;

        // Workspace
        double[] z1, z2, z3;
        double[] y0;
        double[] scal, cont;
        double[] f1, f2, f3;
        double[] fjac, fmas;
        double[] e1, e2r, e2i;
        int[] ip1, ip2;

        public int nstep, ndec, nsol, nrejct, naccpt;

        int nsing;

        // Subroutine
        public int radau5_(int n, S_fp fcn, double x, double[] y,
            double xend, double h, double rtol, double atol,
            J_fp jac, DenseMatrix mas, bool dense)
        {
            // Local variables
            int nit;
            double facl;
            double facr, safe;
            int ijob;
            bool pred;
            double hmax;
            int nmax;
            double thet, expm;
            double quot;
            int nind1, nind2, nind3;
            double quot1, quot2;
            double fnewt;
            double tolst;

            bool implct;
            bool startn;
            double uround;

            /**
             *     INPUT PARAMETERS
             *     ----------------
             *     N           DIMENSION OF THE SYSTEM
             *
             *     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
             *                 VALUE OF F(X,Y):
             *                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),F(N)
             *                    F(1)=...   ETC.
             *                 RPAR, IPAR (SEE BELOW)
             *
             *     X           INITIAL X-VALUE
             *
             *     Y(N)        INITIAL VALUES FOR Y
             *
             *     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
             *
             *     H           INITIAL STEP SIZE GUESS;
             *                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
             *                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
             *                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
             *                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
             *
             *     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
             *                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
             *
             *     ITOL        SWITCH FOR RTOL AND ATOL:
             *                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
             *                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
             *                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
             *                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
             *                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
             *                     RTOL(I)*ABS(Y(I))+ATOL(I).
             *
             *     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
             *                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
             *                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
             *                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
             *                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
             *                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N)
             *                    DFY(1,1)= 0...
             *                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
             *                 FURNISHED BY THE CALLING PROGRAM.
             *                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
             *                    BE FULL AND THE PARTIAL DERIVATIVES ARE
             *                    STORED IN DFY AS
             *                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
             *
             *     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
             *                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
             *                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
             *                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
             *
             *     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
             *     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
             *
             *     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
             *                 MATRIX M.
             *                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
             *                 MATRIX AND NEEDS NOT TO BE DEFINED;
             *                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
             *                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
             *                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
             *                    DOUBLE PRECISION AM(LMAS,N)
             *                    AM(1,1)= 0....
             *                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
             *                    AS FULL MATRIX LIKE
             *                         AM(I,J) = M(I,J)
             *                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
             *                    DIAGONAL-WISE AS
             *                         AM(I-J+MUMAS+1,J) = M(I,J).
             *
             *     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
             *                    IOUT=0: SUBROUTINE IS NEVER CALLED
             *                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
             *
             * ----------------------------------------------------------------------
             *
             *     OUTPUT PARAMETERS
             *     -----------------
             *     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
             *                 (AFTER SUCCESSFUL RETURN X=XEND).
             *
             *     Y(N)        NUMERICAL SOLUTION AT X
             *
             *     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
             *
             *     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
             *                   IDID= 1  COMPUTATION SUCCESSFUL,
             *                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
             *                   IDID=-1  INPUT IS NOT CONSISTENT,
             *                   IDID=-2  LARGER NMAX IS NEEDED,
             *                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
             *                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
             */

            nstep = 0;
            naccpt = 0;
            nrejct = 0;
            ndec = 0;
            nsol = 0;

            nsing = 0;

            this.n = n;

            this.fcn = fcn;
            this.jac = jac;

            // UROUND   Smallest number satisfying 1.0+UROUND>1.0
            //          (the rounding unit, default = 1e-16).
            uround = 1e-16;

            if (uround <= 1e-19 || uround >= 1.0)
            {
                Console.WriteLine(" COEFFICIENTS HAVE 20 DIGITS, UROUND=", uround);
                return -1;
            }

            // Check and change the tolerances
            expm = 0.66666666666666663;

            if (atol <= 0.0 || rtol <= uround * 10.0)
            {
                Console.WriteLine(" TOLERANCES ARE TOO SMALL");
                return -1;
            }
            else
            {
                quot = atol / rtol;
                rtol = Math.Pow(rtol, expm) * 0.1;
                atol = rtol * quot;
            }

            this.rtol = rtol;
            this.atol = atol;

            // NMAX , The maximal number of steps
            nmax = 100000;

            // NIT    The maximum number of Newton iterations for the solution of the implicit system in each step (default = 7).
            nit = 7;

            if (nit <= 0)
            {
                Console.WriteLine(" CURIOUS INPUT IWORK(3)=", nit);
                return -1;
            }

            // STARTN  Switch for starting values of newton iterations
            //    If STARTN = true, the extrapolated collocation solution
            //    is taken as starting value for Newton's method.
            //    If STARTN = false, zero starting values are used.
            //    The latter is recommended if Newton's method has
            //    difficulties with convergence (this is the case when
            //    NSTEP is larger than NACCPT + NREJCT; see output param 0).
            //    default STARTN = false.
            startn = false;

            // The following 3 parameters are important for
            // differential-algebraic systems of INDEX > 1.
            // The function-subroutine should be written such that
            // the index 1,2,3 variables appear in this order.
            // In estimating the error the index 2 variables are
            // multiplied by h, the index 3 variables by h**2.

            // NIND1  dimension of the index 1 variables (must be > 0).
            //        For ODE's this equals the dimension of the system. (default NIND1=N).
            // NIND2  dimension of the index 2 variables. (default NIND2=0).
            // NIND3  dimension of the index 3 variables. (default NIND3=0).
            nind1 = n;
            nind2 = 0;
            nind3 = 0;

            if (nind1 + nind2 + nind3 != n)
            {
                Console.WriteLine(" CURIOUS INPUT FOR IWORK(5,6,7)=", nind1, nind2, nind3);
                return -1;
            }

            // PRED   Switch for step size strategy (default = 1).
            //   IF PRED = 1  mod. predictive controller (Gustafsson)
            //   IF PRED = 2  CLASSICAL STEP SIZE CONTROL
            //  
            //   The choice PRED = 1 seems to produce safer results;
            //   For simple problems, the choice PRED = 2 produces
            //   often slightly faster runs.
            pred = true;

            // SAFE     The safety factor in step size prediction, default 0.9.
            safe = 0.9;

            if (safe <= 0.001 || safe >= 1.0)
            {
                Console.WriteLine(" CURIOUS INPUT FOR WORK(2)=", safe);
                return -1;
            }

            // THET     Decides whether the Jacobian should be recomputed (default 0.001).
            //    Increase THET, to 0.1 say, when Jacobian evaluations are costly. For
            //    small systems THET should be smaller (0.001, say). Negativ THET forces
            //    the code to compute the Jacobian after every accepted step.
            thet = 0.001;

            if (thet >= 1.0)
            {
                Console.WriteLine(" CURIOUS INPUT FOR WORK(3)=", thet);
                return -1;
            }

            // FNEWT   Stopping criterion for Newton's method, usually chosen < 1 (default min(0.03, RTOL^0.5)).
            //    Smaller values of FNEWT make the code slower, but safer.
            tolst = rtol;

            fnewt = Math.Max(uround * 10 / tolst, Math.Min(0.03, Math.Pow(tolst, 0.5)));

            if (fnewt <= uround / tolst)
            {
                Console.WriteLine(" CURIOUS INPUT FOR WORK(4)=", fnewt);
                return -1;
            }

            // QUOT1 and QUOT2: if QUOT1 < HNEW/HOLD < QUOT2, then the step size is not changed (defaults QUOT1 = 1.0, QUOT2 = 1.2).
            //    This saves, together with a large THET, LU-decompositions and computing time for
            //    large systems.
            //    For small systems one may have QUOT1 = 1.0, QUOT2 = 1.2.
            //    For large full systems QUOT1 = 0.99, QUOT2 = 2.0 might be good.
            quot1 = 1.0;
            quot2 = 1.2;

            if (quot1 > 1.0 || quot2 < 1.0)
            {
                Console.WriteLine(" CURIOUS INPUT FOR WORK(5,6)=", quot1, quot2);
                return -1;
            }

            // Maximal step size
            hmax = xend - x;

            // FACL,FACR     Parameters for step size selection

            // THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
            // WORK(8) <= HNEW / HOLD <= WORK(9)
            // DEFAULT VALUES: WORK(8) = 0.2, WORK(9) = 8.0

            facl = 5.0;
            facr = 0.125;
            //facl = 1.0 / work[8];
            //facr = 1.0 / work[9];

            if (facl < 1.0 || facr > 1.0)
            {
                Console.WriteLine(" CURIOUS INPUT WORK(8,9)=", 1.0 / facl, 1.0 / facr);
                return -1;
            }

            // Computation of array entries

            // Implicit, banded or not ?
            implct = mas != null;

            // Mass matrix
            if (implct)
            {
                ijob = 5;
            }
            else
            {
                ijob = 1;
            }

            // Prepare the entry-points for the arrays in work
            z1 = new double[n];
            z2 = new double[n];
            z3 = new double[n];
            y0 = new double[n];
            scal = new double[n];
            f1 = new double[n];
            f2 = new double[n];
            f3 = new double[n];
            cont = new double[n * 4];
            fjac = new double[n * n];
            fmas = new double[n * n]; // TODO: implicit -> 0
            e1 = new double[n * n];
            e2r = new double[n * n];
            e2i = new double[n * n];

            // Entry points for integer workspace
            ip1 = new int[n];
            ip2 = new int[n];

            // Call to core integrator
            int idid = Integrate(n, x, y, xend, hmax, h,
                mas, dense, nmax, uround, safe, thet, fnewt,
                quot1, quot2, nit, ijob, startn, nind1, nind2, nind3,
                pred, facl, facr, implct);

            // Restore tolerances
            expm = 1.0 / expm;

            quot = atol / rtol;
            rtol = Math.Pow(rtol * 10.0, expm);
            atol = rtol * quot;

            if (nstep > nmax)
            {
                Console.WriteLine("(  EXIT OF RADAU5 AT X={0} )", x);
                Console.WriteLine(" MORE THAN NMAX =" + nmax + "STEPS ARE NEEDED");
            }

            return idid;
        }

        // Core integrator for RADAU5
        int Integrate(int n, double x, double[]
            y, double xend, double hmax, double h,
            DenseMatrix mas, bool dense, int nmax,
            double uround, double safe, double thet, double
            fnewt, double quot1, double quot2, int nit, int
            ijob, bool startn, int nind1, int nind2, int nind3,
            bool pred, double facl, double facr, bool implct)
        {
            int i;
            int n2, n3, n4;
            double ak, ak1, ak2, ak3, c1q, c2q, c3q, z1i, z2i, z3i, fac;
            double xph, err = 0, cfac, hacc = 0;
            double hnew, hold;
            bool last;
            double hopt, xold;
            double quot, hhfac, theta = 0, hmaxn;
            bool first;
            double acont3;
            bool index1, index2, index3;
            double faccon;
            double erracc = 0;
            bool reject;
            double facgus;
            double posneg;
            
            // Check the index of the problem
            index1 = nind1 != 0;
            index2 = nind2 != 0;
            index3 = nind3 != 0;

            // Compute mass matrix for implicit case
            if (implct)
            {
                fmas = mas.ToColumnMajorArray();
            }
            
            posneg = d_sign(1.0, xend - x);
            hmaxn = Math.Min(Math.Abs(hmax), Math.Abs(xend - x));
            if (Math.Abs(h) <= uround * 10.0)
            {
                h = 1e-6;
            }
            h = Math.Min(Math.Abs(h), hmaxn);
            h = d_sign(h, posneg);
            hold = h;
            reject = false;
            first = true;
            last = false;
            if ((x + h * 1.0001 - xend) * posneg >= 0.0)
            {
                h = xend - x;
                last = true;
            }
            hopt = h;
            faccon = 1.0;
            cfac = safe * ((nit * 2) + 1);
            nsing = 0;
            xold = x;
            if (dense)
            {
                _xsol = x;
                _hsol = hold;
                for (i = 0; i < n; i++)
                {
                    cont[i] = y[i];
                }
            }

            n2 = n * 2;
            n3 = n * 3;
            n4 = n * 4;

            for (i = 0; i < n; i++)
            {
                scal[i] = atol + rtol * Math.Abs(y[i]);
            }

            hhfac = h;
            fcn(n, x, y, y0);

            ComputeJacobian(n, x, y);

            // Basic integration step
            while (nstep < nmax)
            {
                // Compute the matrices E1 and E2 and their decompositions
                if (!Factorize(h, ijob))
                {
                    h *= 0.5;
                    hhfac = 0.5;

                    reject = true;
                    last = false;

                    if (!caljac)
                    {
                        ComputeJacobian(n, x, y);
                    }

                    continue;
                }

                L30:
                nstep++;

                if (nstep > nmax)
                {
                    return -2;
                }

                if (Math.Abs(h) * 0.1 <= Math.Abs(x) * uround)
                {
                    throw new NumericalBreakdownException("Step size too small, h=" + h);
                }

                if (index2)
                {
                    int end = nind1 + nind2;
                    for (i = nind1; i < end; i++)
                    {
                        scal[i] /= hhfac;
                    }
                }

                if (index3)
                {
                    int end = nind1 + nind2 + nind3;
                    for (i = nind1 + nind2; i < end; i++)
                    {
                        scal[i] /= hhfac * hhfac;
                    }
                }

                xph = x + h;

                // Starting values for newton iteration
                if (first || startn)
                {
                    for (i = 0; i < n; i++)
                    {
                        z1[i] = 0.0;
                        z2[i] = 0.0;
                        z3[i] = 0.0;
                        f1[i] = 0.0;
                        f2[i] = 0.0;
                        f3[i] = 0.0;
                    }
                }
                else
                {
                    c3q = h / hold;
                    c1q = c1 * c3q;
                    c2q = c2 * c3q;

                    for (i = 0; i < n; i++)
                    {
                        ak1 = cont[i + n];
                        ak2 = cont[i + n2];
                        ak3 = cont[i + n3];
                        z1i = c1q * (ak1 + (c1q - c2m1) * (ak2 + (c1q - c1m1) * ak3));
                        z2i = c2q * (ak1 + (c2q - c2m1) * (ak2 + (c2q - c1m1) * ak3));
                        z3i = c3q * (ak1 + (c3q - c2m1) * (ak2 + (c3q - c1m1) * ak3));
                        z1[i] = z1i;
                        z2[i] = z2i;
                        z3[i] = z3i;
                        f1[i] = ti11 * z1i + ti12 * z2i + ti13 * z3i;
                        f2[i] = ti21 * z1i + ti22 * z2i + ti23 * z3i;
                        f3[i] = ti31 * z1i + ti32 * z2i + ti33 * z3i;
                    }
                }

                // Simplified newton iteration.
                int newt = Newton(n, nit, x, y, ref h, ref hhfac, thet, uround, fnewt, ijob, ref faccon, ref theta);

                if (newt < 0)
                {
                    reject = true;
                    last = false;

                    if (!caljac)
                    {
                        ComputeJacobian(n, x, y);
                    }

                    continue;
                }

                // Error estimation
                err = Error(n, fmas, h, y0, y, ijob, x, first, reject);

                // Computation of HNEW (we require 0.2 <= HNEW/H <= 8.0)
                fac = Math.Min(safe, cfac / (newt + (nit << 1)));
                quot = Math.Max(facr, Math.Min(facl, Math.Pow(err, 0.25) / fac));
                hnew = h / quot;

                // Is the error small enough ?
                if (err < 1.0)
                {
                    // Step is accepted
                    first = false;
                    naccpt++;

                    if (pred)
                    {
                        // Predictive controller of gustafsson
                        if (naccpt > 1)
                        {
                            facgus = hacc / h * Math.Pow(err * err / erracc, 0.25) / safe;
                            facgus = Math.Max(facr, Math.Min(facl, facgus));
                            quot = Math.Max(quot, facgus);
                            hnew = h / quot;
                        }
                        hacc = h;
                        erracc = Math.Max(.01, err);
                    }

                    xold = x;
                    hold = h;
                    x = xph;

                    for (i = 0; i < n; i++)
                    {
                        y[i] += z3[i];
                        z2i = z2[i];
                        z1i = z1[i];
                        cont[i + n] = (z2i - z3[i]) / c2m1;
                        ak = (z1i - z2i) / c1mc2;
                        acont3 = z1i / c1;
                        acont3 = (ak - acont3) / c2;
                        cont[i + n2] = (ak - cont[i + n]) / c1m1;
                        cont[i + n3] = cont[i + n2] - acont3;
                    }

                    for (i = 0; i < n; i++)
                    {
                        scal[i] = atol + rtol * Math.Abs(y[i]);
                    }

                    if (dense)
                    {
                        _xsol = x;
                        _hsol = hold;
                        for (i = 0; i < n; i++)
                        {
                            cont[i] = y[i];
                        }
                    }

                    caljac = false;

                    if (last)
                    {
                        h = hopt;
                        return 1;
                    }

                    fcn(n, x, y, y0);

                    hnew = posneg * Math.Min(Math.Abs(hnew), hmaxn);
                    hopt = hnew;
                    hopt = Math.Min(h, hnew);

                    if (reject)
                    {
                        hnew = posneg * Math.Min(Math.Abs(hnew), Math.Abs(h));
                    }

                    reject = false;

                    if ((x + hnew / quot1 - xend) * posneg >= 0.0)
                    {
                        h = xend - x;
                        last = true;
                    }
                    else
                    {
                        double qt = hnew / h;
                        hhfac = h;
                        if (theta <= thet && qt >= quot1 && qt <= quot2)
                        {
                            // No need to change step size (keep Jacobian).
                            goto L30;
                        }
                        h = hnew;
                    }

                    hhfac = h;

                    if (theta > thet)
                    {
                        ComputeJacobian(n, x, y);
                    }
                }
                else
                {
                    // Step is rejected
                    reject = true;
                    last = false;
                    if (first)
                    {
                        h *= 0.1;
                        hhfac = 0.1;
                    }
                    else
                    {
                        hhfac = hnew / h;
                        h = hnew;
                    }
                    if (naccpt >= 1)
                    {
                        ++(nrejct);
                    }
                    if (!caljac)
                    {
                        ComputeJacobian(n, x, y);
                    }
                }
            }

            return -2;
        }

        private bool Factorize(double h, int ijob)
        {
            ndec++;

            int ier = decomr_(n, fjac, fmas, u1 / h, e1, ip1, ijob);

            if (ier != 0)
            {
                ++nsing;
                if (nsing >= 5)
                {
                    throw new Exception("Matrix is repeatedly singular.");
                }

                return false;
            }

            ier = decomc_(n, fjac, fmas, alph / h, beta / h, e2r, e2i, ip2, ijob);

            if (ier != 0)
            {
                ++nsing;
                if (nsing >= 5)
                {
                    throw new Exception("Matrix is repeatedly singular.");
                }

                return false;
            }

            return true;
        }

        private int Newton(int n, int nit, double x, double[] y, ref double h, ref double hhfac, double thet,
            double uround, double fnewt, int ijob, ref double faccon, ref double theta)
        {
            double dynold = 0.0;
            double thqold = 0.0;

            double dyno, dyth, qnewt, thq, denom;
            double d1, d2, d3;
            double fac1 = u1 / h;
            double xph = x + h;

            int n3 = n * 3;
            
            double alphn = alph / h;
            double betan = beta / h;

            faccon = Math.Pow(Math.Max(faccon, uround), 0.8);
            theta = Math.Abs(thet);

            int newt = 0;

            while (newt < nit)
            {
                // Compute the right-hand side
                for (int i = 0; i < n; i++)
                {
                    cont[i] = y[i] + z1[i];
                }
                fcn(n, x + c1 * h, cont, z1);
                for (int i = 0; i < n; i++)
                {
                    cont[i] = y[i] + z2[i];
                }
                fcn(n, x + c2 * h, cont, z2);
                for (int i = 0; i < n; i++)
                {
                    cont[i] = y[i] + z3[i];
                }
                fcn(n, xph, cont, z3);

                // Solve the linear systems
                for (int i = 0; i < n; i++)
                {
                    double a1 = z1[i];
                    double a2 = z2[i];
                    double a3 = z3[i];

                    z1[i] = ti11 * a1 + ti12 * a2 + ti13 * a3;
                    z2[i] = ti21 * a1 + ti22 * a2 + ti23 * a3;
                    z3[i] = ti31 * a1 + ti32 * a2 + ti33 * a3;
                }

                Solve(n, fjac, fmas, fac1, alphn, betan, e1, e2r, e2i, z1, z2, z3, f1, f2, f3, cont, ip1, ip2, ijob);

                newt++;

                dyno = 0.0;

                for (int i = 0; i < n; i++)
                {
                    denom = scal[i];

                    d1 = z1[i] / denom;
                    d2 = z2[i] / denom;
                    d3 = z3[i] / denom;

                    dyno = dyno + d1 * d1 + d2 * d2 + d3 * d3;
                }

                dyno = Math.Sqrt(dyno / n3);

                // Bad convergence or number of iterations to large
                if (newt > 1 && newt < nit)
                {
                    thq = dyno / dynold;
                    if (newt == 2)
                    {
                        theta = thq;
                    }
                    else
                    {
                        theta = Math.Sqrt(thq * thqold);
                    }
                    thqold = thq;
                    if (theta < 0.99)
                    {
                        faccon = theta / (1.0 - theta);
                        dyth = faccon * dyno * Math.Pow(theta, nit - 1 - newt) / fnewt;
                        if (dyth >= 1.0)
                        {
                            qnewt = Math.Max(1e-4, Math.Min(20.0, dyth));
                            hhfac = Math.Pow(qnewt, -1.0 / (nit + 4.0 - 1 - newt)) * 0.8;
                            h = hhfac * h;

                            return -1;
                        }
                    }
                    else
                    {
                        // Unexpected step-rejection
                        h *= 0.5;
                        hhfac = 0.5;

                        return -1;
                    }
                }

                dynold = Math.Max(dyno, uround);

                for (int i = 0; i < n; i++)
                {
                    double f1i = f1[i] + z1[i];
                    double f2i = f2[i] + z2[i];
                    double f3i = f3[i] + z3[i];

                    f1[i] = f1i;
                    f2[i] = f2i;
                    f3[i] = f3i;

                    z1[i] = t11 * f1i + t12 * f2i + t13 * f3i;
                    z2[i] = t21 * f1i + t22 * f2i + t23 * f3i;
                    z3[i] = t31 * f1i + f2i;
                }

                if (faccon * dyno <= fnewt)
                {
                    return newt;
                }
            }

            h *= 0.5;
            hhfac = 0.5;

            return -1;
        }

        private void ComputeJacobian(int n, double x, double[] y)
        {
            // Computation of the Jacobian
            if (jac == null)
            {
                // Compute Jacobian matrix numerically
                double uround = Precision.DoublePrecision;

                // Jacobian is full
                for (int i = 0; i < n; i++)
                {
                    var ysafe = y[i];
                    var delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                    y[i] = ysafe + delt;
                    fcn(n, x, y, cont);

                    for (int j = 0; j < n; j++)
                    {
                        fjac[j - i * n] = (cont[j] - y0[j]) / delt;
                    }
                    y[i] = ysafe;
                }
            }
            else
            {
                // Compute Jacobian matrix analytically
                jac(n, x, y, fjac);
            }

            caljac = true;
        }

        // This function can be used for coninuous output. It provides an
        // approximation to the i-th component of the solution at X.
        // It gives the value of the collocation polynomial, defined for
        // the last successfully computed step (by RADAU5).
        public double Interpolate(int i, double x)
        {
            double s = (x - _xsol) / _hsol;

            return cont[i] + s * (cont[i + n] + (s - c2m1) * (cont[i + n * 2] + (s - c1m1) * cont[i + n * 3]));
        }

        int decomr_(int n, double[] fjac, double[] fmas, double fac1,
            double[] e, int[] ip, int ijob)
        {
            bool dae = ijob == 5;

            int ier = 0;

            if (!dae)
            {
                // B=identity, Jacobian a full matrix
                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        e[i + j * n] = -fjac[i + j * n];
                    }
                    e[j + j * n] += fac1;
                }
                decsol.dec_(n, e, ip, ref ier);
            }
            else
            {
                // B is a full matrix, Jacobian a full matrix
                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        e[i + j * n] = fmas[i + j * n] * fac1 - fjac[i + j * n];
                    }
                }
                decsol.dec_(n, e, ip, ref ier);
            }

            return ier;
        }

        int decomc_(int n, double[] fjac, double[] fmas, double alphn, double betan,
            double[] e2r, double[] e2i, int[] ip2, int ijob)
        {
            bool dae = ijob == 5;

            int ier = 0;

            if (!dae)
            {
                // B=identity, Jacobian a full matrix
                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        e2r[i + j * n] = -fjac[i + j * n];
                        e2i[i + j * n] = 0.0;
                    }
                    e2r[j + j * n] += alphn;
                    e2i[j + j * n] = betan;
                }
                decsol.decc_(n, e2r, e2i, ip2, ref ier);
            }
            else
            {
                // B is a full matrix, Jacobian a full matrix
                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        double bb = fmas[i + j * n];
                        e2r[i + j * n] = bb * alphn - fjac[i + j * n];
                        e2i[i + j * n] = bb * betan;
                    }
                }
                decsol.decc_(n, e2r, e2i, ip2, ref ier);
            }

            return ier;
        }

        int Solve(int n, double[] fjac, double[] fmas,
            double fac1, double alphn, double betan,
            double[] e1, double[] e2r, double[] e2i,
            double[] z1, double[] z2, double[] z3, double[] f1,
            double[] f2, double[] f3, double[] cont, int[] ip1,
            int[] ip2, int ijob)
        {
            nsol++;
            
            bool dae = ijob == 5;

            if (!dae)
            {
                // B=identity, Jacobian a full matrix
                for (int i = 0; i < n; i++)
                {
                    double s2 = -f2[i];
                    double s3 = -f3[i];

                    z1[i] -= f1[i] * fac1;
                    z2[i] += s2 * alphn - s3 * betan;
                    z3[i] += s3 * alphn + s2 * betan;
                }

                decsol.sol_(n, e1, z1, ip1);
                decsol.solc_(n, e2r, e2i, z2, z3, ip2);
            }
            else
            {
                // B is a full matrix, Jacobian a full matrix
                for (int i = 0; i < n; i++)
                {
                    double s1 = 0.0;
                    double s2 = 0.0;
                    double s3 = 0.0;

                    for (int j = 0; j < n; j++)
                    {
                        double bb = fmas[i + j * n];

                        s1 -= bb * f1[j];
                        s2 -= bb * f2[j];
                        s3 -= bb * f3[j];
                    }

                    z1[i] += s1 * fac1;
                    z2[i] += s2 * alphn - s3 * betan;
                    z3[i] += s3 * alphn + s2 * betan;
                }
                decsol.sol_(n, e1, z1, ip1);
                decsol.solc_(n, e2r, e2i, z2, z3, ip2);
            }

            return 0;
        }

        double Error(int n, double[] fmas, double h, double[] y0,
            double[] y, int ijob, double x, bool first, bool reject)
        {
            double d1;
            
            double hee1 = dd1 / h;
            double hee2 = dd2 / h;
            double hee3 = dd3 / h;

            bool dae = ijob == 5;

            if (!dae)
            {
                // B=identity, Jacobian a full matrix
                for (int i = 0; i < n; i++)
                {
                    f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
                    cont[i] = f2[i] + y0[i];
                }
                decsol.sol_(n, e1, cont, ip1);
            }
            else
            {
                // B is a full matrix, Jacobian a full matrix
                for (int i = 0; i < n; i++)
                {
                    f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
                }
                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < n; j++)
                    {
                        sum += fmas[i + j * n] * f1[j];
                    }
                    f2[i] = sum;
                    cont[i] = sum + y0[i];
                }
                decsol.sol_(n, e1, cont, ip1);
            }

            double err = 0.0;

            for (int i = 0; i < n; i++)
            {
                d1 = cont[i] / scal[i];
                err += d1 * d1;
            }

            err = Math.Max(Math.Sqrt(err / n), 1e-10);

            if (err < 1.0)
            {
                return err;
            }

            if (first || reject)
            {
                for (int i = 0; i < n; i++)
                {
                    cont[i] = y[i] + cont[i];
                }
                fcn(n, x, cont, f1);

                for (int i = 0; i < n; i++)
                {
                    cont[i] = f1[i] + f2[i];
                }

                // Full matrix option
                decsol.sol_(n, e1, cont, ip1);

                err = 0.0;

                for (int i = 0; i < n; i++)
                {
                    d1 = cont[i] / scal[i];
                    err += d1 * d1;
                }

                err = Math.Max(Math.Sqrt(err / n), 1e-10);
            }

            return err;
        }
    }
}