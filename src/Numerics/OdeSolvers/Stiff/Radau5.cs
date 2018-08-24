// Based on Fortran code RADAU5
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;

    using S_fp = System.Action<int, double, double[], double[]>;
    using J_fp = System.Action<int, double, double[], double[], int>;

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

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        struct conra5
        {
            public int nn, nn2, nn3, nn4;
            public double xsol, hsol, c2m1, c1m1;
        }

        conra5 conra5_1;

        public int nstep, ndec, nsol, nrejct, naccpt;

        // Subroutine
        public int radau5_(int n, S_fp fcn, double x, double[] y,
            double xend, double h, double rtol, double atol,
            J_fp jac, int ijac, DenseMatrix mas, int iout)
        {
            // Local variables
            int nit, lde1;
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
            int ldjac;
            int ldmas;
            double fnewt;
            double tolst;
            int ldmas2;

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

            // Computation of the row-dimensions of the 2-arrays
            // Jacobian  and  matrices E1, E2
            {
                ldjac = n;
                lde1 = n;
            }
            // Mass matrix
            if (implct)
            {
                ldmas = n;
                ijob = 5;
            }
            else
            {
                ldmas = 0;

                ijob = 1;
            }
            ldmas2 = Math.Max(1, ldmas);

            // Prepare the entry-points for the arrays in work
            var z1 = new double[n];
            var z2 = new double[n];
            var z3 = new double[n];
            var y0 = new double[n];
            var scal = new double[n];
            var f1 = new double[n];
            var f2 = new double[n];
            var f3 = new double[n];
            var con = new double[n << 2];
            var _jac = new double[n * ldjac];
            var _mas = new double[n * ldmas];
            var e1 = new double[n * lde1];
            var e2r = new double[n * lde1];
            var e2i = new double[n * lde1];

            // Entry points for integer workspace
            var ip1 = new int[n];
            var ip2 = new int[n];
            //var iph = new int[n];

            // Call to core integrator
            int idid = radcor_(n, fcn, x, y, xend, hmax, h, rtol, atol,
                jac, ijac, mas, iout, nmax, uround, safe, thet, fnewt,
                quot1, quot2, nit, ijob, startn, nind1, nind2, nind3,
                pred, facl, facr, implct, ldjac,
                lde1, ldmas2, z1, z2, z3, y0,
                 scal, f1, f2, f3,
                _jac, e1, e2r, e2i, _mas,
                ip1, ip2, con);

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

        int radcor_(int n, S_fp fcn, double x, double[]
            y, double xend, double hmax, double h, double rtol, double atol,
            J_fp jac, int ijac, DenseMatrix mas, int iout, int nmax,
            double uround, double safe, double thet, double
            fnewt, double quot1, double quot2, int nit, int
            ijob, bool startn, int nind1, int nind2, int nind3,
            bool pred, double facl, double facr, bool implct, int
            ldjac, int lde1, int ldmas, double[] z1, double[] z2,
            double[] z3, double[] y0, double[] scal, double[] f1,
            double[] f2, double[] f3, double[] fjac, double[] e1,
            double[] e2r, double[] e2i, double[] fmas, int[] ip1,
            int[] ip2, double[] cont)
        {
            double d1, d2, d3;

            int i, j;
            double a1, a2, c1, c2, a3;
            int n2, n3;
            double u1, ak;
            double t11, t12, t13, t21, t22, t23, t31;
            double qt, dd1, dd2, dd3, ak1, ak2, ak3, f1i, f2i, f3i, c1q, c2q, c3q,
                 z1i, z2i, z3i, sq6, fac, ti11, cno;
            int lrc;
            double ti12, ti13, ti21, ti22, ti23, ti31, ti32, ti33;
            int ier = 0;
            double xph, thq, err = 0, fac1, cfac, hacc = 0, c1mc2, beta;
            double alph, hold;
            double delt, hnew;
            bool last;
            double hopt, xold;
            int newt;
            double dyno, dyth, quot, hhfac, betan, alphn, denom, theta, ysafe,
                hmaxn;
            int nsing;
            bool first;
            int irtrn = 0, nrsol, nsolu;
            double qnewt, xosol, acont3;
            bool index1, index2, index3, caljac;
            double faccon;
            double erracc = 0;
            bool reject;
            double facgus;
            double dynold = 0, posneg;
            double thqold = 0;

            // Core integrator for RADAU5
            // Parameters same as in radau5 with workspace added

            // Initialisations

            conra5_1.nn = n;
            conra5_1.nn2 = n << 1;
            conra5_1.nn3 = n * 3;
            lrc = n << 2;

            // Check the index of the problem
            index1 = nind1 != 0;
            index2 = nind2 != 0;
            index3 = nind3 != 0;

            // Compute mass matrix for implicit case
            if (implct)
            {
                fmas = mas.ToColumnMajorArray();
            }
            // Constants
            sq6 = Math.Sqrt(6.0);
            c1 = (4.0 - sq6) / 10.0;
            c2 = (sq6 + 4.0) / 10.0;
            conra5_1.c1m1 = c1 - 1.0;
            conra5_1.c2m1 = c2 - 1.0;
            c1mc2 = c1 - c2;
            dd1 = -(sq6 * 7.0 + 13.0) / 3.0;
            dd2 = (sq6 * 7.0 - 13.0) / 3.0;
            dd3 = -0.33333333333333331;
            u1 = (Math.Pow(81.0, 0.33333333333333331) + 6.0 - Math.Pow(9.0, 0.33333333333333331)) / 30.0;
            alph = (12.0 - Math.Pow(81.0, 0.33333333333333331) + Math.Pow(9.0, 0.33333333333333331)) / 60.0;
            beta = (Math.Pow(81.0, 0.33333333333333331) + Math.Pow(9.0, 0.33333333333333331)) * Math.Sqrt(3.0) / 60.0;
            // Computing 2nd power
            d1 = alph;
            d2 = beta;
            cno = d1 * d1 + d2 * d2;
            u1 = 1.0 / u1;
            alph /= cno;
            beta /= cno;
            t11 = 0.091232394870892942792;
            t12 = -0.14125529502095420843;
            t13 = -0.030029194105147424492;
            t21 = 0.24171793270710701896;
            t22 = 0.20412935229379993199;
            t23 = 0.38294211275726193779;
            t31 = 0.96604818261509293619;
            ti11 = 4.325579890063155351;
            ti12 = 0.33919925181580986954;
            ti13 = 0.54177053993587487119;
            ti21 = -4.1787185915519047273;
            ti22 = -0.32768282076106238708;
            ti23 = 0.47662355450055045196;
            ti31 = -0.50287263494578687595;
            ti32 = 2.5719269498556054292;
            ti33 = -0.59603920482822492497;

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
            cfac = safe * ((nit << 1) + 1);
            nsing = 0;
            xold = x;
            if (iout != 0)
            {
                irtrn = 1;
                nrsol = 1;
                xosol = xold;
                conra5_1.xsol = x;
                for (i = 0; i < n; i++)
                {
                    cont[i] = y[i];
                }
                nsolu = n;
                conra5_1.hsol = hold;
                //solout(&nrsol, &xosol, &conra5_1.xsol, y, cont, &lrc, &nsolu, &irtrn);
                if (irtrn < 0)
                {
                    return 2; // Exit caused by solout
                }
            }

            n2 = n << 1;
            n3 = n * 3;
            
            for (i = 0; i < n; i++)
            {
                scal[i] = atol + rtol * Math.Abs(y[i]);
            }

            hhfac = h;
            fcn(n, x, y, y0);

            // Basic integration step
            L10:
            // Computation of the Jacobian
            if (ijac == 0)
            {
                // Compute Jacobian matrix numerically

                // Jacobian is full
                for (i = 0; i < n; i++)
                {
                    ysafe = y[i];
                    delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                    y[i] = ysafe + delt;
                    fcn(n, x, y, cont);

                    for (j = 0; j < n; j++)
                    {
                        fjac[j - i * ldjac] = (cont[j] - y0[j]) / delt;
                    }
                    y[i] = ysafe;
                }
            }
            else
            {
                // Compute Jacobian matrix analytically
                jac(n, x, y, fjac, ldjac);
            }
            caljac = true;

            L20:
            // Compute the matrices E1 and E2 and their decompositions
            fac1 = u1 / h;
            alphn = alph / h;
            betan = beta / h;
            decomr_(n, fjac, ldjac, fmas, ldmas, fac1, e1, lde1, ip1, ref ier, ijob);
            if (ier != 0)
            {
                goto L78;
            }
            decomc_(n, fjac, ldjac, fmas, ldmas, alphn, betan, e2r, e2i, lde1, ip2, ref ier, ijob);
            if (ier != 0)
            {
                goto L78;
            }
            ++(ndec);
            L30:
            ++(nstep);
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
                    z1i = c1q * (ak1 + (c1q - conra5_1.c2m1) * (ak2 + (c1q - conra5_1.c1m1) * ak3));
                    z2i = c2q * (ak1 + (c2q - conra5_1.c2m1) * (ak2 + (c2q - conra5_1.c1m1) * ak3));
                    z3i = c3q * (ak1 + (c3q - conra5_1.c2m1) * (ak2 + (c3q - conra5_1.c1m1) * ak3));
                    z1[i] = z1i;
                    z2[i] = z2i;
                    z3[i] = z3i;
                    f1[i] = ti11 * z1i + ti12 * z2i + ti13 * z3i;
                    f2[i] = ti21 * z1i + ti22 * z2i + ti23 * z3i;
                    f3[i] = ti31 * z1i + ti32 * z2i + ti33 * z3i;
                }
            }

            // Loop for the simplified newton iteration
            newt = 0;
            faccon = Math.Pow(Math.Max(faccon, uround), 0.8);
            theta = Math.Abs(thet);
            L40:
            if (newt >= nit)
            {
                goto L78;
            }
            // Compute the right-hand side
            for (i = 0; i < n; i++)
            {
                cont[i] = y[i] + z1[i];
            }
            fcn(n, x + c1 * h, cont, z1);
            for (i = 0; i < n; i++)
            {
                cont[i] = y[i] + z2[i];
            }
            fcn(n, x + c2 * h, cont, z2);
            for (i = 0; i < n; i++)
            {
                cont[i] = y[i] + z3[i];
            }
            fcn(n, xph, cont, z3);

            // Solve the linear systems
            for (i = 0; i < n; i++)
            {
                a1 = z1[i];
                a2 = z2[i];
                a3 = z3[i];
                z1[i] = ti11 * a1 + ti12 * a2 + ti13 * a3;
                z2[i] = ti21 * a1 + ti22 * a2 + ti23 * a3;
                z3[i] = ti31 * a1 + ti32 * a2 + ti33 * a3;
            }
            slvrad_(n, fjac, ldjac, fmas, ldmas, fac1, alphn, betan, e1,
                 e2r, e2i, lde1, z1, z2, z3, f1, f2, f3, cont, ip1, ip2, ref ier, ijob);
            ++(nsol);
            ++newt;
            dyno = 0.0;
            for (i = 0; i < n; i++)
            {
                denom = scal[i];
                // Computing 2nd power
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
                    dyth = faccon * dyno * Math.Pow(theta, nit - 1 - newt) / fnewt; // TODO: pow_di
                    if (dyth >= 1.0)
                    {
                        qnewt = Math.Max(1e-4, Math.Min(20.0, dyth));
                        hhfac = Math.Pow(qnewt, -1.0 / (nit + 4.0 - 1 - newt)) * 0.8;
                        h = hhfac * h;
                        reject = true;
                        last = false;
                        if (caljac)
                        {
                            goto L20;
                        }
                        goto L10;
                    }
                }
                else
                {
                    goto L78;
                }
            }
            dynold = Math.Max(dyno, uround);
            for (i = 0; i < n; i++)
            {
                f1i = f1[i] + z1[i];
                f2i = f2[i] + z2[i];
                f3i = f3[i] + z3[i];
                f1[i] = f1i;
                f2[i] = f2i;
                f3[i] = f3i;
                z1[i] = t11 * f1i + t12 * f2i + t13 * f3i;
                z2[i] = t21 * f1i + t22 * f2i + t23 * f3i;
                z3[i] = t31 * f1i + f2i;
            }
            if (faccon * dyno > fnewt)
            {
                goto L40;
            }
            // Error estimation
            estrad_(n, fjac, ldjac, h, dd1, dd2, dd3, fcn, y0,
                y, ijob, x, e1, lde1, z1, z2, z3,
                cont, f1, f2, ip1, scal, ref err, first, reject);

            // Computation of HNEW
            // We require .2<=HNEW/H<=8.
            fac = Math.Min(safe, cfac / (newt + (nit << 1)));
            quot = Math.Max(facr, Math.Min(facl, Math.Pow(err, 0.25) / fac));
            hnew = h / quot;

            // Is the error small enough ?
            if (err < 1.0)
            {
                // Step is accepted
                first = false;
                ++(naccpt);
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
                    cont[i + n] = (z2i - z3[i]) / conra5_1.c2m1;
                    ak = (z1i - z2i) / c1mc2;
                    acont3 = z1i / c1;
                    acont3 = (ak - acont3) / c2;
                    cont[i + n2] = (ak - cont[i + n]) / conra5_1.c1m1;
                    cont[i + n3] = cont[i + n2] - acont3;
                }
                
                for (i = 0; i < n; i++)
                {
                    scal[i] = atol + rtol * Math.Abs(y[i]);
                }

                if (iout != 0)
                {
                    nrsol = naccpt + 1;
                    conra5_1.xsol = x;
                    xosol = xold;
                    for (i = 0; i < n; i++)
                    {
                        cont[i] = y[i];
                    }
                    nsolu = n;
                    conra5_1.hsol = hold;
                    //solout(&nrsol, &xosol, &conra5_1.xsol, y, cont, &lrc, &nsolu, &irtrn);
                    if (irtrn < 0)
                    {
                        return 2; // Exit caused by solout
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
                    qt = hnew / h;
                    hhfac = h;
                    if (theta <= thet && qt >= quot1 && qt <= quot2)
                    {
                        goto L30;
                    }
                    h = hnew;
                }
                hhfac = h;
                if (theta <= thet)
                {
                    goto L20;
                }
                goto L10;
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
                if (caljac)
                {
                    goto L20;
                }
                goto L10;
            }
            // Unexpected step-rejection
            L78:
            if (ier != 0)
            {
                ++nsing;
                if (nsing >= 5)
                {
                    throw new Exception(" MATRIX IS REPEATEDLY SINGULAR, IER=" + ier);
                }
            }
            h *= 0.5;
            hhfac = 0.5;
            reject = true;
            last = false;
            if (caljac)
            {
                goto L20;
            }
            goto L10;
        }

        // This function can be used for coninuous output. it provides an
        // approximation to the i-th component of the solution at X.
        // It gives the value of the collocation polynomial, defined for
        // the last successfully computed step (by RADAU5).
        double contr5_(int i, double x, double[] cont, int lrc)
        {
            double s = (x - conra5_1.xsol) / conra5_1.hsol;

            return cont[i] + s * (cont[i + conra5_1.nn] + (s - conra5_1.c2m1)
                * (cont[i + conra5_1.nn2] + (s - conra5_1.c1m1) * cont[i + conra5_1.nn3]));
        }
        
        int decomr_(int n, double[] fjac, int ldjac,
            double[] fmas, int ldmas, double fac1, double[] e,
            int lde, int[] ip, ref int ier, int ijob)
        {
            switch (ijob)
            {
                case 1: goto L1;
                case 5: goto L5;
            }

            L1:
            /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    e[i + j * lde] = -fjac[i + j * ldjac];
                }
                e[j + j * lde] += fac1;
            }
            decsol.dec_(n, lde, e, ip, ref ier);
            return 0;

            L5:
            /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    e[i + j * lde] = fmas[i + j * ldmas] * fac1 - fjac[i + j * ldjac];
                }
            }
            decsol.dec_(n, lde, e, ip, ref ier);
            return 0;
        }

        int decomc_(int n, double[] fjac, int ldjac,
            double[] fmas, int ldmas, double alphn, double betan,
            double[] e2r, double[] e2i, int lde1, int[] ip2, ref int ier, int ijob)
        {
            switch (ijob)
            {
                case 1: goto L1;
                case 5: goto L5;
            }

            L1:
            /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    e2r[i + j * lde1] = -fjac[i + j * ldjac];
                    e2i[i + j * lde1] = 0.0;
                }
                e2r[j + j * lde1] += alphn;
                e2i[j + j * lde1] = betan;
            }
            decsol.decc_(n, lde1, e2r, e2i, ip2, ref ier);
            return 0;

            L5:
            /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    double bb = fmas[i + j * ldmas];
                    e2r[i + j * lde1] = bb * alphn - fjac[i + j * ldjac];
                    e2i[i + j * lde1] = bb * betan;
                }
            }
            decsol.decc_(n, lde1, e2r, e2i, ip2, ref ier);
            return 0;
        }

        int slvrad_(int n, double[] fjac, int ldjac,
            double[] fmas, int ldmas,
            double fac1, double alphn, double betan,
            double[] e1, double[] e2r, double[] e2i, int lde1,
            double[] z1, double[] z2, double[] z3, double[] f1,
            double[] f2, double[] f3, double[] cont, int[] ip1,
            int[] ip2, ref int ier, int ijob)
        {
            /* Local variables */
            int i, j;
            double s1, s2, s3, bb;

            switch (ijob)
            {
                case 1: goto L1;
                case 5: goto L5;
            }

            /* ----------------------------------------------------------- */

            L1:
            /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
            for (i = 0; i < n; i++)
            {
                s2 = -f2[i];
                s3 = -f3[i];
                z1[i] -= f1[i] * fac1;
                z2[i] = z2[i] + s2 * alphn - s3 * betan;
                z3[i] = z3[i] + s3 * alphn + s2 * betan;
            }
            decsol.sol_(n, lde1, e1, z1, ip1);
            decsol.solc_(n, lde1, e2r, e2i, z2, z3, ip2);
            return 0;

            /* ----------------------------------------------------------- */

            L5:
            /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
            for (i = 0; i < n; i++)
            {
                s1 = 0.0;
                s2 = 0.0;
                s3 = 0.0;
                for (j = 0; j < n; j++)
                {
                    bb = fmas[i + j * ldmas];
                    s1 -= bb * f1[j];
                    s2 -= bb * f2[j];
                    s3 -= bb * f3[j];
                }
                z1[i] += s1 * fac1;
                z2[i] = z2[i] + s2 * alphn - s3 * betan;
                z3[i] = z3[i] + s3 * alphn + s2 * betan;
            }
            decsol.sol_(n, lde1, e1, z1, ip1);
            decsol.solc_(n, lde1, e2r, e2i, z2, z3, ip2);
            return 0;
        }

        int estrad_(int n,
            double[] fmas, int ldmas, double h, double dd1,
            double dd2, double dd3, S_fp fcn, double[] y0,
            double[] y, int ijob, double x, double[] e1, int lde1, double[] z1,
            double[] z2, double[] z3, double[] cont, double[] f1,
            double[] f2, int[] ip1, double[] scal,
            ref double err, bool first, bool reject)
        {
            double d1;

            double sum, hee1, hee2, hee3;

            hee1 = dd1 / h;
            hee2 = dd2 / h;
            hee3 = dd3 / h;

            switch (ijob)
            {
                case 1: goto L1;
                case 5: goto L5;
            }

            L1:
            /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
            for (int i = 0; i < n; i++)
            {
                f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
                cont[i] = f2[i] + y0[i];
            }
            decsol.sol_(n, lde1, e1, cont, ip1);
            goto L77;

            L5:
            /* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
            for (int i = 0; i < n; i++)
            {
                f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
            }
            for (int i = 0; i < n; i++)
            {
                sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    sum += fmas[i + j * ldmas] * f1[j];
                }
                f2[i] = sum;
                cont[i] = sum + y0[i];
            }
            decsol.sol_(n, lde1, e1, cont, ip1);
            goto L77;

            L77:
            err = 0.0;
            for (int i = 0; i < n; i++)
            {
                /* Computing 2nd power */
                d1 = cont[i] / scal[i];
                err += d1 * d1;
            }
            /* Computing MAX */
            d1 = Math.Sqrt(err / n);
            err = Math.Max(d1, 1e-10);

            if (err < 1.0)
            {
                return 0;
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
                switch (ijob)
                {
                    case 1: goto L31;
                    case 5: goto L31;
                }
                /* ------ FULL MATRIX OPTION */
                L31:
                decsol.sol_(n, lde1, e1, cont, ip1);
                goto L88;

                L88:
                err = 0.0;

                for (int i = 0; i < n; i++)
                {
                    /* Computing 2nd power */
                    d1 = cont[i] / scal[i];
                    err += d1 * d1;
                }
                /* Computing MAX */
                d1 = Math.Sqrt(err / n);
                err = Math.Max(d1, 1e-10);
            }
            return 0;
        }
    }
}