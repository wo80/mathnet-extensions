
namespace MathNet.Numerics.OdeSolvers
{
    using System;
    using System.Diagnostics;
    
    public class DormandPrince5
    {
        #region Runge-Kutta coefficients

        private const double c2 = 0.2;
        private const double c3 = 0.3;
        private const double c4 = 0.8;
        private const double c5 = 0.88888888888888884;

        private const double a21 = 0.2;
        private const double a31 = 0.074999999999999997;
        private const double a32 = 0.22500000000000001;
        private const double a41 = 0.97777777777777775;
        private const double a42 = -3.7333333333333334;
        private const double a43 = 3.5555555555555554;
        private const double a51 = 2.9525986892242035;
        private const double a52 = -11.595793324188385;
        private const double a53 = 9.8228928516994358;
        private const double a54 = -0.29080932784636487;
        private const double a61 = 2.8462752525252526;
        private const double a62 = -10.757575757575758;
        private const double a63 = 8.9064227177434727;
        private const double a64 = 0.27840909090909088;
        private const double a65 = -0.2735313036020583;
        private const double a71 = 0.091145833333333329;
        private const double a73 = 0.44923629829290207;
        private const double a74 = 0.65104166666666663;
        private const double a75 = -0.322376179245283;
        private const double a76 = 0.13095238095238096;

        private const double e1 = 0.0012326388888888888;
        private const double e3 = -0.0042527702905061394;
        private const double e4 = 0.036979166666666667;
        private const double e5 = -0.05086379716981132;
        private const double e6 = 0.041904761904761903;
        private const double e7 = -0.025000000000000001;

        // Dense output of Shampine (1986)
        private const double d1 = -1.1270175653862835;
        private const double d3 = 2.675424484351598;
        private const double d4 = -5.6855269615885042;
        private const double d5 = 3.5219323679207912;
        private const double d6 = -1.7672812570757455;
        private const double d7 = 2.3824689317781438;

        #endregion

        int n;

        double[] dxdt, xtemp, xerr;
        double[] k2, k3, k4, k5, k6;
        double[] r;

        double[] rtol, atol;

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        class condo_
        {
            public double xold, hout;
        }

        condo_ condo_1 = new condo_();

        /* ----------------------------------------------------------
         *     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
         *     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
         *     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER (4)5
         *     DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND
         *     DENSE OUTPUT).
         *
         *     AUTHORS: E. HAIRER AND G. WANNER
         *              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
         *              CH-1211 GENEVE 24, SWITZERLAND
         *              E-MAIL:  Ernst.Hairer@math.unige.ch
         *                       Gerhard.Wanner@math.unige.ch
         *
         *     THIS CODE IS DESCRIBED IN:
         *         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
         *         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
         *         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
         *         SPRINGER-VERLAG (1993)
         *
         *     VERSION OF APRIL 25, 1996
         *     (latest correction of a small bug: August 8, 2005)
         *
         *     INPUT PARAMETERS
         *     ----------------
         *     N           DIMENSION OF THE SYSTEM
         *
         *     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
         *                 VALUE OF F(X,Y):
         *                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
         *                    DOUBLE PRECISION X,Y(N),F(N)
         *                    F(1)=...   ETC.
         *
         *     X           INITIAL X-VALUE
         *
         *     Y(N)        INITIAL VALUES FOR Y
         *
         *     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
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
         *     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
         *                 NUMERICAL SOLUTION DURING INTEGRATION.
         *                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
         *                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
         *                 IT MUST HAVE THE FORM
         *                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
         *                                       RPAR,IPAR,IRTRN)
         *                    DIMENSION Y(N),CON(5*ND),ICOMP(ND)
         *                    ....
         *                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
         *                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
         *                    THE FIRST GRID-POINT).
         *                 "XOLD" IS THE PRECEEDING GRID-POINT.
         *                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
         *                    IS SET <0, DOPRI5 WILL RETURN TO THE CALLING PROGRAM.
         *                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
         *                    SET  IRTRN = 2
         *
         *          -----  CONTINUOUS OUTPUT: -----
         *                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
         *                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
         *                 THE FUNCTION
         *                        >>>   CONTD5(I,S,CON,ICOMP,ND)   <<<
         *                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
         *                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
         *                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
         *
         *     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
         *                    IOUT=0: SUBROUTINE IS NEVER CALLED
         *                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT.
         *                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
         *                            (IN THIS CASE WORK(5) MUST BE SPECIFIED)
         *
         *     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
         *                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE.
         *                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
         *                 "LWORK" MUST BE AT LEAST  8*N+5*NRDENS+21
         *                 WHERE  NRDENS = IWORK(5)
         *
         *     LWORK       DECLARED LENGHT OF ARRAY "WORK".
         *
         *     IWORK       INTEGER WORKING SPACE OF LENGHT "LIWORK".
         *                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS FOR THE CODE.
         *                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
         *                 "LIWORK" MUST BE AT LEAST NRDENS+21 .
         *
         *     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
         *
         *     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
         *                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
         *                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.
         *
         * -----------------------------------------------------------------------
         *
         *     SOPHISTICATED SETTING OF PARAMETERS
         *     -----------------------------------
         *              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...0) ALLOW
         *              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
         *              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES.
         *
         *    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.
         *
         *    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
         *              DEFAULT 0.9D0.
         *
         *    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
         *              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
         *                 WORK(3) <= HNEW/HOLD <= WORK(4)
         *              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=10.D0
         *
         *    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL
         *              (SEE SECTION IV.2). LARGER VALUES OF BETA ( <= 0.1 )
         *              MAKE THE STEP SIZE CONTROL MORE STABLE. DOPRI5 NEEDS
         *              A LARGER BETA THAN HIGHAM & HALL. NEGATIVE WORK(5)
         *              PROVOKE BETA=0.
         *              DEFAULT 0.04D0.
         *
         *    WORK(6)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
         *
         *    WORK(7)   INITIAL STEP SIZE, FOR WORK(7)=0.D0 AN INITIAL GUESS
         *              IS COMPUTED WITH HELP OF THE FUNCTION HINIT
         *
         *    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
         *              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
         *
         *    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
         *              IF IWORK(2).EQ.1  METHOD DOPRI5 OF DORMAND AND PRINCE
         *              (TABLE 5.2 OF SECTION II.5).
         *              AT THE MOMENT THIS IS THE ONLY POSSIBLE CHOICE.
         *              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1.
         *
         *    IWORK(3)  SWITCH FOR PRINTING ERROR MESSAGES
         *              IF IWORK(3).LT.0 NO MESSAGES ARE BEING PRINTED
         *              IF IWORK(3).GT.0 MESSAGES ARE PRINTED WITH
         *              WRITE (IWORK(3),*) ...
         *              DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=6
         *
         *    IWORK(4)  TEST FOR STIFFNESS IS ACTIVATED AFTER STEP NUMBER
         *              J*IWORK(4) (J INTEGER), PROVIDED IWORK(4).GT.0.
         *              FOR NEGATIVE IWORK(4) THE STIFFNESS TEST IS
         *              NEVER ACTIVATED; DEFAULT VALUE IS IWORK(4)=1000
         *
         *    IWORK(5)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
         *              IS REQUIRED; DEFAULT VALUE IS IWORK(5)=0;
         *              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE
         *              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN
         *              IWORK(21),...,IWORK(NRDENS+20);
         *              FOR  NRDENS=N  THIS IS DONE BY THE CODE.
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
         *                   IDID=-3  STEP SIZE BECOMES TOO SMALL.
         *                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED).
         *
         *   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
         *   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
         *   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
         *   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
         *                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
         * ----------------------------------------------------------------------- */
        public int dopri5_(Action<double, double[], double[]> fcn, double x, double[] y, double xend, double[] rtol, double[] atol,
                int itol, double[] work, int[] iwork)
        {
            /* Local variables */
            double h;
            double fac1, fac2;
            double beta, safe;
            double hmax;
            int nmax;
            
            int nstiff, nrdens, iprint;
            double uround;

            /* Function Body */
            int nfcn = 0;
            int nstep = 0;
            int naccpt = 0;
            int nrejct = 0;

            bool arret = false;

            this.n = y.Length;

            // IPRINT for monitoring the printing
            if (iwork[2] == 0)
            {
                iprint = 6;
            }
            else
            {
                iprint = iwork[2];
            }

            // NMAX the maximal number of steps
            if (iwork[0] == 0)
            {
                nmax = 100000;
            }
            else
            {
                nmax = iwork[0];
                if (nmax <= 0)
                {
                    if (iprint > 0)
                    {
                        Console.WriteLine(" Wrong input IWORK(1)=" + iwork[0]);
                    }
                    arret = true;
                }
            }

            // NSTIFF parameter for stiffness detection
            nstiff = iwork[3];
            if (nstiff == 0)
            {
                nstiff = 1000;
            }
            if (nstiff < 0)
            {
                nstiff = nmax + 10;
            }

            // NRDENS number of dense output components
            nrdens = iwork[4];
            if (nrdens < 0 || nrdens > n)
            {
                if (iprint > 0)
                {
                    Console.WriteLine(" Curious input IWORK(5)=" + iwork[4]);
                }
                arret = true;
            }
            else
            {
                if (nrdens == n)
                {
                    for (int i = 0; i < nrdens; ++i)
                    {
                        iwork[i + 20] = i;
                    }
                }
            }

            // UROUND smallest number satisfying 1.0+UROUND>1.0
            if (work[0] == 0.0)
            {
                uround = 2.3e-16;
            }
            else
            {
                uround = work[0];
                if (uround <= 1e-35 || uround >= 1.0)
                {
                    if (iprint > 0)
                    {
                        Console.WriteLine(" Which machine do you have? Your UROUND was:" + work[0]);
                    }
                    arret = true;
                }
            }

            // Safety factor
            if (work[1] == 0.0)
            {
                safe = 0.9;
            }
            else
            {
                safe = work[1];
                if (safe >= 1.0 || safe <= 1e-4)
                {
                    if (iprint > 0)
                    {
                        Console.WriteLine(" Curious input for safety factor WORK(2)=" + work[1]);
                    }
                    arret = true;
                }
            }

            // FAC1,FAC2 parameters for step size selection
            if (work[2] == 0.0)
            {
                fac1 = 0.2;
            }
            else
            {
                fac1 = work[2];
            }
            if (work[3] == 0.0)
            {
                fac2 = 10.0;
            }
            else
            {
                fac2 = work[3];
            }

            // BETA for step control stabilization
            if (work[4] == 0.0)
            {
                beta = 0.04;
            }
            else
            {
                if (work[4] < 0.0)
                {
                    beta = 0.0;
                }
                else
                {
                    beta = work[4];
                    if (beta > 0.2)
                    {
                        if (iprint > 0)
                        {
                            Console.WriteLine(" Curious input for BETA: WORK(5)=" + work[4]);
                        }
                        arret = true;
                    }
                }
            }

            // Maximal step size
            if (work[5] == 0.0)
            {
                hmax = xend - x;
            }
            else
            {
                hmax = work[5];
            }

            // Initial step size
            h = work[6];

            // Prepare the entry-points for the arrays in work
            xtemp = new double[n];
            xerr = new double[n];
            dxdt = new double[n];
            k2 = new double[n];
            k3 = new double[n];
            k4 = new double[n];
            k5 = new double[n];
            k6 = new double[n];
            double[] ysti = new double[n];
            r = new double[5 * nrdens];
            int[] icomp = new int[nrdens];
            for (int i = 0; i < nrdens; i++)
            {
                icomp[i] = iwork[20 + i];
            }

            // When a fail has occured, we return with idid=-1
            if (arret)
            {
                return -1;
            }

            this.rtol = rtol;
            this.atol = atol;

            // Call to core integrator
            int idid = dopcor_(fcn, x, y, xend, hmax, ref h,
                itol, iprint, nmax, uround,
                nstiff, safe, beta, fac1, fac2,
                ysti, icomp, nrdens,
                ref nfcn, ref nstep, ref naccpt, ref nrejct);

            work[6] = h;
            iwork[16] = nfcn;
            iwork[17] = nstep;
            iwork[18] = naccpt;
            iwork[19] = nrejct;

            return idid;
        }

        /* ---------------------------------------------------------- */
        /*     Core integrator for DOPRI5 */
        /*     Parameters same as in DOPRI5 with workspace added */
        /* ---------------------------------------------------------- */
        int dopcor_(Action<double, double[], double[]> fcn, double t, double[] x,
            double xend, double hmax, ref double dt, int itol, int iprint,
            int nmax, double uround,
            int nstiff, double safe, double beta,
            double fac1, double fac2, double[] ysti, int[] icomp,
            int nrd, ref int nfcn, ref int nstep, ref int naccpt, ref int nrejct)
        {
            /* System generated locals */
            double d__1;

            /* Local variables */
            int i, j;
            double yd0, fac, err, fac11;
            int iord;
            bool last;
            double hnew, bspl, facc1, facc2, expo1, hlamb, ydiff;
            int iasti;

            double stden;
            int irtrn = 0;
            double stnum, facold;
            bool reject;

            double posneg;
            int nonsti = 0;


            /* Function Body */
            facold = 1e-4;
            expo1 = 0.2 - beta * 0.75;
            facc1 = 1.0 / fac1;
            facc2 = 1.0 / fac2;
            posneg = d_sign(1.0, xend - t);

            // Initial preparations
            last = false;
            hlamb = 0.0;
            iasti = 0;
            fcn(t, x, dxdt);
            hmax = Math.Abs(hmax);
            iord = 5;
            if (dt == 0.0)
            {
                dt = hinit_(fcn, t, x, xend, posneg, dxdt, k2, k3, iord, hmax, itol);
            }
            nfcn += 2;
            reject = false;
            condo_1.xold = t;

            // Basic integration step
            L1:
            if (nstep > nmax)
            {
                if (iprint > 0)
                {
                    Debug.WriteLine("Exit of DOPRI5 at X=" + t);
                    Console.WriteLine(" More than NMAX =" + nmax + " steps are needed");
                }
                return -2;
            }

            if (Math.Abs(dt) * 0.1 <= Math.Abs(t) * uround)
            {
                if (iprint > 0)
                {
                    Debug.WriteLine("Exit of DOPRI5 at X=" + t);
                    Console.WriteLine(" Step size too small, H=" + dt);
                }
                return -3;
            }

            if ((t + dt * 1.01 - xend) * posneg > 0.0)
            {
                dt = xend - t;
                last = true;
            }
            ++(nstep);

            // The first 6 stages
            if (irtrn >= 2)
            {
                fcn(t, x, dxdt);
            }

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * a21 * dxdt[i];
            }

            fcn(t + c2 * dt, xtemp, k2);

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (a31 * dxdt[i] + a32 * k2[i]);
            }

            fcn(t + c3 * dt, xtemp, k3);

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (a41 * dxdt[i] + a42 * k2[i] + a43 * k3[i]);
            }

            fcn(t + c4 * dt, xtemp, k4);

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (a51 * dxdt[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
            }

            fcn(t + c5 * dt, xtemp, k5);

            for (i = 0; i < n; ++i)
            {
                ysti[i] = x[i] + dt * (a61 * dxdt[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            }

            double th = t + dt;
            fcn(th, ysti, k6);

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (a71 * dxdt[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            }

            fcn(th, xtemp, k2);

            //if (iout >= 2)
            {
                for (j = 0; j < nrd; ++j)
                {
                    i = icomp[j];
                    r[nrd * 4 + j] = dt * (d1 * dxdt[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * k2[i]);
                }
            }

            for (i = 0; i < n; ++i)
            {
                xerr[i] = dt * (e1 * dxdt[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k2[i]);
            }

            nfcn += 6;

            // Error estimation
            err = Error(dt, x, itol);

            // Computation of hnew
            fac11 = Math.Pow(err, expo1);

            // LUND-stabilization
            fac = fac11 / Math.Pow(facold, beta);

            // We require  FAC1 <= HNEW/H <= FAC2
            fac = Math.Max(facc2, Math.Min(facc1, fac / safe));
            hnew = dt / fac;
            if (err <= 1.0)
            {
                // Step is accepted
                facold = Math.Max(err, 1e-4);
                ++(naccpt);

                // Stiffness detection
                if (naccpt % nstiff == 0 || iasti > 0)
                {
                    stnum = 0.0;
                    stden = 0.0;
                    for (i = 0; i < n; ++i)
                    {
                        /* Computing 2nd power */
                        d__1 = k2[i] - k6[i];
                        stnum += d__1 * d__1;
                        /* Computing 2nd power */
                        d__1 = xtemp[i] - ysti[i];
                        stden += d__1 * d__1;
                    }
                    if (stden > 0.0)
                    {
                        hlamb = dt * Math.Sqrt(stnum / stden);
                    }
                    if (hlamb > 3.25)
                    {
                        nonsti = 0;
                        ++iasti;
                        if (iasti == 15)
                        {
                            if (iprint > 0)
                            {
                                Console.WriteLine(" The problem seems to become stiff at X = " + (t));
                            }
                            if (iprint <= 0)
                            {
                                return -4;
                            }
                        }
                    }
                    else
                    {
                        ++nonsti;
                        if (nonsti == 6)
                        {
                            iasti = 0;
                        }
                    }
                }
                
                //if (iout >= 2)
                {
                    for (j = 0; j < nrd; ++j)
                    {
                        i = icomp[j];
                        yd0 = x[i];
                        ydiff = xtemp[i] - yd0;
                        bspl = dt * dxdt[i] - ydiff;
                        r[j] = x[i];
                        r[nrd + j] = ydiff;
                        r[nrd * 2 + j] = bspl;
                        r[nrd * 3 + j] = -(dt) * k2[i] + ydiff - bspl;
                        /* L43: */
                    }
                }

                for (i = 0; i < n; ++i)
                {
                    dxdt[i] = k2[i];
                    x[i] = xtemp[i];
                }
                condo_1.xold = t;
                t = th;
                
                // Normal exit
                if (last)
                {
                    dt = hnew;
                    return 1;
                }
                if (Math.Abs(hnew) > hmax)
                {
                    hnew = posneg * hmax;
                }
                if (reject)
                {
                    hnew = posneg * Math.Min(Math.Abs(hnew), Math.Abs(dt));
                }
                reject = false;
            }
            else
            {
                // Step is rejected
                hnew = dt / Math.Min(facc1, fac11 / safe);
                reject = true;
                if (naccpt >= 1)
                {
                    ++(nrejct);
                }
                last = false;
            }
            dt = hnew;
            goto L1;
        }

        double Error(double dt, double[] x, int itol)
        {
            double err = 0.0, sk, temp;
            if (itol == 0)
            {
                double atoli = atol[0];
                double rtoli = rtol[0];

                for (int i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Max(Math.Abs(x[i]), Math.Abs(xtemp[i]));
                    /* Computing 2nd power */
                    temp = xerr[i] / sk;
                    err += temp * temp;
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Max(Math.Abs(x[i]), Math.Abs(xtemp[i]));
                    /* Computing 2nd power */
                    temp = xerr[i] / sk;
                    err += temp * temp;
                }
            }

            return Math.Sqrt(err / n);
        }

        /* ---------------------------------------------------------- */
        /* ----  Computation of an initial step size guess */
        /* ---------------------------------------------------------- */
        double hinit_(Action<double, double[], double[]> fcn, double x, double[] y,
            double xend, double posneg, double[] f0, double[] f1,
            double[] y1, int iord, double hmax, int itol)
        {
            /* System generated locals */
            double d__1;

            /* Local variables */
            double h;
            int i;
            double h1, sk, dnf, dny, der2, der12, atoli, rtoli;

            // Compute a first guess for explicit euler as
            //   H = 0.01 * NORM (Y0) / NORM (F0)
            // The increment for explicit euler is small
            // compared to the solution

            /* Function Body */
            dnf = 0.0;
            dny = 0.0;
            atoli = atol[0];
            rtoli = rtol[0];
            if (itol == 0)
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Abs(y[i]);
                    /* Computing 2nd power */
                    d__1 = f0[i] / sk;
                    dnf += d__1 * d__1;
                    /* Computing 2nd power */
                    d__1 = y[i] / sk;
                    dny += d__1 * d__1;
                }
            }
            else
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Abs(y[i]);
                    /* Computing 2nd power */
                    d__1 = f0[i] / sk;
                    dnf += d__1 * d__1;
                    /* Computing 2nd power */
                    d__1 = y[i] / sk;
                    dny += d__1 * d__1;
                }
            }
            if (dnf <= 1e-10 || dny <= 1e-10)
            {
                h = 1e-6;
            }
            else
            {
                h = Math.Sqrt(dny / dnf) * 0.01;
            }
            h = Math.Min(h, hmax);
            h = d_sign(h, posneg);

            // Perform an explicit euler step
            for (i = 0; i < n; ++i)
            {
                y1[i] = y[i] + h * f0[i];
            }

            fcn(x + h, y1, f1);

            // Estimate the second derivative of the solution
            der2 = 0.0;
            if (itol == 0)
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Abs(y[i]);
                    /* Computing 2nd power */
                    d__1 = (f1[i] - f0[i]) / sk;
                    der2 += d__1 * d__1;
                }
            }
            else
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Abs(y[i]);
                    /* Computing 2nd power */
                    d__1 = (f1[i] - f0[i]) / sk;
                    der2 += d__1 * d__1;
                }
            }
            der2 = Math.Sqrt(der2) / h;

            // Step size is computed such that
            //  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
            der12 = Math.Max(Math.Abs(der2), Math.Sqrt(dnf));
            if (der12 <= 1e-15)
            {
                h1 = Math.Max(1e-6, Math.Abs(h) * 0.001);
            }
            else
            {
                h1 = Math.Pow(0.01 / der12, 1.0 / iord);
            }

            h = Math.Min(Math.Min(Math.Abs(h) * 100, h1), hmax);

            return d_sign(h, posneg);
        }

        /* ---------------------------------------------------------- */
        /*     This function can be used for continuous output in connection */
        /*     with the output-subroutine for DOPRI5. It provides an */
        /*     approximation to the II-th component of the solution at X. */
        /* ---------------------------------------------------------- */
        public double contd5_(int ii, double x, double[] con, int[] icomp, int nd)
        {
            // Compute place of II-th component

            int i = -1;

            for (int j = 0; j < nd; ++j)
            {
                if (icomp[j] == ii)
                {
                    i = j;
                }
            }

            if (i < 0)
            {
                Console.WriteLine(" No dense output available for comp. " + ii);
                return 0.0;
            }

            double s = (x - condo_1.xold) / condo_1.hout;
            double s1 = 1.0 - s;

            return con[i] + s * (con[nd + i] + s1 * (con[2 * nd + i] + s * (con[3 * nd + i] + s1 * con[4 * nd + i])));
        }
    }
}
