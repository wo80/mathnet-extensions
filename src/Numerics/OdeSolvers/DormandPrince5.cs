
namespace MathNet.Numerics.OdeSolvers
{
    using System;
    using System.Diagnostics;

    public class DormandPrince5
    {
        public delegate void U_fp(int n, double x, double[] y, double[] y_out, double[] rpar, int[] ipar);
        public delegate void S_fp(int i, double x_old, double x, double[] y, int n, double[] cont, int[] icomp, int nrd, double[] rpar, int[] ipar, int irtrn);
        
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

        // ---- DENSE OUTPUT OF SHAMPINE (1986) private const double /
        private const double d1 = -1.1270175653862835;
        private const double d3 = 2.675424484351598;
        private const double d4 = -5.6855269615885042;
        private const double d5 = 3.5219323679207912;
        private const double d6 = -1.7672812570757455;
        private const double d7 = 2.3824689317781438;

        #endregion

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        class condo5_
        {
            public double xold, hout;
        }
        
        condo5_ condo5_1 = new condo5_();//		double xold, hout;

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
        public int dopri5_(int n, U_fp fcn, double x, double[] y, double xend, double[] rtol, double[] atol,
                int itol, S_fp solout, int iout, double[] work, int lwork,
                int[] iwork, int liwork, double[] rpar, int[] ipar, out int idid)
        {
            /* Local variables */
            double h;
            int i;
            double fac1, fac2;
            double beta, safe;
            int nfcn, meth;
            double hmax;
            int nmax;
            bool arret;
            int nstep, naccpt, nrejct;

            int nstiff, nrdens, iprint, istore;
            double uround;
            
            /* Function Body */
            nfcn = 0;
            nstep = 0;
            naccpt = 0;
            nrejct = 0;
            arret = false;
            /* -------- IPRINT FOR MONITORING THE PRINTING */
            if (iwork[2] == 0)
            {
                iprint = 6;
            }
            else
            {
                iprint = iwork[2];
            }
            /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
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
                        Console.WriteLine(" WRONG INPUT IWORK(1)=" + iwork[0]);
                    }
                    arret = true;
                }
            }
            /* -------- METH   COEFFICIENTS OF THE METHOD */
            if (iwork[1] == 0)
            {
                meth = 1;
            }
            else
            {
                meth = iwork[1];
                if (meth <= 0 || meth >= 4)
                {
                    if (iprint > 0)
                    {
                        Console.WriteLine(" CURIOUS INPUT IWORK(2)=" + iwork[1]);
                    }
                    arret = true;
                }
            }
            /* -------- NSTIFF   PARAMETER FOR STIFFNESS DETECTION */
            nstiff = iwork[3];
            if (nstiff == 0)
            {
                nstiff = 1000;
            }
            if (nstiff < 0)
            {
                nstiff = nmax + 10;
            }
            /* -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS */
            nrdens = iwork[4];
            if (nrdens < 0 || nrdens > n)
            {
                if (iprint > 0)
                {
                    Console.WriteLine(" CURIOUS INPUT IWORK(5)=" + iwork[4]);
                }
                arret = true;
            }
            else
            {
                if (nrdens > 0 && iout < 2)
                {
                    if (iprint > 0)
                    {
                        Console.WriteLine(" WARNING: PUT IOUT=2 FOR DENSE OUTPUT ");
                    }
                }
                if (nrdens == n)
                {
                    for (i = 0; i < nrdens; ++i)
                    {
                        /* L16: */
                        iwork[i + 20] = i;
                    }
                }
            }
            /* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
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
                        Console.WriteLine(" WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:" + work[0]);
                    }
                    arret = true;
                }
            }
            /* -------  SAFETY FACTOR ------------- */
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
                        Console.WriteLine(" CURIOUS INPUT FOR SAFETY FACTOR WORK(2)=" + work[1]);
                    }
                    arret = true;
                }
            }
            /* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
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
            /* --------- BETA FOR STEP CONTROL STABILIZATION ----------- */
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
                            Console.WriteLine(" CURIOUS INPUT FOR BETA: WORK(5)=" + work[4]);
                        }
                        arret = true;
                    }
                }
            }
            /* -------- MAXIMAL STEP SIZE */
            if (work[5] == 0.0)
            {
                hmax = xend - x;
            }
            else
            {
                hmax = work[5];
            }
            /* -------- INITIAL STEP SIZE */
            h = work[6];
            /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
            double[] y1 = new double[n];
            double[] k1 = new double[n];
            double[] k2 = new double[n];
            double[] k3 = new double[n];
            double[] k4 = new double[n];
            double[] k5 = new double[n];
            double[] k6 = new double[n];
            double[] ysti = new double[n];
            double[] cont = new double[5 * nrdens];
            int[] icomp = new int[nrdens];
            for (i = 0; i < nrdens; i++)
            {
                icomp[i] = iwork[20 + i];
            }

            /* ------ TOTAL STORAGE REQUIREMENT ----------- */
            istore = 21 + 8 * n + nrdens * 5 - 1;

            if (istore > lwork)
            {
                if (iprint > 0)
                {
                    Console.WriteLine(" INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=" + istore);
                }
                arret = true;
            }

            istore = 21 + nrdens - 1;
            if (istore > liwork)
            {
                if (iprint > 0)
                {
                    Console.WriteLine(" INSUFFICIENT STORAGE FOR IWORK, MIN. LIWORK=" + istore);
                }
                arret = true;
            }

            /* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
            if (arret)
            {
                idid = -1;
                return 0;
            }

            /* -------- CALL TO CORE INTEGRATOR ------------ */
            dopcor_(n, fcn, x, y, xend, hmax, ref h, rtol, atol,
                itol, iprint, solout, iout, out idid, nmax, uround, meth,
                nstiff, safe, beta, fac1, fac2, y1, k1,
                k2, k3, k4, k5, k6,
                ysti, cont, icomp, nrdens, rpar, ipar,
                ref nfcn, ref nstep, ref naccpt, ref nrejct);

            work[6] = h;
            iwork[16] = nfcn;
            iwork[17] = nstep;
            iwork[18] = naccpt;
            iwork[19] = nrejct;
            /* ----------- RETURN ----------- */
            return 0;
        }
        
        /*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

        /* Subroutine */
        int dopcor_(int n, U_fp fcn, double x, double[] y,
            double xend, double hmax, ref double h, double[] rtol, double[] atol, int itol, int iprint, S_fp solout,
            int iout, out int idid, int nmax, double uround,
            int meth, int nstiff, double safe, double beta,
            double fac1, double fac2, double[] y1, double[] k1,
            double[] k2, double[] k3, double[] k4, double[] k5,
            double[] k6, double[] ysti, double[] cont, int[] icomp,
            int nrd, double[] rpar, int[] ipar, ref int nfcn, ref int nstep, ref int naccpt, ref int nrejct)
        {
            /* System generated locals */
            double d__1;

            /* Local variables */
            int i, j;
            double sk, yd0, fac, err, xph, fac11;
            int iord;
            bool last;
            double hnew, bspl, facc1, facc2, expo1, hlamb, ydiff, atoli;
            int iasti;

            double stden, rtoli;
            int irtrn;
            double stnum, facold;
            bool reject;

            double posneg;
            int nonsti = 0;

            /* ---------------------------------------------------------- */
            /*     CORE INTEGRATOR FOR DOPRI5 */
            /*     PARAMETERS SAME AS IN DOPRI5 WITH WORKSPACE ADDED */
            /* ---------------------------------------------------------- */

            /* Function Body */
            if (meth == 1)
            {
                //cdopri_(...);
            }
            facold = 1e-4;
            expo1 = 0.2 - beta * 0.75;
            facc1 = 1.0 / fac1;
            facc2 = 1.0 / fac2;
            posneg = d_sign(1.0, xend - x);
            /* --- INITIAL PREPARATIONS */
            atoli = atol[0];
            rtoli = rtol[0];
            last = false;
            hlamb = 0.0;
            iasti = 0;
            fcn(n, x, y, k1, rpar, ipar);
            hmax = Math.Abs(hmax);
            iord = 5;
            if (h == 0.0)
            {
                h = hinit_(n, fcn, x, y, xend, posneg, k1, k2, k3, iord, hmax, atol, rtol, itol, rpar, ipar);
            }
            nfcn += 2;
            reject = false;
            condo5_1.xold = x;
            if (iout != 0)
            {
                irtrn = 1;
                condo5_1.hout = h;

                solout(naccpt + 1, condo5_1.xold, x, y, n, cont, icomp, nrd, rpar, ipar, irtrn);
                if (irtrn < 0)
                {
                    goto L79;
                }
            }
            else
            {
                irtrn = 0;
            }
            /* --- BASIC INTEGRATION STEP */
            L1:
            if (nstep > nmax)
            {
                goto L78;
            }
            if (Math.Abs(h) * 0.1 <= Math.Abs(x) * uround)
            {
                goto L77;
            }
            if ((x + h * 1.01 - xend) * posneg > 0.0)
            {
                h = xend - x;
                last = true;
            }
            ++(nstep);
            /* --- THE FIRST 6 STAGES */
            if (irtrn >= 2)
            {
                fcn(n, x, y, k1, rpar, ipar);
            }
            for (i = 0; i < n; ++i)
            {
                /* L22: */
                y1[i] = y[i] + h * a21 * k1[i];
            }

            fcn(n, x + c2 * h, y1, k2, rpar, ipar);

            for (i = 0; i < n; ++i)
            {
                /* L23: */
                y1[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
            }

            fcn(n, x + c3 * h, y1, k3, rpar, ipar);

            for (i = 0; i < n; ++i)
            {
                /* L24: */
                y1[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
            }

            fcn(n, x + c4 * h, y1, k4, rpar, ipar);

            for (i = 0; i < n; ++i)
            {
                /* L25: */
                y1[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
            }

            fcn(n, x + c5 * h, y1, k5, rpar, ipar);

            for (i = 0; i < n; ++i)
            {
                /* L26: */
                ysti[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            }

            xph = x + h;
            fcn(n, xph, ysti, k6, rpar, ipar);

            for (i = 0; i < n; ++i)
            {
                /* L27: */
                y1[i] = y[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            }

            fcn(n, xph, y1, k2, rpar, ipar);

            if (iout >= 2)
            {
                for (j = 0; j < nrd; ++j)
                {
                    i = icomp[j];
                    cont[nrd * 4 + j] = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * k2[i]);
                }
            }

            for (i = 0; i < n; ++i)
            {
                k4[i] = (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k2[i]) * h;
            }

            nfcn += 6;
            /* --- ERROR ESTIMATION */
            err = 0.0;
            if (itol == 0)
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Max(Math.Abs(y[i]), Math.Abs(y1[i]));
                    /* L41: */
                    /* Computing 2nd power */
                    d__1 = k4[i] / sk;
                    err += d__1 * d__1;
                }
            }
            else
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Max(Math.Abs(y[i]), Math.Abs(y1[i]));
                    /* L42: */
                    /* Computing 2nd power */
                    d__1 = k4[i] / sk;
                    err += d__1 * d__1;
                }
            }
            err = Math.Sqrt(err / n);
            /* --- COMPUTATION OF HNEW */
            fac11 = Math.Pow(err, expo1);
            /* --- LUND-STABILIZATION */
            fac = fac11 / Math.Pow(facold, beta);
            /* --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2 */
            fac = Math.Max(facc2, Math.Min(facc1, fac / safe));
            hnew = h / fac;
            if (err <= 1.0)
            {
                /* --- STEP IS ACCEPTED */
                facold = Math.Max(err, 1e-4);
                ++(naccpt);
                /* ------- STIFFNESS DETECTION */
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
                        d__1 = y1[i] - ysti[i];
                        stden += d__1 * d__1;
                        /* L64: */
                    }
                    if (stden > 0.0)
                    {
                        hlamb = h * Math.Sqrt(stnum / stden);
                    }
                    if (hlamb > 3.25)
                    {
                        nonsti = 0;
                        ++iasti;
                        if (iasti == 15)
                        {
                            if (iprint > 0)
                            {
                                Console.WriteLine(" THE PROBLEM SEEMS TO BECOME STIFF AT X = " + (x));
                            }
                            if (iprint <= 0)
                            {
                                goto L76;
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
                if (iout >= 2)
                {
                    for (j = 0; j < nrd; ++j)
                    {
                        i = icomp[j];
                        yd0 = y[i];
                        ydiff = y1[i] - yd0;
                        bspl = h * k1[i] - ydiff;
                        cont[j] = y[i];
                        cont[nrd + j] = ydiff;
                        cont[nrd * 2 + j] = bspl;
                        cont[nrd * 3 + j] = -(h) * k2[i] + ydiff - bspl;
                        /* L43: */
                    }
                }
                for (i = 0; i < n; ++i)
                {
                    k1[i] = k2[i];
                    /* L44: */
                    y[i] = y1[i];
                }
                condo5_1.xold = x;
                x = xph;
                if (iout != 0)
                {
                    condo5_1.hout = h;
                    solout(naccpt + 1, condo5_1.xold, x, y, n, cont, icomp, nrd, rpar, ipar, irtrn);
                    if (irtrn < 0)
                    {
                        goto L79;
                    }
                }
                /* ------- NORMAL EXIT */
                if (last)
                {
                    h = hnew;
                    idid = 1;
                    return 0;
                }
                if (Math.Abs(hnew) > hmax)
                {
                    hnew = posneg * hmax;
                }
                if (reject)
                {
                    hnew = posneg * Math.Min(Math.Abs(hnew), Math.Abs(h));
                }
                reject = false;
            }
            else
            {
                /* --- STEP IS REJECTED */
                hnew = h / Math.Min(facc1, fac11 / safe);
                reject = true;
                if (naccpt >= 1)
                {
                    ++(nrejct);
                }
                last = false;
            }
            h = hnew;
            goto L1;
            /* --- FAIL EXIT */
            L76:
            idid = -4;
            return 0;
            L77:
            if (iprint > 0)
            {
                Debug.WriteLine("EXIT OF DOPRI5 AT X=" + x);
                Console.WriteLine(" STEP SIZE T0O SMALL, H=" + h);
            }
            idid = -3;
            return 0;
            L78:
            if (iprint > 0)
            {
                Debug.WriteLine("EXIT OF DOPRI5 AT X=" + x);
                Console.WriteLine(" MORE THAN NMAX =" + nmax + " STEPS ARE NEEDED");
            }
            idid = -2;
            return 0;
            L79:
            if (iprint > 0)
            {
                Debug.WriteLine("EXIT OF DOPRI5 AT X=" + x);
            }
            idid = 2;
            return 0;
        }

        /* ---------------------------------------------------------- */
        /* ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS */
        /* ---------------------------------------------------------- */
        double hinit_(int n, U_fp fcn, double x, double[] y,
            double xend, double posneg, double[] f0, double[] f1,
            double[] y1, int iord, double hmax, double[] atol,
            double[] rtol, int itol, double[] rpar, int[] ipar)
        {
            /* System generated locals */
            double d__1;

            /* Local variables */
            double h;
            int i;
            double h1, sk, dnf, dny, der2, der12, atoli, rtoli;

            /* ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS */
            /* ----   H = 0.01 * NORM (Y0) / NORM (F0) */
            /* ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL */
            /* ---- COMPARED TO THE SOLUTION */

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
                    /* L10: */
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
                    /* L11: */
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
                h = Math.Sqrt(dny / dnf) * .01;
            }
            h = Math.Min(h, hmax);
            h = d_sign(h, posneg);
            /* ---- PERFORM AN EXPLICIT EULER STEP */
            for (i = 0; i < n; ++i)
            {
                /* L12: */
                y1[i] = y[i] + h * f0[i];
            }
            fcn(n, x + h, y1, f1, rpar, ipar);
            /* ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION */
            der2 = 0.0;
            if (itol == 0)
            {
                for (i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Abs(y[i]);
                    /* L15: */
                    /* Computing 2nd power */
                    d__1 = (f1[i] - f0[i]) / sk;
                    der2 += d__1 * d__1;
                }
            }
            else
            {
                for (i = 1; i <= n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Abs(y[i]);
                    /* L16: */
                    /* Computing 2nd power */
                    d__1 = (f1[i] - f0[i]) / sk;
                    der2 += d__1 * d__1;
                }
            }
            der2 = Math.Sqrt(der2) / h;
            /* ---- STEP SIZE IS COMPUTED SUCH THAT */
            /* ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01 */
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
        /*     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION */
        /*     WITH THE OUTPUT-SUBROUTINE FOR DOPRI5. IT PROVIDES AN */
        /*     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X. */
        /* ---------------------------------------------------------- */
        public double contd5_(int ii, double x, double[] con, int[] icomp, int nd)
        {
            /* ----- COMPUTE PLACE OF II-TH COMPONENT */
            
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
                Console.WriteLine(" NO DENSE OUTPUT AVAILABLE FOR COMP. " + ii);
                return 0.0;
            }

            double theta = (x - condo5_1.xold) / condo5_1.hout; // condo5_2.h
            double theta1 = 1.0 - theta;

            return con[i] + theta * (con[nd + i] + theta1 * (con[(nd << 1) + i] + theta * (con[nd * 3 + i] + theta1 * con[(nd << 2) + i])));
        }
    }
}
