
namespace MathNet.Numerics.OdeSolvers
{
    using System;
    using System.Diagnostics;

    public class DormandPrince853
    {
        public delegate void U_fp(int n, double x, double[] y, double[] y_out, double[] rpar, int[] ipar);
        public delegate void S_fp(int i, double xold, double x, double[] y, int n, double[] con, int[] icomp, int nd, double[] rpar, int[] ipar, int irtrn, ref double xout);
        
        double d_sign(double a, double b)
        {

            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        class condo8_
        {
            public double xold, hout;
        }

        condo8_ condo8_1 = new condo8_();//		double xold, hout;

        /* ----------------------------------------------------------
         *     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
         *     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
         *     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER 8(5,3)
         *     DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND
         *     DENSE OUTPUT)
         *
         *     AUTHORS: E. HAIRER AND G. WANNER
         *              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
         *              CH-1211 GENEVE 24, SWITZERLAND
         *              E-MAIL:  Ernst.Hairer@unige.ch
         *                       Gerhard.Wanner@unige.ch
         *
         *     THIS CODE IS DESCRIBED IN:
         *         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
         *         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
         *         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
         *         SPRINGER-VERLAG (1993)
         *
         *     VERSION OF OCTOBER 11, 2009
         *      (new option IOUT=3 for sparse dense output)
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
         *                 ATOL SHOULD BE STRICTLY POSITIVE (POSSIBLY VERY SMALL)
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
         *                 IF IOUT.GE.1, IT IS CALLED DURING INTEGRATION.
         *                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
         *                 IT MUST HAVE THE FORM
         *                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
         *                                       RPAR,IPAR,IRTRN,XOUT)
         *                    DIMENSION Y(N),CON(8*ND),ICOMP(ND)
         *                    ....
         *                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
         *                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
         *                    THE FIRST GRID-POINT).
         *                 "XOLD" IS THE PRECEEDING GRID-POINT.
         *                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
         *                    IS SET <0, DOP853 WILL RETURN TO THE CALLING PROGRAM.
         *                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
         *                    SET  IRTRN = 2
         *                 "XOUT" CAN BE USED FOR EFFICIENT INTERMEDIATE OUTPUT
         *                    IF ONE PUTS IOUT=3
         *                    WHEN NR=1 DEFINE THE FIRST OUTPUT POINT XOUT IN SOLOUT.
         *                      THE SUBROUTINE SOLOUT WILL BE CALLED ONLY WHEN
         *                      XOUT IS IN THE INTERVAL [XOLD,X]; DURING THIS CALL
         *                      A NEW VALUE FOR XOUT CAN BE DEFINED, ETC.
         *
         *          -----  CONTINUOUS OUTPUT: -----
         *                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
         *                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
         *                 THE FUNCTION
         *                        >>>   CONTD8(I,S,CON,ICOMP,ND)   <<<
         *                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
         *                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
         *                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
         *
         *     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
         *                    IOUT=0: SUBROUTINE IS NEVER CALLED
         *                    IOUT=1: SUBROUTINE IS CALLED AFTER EVERY SUCCESSFUL STEP
         *                    IOUT=2: DENSE OUTPUT IS PERFORMED AFTER EVERY SUCCESSFUL STEP
         *                            (IN THIS CASE IWORK(5) MUST BE SPECIFIED)
         *                    IOUT=3: DENSE OUTPUT IS PERFORMED IN STEPS DEFINED BY THE USER
         *                            (SEE "XOUT" ABOVE)
         *
         *     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
         *                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE.
         *                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
         *                 "LWORK" MUST BE AT LEAST  11*N+8*NRDENS+21
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
         *              DEFAULT VALUES: WORK(3)=0.333D0, WORK(4)=6.D0
         *
         *    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL
         *              (SEE SECTION IV.2). POSITIVE VALUES OF BETA ( <= 0.04 )
         *              MAKE THE STEP SIZE CONTROL MORE STABLE.
         *              NEGATIVE WORK(5) PROVOKE BETA=0.
         *              DEFAULT 0.0D0.
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
         *              IF IWORK(2).EQ.1  METHOD DOP853 OF DORMAND AND PRINCE
         *              (SECTION II.6).
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
        public int dop853_(int n, U_fp fcn, double x, double[] y, double xend,
            double[] rtol, double[] atol, int itol, S_fp solout, int iout, double[] work, int lwork,
            int[] iwork, int liwork, double[] rpar, int[] ipar, out int idid)
        {
            /* System generated locals */
            int i1;

            /* Local variables */
            double h;
            int i;
            double fac1, fac2;
            double beta, safe;
            int nfcn, meth;
            double hmax;
            int nmax;
            
            bool arret;
            int nstep, naccpt, nrejct, nstiff, nrdens, iprint, istore;
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
                        Console.WriteLine(" WARNING: PUT IOUT=2 OR IOUT=3 FOR DENSE OUTPUT ");
                    }
                }
                if (nrdens == n)
                {
                    for (i = 0; i < nrdens; ++i)
                    {
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
                fac1 = 0.333;
            }
            else
            {
                fac1 = work[2];
            }
            if (work[3] == 0.0)
            {
                fac2 = 6.0;
            }
            else
            {
                fac2 = work[3];
            }
            /* --------- BETA FOR STEP CONTROL STABILIZATION ----------- */
            if (work[4] == 0.0)
            {
                beta = 0.0;
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
            double[] k7 = new double[n];
            double[] k8 = new double[n];
            double[] k9 = new double[n];
            double[] k10 = new double[n];
            double[] cont = new double[nrdens << 3];
            int[] icomp = new int[nrdens];
            for (i = 0; i < nrdens; i++)
            {
                icomp[i] = iwork[20 + i];
            }

            /* ------ TOTAL STORAGE REQUIREMENT ----------- */
            istore = 21 + 12 * n + (nrdens << 3) - 1;
            if (istore > lwork)
            {
                if (iprint > 0)
                {
                    Console.WriteLine(" INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=" + istore);
                }
                //arret = true;
            }

            istore = 21 + nrdens - 1;
            if (istore > liwork)
            {
                if (iprint > 0)
                {
                    Console.WriteLine(" INSUFFICIENT STORAGE FOR IWORK, MIN. LIWORK=" + istore);
                }
                //arret = true;
            }

            /* -------- WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
            if (arret)
            {
                idid = -1;
                return 0;
            }

            /* -------- CALL TO CORE INTEGRATOR ------------ */
            dp86co_(n, fcn, x, y, xend, hmax, ref h, rtol, atol,
                itol, iprint, solout, iout, out idid, nmax, uround, meth,
                nstiff, safe, beta, fac1, fac2, k1, k2,
                k3, k4, k5, k6, k7,
                k8, k9, k10, y1, cont,
                icomp, nrdens, rpar, ipar,
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
        int dp86co_(int n, U_fp fcn, double x, double[] y, double xend, double hmax, ref double h,
            double[] rtol, double[] atol, int itol, int iprint, S_fp solout,
            int iout, out int idid, int nmax, double uround,
            int meth, int nstiff, double safe, double beta,
            double fac1, double fac2, double[] k1, double[] k2,
            double[] k3, double[] k4, double[] k5, double[] k6,
            double[] k7, double[] k8, double[] k9, double[] k10,
            double[] y1, double[] cont, int[] icomp, int nrd,
            double[] rpar, int[] ipar, ref int nfcn, ref int nstep,
            ref int naccpt, ref int nrejct)
        {
            /* System generated locals */
            double d__1;

            /* Local variables */
            int i, j;
            double xph, err, sk, fac, err2, fac11, deno;
            int iord;
            bool last;
            double erri, hnew, bspl, facc1, facc2, xout = 0.0, expo1, hlamb, ydiff,
                atoli;
            int iasti;

            double stden;
            bool _event;
            double rtoli;
            int irtrn = 0;
            double stnum, facold;
            bool reject;
            double posneg;
            int nonsti = 0;

            /* ---------------------------------------------------------- */
            /*     CORE INTEGRATOR FOR DOP853 */
            /*     PARAMETERS SAME AS IN DOP853 WITH WORKSPACE ADDED */
            /* ---------------------------------------------------------- */

            /* Function Body */
            facold = 1e-4;
            expo1 = 0.125 - beta * 0.2;
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
            iord = 8;
            if (h == 0.0)
            {
                h = hinit_(n, fcn, x, y, xend, posneg, k1, k2,
                    k3, iord, hmax, atol, rtol, itol, rpar, ipar);
            }
            nfcn += 2;
            reject = false;
            condo8_1.xold = x;
            if (iout != 0)
            {
                irtrn = 1;
                condo8_1.hout = 1.0;
                solout(naccpt + 1, condo8_1.xold, x, y, n, cont, icomp, nrd, rpar, ipar, irtrn, ref xout);
                if (irtrn < 0)
                {
                    goto L79;
                }
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
            /* --- THE TWELVE STAGES */
            if (irtrn >= 2)
            {
                fcn(n, x, y, k1, rpar, ipar);
            }
            
            for (i = 0; i < n; ++i)
            {
                /* L22: */
                y1[i] = y[i] + h * 0.0526001519587677318785587544488 * k1[i];
            }

            fcn(n, x + h * 0.0526001519587677318785587544488, y1, k2, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L23: */
                y1[i] = y[i] + h * (k1[i] * 0.0197250569845378994544595329183
                    + k2[i] * 0.0591751709536136983633785987549);
            }

            fcn(n, x + h * 0.0789002279381515978178381316732, y1, k3, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L24: */
                y1[i] = y[i] + h * (k1[i] * 0.0295875854768068491816892993775
                    + k3[i] * 0.0887627564304205475450678981324);
            }

            fcn(n, x + h * 0.11835034190722739672675719751, y1, k4, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L25: */
                y1[i] = y[i] + h * (k1[i] * 0.241365134159266685502369798665
                    + k3[i] * -0.884549479328286085344864962717
                    + k4[i] * 0.924834003261792003115737966543);
            }

            fcn(n, x + h * 0.28164965809277260327324280249, y1, k5, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L26: */
                y1[i] = y[i] + h * (k1[i] * 0.037037037037037037037037037037
                    + k4[i] * 0.170828608729473871279604482173
                    + k5[i] * 0.125467687566822425016691814123);
            }

            fcn(n, x + h * 0.333333333333333333333333333333, y1, k6, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L27: */
                y1[i] = y[i] + h * (k1[i] * 0.037109375
                    + k4[i] * 0.170252211019544039314978060272
                    + k5[i] * 0.0602165389804559606850219397283
                    + k6[i] * -0.017578125);
            }

            fcn(n, x + h * 0.25, y1, k7, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L28: */
                y1[i] = y[i] + h * (k1[i] * 0.0370920001185047927108779319836
                    + k4[i] * 0.170383925712239993810214054705
                    + k5[i] * 0.107262030446373284651809199168
                    + k6[i] * -0.0153194377486244017527936158236
                    + k7[i] * 0.00827378916381402288758473766002);
            }

            fcn(n, x + h * 0.307692307692307692307692307692, y1, k8, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L29: */
                y1[i] = y[i] + h * (k1[i] * 0.624110958716075717114429577812
                    + k4[i] * -3.36089262944694129406857109825
                    + k5[i] * -0.868219346841726006818189891453
                    + k6[i] * 27.5920996994467083049415600797
                    + k7[i] * 20.1540675504778934086186788979
                    + k8[i] * -43.4898841810699588477366255144);
            }

            fcn(n, x + h * 0.651282051282051282051282051282, y1, k9, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L30: */
                y1[i] = y[i] + h * (k1[i] * 0.477662536438264365890433908527
                    + k4[i] * -2.48811461997166764192642586468
                    + k5[i] * -0.590290826836842996371446475743
                    + k6[i] * 21.2300514481811942347288949897
                    + k7[i] * 15.2792336328824235832596922938
                    + k8[i] * -33.2882109689848629194453265587
                    + k9[i] * -0.0203312017085086261358222928593);
            }

            fcn(n, x + h * 0.6, y1, k10, rpar, ipar);
            
            for (i = 0; i < n; ++i)
            {
                /* L31: */
                y1[i] = y[i] + h * (k1[i] * -0.93714243008598732571704021658
                    + k4[i] * 5.18637242884406370830023853209
                    + k5[i] * 1.09143734899672957818500254654
                    + k6[i] * -8.14978701074692612513997267357
                    + k7[i] * -18.5200656599969598641566180701
                    + k8[i] * 22.7394870993505042818970056734
                    + k9[i] * 2.49360555267965238987089396762
                    + k10[i] * -3.0467644718982195003823669022);
            }

            fcn(n, x + h * 0.857142857142857142857142857142, y1, k2, rpar, ipar);
            xph = x + h;
            
            for (i = 0; i < n; ++i)
            {
                /* L32: */
                y1[i] = y[i] + h * (k1[i] * 2.27331014751653820792359768449
                    + k4[i] * -10.5344954667372501984066689879
                    + k5[i] * -2.00087205822486249909675718444
                    + k6[i] * -17.9589318631187989172765950534
                    + k7[i] * 27.9488845294199600508499808837
                    + k8[i] * -2.85899827713502369474065508674
                    + k9[i] * -8.87285693353062954433549289258
                    + k10[i] * 12.3605671757943030647266201528
                    + k2[i] * 0.643392746015763530355970484046);
            }

            fcn(n, xph, y1, k3, rpar, ipar);
            nfcn += 11;
            
            for (i = 0; i < n; ++i)
            {
                k4[i] = k1[i] * 0.0542937341165687622380535766363
                    + k6[i] * 4.45031289275240888144113950566
                    + k7[i] * 1.89151789931450038304281599044
                    + k8[i] * -5.8012039600105847814672114227
                    + k9[i] * 0.31116436695781989440891606237
                    + k10[i] * -0.152160949662516078556178806805
                    + k2[i] * 0.201365400804030348374776537501
                    + k3[i] * 0.0447106157277725905176885569043;
                /* L35: */
                k5[i] = y[i] + h * k4[i];
            }
            /* --- ERROR ESTIMATION */
            err = 0.0;
            err2 = 0.0;
            if (itol == 0)
            {
                
                for (i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Max(Math.Abs(y[i]), Math.Abs(k5[i]));
                    erri = k4[i] - k1[i] * 0.244094488188976377952755905512
                        - k9[i] * 0.733846688281611857341361741547
                        - k3[i] * 0.0220588235294117647058823529412;
                    /* Computing 2nd power */
                    d__1 = erri / sk;
                    err2 += d__1 * d__1;
                    erri = k1[i] * 0.01312004499419488073250102996
                        + k6[i] * -1.225156446376204440720569753
                        + k7[i] * -0.4957589496572501915214079952
                        + k8[i] * 1.664377182454986536961530415
                        + k9[i] * -0.350328848749973681688648729
                        + k10[i] * 0.3341791187130174790297318841
                        + k2[i] * 0.08192320648511571246570742613
                        + k3[i] * -0.02235530786388629525884427845;
                    /* L41: */
                    /* Computing 2nd power */
                    d__1 = erri / sk;
                    err += d__1 * d__1;
                }
            }
            else
            {
                
                for (i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Max(Math.Abs(y[i]), Math.Abs(k5[i]));
                    erri = k4[i] - k1[i] * 0.244094488188976377952755905512
                        - k9[i] * 0.733846688281611857341361741547
                        - k3[i] * 0.0220588235294117647058823529412;
                    /* Computing 2nd power */
                    d__1 = erri / sk;
                    err2 += d__1 * d__1;
                    erri = k1[i] * 0.01312004499419488073250102996
                        + k6[i] * -1.225156446376204440720569753
                        + k7[i] * -0.4957589496572501915214079952
                        + k8[i] * 1.664377182454986536961530415
                        + k9[i] * -0.350328848749973681688648729
                        + k10[i] * 0.3341791187130174790297318841
                        + k2[i] * 0.08192320648511571246570742613
                        + k3[i] * -0.02235530786388629525884427845;
                    /* L42: */
                    /* Computing 2nd power */
                    d__1 = erri / sk;
                    err += d__1 * d__1;
                }
            }
            deno = err + err2 * 0.01;
            if (deno <= 0.0)
            {
                deno = 1.0;
            }
            err = Math.Abs(h) * err * Math.Sqrt(1.0 / (n * deno));
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
                fcn(n, xph, k5, k4, rpar, ipar);
                ++(nfcn);
                /* ------- STIFFNESS DETECTION */
                if (naccpt % nstiff == 0 || iasti > 0)
                {
                    stnum = 0.0;
                    stden = 0.0;
                    
                    for (i = 0; i < n; ++i)
                    {
                        /* Computing 2nd power */
                        d__1 = k4[i] - k3[i];
                        stnum += d__1 * d__1;
                        /* Computing 2nd power */
                        d__1 = k5[i] - y1[i];
                        stden += d__1 * d__1;
                        /* L64: */
                    }
                    if (stden > 0.0)
                    {
                        hlamb = Math.Abs(h) * Math.Sqrt(stnum / stden);
                    }
                    if (hlamb > 6.1)
                    {
                        nonsti = 0;
                        ++iasti;
                        if (iasti == 15)
                        {
                            if (iprint > 0)
                            {
                                Console.WriteLine(" THE PROBLEM SEEMS TO BECOME STIFF AT X = " + x);
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
                /* ------- FINAL PREPARATION FOR DENSE OUTPUT */
                _event = iout == 3 && xout <= xph;
                if (iout == 2 || _event)
                {
                    /* ----    SAVE THE FIRST FUNCTION EVALUATIONS */
                    for (j = 0; j < nrd; ++j)
                    {
                        i = icomp[j];
                        cont[j] = y[i];
                        ydiff = k5[i] - y[i];
                        cont[j + nrd] = ydiff;
                        bspl = h * k1[i] - ydiff;
                        cont[j + (nrd << 1)] = bspl;
                        cont[j + nrd * 3] = ydiff - h * k4[i] - bspl;
                        cont[j + (nrd << 2)] = k1[i] *
                            -8.4289382761090128651353491142 + k6[i] *
                            .5667149535193777696253178359 + k7[i] *
                            -3.0689499459498916912797304727 + k8[i] *
                            2.384667656512069828772814968 + k9[i] *
                            2.1170345824450282767155149946 + k10[i] *
                            -0.8713915837779729920678990749 + k2[i] *
                            2.240437430260788275854177165 + k3[i] *
                            .6315787787694688181557024929;
                        cont[j + nrd * 5] = k1[i] * 10.427508642579134603413151009
                            + k6[i] * 242.28349177525818288430175319 + k7[i] *
                            165.20045171727028198505394887 + k8[i] *
                            -374.54675472269020279518312152 + k9[i] *
                            -22.113666853125306036270938578 + k10[i] *
                            7.7334326684722638389603898808 + k2[i] *
                            -30.674084731089398182061213626 + k3[i] *
                            -9.3321305264302278729567221706;
                        cont[j + nrd * 6] = k1[i] * 19.985053242002433820987653617
                            + k6[i] * -387.03730874935176555105901742 + k7[i]
                            * -189.17813819516756882830838328 + k8[i] *
                            527.80815920542364900561016686 + k9[i] *
                            -11.573902539959630126141871134 + k10[i] *
                            6.8812326946963000169666922661 + k2[i] *
                            -1.000605096691083840318386098 + k3[i] *
                            .7777137798053443209286926574;
                        cont[j + nrd * 7] = k1[i] *
                            -25.693933462703749003312586129 + k6[i] *
                            -154.18974869023643374053993627 + k7[i] *
                            -231.52937917604549567536039109 + k8[i] *
                            357.6391179106141237828534991 + k9[i] *
                            93.405324183624310003907691704 + k10[i] *
                            -37.458323136451633156875139351 + k2[i] *
                            104.09964950896230045147246184 + k3[i] *
                            29.840293426660503123344363579;
                        /* L62: */
                    }
                    /* ---     THE NEXT THREE FUNCTION EVALUATIONS */
                    
                    for (i = 0; i < n; ++i)
                    {
                        /* L51: */
                        y1[i] = y[i] + h * (k1[i] *
                            .0561675022830479523392909219681 + k7[i] *
                            .253500210216624811088794765333 + k8[i] *
                            -0.246239037470802489917441475441 + k9[i] *
                            -0.124191423263816360469010140626 + k10[i] *
                            .15329179827876569731206322685 + k2[i] *
                            .00820105229563468988491666602057 + k3[i] *
                            .00756789766054569976138603589584 + k4[i] *
                            -0.008298);
                    }

                    fcn(n, x + h * 0.1, y1, k10, rpar, ipar);
                    
                    for (i = 0; i < n; ++i)
                    {
                        /* L52: */
                        y1[i] = y[i] + h * (k1[i] *
                            .0318346481635021405060768473261 + k6[i] *
                            .0283009096723667755288322961402 + k7[i] *
                            .0535419883074385676223797384372 + k8[i] *
                            -0.0549237485713909884646569340306 + k2[i] *
                            -1.08347328697249322858509316994e-4 + k3[i] *
                            3.82571090835658412954920192323e-4 + k4[i] *
                            -3.40465008687404560802977114492e-4 + k10[i] *
                            .141312443674632500278074618366);
                    }

                    fcn(n, x + h * 0.2, y1, k2, rpar, ipar);
                    
                    for (i = 0; i < n; ++i)
                    {
                        /* L53: */
                        y1[i] = y[i] + h * (k1[i] *
                            -0.428896301583791923408573538692 + k6[i] *
                            -4.69762141536116384314449447206 + k7[i] *
                            7.68342119606259904184240953878 + k8[i] *
                            4.06898981839711007970213554331 + k9[i] *
                            .356727187455281109270669543021 + k4[i] *
                            -0.00139902416515901462129418009734 + k10[i] *
                            2.9475147891527723389556272149 + k2[i] *
                            -9.15095847217987001081870187138);
                    }

                    fcn(n, x + h * 0.777777777777777777777777777778, y1, k3, rpar, ipar);
                    nfcn += 3;
                    /* ---     FINAL PREPARATION */
                    for (j = 0; j < nrd; ++j)
                    {
                        i = icomp[j];
                        cont[j + (nrd << 2)] = h * (cont[j + (nrd << 2)]
                            + k4[i] * -0.0889903364513333108206981174
                            + k10[i] * 18.148505520854727256656404962
                            + k2[i] * -9.1946323924783554000451984436
                            + k3[i] * -4.4360363875948939664310572);
                        cont[j + nrd * 5] = h * (cont[j + nrd * 5]
                            + k4[i] * 15.697238121770843886131091075
                            + k10[i] * -31.139403219565177677282850411
                            + k2[i] * -9.3529243588444783865713862664
                            + k3[i] * 35.81684148639408375246589854);
                        cont[j + nrd * 6] = h * (cont[j + nrd * 6]
                            + k4[i] * -2.7782057523535084065932004339
                            + k10[i] * -60.196695231264120758267380846
                            + k2[i] * 84.320405506677161018159903784
                            + k3[i] * 11.99229113618278932803513003);
                        cont[j + nrd * 7] = h * (cont[j + nrd * 7]
                            + k4[i] * -43.533456590011143754432175058
                            + k10[i] * 96.3245539591882829483949506
                            + k2[i] * -39.177261675615439165231486172
                            + k3[i] * -149.72683625798562581422125276);
                        /* L63: */
                    }
                    condo8_1.hout = h;
                }
                
                for (i = 0; i < n; ++i)
                {
                    k1[i] = k4[i];
                    /* L67: */
                    y[i] = k5[i];
                }
                condo8_1.xold = x;
                x = xph;
                if (iout == 1 || iout == 2 || _event)
                {
                    solout(naccpt + 1, condo8_1.xold, x, y, n, cont, icomp, nrd, rpar, ipar, irtrn, ref xout);
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
                    /* Computing MIN */
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
                Debug.WriteLine("EXIT OF DOP853 AT X=" + x);
                Console.WriteLine(" STEP SIZE TOO SMALL, H=" + h);
            }
            idid = -3;
            return 0;

            L78:
            if (iprint > 0)
            {
                Debug.WriteLine("EXIT OF DOP853 AT X=" + x);
                Console.WriteLine(" MORE THAN NMAX =" + nmax + "STEPS ARE NEEDED");
            }
            idid = -2;
            return 0;

            L79:
            if (iprint > 0)
            {
                Debug.WriteLine("EXIT OF DOP853 AT X=" + x);
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
                h = Math.Sqrt(dny / dnf) * 0.01;
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
                for (i = 0; i < n; ++i)
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
        /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONNECTION */
        /*     WITH THE OUTPUT-SUBROUTINE FOR DOP853. IT PROVIDES AN */
        /*     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X. */
        /* ---------------------------------------------------------- */
        public double contd8_(int ii, double x, double[] con, int[] icomp, int nd)
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

            double s = (x - condo8_1.xold) / condo8_1.hout; // condo8_2.h
            double s1 = 1.0 - s;
            double conpar = con[i + (nd << 2)] + s * (con[i + nd * 5] + s1 * (con[i + nd * 6] + s * con[i + nd * 7]));

            return con[i] + s * (con[i + nd] + s1 * (con[i + (nd << 1)] + s * (con[i + nd * 3] + s1 * conpar)));
        }
    }
}
