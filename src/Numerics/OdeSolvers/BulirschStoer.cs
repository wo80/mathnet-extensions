// Based on Fortran code ODEX
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers
{
    using System;
    using S_fp = System.Action<int, double, double[], double[]>;

    /// <summary>
    /// Numerical solution of a system of first order ordinary differential equations  y'=f(x,y).
    /// 
    /// This is an extrapolation-algorithm (GBS / Gragg–Bulirsch–Stoer), based on the explicit
    /// midpoint rule (with stepsize control, order selection and dense output).
    ///
    /// Authors: E. Hairer and G. Wanner (dense output written by E. Hairer and A. Ostermann)
    ///
    /// This code is described in:
    ///         E. Hairer, S.P. Norsett and G. Wanner
    ///         Solving Ordinary Differential Equations I. Nonstiff Problems (2nd edition)
    ///         Springer-Verlag (1993)
    /// </summary>
    public class BulirschStoer
    {
        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        /* Common Block Declarations */

        struct conodx1
        {
            public double xoldd, hhh;
            public int kmit;
        }
        struct conodx2
        {
            public double xold, h;
            public int imit;
        }

        conodx1 conodx_1;
        conodx2 conodx_2;
        
        public int odex_(int n, S_fp fcn, double x, double[] y,
             double xend, double h, double[] rtol, double[]
            atol, int itol, S_fp solout, int iout, double[] work,
            int lwork, int[] iwork, int liwork, double[] rpar,
            int[] ipar, int idid)
        {
            /* Local variables */
            int i, km, nrd;
            double fac1, fac2, fac3, fac4;
            int nfcn;
            double hmax;
            int ncom = 0, nmax;
            double safe1, safe2, safe3;
            int jstab, mudif, iderr = 0, mstab;
            bool arret;
            int nstep, nsequ, lfsafe, naccpt, nrejct, nrdens;
            int istore;
            double uround;

            /**
             *     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
             *     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
             *     THIS IS AN EXTRAPOLATION-ALGORITHM (GBS), BASED ON THE
             *     EXPLICIT MIDPOINT RULE (WITH STEPSIZE CONTROL,
             *     ORDER SELECTION AND DENSE OUTPUT).
             *
             *     AUTHORS: E. HAIRER AND G. WANNER
             *              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
             *              CH-1211 GENEVE 24, SWITZERLAND
             *              E-MAIL:  Ernst.Hairer@unige.ch
             *                       Gerhard.Wanner@unige.ch
             *              DENSE OUTPUT WRITTEN BY E. HAIRER AND A. OSTERMANN
             *
             *     THIS CODE IS DESCRIBED IN SECTION II.9 OF THE BOOK:
             *         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
             *         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
             *         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
             *         SPRINGER-VERLAG (1993)
             *
             *     BASED ON VERSION SEPTEMBER 30, 1995
             *         SMALL CORRECTIONS ON OCTOBER 11, 2009
             *         SMALL CHANGE FOR FEATURE IDID=2, OCTOBER 4, 2015
             *               INTERRUPTED BY SOLOUT (BY C. LUDWIG)
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
             *     H           INITIAL STEP SIZE GUESS;
             *                 H=1.D0/(NORM OF F'), USUALLY 1.D-1 OR 1.D-3, IS GOOD.
             *                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
             *                 ADAPTS ITS STEP SIZE. WHEN YOU ARE NOT SURE, THEN
             *                 STUDY THE CHOSEN VALUES FOR A FEW
             *                 STEPS IN SUBROUTINE "SOLOUT".
             *                 (IF H=0.D0, THE CODE PUTS H=1.D-4).
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
             *                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND,
             *                                       RPAR,IPAR,IRTRN)
             *                    DIMENSION X,Y(N),CON(NCON),ICOMP(ND)
             *                    ....
             *                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
             *                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
             *                    THE FIRST GRID-POINT).
             *                 "XOLD" IS THE PRECEEDING GRID-POINT.
             *                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
             *                    IS SET <0, ODEX WILL RETURN TO THE CALLING PROGRAM.
             *
             *          -----  CONTINUOUS OUTPUT (IF IOUT=2): -----
             *                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
             *                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
             *                 THE DOUBLE PRECISION FUNCTION
             *                    >>>   CONTEX(I,S,CON,NCON,ICOMP,ND)
             *                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
             *                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
             *                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
             *
             *     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
             *                    IOUT=0: SUBROUTINE IS NEVER CALLED
             *                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
             *                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
             *
             *     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
             *                 SERVES AS WORKING SPACE FOR ALL VECTORS.
             *                 "LWORK" MUST BE AT LEAST
             *                    N*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS
             *                 WHERE NRDENS=IWORK(8) (SEE BELOW) AND
             *                        KM=9                IF IWORK(2)=0
             *                        KM=IWORK(2)         IF IWORK(2).GT.0
             *                 WORK(1),...,WORK(20) SERVE AS PARAMETERS
             *                 FOR THE CODE. FOR STANDARD USE, SET THESE
             *                 PARAMETERS TO ZERO BEFORE CALLING.
             *
             *     LWORK       DECLARED LENGTH OF ARRAY "WORK".
             *
             *     IWORK       int WORKING SPACE OF LENGTH "LIWORK".
             *                 "LIWORK" MUST BE AT LEAST
             *                               2*KM+21+NRDENS
             *                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS
             *                 FOR THE CODE. FOR STANDARD USE, SET THESE
             *                 PARAMETERS TO ZERO BEFORE CALLING.
             *
             *     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
             *
             *     RPAR, IPAR  REAL AND int PARAMETERS (OR PARAMETER ARRAYS) WHICH
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
             *    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
             *
             *    WORK(3)   STEP SIZE IS REDUCED BY FACTOR WORK(3), IF THE
             *              STABILITY CHECK IS NEGATIVE, DEFAULT 0.5.
             *
             *    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION
             *              THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
             *              CHOSEN SUBJECT TO THE RESTRICTION
             *                 FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN
             *              WHERE FACMIN=WORK(4)**(1/(2*J-1))
             *              DEFAULT VALUES: WORK(4)=0.02D0, WORK(5)=4.D0
             *
             *    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION
             *              STEP SIZE IS DECREASED IF    W(K-1) <= W(K)*WORK(6)
             *              STEP SIZE IS INCREASED IF    W(K) <= W(K-1)*WORK(7)
             *              DEFAULT VALUES: WORK(6)=0.8D0, WORK(7)=0.9D0
             *
             *    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM
             *             HNEW=H*WORK(9)*(WORK(8)*TOL/ERR)**(1/(J-1))
             *             DEFAULT VALUES: WORK(8)=0.65D0,
             *                        WORK(9)=0.94D0  IF "HOPE FOR CONVERGENCE"
             *                        WORK(9)=0.90D0  IF "NO HOPE FOR CONVERGENCE"
             *
             *    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
             *              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 10000.
             *
             *    IWORK(2)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
             *              TABLE. THE DEFAULT VALUE (FOR IWORK(2)=0) IS 9.
             *              IF IWORK(2).NE.0 THEN IWORK(2) SHOULD BE .GE.3.
             *
             *    IWORK(3)  SWITCH FOR THE STEP SIZE SEQUENCE (EVEN NUMBERS ONLY)
             *              IF IWORK(3).EQ.1 THEN 2,4,6,8,10,12,14,16,...
             *              IF IWORK(3).EQ.2 THEN 2,4,8,12,16,20,24,28,...
             *              IF IWORK(3).EQ.3 THEN 2,4,6,8,12,16,24,32,...
             *              IF IWORK(3).EQ.4 THEN 2,6,10,14,18,22,26,30,...
             *              IF IWORK(3).EQ.5 THEN 4,8,12,16,20,24,28,32,...
             *              THE DEFAULT VALUE IS IWORK(3)=1 IF IOUT.LE.1;
             *              THE DEFAULT VALUE IS IWORK(3)=4 IF IOUT.GE.2.
             *
             *    IWORK(4)  STABILITY CHECK IS ACTIVATED AT MOST IWORK(4) TIMES IN
             *              ONE LINE OF THE EXTRAP. TABLE, DEFAULT IWORK(4)=1.
             *
             *    IWORK(5)  STABILITY CHECK IS ACTIVATED ONLY IN THE LINES
             *              1 TO IWORK(5) OF THE EXTRAP. TABLE, DEFAULT IWORK(5)=1.
             *
             *    IWORK(6)  IF  IWORK(6)=0  ERROR ESTIMATOR IN THE DENSE
             *              OUTPUT FORMULA IS ACTIVATED. IT CAN BE SUPPRESSED
             *              BY PUTTING IWORK(6)=1.
             *              DEFAULT IWORK(6)=0  (IF IOUT.GE.2).
             *
             *    IWORK(7)  DETERMINES THE DEGREE OF INTERPOLATION FORMULA
             *              MU = 2 * KAPPA - IWORK(7) + 1
             *              IWORK(7) SHOULD LIE BETWEEN 1 AND 6
             *              DEFAULT IWORK(7)=4  (IF IWORK(7)=0).
             *
             *    IWORK(8)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
             *              IS REQUIRED
             *
             *    IWORK(21),...,IWORK(NRDENS+20) INDICATE THE COMPONENTS, FOR WHICH
             *              DENSE OUTPUT IS REQUIRED
             *
             * ----------------------------------------------------------------------C
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
             *                   IDID=1  COMPUTATION SUCCESSFUL,
             *                   IDID=2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
             *                   IDID=-1 COMPUTATION UNSUCCESSFUL.
             *
             *   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
             *   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
             *   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
             *   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
             *                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
             */

            /* Parameter adjustments */
            //--y;
            //--rtol;
            //--atol;
            //--work;
            //--iwork;
            //--rpar;
            //--ipar;

            /* Function Body */
            nfcn = 0;
            nstep = 0;
            naccpt = 0;
            nrejct = 0;
            arret = false;
            /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
            if (iwork[0] == 0)
            {
                nmax = 10000;
            }
            else
            {
                nmax = iwork[0];
                if (nmax <= 0)
                {

                    Console.WriteLine(" WRONG INPUT IWORK(1)=", iwork[0]);

                    arret = true;
                }
            }
            /* -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION */
            if (iwork[2] == 0)
            {
                km = 9;
            }
            else
            {
                km = iwork[2];
                if (km <= 2)
                {

                    Console.WriteLine(" CURIOUS INPUT IWORK(2)=", iwork[2]);

                    arret = true;
                }
            }
            /* -------- NSEQU     CHOICE OF STEP SIZE SEQUENCE */
            nsequ = iwork[3];
            if (iwork[3] == 0 && iout <= 1)
            {
                nsequ = 1;
            }
            if (iwork[3] == 0 && iout >= 2)
            {
                nsequ = 4;
            }
            if (nsequ <= 0 || nsequ >= 6)
            {

                Console.WriteLine(" CURIOUS INPUT IWORK(3)=", iwork[3]);

                arret = true;
            }
            if (nsequ <= 3 && iout >= 2)
            {

                Console.WriteLine(" IWORK(3) NOT COMPATIBLE WITH IOUT");

                arret = true;
            }
            /* -------- MSTAB     PARAMETER FOR STABILITY CHECK */
            if (iwork[4] == 0)
            {
                mstab = 1;
            }
            else
            {
                mstab = iwork[4];
            }
            /* -------- JSTAB     PARAMETER FOR STABILITY CHECK */
            if (iwork[5] == 0)
            {
                jstab = 2;
            }
            else
            {
                jstab = iwork[5];
            }
            /* -------- IDERR  PARAMETER FOR ERROR ESTIMATION IN DENSE OUTPUT */
            if (iwork[6] == 0)
            {
                if (iout <= 1)
                {
                    iderr = 1;
                }
                if (iout >= 2)
                {
                    iderr = 0;
                }
            }
            else
            {
                iderr = iwork[6];
                if (iout <= 1)
                {

                    Console.WriteLine(" ERROR ESTIMATION IN DENSE OUTPUT NOT POSSIBLE, WRONG IWORK(6)=", iwork[6]);

                    arret = true;
                }
            }
            /* -------- MUDIF */
            if (iwork[7] == 0)
            {
                mudif = 4;
            }
            else
            {
                mudif = iwork[7];
                if (mudif <= 0 || mudif >= 7)
                {

                    Console.WriteLine(" WRONG INPUT IWORK(7)=", iwork[7]);

                    arret = true;
                }
            }
            /* -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS */
            nrdens = iwork[8];
            if (nrdens < 0 || nrdens > n)
            {

                Console.WriteLine(" CURIOUS INPUT IWORK(8)=", iwork[8]);

                arret = true;
            }
            if (nrdens == n)
            {
                for (i = 0; i < nrdens; ++i)
                {
                    /* L17: */
                    iwork[i + 20] = i;
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

                    Console.WriteLine(" WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:", work[0]);

                    arret = true;
                }
            }
            /* -------- MAXIMAL STEP SIZE */
            if (work[2] == 0.0)
            {
                hmax = xend - x;
            }
            else
            {
                hmax = Math.Abs(work[2]);
            }
            /* -------- STEP SIZE REDUCTION FACTOR */
            if (work[3] == 0.0)
            {
                safe3 = .5;
            }
            else
            {
                safe3 = work[3];
                if (safe3 <= uround || safe3 >= 1.0)
                {

                    Console.WriteLine(" CURIOUS INPUT WORK(3)=", work[3]);

                    arret = true;
                }
            }
            /* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
            if (work[4] == 0.0)
            {
                fac1 = .02;
            }
            else
            {
                fac1 = work[4];
            }
            if (work[5] == 0.0)
            {
                fac2 = 4.0;
            }
            else
            {
                fac2 = work[5];
            }
            /* -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION */
            if (work[6] == 0.0)
            {
                fac3 = .8;
            }
            else
            {
                fac3 = work[6];
            }
            if (work[7] == 0.0)
            {
                fac4 = .9;
            }
            else
            {
                fac4 = work[7];
            }
            /* ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION */
            if (work[8] == 0.0)
            {
                safe1 = .65;
            }
            else
            {
                safe1 = work[8];
            }
            if (work[9] == 0.0)
            {
                safe2 = .94;
            }
            else
            {
                safe2 = work[9];
            }
            /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
            lfsafe = (km << 1) * km + km;
            var dy = new double[n];// 21;
            var yh1 = new double[n];// dy + n;
            var yh2 = new double[n];// yh1 + n;
            var dz = new double[n];// yh2 + n;
            var scal = new double[n];// dz + n;
            var t = new double[km * n];// scal + n;
            var fs = new double[lfsafe * nrdens];// t + km * n;
            var ys = new double[km * nrdens];// fs + lfsafe * nrdens;
            var hh = new double[km];// ys + km * nrdens;
            var w = new double[km];// hh + km;
            var a = new double[km];// w + km;
            var fac = new double[km << 1];// a + km;
            /* ------ TOTAL STORAGE REQUIREMENT ----------- */
            var co = new double[((km << 1) + 5) * nrdens];// fac + (km << 1);

            //istore = co + ((km << 1) + 5) * nrdens - 1;
            //if (istore > *lwork)
            //{
            //    Console.WriteLine(" INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", istore);
            //    arret = true;
            //}

            /* ------- ENTRY POINTS FOR int WORKSPACE ----- */
            var icom = new int[nrdens];// 21;
            var nj = new int[km];// icom + nrdens;

            /* --------- TOTAL REQUIREMENT --------------- */
            var ip = new int[km];// nj + km;

            //istore = ip + km;
            //if (istore > *liwork)
            //{
            //    Console.WriteLine(" INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", istore);
            //    arret = true;
            //}

            /* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
            if (arret)
            {
                idid = -1;
                return 0;
            }
            /* -------- CALL TO CORE INTEGRATOR ------------ */
            nrd = Math.Max(1, nrdens);
            /* Computing MAX */
            ncom = Math.Max(1, ((km << 1) + 5) * nrdens);
            odxcor_(n, fcn, x, y, xend, hmax, h, rtol, atol,
                itol, km, solout, iout, idid, nmax, uround, dy,
                yh1, yh2, dz, scal, fs, ys, t, hh, w, a, co, ncom, icom, nj, ip
                , nsequ, mstab, jstab, lfsafe, safe1, safe2, safe3, fac1,
                fac2, fac3, fac4, iderr, fac, mudif, nrd, rpar,
                 ipar, ref nfcn, ref nstep, ref naccpt, ref nrejct);

            iwork[17] = nfcn;
            iwork[18] = nstep;
            iwork[19] = naccpt;
            iwork[20] = nrejct;
            /* ----------- RETURN ----------- */
            return 0;
        } /* odex_ */




        /*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

        int odxcor_(int n, S_fp fcn, double x, double[]
            y, double xend, double hmax, double h, double[]
            rtol, double[] atol, int itol, int km, S_fp solout,
            int iout, int idid, int nmax, double uround,
            double[] dy, double[] yh1, double[] yh2, double[] dz,
            double[] scal, double[] fsafe, double[] ysafe, double[] t,
             double[] hh, double[] w, double[] a, double[] dens,
            int ncom, int[] icomp, int[] nj, int[] ipoint, int
            nsequ, int mstab, int jstab, int lfsafe, double
            safe1, double safe2, double safe3, double fac1,
            double fac2, double fac3, double fac4, int iderr,
            double[] errfac, int mudif, int nrd, double[] rpar,
            int[] ipar, ref int nfcn, ref int nstep, ref int naccpt,
            ref int nrejct)
        {
            /* Format strings */
            //static char fmt_979[] = "(  EXIT OF ODEX AT X= ,d14.7,    H= ,d14.7)";

            /* System generated locals */
            int i1, i2, i3, i4;
            double d1;

            /* Local variables */
            int i, j, k, l, kc = 0, kk, mu;
            double fac = 0;
            int kmi, kln;
            double err;
            int krn, ipt, kbeg, lbeg, lend, ncon;
            bool last;
            double prod;
            bool atov;
            double xold;
            int kopt;
            double errx;
            int njadd;
            double facnj;
            int irtrn = 0;
            double dblenj;
            bool reject;
            double factor, hoptde, errold, posneg;
            double errint;


            /* ---------------------------------------------------------- */
            /*     CORE INTEGRATOR FOR ODEX */
            /*     PARAMETERS SAME AS IN ODEX WITH WORKSPACE ADDED */
            /* ---------------------------------------------------------- */
            /*         DECLARATIONS */
            /* ---------------------------------------------------------- */
            /* --- DEFINE THE STEP SIZE SEQUENCE */
            /* Parameter adjustments */
            //--scal;
            //--dz;
            //--yh2;
            //--yh1;
            //--dy;
            //--y;
            //--rtol;
            //--atol;
            //--errfac;
            //--ipoint;
            //--nj;
            //--a;
            //--w;
            //--hh;
            //t_dim1 = km;
            //t_offset = 1 + km;
            //t -= t_offset;
            //--dens;
            //--icomp;
            //km = km;
            //ysafe_offset = 1 + km;
            //ysafe -= ysafe_offset;
            //fsafe_dim1 = *lfsafe;
            //fsafe_offset = 1 + lfsafe;
            //fsafe -= fsafe_offset;
            //--rpar;
            //--ipar;

            /* Function Body */
            if (nsequ == 1)
            {

                for (i = 0; i < km; ++i)
                {
                    /* L1: */
                    nj[i] = i << 1;
                }
            }
            if (nsequ == 2)
            {
                nj[0] = 2;

                for (i = 1; i < km; ++i)
                {
                    /* L2: */
                    nj[i] = (i << 2) - 4;
                }
            }
            if (nsequ == 3)
            {
                nj[0] = 2;
                nj[1] = 4;
                nj[2] = 6;

                for (i = 3; i < km; ++i)
                {
                    /* L11: */
                    nj[i] = nj[i - 2] << 1;
                }
            }
            if (nsequ == 4)
            {
                for (i = 0; i < km; ++i)
                {
                    /* L3: */
                    nj[i] = (i << 2) - 2;
                }
            }
            if (nsequ == 5)
            {
                for (i = 0; i < km; ++i)
                {
                    /* L6: */
                    nj[i] = i << 2;
                }
            }
            /* --- DEFINE THE A(I) FOR ORDER SELECTION */
            a[0] = nj[0] + 1.0;
            for (i = 1; i < km; ++i)
            {
                /* L4: */
                a[i] = a[i - 1] + nj[i];
            }
            /* --- INITIAL SCALING */
            for (i = 0; i < n; ++i)
            {
                if (itol == 0)
                {
                    scal[i] = atol[0] + rtol[0] * Math.Abs(y[i]);
                }
                else
                {
                    scal[i] = atol[i] + rtol[i] * Math.Abs(y[i]);
                }
                /* L8: */
            }
            /* --- INITIAL PREPARATIONS */
            posneg = d_sign(1.0, xend - x);
            k = Math.Max(2, Math.Min(km - 1, (int)(-Math.Log10(rtol[0] + 1e-40) * 0.6 + 1.5)));
            hmax = Math.Abs(hmax);
            h = Math.Max(Math.Abs(h), 1e-4);
            h = posneg * Math.Min(Math.Min(h, hmax), Math.Abs(xend - x) / 2.0);
            if (iout >= 1)
            {
                if (iout >= 2)
                {
                    ipoint[0] = 0;
                    for (i = 0; i < km; ++i)
                    {
                        njadd = (i << 2) - 2;
                        if (nj[i] > njadd)
                        {
                            ++njadd;
                        }
                        /* L5: */
                        ipoint[i + 1] = ipoint[i] + njadd;
                    }
                    i1 = km << 1;
                    for (mu = 0; mu < i1; ++mu)
                    {
                        errx = Math.Sqrt(mu / (mu + 4.0)) * 0.5;
                        /* Computing 2nd power */
                        d1 = mu + 4.0;
                        prod = 1.0 / (d1 * d1);
                        for (j = 0; j < mu; ++j)
                        {
                            /* L7: */
                            prod = prod * errx / j;
                        }
                        /* L9: */
                        errfac[mu] = prod;
                    }
                    ipt = 0;
                }
                irtrn = 0;
                xold = x;
                i1 = naccpt + 1;
                //solout(&i1, &xold, x, y, n, dens, &ncon, icomp, nrd, &irtrn);
                if (irtrn < 0)
                {
                    goto L130;
                }
            }
            err = 0.0;
            errold = 1e10;
            hoptde = posneg * hmax;
            w[0] = 0.0;
            reject = false;
            last = false;
            L10:
            atov = false;
            /* --- IS XEND REACHED IN THE NEXT STEP? */
            if (Math.Abs(xend - x) * 0.1 <= Math.Abs(x) * uround)
            {
                goto L110;
            }
            h = posneg * Math.Min(Math.Min(Math.Min(Math.Abs(h), Math.Abs(xend - x)), hmax), Math.Abs(hoptde));
            if ((x + h * 1.01 - xend) * posneg > 0.0)
            {
                h = xend - x;
                last = true;
            }
            if (nstep == 0 || iout != 2)
            {
                fcn(n, x, y, dz);
            }
            ++(nfcn);
            /* --- THE FIRST AND LAST STEP */
            if (nstep == 0 || last)
            {
                ipt = 0;
                ++(nstep);
                for (j = 0; j < k; ++j)
                {
                    kc = j;
                    midex_(j, x, y, h, hmax, n, (S_fp)fcn, dy, yh1, yh2, dz, t, nj, hh, w, ref err,
                         ref fac, a, safe1, uround, fac1, fac2, safe2, scal,
                        ref atov, safe3, ref reject, km, rtol, atol, itol,
                        mstab, jstab, errold, fsafe, lfsafe, iout,
                         ipt, ysafe, icomp, nrd, ref nfcn);
                    if (atov)
                    {
                        goto L10;
                    }
                    /* L20: */
                    if (j > 1 && err <= 1.0)
                    {
                        goto L60;
                    }
                }
                goto L55;
            }
            /* --- BASIC INTEGRATION STEP */
            L30:
            ipt = 0;
            ++(nstep);
            if (nstep >= nmax)
            {
                goto L120;
            }
            kc = k - 1;
            for (j = 0; j < kc; ++j)
            {
                midex_(j, x, y, h, hmax, n, (S_fp)fcn, dy, yh1, yh2
                    , dz, t, nj, hh, w, ref err, ref fac, a,
                     safe1, uround, fac1, fac2, safe2, scal, ref atov, safe3,
                    ref reject, km, rtol, atol, itol, mstab, jstab, errold,
                    fsafe, lfsafe, iout, ipt, ysafe
                    , icomp, nrd, ref nfcn);
                if (atov)
                {
                    goto L10;
                }
                /* L40: */
            }
            /* --- CONVERGENCE MONITOR */
            if (k == 2 || reject)
            {
                goto L50;
            }
            if (err <= 1.0)
            {
                goto L60;
            }
            /* Computing 2nd power */
            d1 = nj[k + 1] * nj[k] / 4.0;
            if (err > d1 * d1)
            {
                goto L100;
            }
            L50:
            midex_(k, x, y, h, hmax, n, fcn, dy, yh1, yh2, dz,
                 t, nj, hh, w, ref err, ref fac, a,
                safe1, uround, fac1, fac2, safe2, scal, ref atov, safe3, ref reject,
                 km, rtol, atol, itol, mstab, jstab, errold, fsafe,
                 lfsafe, iout, ipt, ysafe, icomp, nrd, ref nfcn);
            if (atov)
            {
                goto L10;
            }
            kc = k;
            if (err <= 1.0)
            {
                goto L60;
            }
            /* --- HOPE FOR CONVERGENCE IN LINE K+1 */
            L55:
            /* Computing 2nd power */
            d1 = nj[k + 1] / 2.0;
            if (err > d1 * d1)
            {
                goto L100;
            }
            kc = k + 1;
            midex_(kc, x, y, h, hmax, n, (S_fp)fcn, dy, yh1, yh2, dz,
                 t, nj, hh, w, ref err, ref fac, a,
                safe1, uround, fac1, fac2, safe2, scal, ref atov, safe3, ref reject,
                 km, rtol, atol, itol, mstab, jstab, errold, fsafe,
                 lfsafe, iout, ipt, ysafe, icomp, nrd, ref nfcn);
            if (atov)
            {
                goto L10;
            }
            if (err > 1.0)
            {
                goto L100;
            }
            /* --- STEP IS ACCEPTED */
            L60:
            xold = x;
            x += h;
            if (iout >= 2)
            {
                /* ---  KMIT = MU OF THE PAPER */
                conodx_1.kmit = (kc << 1) - mudif + 1;
                for (i = 0; i < nrd; ++i)
                {
                    /* L69: */
                    dens[i] = y[icomp[i]];
                }
                conodx_1.xoldd = xold;
                conodx_1.hhh = h;
                for (i = 0; i < nrd; ++i)
                {
                    /* L76: */
                    dens[nrd + i] = h * dz[icomp[i]];
                }
                kln = nrd << 1;
                for (i = 0; i < nrd; ++i)
                {
                    /* L176: */
                    dens[kln + i] = t[icomp[i] * km + 1];
                }
                /* --- COMPUTE SOLUTION AT MID-POINT ---- */
                for (j = 1; j < kc; ++j)
                {
                    dblenj = nj[j];
                    for (l = j; l >= 2; --l)
                    {
                        /* Computing 2nd power */
                        d1 = dblenj / nj[l - 1];
                        factor = d1 * d1 - 1.0;
                        for (i = 0; i < nrd; ++i)
                        {
                            ysafe[l - 1 + i * km] = ysafe[l + i * km] + (ysafe[l + i * km] - ysafe[l - 1 + i * km]) / factor;
                            /* L473: */
                        }
                    }
                }
                krn = nrd << 2;
                for (i = 0; i < nrd; ++i)
                {
                    /* L474: */
                    dens[krn + i] = ysafe[i * km + 1];
                }
                /* --- COMPUTE FIRST DERIVATIVE AT RIGHT END ---- */
                for (i = 0; i < n; ++i)
                {
                    /* L478: */
                    yh1[i] = t[i * km + 1];
                }
                fcn(n, x, yh1, yh2);
                krn = nrd * 3;
                for (i = 0; i < nrd; ++i)
                {
                    /* L274: */
                    dens[krn + i] = yh2[icomp[i]] * h;
                }
                /* --- THE LOOP --- */
                i2 = conodx_1.kmit;
                for (kmi = 0; kmi < i2; ++kmi)
                {
                    /* --- COMPUTE KMI-TH DERIVATIVE AT MID-POINT ---- */
                    kbeg = (kmi + 1) / 2;
                    i1 = kc;
                    for (kk = kbeg; kk <= i1; ++kk)
                    {
                        facnj = Math.Pow(nj[kk] / 2.0, kmi - 1); // TODO: pow_di
                        ipt = ipoint[kk + 1] - (kk << 1) + kmi;
                        for (i = 0; i < nrd; ++i)
                        {
                            /* L371: */
                            ysafe[kk + i * km] = fsafe[ipt + i * lfsafe] * facnj;
                        }
                        /* L375: */
                    }
                    for (j = kbeg; j < kc; ++j)
                    {
                        dblenj = nj[j];
                        i3 = kbeg + 1;
                        for (l = j; l >= i3; --l)
                        {
                            /* Computing 2nd power */
                            d1 = dblenj / nj[l - 1];
                            factor = d1 * d1 - 1.0;
                            for (i = 0; i < nrd; ++i)
                            {
                                ysafe[l - 1 + i * km] = ysafe[l + i * km] + (ysafe[l + i * km] - ysafe[l - 1 + i * km]) / factor;
                                /* L373: */
                            }
                        }
                    }
                    krn = (kmi + 4) * nrd;
                    for (i = 0; i < nrd; ++i)
                    {
                        /* L374: */
                        dens[krn + i] = ysafe[kbeg + i * km] * h;
                    }
                    if (kmi == conodx_1.kmit)
                    {
                        goto L180;
                    }
                    /* --- COMPUTE DIFFERENCES */
                    for (kk = (kmi + 2) / 2 - 1; kk < kc; ++kk)
                    {
                        lbeg = ipoint[kk + 1];
                        lend = ipoint[kk] + kmi + 1;
                        if (kmi == 1 && nsequ == 4)
                        {
                            lend += 2;
                        }
                        for (l = lbeg; l >= lend; l += -2)
                        {
                            for (i = 0; i < nrd; ++i)
                            {
                                /* L64: */
                                fsafe[l + i * lfsafe] -= fsafe[l - 2 + i * lfsafe];
                            }
                        }
                        if (kmi == 1 && nsequ == 4)
                        {
                            l = lend - 2;
                            for (i = 0; i < nrd; ++i)
                            {
                                /* L65: */
                                fsafe[l + i * lfsafe] -= dz[icomp[i]];
                            }
                        }
                        /* L66: */
                    }
                    /* --- COMPUTE DIFFERENCES */
                    i4 = kc;
                    for (kk = (kmi + 2) / 2; kk <= i4; ++kk)
                    {
                        lbeg = ipoint[kk + 1] - 1;
                        lend = ipoint[kk] + kmi + 2;
                        for (l = lbeg; l >= lend; l += -2)
                        {
                            for (i = 0; i < nrd; ++i)
                            {
                                /* L164: */
                                fsafe[l + i * lfsafe] -= fsafe[l - 2 + i * lfsafe];
                            }
                        }
                        /* L166: */
                    }
                    L180:
                    ;
                }
                interp_(nrd, dens, conodx_1.kmit);
                /* --- ESTIMATION OF INTERPOLATION ERROR */
                if (iderr == 0 && conodx_1.kmit >= 1)
                {
                    errint = 0.0;
                    for (i = 0; i < nrd; ++i)
                    {
                        /* L187: */
                        /* Computing 2nd power */
                        d1 = dens[(conodx_1.kmit + 4) * nrd + i] / scal[icomp[i]];
                        errint += d1 * d1;
                    }
                    errint = Math.Sqrt(errint / nrd) * errfac[conodx_1.kmit];
                    hoptde = h / Math.Max(Math.Pow(errint, 1.0 / (conodx_1.kmit + 4)), 0.01);
                    if (errint > 10.0)
                    {
                        h = hoptde;
                        x = xold;
                        ++(nrejct);
                        reject = true;
                        goto L10;
                    }
                }
                for (i = 0; i < n; ++i)
                {
                    /* L189: */
                    dz[i] = yh2[i];
                }
            }
            for (i = 0; i < n; ++i)
            {
                /* L70: */
                y[i] = t[i * km + 1];
            }
            ++(naccpt);
            if (iout >= 1)
            {
                i2 = naccpt + 1;
                //solout(&i2, &xold, x, y, n, dens, ncom, icomp, nrd, &irtrn);
                if (irtrn < 0)
                {
                    goto L130;
                }
            }
            /* --- COMPUTE OPTIMAL ORDER */
            if (kc == 2)
            {
                kopt = Math.Min(3, km - 1);
                if (reject)
                {
                    kopt = 2;
                }
                goto L80;
            }
            if (kc <= k)
            {
                kopt = kc;
                if (w[kc - 1] < w[kc] * fac3)
                {
                    kopt = kc - 1;
                }
                if (w[kc] < w[kc - 1] * fac4)
                {
                    kopt = Math.Min(kc + 1, km - 1);
                }
            }
            else
            {
                kopt = kc - 1;
                if (kc > 3 && w[kc - 2] < w[kc - 1] * fac3)
                {
                    kopt = kc - 2;
                }
                if (w[kc] < w[kopt] * fac4)
                {
                    kopt = Math.Min(kc, km - 1);
                }
            }
            /* --- AFTER A REJECTED STEP */
            L80:
            if (reject)
            {
                k = Math.Min(kopt, kc);
                h = posneg * Math.Min(Math.Abs(h), Math.Abs(hh[k]));
                reject = false;
                goto L10;
            }
            /* --- COMPUTE STEPSIZE FOR NEXT STEP */
            if (kopt <= kc)
            {
                h = hh[kopt];
            }
            else
            {
                if (kc < k && w[kc] < w[kc - 1] * fac4)
                {
                    h = hh[kc] * a[kopt + 1] / a[kc];
                }
                else
                {
                    h = hh[kc] * a[kopt] / a[kc];
                }
            }
            k = kopt;
            h = posneg * Math.Abs(h);
            goto L10;
            /* --- STEP IS REJECTED */
            L100:
            k = Math.Min(Math.Min(k, kc), km - 1);
            if (k > 2 && w[k - 1] < w[k] * fac3)
            {
                --k;
            }
            ++(nrejct);
            h = posneg * hh[k];
            reject = true;
            goto L30;
            /* --- SOLUTION EXIT */
            L110:
            idid = 1;
            return 0;
            /* --- INTERRUPTED BY SOLOUT */
            L130:
            idid = 2;
            return 0;
            /* --- FAIL EXIT */
            L120:

            Console.WriteLine((x));
            Console.WriteLine((h));

            idid = -1;
            return 0;
        } /* odxcor_ */


        int midex_(int j, double x, double[] y,
            double h, double hmax, int n, S_fp fcn, double[]
            dy, double[] yh1, double[] yh2, double[] dz, double[] t,
            int[] nj, double[] hh, double[] w, ref double err,
            ref double fac, double[] a, double safe1, double uround,
             double fac1, double fac2, double safe2, double[]
            scal, ref bool atov, double safe3, ref bool reject, int km,
            double[] rtol, double[] atol, int itol, int mstab,
            int jstab, double errold, double[] fsafe, int
            lfsafe, int iout, int ipt, double[] ysafe, int[]
            icomp, int nrd, ref int nfcn)
        {
            /* System generated locals */
            double d1;

            /* Local variables */
            int i, l, m;
            double hj;
            int mm;
            double ys, t1i, del1, del2, expo, quot;
            int njmid;
            double facmin, dblenj;

            /* --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE */
            /* --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATION */
            /* --- OF THE OPTIMAL STEPSIZE */
            /* Parameter adjustments */
            //--scal;
            //--dz;
            //--yh2;
            //--yh1;
            //--dy;
            //--y;
            //--a;
            //--w;
            //--hh;
            //--nj;
            //km = km;
            //t_offset = 1 + km;
            //t -= t_offset;
            //--rtol;
            //--atol;
            //--icomp;
            //km = km;
            //ysafe_offset = 1 + km;
            //ysafe -= ysafe_offset;
            //fsafe_dim1 = *lfsafe;
            //fsafe_offset = 1 + lfsafe;
            //fsafe -= fsafe_offset;
            //--rpar;
            //--ipar;

            /* Function Body */
            hj = h / nj[j];
            /* --- EULER STARTING STEP */
            for (i = 0; i < n; ++i)
            {
                yh1[i] = y[i];
                /* L30: */
                yh2[i] = y[i] + hj * dz[i];
            }
            /* --- EXPLICIT MIDPOINT RULE */
            m = nj[j] - 1;
            njmid = nj[j] / 2;
            for (mm = 0; mm < m; ++mm)
            {
                if (iout >= 2 && mm == njmid)
                {
                    for (i = 0; i < nrd; ++i)
                    {
                        /* L31: */
                        ysafe[j + i * km] = yh2[icomp[i]];
                    }
                }
                fcn(n, x + hj * mm, yh2, dy);
                if (iout >= 2 && Math.Abs(mm - njmid) <= (j << 1) - 1)
                {
                    ++(ipt);
                    for (i = 0; i < nrd; ++i)
                    {
                        /* L32: */
                        fsafe[ipt + i * lfsafe] = dy[icomp[i]];
                    }
                }
                for (i = 0; i < n; ++i)
                {
                    ys = yh1[i];
                    yh1[i] = yh2[i];
                    /* L34: */
                    yh2[i] = ys + hj * 2.0 * dy[i];
                }
                if (mm <= mstab && j <= jstab)
                {
                    /* --- STABILITY CHECK */
                    del1 = 0.0;
                    for (i = 0; i < n; ++i)
                    {
                        /* L21: */
                        /* Computing 2nd power */
                        d1 = dz[i] / scal[i];
                        del1 += d1 * d1;
                    }
                    del2 = 0.0;
                    for (i = 0; i < n; ++i)
                    {
                        /* L26: */
                        /* Computing 2nd power */
                        d1 = (dy[i] - dz[i]) / scal[i];
                        del2 += d1 * d1;
                    }
                    quot = del2 / Math.Max(uround, del1);
                    if (quot > 4.0)
                    {
                        ++(nfcn);
                        goto L79;
                    }
                }
                /* L35: */
            }
            /* --- FINAL SMOOTHING STEP */
            fcn(n, x + h, yh2, dy);
            if (iout >= 2 && njmid <= (j << 1) - 1)
            {
                ++(ipt);
                for (i = 0; i < nrd; ++i)
                {
                    /* L39: */
                    fsafe[ipt + i * lfsafe] = dy[icomp[i]];
                }
            }
            for (i = 0; i < n; ++i)
            {
                /* L40: */
                t[j + i * km] = (yh1[i] + yh2[i] + hj * dy[i]) / 2.0;
            }
            nfcn += nj[j];
            /* --- POLYNOMIAL EXTRAPOLATION */
            if (j == 1)
            {
                return 0;
            }
            dblenj = nj[j];
            for (l = j; l >= 2; --l)
            {
                /* Computing 2nd power */
                d1 = dblenj / nj[l - 1];
                fac = d1 * d1 - 1.0;
                for (i = 0; i < n; ++i)
                {
                    t[l - 1 + i * km] = t[l + i * km] + (t[l + i * km] - t[l - 1 + i * km]) / fac;
                    /* L60: */
                }
            }
            err = 0.0;
            /* --- SCALING */
            for (i = 0; i < n; ++i)
            {
                t1i = Math.Max(Math.Abs(y[i]), Math.Abs(t[i * km + 1]));
                if (itol == 0)
                {
                    scal[i] = atol[0] + rtol[0] * t1i;
                }
                else
                {
                    scal[i] = atol[i] + rtol[i] * t1i;
                }
                /* L65: */
                /* Computing 2nd power */
                d1 = (t[i * km + 1] - t[i * km + 2]) / scal[i];
                err += d1 * d1;
            }
            err = Math.Sqrt(err / n);
            if (err * uround >= 1.0)
            {
                goto L79;
            }
            if (j > 2 && err >= errold)
            {
                goto L79;
            }
            errold = Math.Max(err * 4, 1.0);
            /* --- COMPUTE OPTIMAL STEPSIZES */
            expo = 1.0 / ((j << 1) - 1);
            facmin = Math.Pow(fac1, expo);
            fac = Math.Min(fac2 / facmin, Math.Max(facmin, Math.Pow(err / safe1, expo) / safe2));
            fac = 1.0 / fac;
            hh[j] = Math.Min(Math.Abs(h) * fac, hmax);
            w[j] = a[j] / hh[j];
            return 0;
            L79:
            atov = true;
            h *= safe3;
            reject = true;
            return 0;
        } /* midex_ */


        int interp_(int n, double[] y, int imit)
        {
            /* Local variables */
            double[] a = new double[31];
            int i;
            double y0, y1;
            int im;
            double ph0, ph1, ph2, ph3, yp0, yp1, fac1, fac2, aspl, bspl, ydiff;

            /* --- COMPUTES THE COEFFICIENTS OF THE INTERPOLATION FORMULA */
            /* --- BEGIN WITH HERMITE INTERPOLATION */
            /* Parameter adjustments */
            //--y;

            /* Function Body */
            for (i = 0; i < n; ++i)
            {
                y0 = y[i];
                y1 = y[(n << 1) + i];
                yp0 = y[n + i];
                yp1 = y[n * 3 + i];
                ydiff = y1 - y0;
                aspl = -yp1 + ydiff;
                bspl = yp0 - ydiff;
                y[n + i] = ydiff;
                y[(n << 1) + i] = aspl;
                y[n * 3 + i] = bspl;
                if (imit < 0)
                {
                    goto L100;
                }
                /* --- COMPUTE THE DERIVATIVES OF HERMITE AT MIDPOINT */
                ph0 = (y0 + y1) * 0.5 + (aspl + bspl) * 0.125;
                ph1 = ydiff + (aspl - bspl) * 0.25;
                ph2 = -(yp0 - yp1);
                ph3 = (bspl - aspl) * 6.0;
                /* --- COMPUTE THE FURTHER COEFFICIENTS */
                if (imit < 1)
                {
                    goto L20;
                }
                a[0] = (y[n * 5 + i] - ph1) * 16.0;
                if (imit < 3)
                {
                    goto L20;
                }
                a[3] = (y[n * 7 + i] - ph3 + a[0] * 3) * 16.0;
                if (imit < 5)
                {
                    goto L20;
                }

                for (im = 4; im < imit; im += 2)
                {
                    fac1 = im * (im - 1) / 2.0;
                    fac2 = fac1 * (im - 2) * (im - 3) * 2.0;
                    /* L10: */
                    a[im] = (y[(im + 4) * n + i] + fac1 * a[im - 2] - fac2 * a[im - 4]) * 16.0;
                }
                L20:
                a[0] = (y[(n << 2) + i] - ph0) * 16.0;
                if (imit < 2)
                {
                    goto L60;
                }
                a[2] = (y[n * 6 + i] - ph2 + a[0]) * 16.0;
                if (imit < 4)
                {
                    goto L60;
                }

                for (im = 3; im < imit; im += 2)
                {
                    fac1 = im * (im - 1) / 2.0;
                    fac2 = (double)(im * (im - 1) * (im - 2) * (im - 3));
                    /* L30: */
                    a[im] = (y[n * (im + 4) + i] + a[im - 2] * fac1 - a[im - 4] * fac2) * 16.0;
                }
                L60:
                for (im = -1; im < imit; ++im)
                {
                    /* L70: */
                    y[n * (im + 4) + i] = a[im];
                }
                L100:
                ;
            }
            return 0;
        } /* interp_ */


        double contex_(int ii, double x, double[] y, int ncon, int[] icomp, int n)
        {
            /* System generated locals */
            double ret_val, d1;

            /* Local variables */
            int i, j, im;
            double theta, theta1, thetah, phthet;

            /* ---------------------------------------------------------- */
            /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONECTION */
            /*     WITH THE OUTPUT-SUBROUTINE FOR ODEX. IT PROVIDES AN */
            /*     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.0 */
            /* ---------------------------------------------------------- */
            /* ----- COMPUTE PLACE OF II-TH COMPONENT */
            /* Parameter adjustments */
            //--y;
            //--icomp;

            /* Function Body */
            i = 0;
            for (j = 0; j < n; ++j)
            {
                if (icomp[j] == ii)
                {
                    i = j;
                }
                /* L5: */
            }
            if (i == 0)
            {

                Console.WriteLine(" NO DENSE OUTPUT AVAILABLE FOR COMP.", ii);
                return 0;
            }

            /* ----- COMPUTE THE INTERPOLATED VALUE */
            theta = (x - conodx_2.xold) / conodx_2.h;
            theta1 = 1.0 - theta;
            phthet = y[i] + theta * (y[n + i] + theta1 * (y[(n << 1) + i] * theta + y[n * 3 + i] * theta1));
            if (conodx_2.imit < 0)
            {
                ret_val = phthet;
                return ret_val;
            }
            thetah = theta - .5;
            ret_val = y[n * (conodx_2.imit + 4) + i];
            for (im = conodx_2.imit; im >= 1; --im)
            {
                /* L70: */
                ret_val = y[n * (im + 3) + i] + ret_val * thetah / im;
            }
            /* Computing 2nd power */
            d1 = theta * theta1;
            ret_val = phthet + d1 * d1 * ret_val;
            return ret_val;
        } /* contex_ */
    }
}