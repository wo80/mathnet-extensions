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
    /// E. Hairer, S.P. Norsett and G. Wanner
    /// Solving Ordinary Differential Equations I. Nonstiff Problems (2nd edition)
    /// Springer-Verlag (1993)
    /// </summary>
    public class BulirschStoer
    {
        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

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

        double rtol, atol;

        public int nstep, naccpt, nrejct; // TODO: remove

        public int odex_(int n, S_fp fcn, double x, double[] y, double xend, double h, double rtol, double atol, int iout)
        {
            int km;
            double fac1, fac2, fac3, fac4;
            double hmax;
            int nmax;
            double safe1, safe2, safe3;
            int jstab, mudif, iderr = 0, mstab;
            int nsequ, lfsafe;
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
             *     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES.
             *                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
             *                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
             *
             *     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
             *                    IOUT=0: SUBROUTINE IS NEVER CALLED
             *                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
             *                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
             *
             * ----------------------------------------------------------------------
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
             */

            this.rtol = rtol;
            this.atol = atol;

            nstep = 0;
            naccpt = 0;
            nrejct = 0;

            // NMAX , THE MAXIMAL NUMBER OF STEPS
            nmax = 10000;

            // KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION TABLE.
            // THE DEFAULT VALUE(FOR IWORK(2)= 0) IS 9.
            // IF IWORK(2).NE.0 THEN IWORK(2) SHOULD BE .GE.3.
            km = 9;

            if (km <= 2)
            {
                Console.WriteLine(" CURIOUS INPUT KM=", km);
                return -1;
            }

            // NSEQU     CHOICE OF STEP SIZE SEQUENCE
            // SWITCH FOR THE STEP SIZE SEQUENCE(EVEN NUMBERS ONLY)
            // IF NSEQU = 1 THEN 2,4,6,8,10,12,14,16,...
            // IF NSEQU = 2 THEN 2,4,8,12,16,20,24,28,...
            // IF NSEQU = 3 THEN 2,4,6,8,12,16,24,32,...
            // IF NSEQU = 4 THEN 2,6,10,14,18,22,26,30,...
            // IF NSEQU = 5 THEN 4,8,12,16,20,24,28,32,...
            // THE DEFAULT VALUE IS NSEQU = 1 IF IOUT <= 1;
            // THE DEFAULT VALUE IS NSEQU = 4 IF IOUT >= 2.
            nsequ = iout <= 1 ? 1 : 4;

            if (nsequ <= 0 || nsequ >= 6)
            {
                Console.WriteLine(" CURIOUS INPUT NSEQU=", nsequ);
                return -1;
            }

            if (nsequ <= 3 && iout >= 2)
            {
                Console.WriteLine(" NSEQU NOT COMPATIBLE WITH IOUT");
                return -1;
            }
            
            // MSTAB     PARAMETER FOR STABILITY CHECK
            // STABILITY CHECK IS ACTIVATED AT MOST IWORK(4) TIMES IN
            // ONE LINE OF THE EXTRAP. TABLE, DEFAULT IWORK(4)= 1.
            mstab = 1;

            // JSTAB     PARAMETER FOR STABILITY CHECK
            // STABILITY CHECK IS ACTIVATED ONLY IN THE LINES
            // 1 TO IWORK(5) OF THE EXTRAP.TABLE, DEFAULT IWORK(5)= 1.
            jstab = 2;

            // IDERR  PARAMETER FOR ERROR ESTIMATION IN DENSE OUTPUT
            // IF IWORK(6)= 0  ERROR ESTIMATOR IN THE DENSE
            // OUTPUT FORMULA IS ACTIVATED.IT CAN BE SUPPRESSED
            // BY PUTTING IWORK(6)= 1.
            // DEFAULT IWORK(6) = 0 (IF IOUT.GE.2).
            iderr = (iout <= 1) ? 1 : 0;

            if (iderr == 1 && iout <= 1)
            {
                Console.WriteLine(" ERROR ESTIMATION IN DENSE OUTPUT NOT POSSIBLE, WRONG IDERR=", iderr);
                return -1;
            }

            // MUDIF      DETERMINES THE DEGREE OF INTERPOLATION FORMULA
            // MU = 2 * KAPPA - MUDIF + 1
            // MUDIF SHOULD LIE BETWEEN 1 AND 6
            // DEFAULT MUDIF = 4.
            mudif = 4;

            if (mudif <= 0 || mudif >= 7)
            {
                Console.WriteLine(" WRONG INPUT MUDIF=", mudif);
                return -1;
            }

            // UROUND   SMALLEST NUMBER SATISFYING 1.0+UROUND > 1.0
            uround = 2.3e-16;

            // MAXIMAL STEP SIZE
            hmax = xend - x;

            // STEP SIZE REDUCTION FACTOR
            // STEP SIZE IS REDUCED BY FACTOR WORK(3), IF THE
            // STABILITY CHECK IS NEGATIVE, DEFAULT 0.5.
            safe3 = 0.5;

            if (safe3 <= uround || safe3 >= 1.0)
            {
                Console.WriteLine(" CURIOUS INPUT safe3=", safe3);
                return -1;
            }

            // FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
            // THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
            // CHOSEN SUBJECT TO THE RESTRICTION
            //    FACMIN / WORK(5) <= HNEW(J) / HOLD <= 1 / FACMIN
            // WHERE FACMIN = WORK(4) * *(1 / (2 * J - 1))
            // DEFAULT VALUES: WORK(4) = 0.02, WORK(5) = 4.0
            fac1 = 0.02;
            fac2 = 4.0;

            // FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION
            // STEP SIZE IS DECREASED IF    W(K - 1) <= W(K) * WORK(6)
            // STEP SIZE IS INCREASED IF    W(K) <= W(K - 1) * WORK(7)
            // DEFAULT VALUES: WORK(6) = 0.8, WORK(7) = 0.9
            fac3 = 0.8;
            fac4 = 0.9;

            // SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION
            //              SAFETY FACTORS FOR STEP CONTROL ALGORITHM
            //    HNEW = H * WORK(9) * (WORK(8) * TOL / ERR) * *(1 / (J - 1))
            // DEFAULT VALUES: WORK(8) = 0.65,
            //    WORK(9) = 0.94  IF "HOPE FOR CONVERGENCE"
            //    WORK(9) = 0.90  IF "NO HOPE FOR CONVERGENCE"
            safe1 = 0.65;
            safe2 = 0.94;

            // PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK
            lfsafe = (2 * km) * km + km;
            var dy = new double[n];
            var yh1 = new double[n];
            var yh2 = new double[n];
            var dz = new double[n];
            var scal = new double[n];
            var t = new double[km * n];
            var hh = new double[km];
            var w = new double[km];
            var a = new double[km];
            var fac = new double[2 * km];
            var fs = new double[lfsafe * n]; // For interp
            var ys = new double[km * n]; // For interp
            var co = new double[((2 * km) + 5) * n]; // For interp

            // ENTRY POINTS FOR int WORKSPACE
            var nj = new int[km];
            var ip = new int[km + 1];

            // CALL TO CORE INTEGRATOR
            odxcor_(n, fcn, x, y, xend, hmax, h, km, iout, nmax, uround, dy,
                yh1, yh2, dz, scal, fs, ys, t, hh, w, a, co, nj, ip,
                nsequ, mstab, jstab, lfsafe, safe1, safe2, safe3, fac1,
                fac2, fac3, fac4, iderr, fac, mudif);

            return 0;
        }

        // ... AND HERE IS THE CORE INTEGRATOR

        int odxcor_(int n, S_fp fcn, double x, double[]
            y, double xend, double hmax, double h, int km,
            int iout, int nmax, double uround,
            double[] dy, double[] yh1, double[] yh2, double[] dz,
            double[] scal, double[] fsafe, double[] ysafe, double[] t,
            double[] hh, double[] w, double[] a, double[] dens,
            int[] nj, int[] ipoint, int nsequ,
            int mstab, int jstab, int lfsafe,
            double safe1, double safe2, double safe3, double fac1,
            double fac2, double fac3, double fac4, int iderr,
            double[] errfac, int mudif)
        {
            int i1, i2;
            double d1;

            int i, j, k, kc = 0, mu;
            double fac = 0;
            double err;
            int ipt;
            bool last;
            double prod;
            bool atov;
            double xold;
            int kopt;
            double errx;
            int njadd;
            int irtrn = 0;
            bool reject;
            double hoptde, errold, posneg;

            // CORE INTEGRATOR FOR ODEX
            // PARAMETERS SAME AS IN ODEX WITH WORKSPACE ADDED

            // DEFINE THE STEP SIZE SEQUENCE
            if (nsequ == 1)
            {
                for (i = 0; i < km; ++i)
                {
                    nj[i] = 2 * (i + 1);
                }
            }
            if (nsequ == 2)
            {
                nj[0] = 2;

                for (i = 1; i < km; ++i)
                {
                    nj[i] = 4 * (i + 1) - 4;
                }
            }
            if (nsequ == 3)
            {
                nj[0] = 2;
                nj[1] = 4;
                nj[2] = 6;

                for (i = 3; i < km; ++i)
                {
                    nj[i] = 2 * nj[i - 2];
                }
            }
            if (nsequ == 4)
            {
                for (i = 0; i < km; ++i)
                {
                    nj[i] = 4 * (i + 1) - 2;
                }
            }
            if (nsequ == 5)
            {
                for (i = 0; i < km; ++i)
                {
                    nj[i] = 4 * (i + 1);
                }
            }
            // DEFINE THE A(I) FOR ORDER SELECTION
            a[0] = nj[0] + 1.0;
            for (i = 1; i < km; ++i)
            {
                a[i] = a[i - 1] + nj[i];
            }
            // INITIAL SCALING
            for (i = 0; i < n; ++i)
            {
                scal[i] = atol + rtol * Math.Abs(y[i]);
            }
            // INITIAL PREPARATIONS
            posneg = d_sign(1.0, xend - x);
            k = Math.Max(2, Math.Min(km - 1, (int)(-Math.Log10(rtol + 1e-40) * 0.6 + 1.5)));
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
                        njadd = ((i + 1) << 2) - 2;
                        if (nj[i] > njadd)
                        {
                            ++njadd;
                        }
                        ipoint[i + 1] = ipoint[i] + njadd;
                    }
                    i1 = 2 * km + 1;
                    for (mu = 1; mu < i1; ++mu)
                    {
                        errx = Math.Sqrt(mu / (mu + 4.0)) * 0.5;
                        // Computing 2nd power
                        d1 = mu + 4.0;
                        prod = 1.0 / (d1 * d1);
                        for (j = 1; j <= mu; ++j)
                        {
                            prod = prod * errx / j;
                        }
                        errfac[mu - 1] = prod;
                    }
                    ipt = 0;
                }
                irtrn = 0;
                xold = x;
                i1 = naccpt + 1;
                //solout(&i1, &xold, x, y, n, dens, &ncon, icomp, nrd, &irtrn);
                if (irtrn < 0)
                {
                    // INTERRUPTED BY SOLOUT
                    return 2;
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
            // IS XEND REACHED IN THE NEXT STEP?
            if (Math.Abs(xend - x) * 0.1 <= Math.Abs(x) * uround)
            {
                // SOLUTION EXIT
                return 1;
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

            // THE FIRST AND LAST STEP
            if (nstep == 0 || last)
            {
                ipt = 0;
                ++(nstep);
                for (j = 0; j < k; ++j)
                {
                    kc = j + 1;
                    midex_(j, x, y, ref h, hmax, n, fcn, dy, yh1, yh2, dz, t, nj, hh, w, ref err,
                        ref fac, a, safe1, uround, fac1, fac2, safe2, scal,
                        ref atov, safe3, ref reject, km,
                        mstab, jstab, ref errold, fsafe, lfsafe, iout,
                        ref ipt, ysafe);
                    if (atov)
                    {
                        goto L10;
                    }
                    if (j > 0 && err <= 1.0)
                    {
                        goto L60;
                    }
                }
                goto L55;
            }
            // BASIC INTEGRATION STEP
            L30:
            ipt = 0;
            ++(nstep);
            if (nstep >= nmax)
            {
                // FAIL EXIT
                //Console.WriteLine("(  EXIT OF ODEX AT X={0}, H={1})", x, h);
                return -1;
            }
            kc = k - 1;
            for (j = 0; j < kc; ++j)
            {
                midex_(j, x, y, ref h, hmax, n, (S_fp)fcn, dy, yh1, yh2
                    , dz, t, nj, hh, w, ref err, ref fac, a,
                     safe1, uround, fac1, fac2, safe2, scal, ref atov, safe3,
                    ref reject, km, mstab, jstab, ref errold,
                    fsafe, lfsafe, iout, ref ipt, ysafe);
                if (atov)
                {
                    goto L10;
                }
            }
            // CONVERGENCE MONITOR
            if (k == 2 || reject)
            {
                goto L50;
            }
            if (err <= 1.0)
            {
                goto L60;
            }
            // Computing 2nd power
            d1 = nj[k] * nj[k - 1] / 4.0;
            if (err > d1 * d1)
            {
                goto L100;
            }
            L50:
            midex_(k - 1, x, y, ref h, hmax, n, fcn, dy, yh1, yh2, dz,
                 t, nj, hh, w, ref err, ref fac, a,
                safe1, uround, fac1, fac2, safe2, scal, ref atov, safe3, ref reject,
                 km, mstab, jstab, ref errold, fsafe,
                 lfsafe, iout, ref ipt, ysafe);
            if (atov)
            {
                goto L10;
            }
            kc = k;
            if (err <= 1.0)
            {
                goto L60;
            }
            // HOPE FOR CONVERGENCE IN LINE K+1
            L55:
            // Computing 2nd power
            d1 = nj[k] / 2.0;
            if (err > d1 * d1)
            {
                goto L100;
            }
            kc = k + 1;
            midex_(kc - 1, x, y, ref h, hmax, n, (S_fp)fcn, dy, yh1, yh2, dz,
                 t, nj, hh, w, ref err, ref fac, a,
                safe1, uround, fac1, fac2, safe2, scal, ref atov, safe3, ref reject,
                 km, mstab, jstab, ref errold, fsafe,
                 lfsafe, iout, ref ipt, ysafe);
            if (atov)
            {
                goto L10;
            }
            if (err > 1.0)
            {
                goto L100;
            }
            // STEP IS ACCEPTED
            L60:
            xold = x;
            x += h;
            if (iout >= 2)
            {
                PrepareInterpolation(kc, n, km, nsequ, mudif, fcn, h, x, xold,
                    y, dens, t, dz, nj, yh1, yh2, ysafe, fsafe, lfsafe, ipoint);

                interp_(n, dens, conodx_1.kmit);

                // ESTIMATION OF INTERPOLATION ERROR
                if (iderr == 0 && conodx_1.kmit >= 1)
                {
                    double errint = 0.0;
                    for (i = 0; i < n; ++i)
                    {
                        // Computing 2nd power
                        d1 = dens[(conodx_1.kmit + 4) * n + i] / scal[i];
                        errint += d1 * d1;
                    }
                    errint = Math.Sqrt(errint / n) * errfac[conodx_1.kmit - 1];
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
                    dz[i] = yh2[i];
                }
            }
            for (i = 0; i < n; ++i)
            {
                y[i] = t[i * km];
            }
            ++(naccpt);
            if (iout >= 1)
            {
                i2 = naccpt + 1;
                //solout(&i2, &xold, x, y, n, dens, ncom, icomp, nrd, &irtrn);
                if (irtrn < 0)
                {
                    // INTERRUPTED BY SOLOUT
                    return 2;
                }
            }
            // COMPUTE OPTIMAL ORDER
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
                if (w[kc - 2] < w[kc - 1] * fac3)
                {
                    kopt = kc - 1;
                }
                if (w[kc - 1] < w[kc - 2] * fac4)
                {
                    kopt = Math.Min(kc + 1, km - 1);
                }
            }
            else
            {
                kopt = kc - 1;
                if (kc > 3 && w[kc - 3] < w[kc - 2] * fac3)
                {
                    kopt = kc - 2;
                }
                if (w[kc - 1] < w[kopt - 1] * fac4)
                {
                    kopt = Math.Min(kc, km - 1);
                }
            }
            // AFTER A REJECTED STEP
            L80:
            if (reject)
            {
                k = Math.Min(kopt, kc);
                h = posneg * Math.Min(Math.Abs(h), Math.Abs(hh[k - 1]));
                reject = false;
                goto L10;
            }
            // COMPUTE STEPSIZE FOR NEXT STEP
            if (kopt <= kc)
            {
                h = hh[kopt - 1];
            }
            else
            {
                if (kc < k && w[kc - 1] < w[kc - 2] * fac4)
                {
                    h = hh[kc - 1] * a[kopt] / a[kc - 1];
                }
                else
                {
                    h = hh[kc - 1] * a[kopt - 1] / a[kc - 1];
                }
            }
            k = kopt;
            h = posneg * Math.Abs(h);
            goto L10;
            // STEP IS REJECTED
            L100:
            k = Math.Min(Math.Min(k, kc), km - 1);
            if (k > 2 && w[k - 2] < w[k - 1] * fac3)
            {
                --k;
            }
            ++(nrejct);
            h = posneg * hh[k - 1];
            reject = true;
            goto L30;
        }

        private void PrepareInterpolation(int kc, int n, int km, int nsequ, int mudif, S_fp fcn, double h, double x, double xold,
            double[] y, double[] dens, double[] t, double[] dz, int[] nj, double[] yh1, double[] yh2,
            double[] ysafe, double[] fsafe, int lfsafe, int[] ipoint)
        {
            int i, j, k, l, i2;
            double d1;

            // KMIT = MU OF THE PAPER
            conodx_1.kmit = 2 * kc - mudif + 1;
            for (i = 0; i < n; ++i)
            {
                dens[i] = y[i];
            }

            conodx_1.xoldd = xold;
            conodx_1.hhh = h;
            for (i = 0; i < n; ++i)
            {
                dens[n + i] = h * dz[i];
            }

            int kln = 2 * n;
            for (i = 0; i < n; ++i)
            {
                dens[kln + i] = t[i * km];
            }

            // COMPUTE SOLUTION AT MID-POINT
            for (j = 1; j < kc; ++j)
            {
                double dblenj = nj[j];
                for (l = j; l >= 1; --l)
                {
                    // Computing 2nd power
                    d1 = dblenj / nj[l - 1];
                    double factor = d1 * d1 - 1.0;
                    for (i = 0; i < n; ++i)
                    {
                        ysafe[l - 1 + i * km] = ysafe[l + i * km] + (ysafe[l + i * km] - ysafe[l - 1 + i * km]) / factor;
                    }
                }
            }
            int krn = n << 2;
            for (i = 0; i < n; ++i)
            {
                dens[krn + i] = ysafe[i * km];
            }

            // COMPUTE FIRST DERIVATIVE AT RIGHT END
            for (i = 0; i < n; ++i)
            {
                yh1[i] = t[i * km];
            }
            fcn(n, x, yh1, yh2);
            krn = 3 * n;
            for (i = 0; i < n; ++i)
            {
                dens[krn + i] = yh2[i] * h;
            }

            // THE LOOP
            i2 = conodx_1.kmit;
            for (int kmi = 1; kmi <= i2; ++kmi)
            {
                // COMPUTE KMI-TH DERIVATIVE AT MID-POINT
                int kbeg = (kmi + 1) / 2;
                for (int kk = kbeg - 1; kk < kc; ++kk)
                {
                    double facnj = Math.Pow(nj[kk] / 2.0, kmi - 1); // TODO: pow_di
                    int ipt = ipoint[kk + 1] - 2 * (kk + 1) + kmi - 1;
                    for (i = 0; i < n; ++i)
                    {
                        ysafe[kk + i * km] = fsafe[ipt + i * lfsafe] * facnj;
                    }
                }

                for (j = kbeg; j < kc; ++j)
                {
                    double dblenj = nj[j];
                    for (l = j; l >= kbeg; --l)
                    {
                        // Computing 2nd power
                        d1 = dblenj / nj[l - 1];
                        double factor = d1 * d1 - 1.0;
                        for (i = 0; i < n; ++i)
                        {
                            ysafe[l - 1 + i * km] = ysafe[l + i * km] + (ysafe[l + i * km] - ysafe[l - 1 + i * km]) / factor;
                        }
                    }
                }
                krn = (kmi + 4) * n;
                for (i = 0; i < n; ++i)
                {
                    dens[krn + i] = ysafe[kbeg - 1 + i * km] * h;
                }

                if (kmi == conodx_1.kmit)
                {
                    continue;
                }

                // COMPUTE DIFFERENCES
                for (int kk = (kmi + 2) / 2 - 1; kk < kc; ++kk)
                {
                    int lbeg = ipoint[kk + 1] - 1;
                    int lend = ipoint[kk] + kmi;
                    if (kmi == 1 && nsequ == 4)
                    {
                        lend += 2;
                    }
                    for (l = lbeg; l >= lend; l += -2)
                    {
                        for (i = 0; i < n; ++i)
                        {
                            fsafe[l + i * lfsafe] -= fsafe[l - 2 + i * lfsafe];
                        }
                    }
                    if (kmi == 1 && nsequ == 4)
                    {
                        l = lend - 2;
                        for (i = 0; i < n; ++i)
                        {
                            fsafe[l + i * lfsafe] -= dz[i];
                        }
                    }
                }

                // COMPUTE DIFFERENCES
                for (int kk = (kmi + 2) / 2 - 1; kk < kc; ++kk)
                {
                    int lbeg = ipoint[kk + 1] - 2;
                    int lend = ipoint[kk] + kmi + 1;
                    for (l = lbeg; l >= lend; l += -2)
                    {
                        for (i = 0; i < n; ++i)
                        {
                            fsafe[l + i * lfsafe] -= fsafe[l - 2 + i * lfsafe];
                        }
                    }
                }
            }
        }

        int midex_(int j, double x, double[] y,
            ref double h, double hmax, int n, S_fp fcn, double[] dy,
            double[] yh1, double[] yh2, double[] dz, double[] t,
            int[] nj, double[] hh, double[] w, ref double err,
            ref double fac, double[] a, double safe1, double uround,
            double fac1, double fac2, double safe2, double[] scal,
            ref bool atov, double safe3, ref bool reject, int km,
            int mstab, int jstab, ref double errold, double[] fsafe,
            int lfsafe, int iout, ref int ipt, double[] ysafe)
        {
            double d1;

            int i, l, m;
            double hj;
            int mm;
            double ys, t1i, del1, del2, expo, quot;
            int njmid;
            double facmin, dblenj;

            // THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
            // EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATION
            // OF THE OPTIMAL STEPSIZE

            hj = h / nj[j];
            // EULER STARTING STEP
            for (i = 0; i < n; ++i)
            {
                yh1[i] = y[i];
                yh2[i] = y[i] + hj * dz[i];
            }
            // EXPLICIT MIDPOINT RULE
            m = nj[j];
            njmid = nj[j] / 2;
            for (mm = 1; mm < m; ++mm)
            {
                if (iout >= 2 && mm == njmid)
                {
                    for (i = 0; i < n; ++i)
                    {
                        ysafe[j + i * km] = yh2[i];
                    }
                }
                fcn(n, x + hj * mm, yh2, dy);
                if (iout >= 2 && Math.Abs(mm - njmid) <= 2 * j + 1)
                {
                    for (i = 0; i < n; ++i)
                    {
                        fsafe[ipt + i * lfsafe] = dy[i];
                    }
                    ++(ipt);
                }
                for (i = 0; i < n; ++i)
                {
                    ys = yh1[i];
                    yh1[i] = yh2[i];
                    yh2[i] = ys + hj * 2.0 * dy[i];
                }
                if (mm <= mstab && j + 1 <= jstab)
                {
                    // STABILITY CHECK
                    del1 = 0.0;
                    for (i = 0; i < n; ++i)
                    {
                        // Computing 2nd power
                        d1 = dz[i] / scal[i];
                        del1 += d1 * d1;
                    }
                    del2 = 0.0;
                    for (i = 0; i < n; ++i)
                    {
                        // Computing 2nd power
                        d1 = (dy[i] - dz[i]) / scal[i];
                        del2 += d1 * d1;
                    }
                    quot = del2 / Math.Max(uround, del1);
                    if (quot > 4.0)
                    {
                        goto L79;
                    }
                }
            }
            // FINAL SMOOTHING STEP
            fcn(n, x + h, yh2, dy);
            if (iout >= 2 && njmid <= 2 * j + 1)
            {
                for (i = 0; i < n; ++i)
                {
                    fsafe[ipt + i * lfsafe] = dy[i];
                }
                ++(ipt);
            }
            for (i = 0; i < n; ++i)
            {
                t[j + i * km] = (yh1[i] + yh2[i] + hj * dy[i]) / 2.0;
            }

            // POLYNOMIAL EXTRAPOLATION
            if (j == 0)
            {
                return 0;
            }
            dblenj = nj[j];
            for (l = j; l >= 1; --l)
            {
                // Computing 2nd power
                d1 = dblenj / nj[l - 1];
                fac = d1 * d1 - 1.0;
                for (i = 0; i < n; ++i)
                {
                    t[l - 1 + i * km] = t[l + i * km] + (t[l + i * km] - t[l - 1 + i * km]) / fac;
                }
            }
            err = 0.0;
            // SCALING
            for (i = 0; i < n; ++i)
            {
                t1i = Math.Max(Math.Abs(y[i]), Math.Abs(t[i * km]));
                scal[i] = atol + rtol * t1i;
                // Computing 2nd power
                d1 = (t[i * km] - t[i * km + 1]) / scal[i];
                err += d1 * d1;
            }
            err = Math.Sqrt(err / n);
            if (err * uround >= 1.0)
            {
                goto L79;
            }
            if (j > 1 && err >= errold)
            {
                goto L79;
            }
            errold = Math.Max(err * 4, 1.0);

            // COMPUTE OPTIMAL STEPSIZES
            expo = 1.0 / (2 * j + 1);
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
        }

        int interp_(int n, double[] y, int imit)
        {
            double[] a = new double[31];
            int i;
            double y0, y1;
            int im;
            double ph0, ph1, ph2, ph3, yp0, yp1, fac1, fac2, aspl, bspl, ydiff;

            // COMPUTES THE COEFFICIENTS OF THE INTERPOLATION FORMULA
            // BEGIN WITH HERMITE INTERPOLATION

            for (i = 0; i < n; ++i)
            {
                y0 = y[i];
                y1 = y[2 * n + i];
                yp0 = y[n + i];
                yp1 = y[3 * n + i];
                ydiff = y1 - y0;
                aspl = -yp1 + ydiff;
                bspl = yp0 - ydiff;
                y[n + i] = ydiff;
                y[2 * n + i] = aspl;
                y[3 * n + i] = bspl;
                if (imit < 0)
                {
                    goto L100;
                }
                // COMPUTE THE DERIVATIVES OF HERMITE AT MIDPOINT
                ph0 = (y0 + y1) * 0.5 + (aspl + bspl) * 0.125;
                ph1 = ydiff + (aspl - bspl) * 0.25;
                ph2 = -(yp0 - yp1);
                ph3 = (bspl - aspl) * 6.0;
                // COMPUTE THE FURTHER COEFFICIENTS
                if (imit < 1)
                {
                    goto L20;
                }
                a[1] = (y[n * 5 + i] - ph1) * 16.0;
                if (imit < 3)
                {
                    goto L20;
                }
                a[3] = (y[n * 7 + i] - ph3 + a[1] * 3) * 16.0;
                if (imit < 5)
                {
                    goto L20;
                }

                for (im = 4; im < imit; im += 2)
                {
                    fac1 = im * (im - 1) / 2.0;
                    fac2 = fac1 * (im - 2) * (im - 3) * 2.0;
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

                for (im = 4; im <= imit; im += 2)
                {
                    fac1 = im * (im - 1) / 2.0;
                    fac2 = (double)(im * (im - 1) * (im - 2) * (im - 3));
                    a[im] = (y[n * (im + 4) + i] + a[im - 2] * fac1 - a[im - 4] * fac2) * 16.0;
                }
                L60:
                for (im = 0; im <= imit; ++im)
                {
                    y[n * (im + 4) + i] = a[im];
                }
                L100:
                ;
            }
            return 0;
        }

        double contex_(int ii, double x, double[] y, int ncon, int[] icomp, int n)
        {
            double ret_val, d1;

            int i, j, im;
            double theta, theta1, thetah, phthet;

            // THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONECTION
            // WITH THE OUTPUT-SUBROUTINE FOR ODEX. IT PROVIDES AN
            // APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.0

            // COMPUTE PLACE OF II-TH COMPONENT

            i = 0;
            for (j = 0; j < n; ++j)
            {
                if (icomp[j] == ii)
                {
                    i = j;
                }
            }

            if (i == 0)
            {
                Console.WriteLine(" NO DENSE OUTPUT AVAILABLE FOR COMP.", ii);
                return 0;
            }

            // COMPUTE THE INTERPOLATED VALUE
            theta = (x - conodx_2.xold) / conodx_2.h;
            theta1 = 1.0 - theta;
            phthet = y[i] + theta * (y[n + i] + theta1 * (y[2 * n + i] * theta + y[3 * n + i] * theta1));

            if (conodx_2.imit < 0)
            {
                ret_val = phthet;
                return ret_val;
            }

            thetah = theta - .5;
            ret_val = y[n * (conodx_2.imit + 4) + i];
            for (im = conodx_2.imit; im >= 1; --im)
            {
                ret_val = y[n * (im + 3) + i] + ret_val * thetah / im;
            }

            // Computing 2nd power
            d1 = theta * theta1;
            return phthet + d1 * d1 * ret_val;
        }
    }
}