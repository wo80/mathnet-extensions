// Based on Fortran code SEULEX
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;
    using S_fp = System.Action<int, double, double[], double[]>;
    using J_fp = System.Action<int, double, double[], double[], int>;
    using M_fp = System.Action<int, double[], int>;

    /// <summary>
    /// Numerical solution of a stiff (or differential algebraic) system of first order
    /// ordinary differential equations My'=f(x,y).
    ///
    /// This is an extrapolation-algorithm, based on the linearly implicit euler method
    /// (with step size control and order selection).
    ///
    /// Authors: E. Hairer and G. Wanner (inclusion of dense output by E. Hairer and A. Ostermann)
    ///
    /// This code is part of the book:
    ///      E. Hairer and G. Wanner
    ///      Solving Ordinary Differential Equations II.
    ///      Stiff and differential-algebraic problems. (2nd edition)
    ///      Springer-Verlag (1996)
    /// </summary>
    public class SemiImplicitEuler
    {

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        struct coseu
        {
            public double xold, h;
            public int nrd, kright;
        }

        coseu coseu_1;

        double fac1, fac2, fac3, fac4, safe1, safe2;

        double rtol, atol;

        public int nstep, naccpt, nrejct, nsol, ndec; // TODO: remove

        /**
         *     N           DIMENSION OF THE SYSTEM
         *
         *     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
         *                 VALUE OF F(X,Y)
         *
         *     IFCN        GIVES INFORMATION ON FCN:
         *                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
         *                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
         *
         *     X           INITIAL X-VALUE
         *
         *     Y(N)        INITIAL VALUES FOR Y
         *
         *     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
         *
         *     H           INITIAL STEP SIZE GUESS;
         *                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
         *                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
         *                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
         *                 ADAPTS ITS STEP SIZE (IF H=0.D0, THE CODE PUTS H=1.D-6).
         *
         *     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES.
         *
         *     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
         *                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
         *                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
         *                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
         *                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
         *                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
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
         *     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
         *                 MATRIX M.
         *                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
         *                 MATRIX AND NEEDS NOT TO BE DEFINED;
         *                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
         *                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
         *                 IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
         *                    AS FULL MATRIX LIKE
         *                         AM(I,J) = M(I,J)
         *
         *     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
         *                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
         *                       MATRIX, MAS IS NEVER CALLED.
         *                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
         *
         *     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT:
         *                    IOUT=0: SUBROUTINE IS NEVER CALLED
         *                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
         *                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
         *
         *     OUTPUT PARAMETERS
         *     -----------------
         *     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
         *                 (AFTER SUCCESSFUL RETURN X=XEND)
         *
         *     Y(N)        SOLUTION AT X
         *
         *     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
         *
         *     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
         *                   IDID=1  COMPUTATION SUCCESSFUL,
         *                   IDID=-1 COMPUTATION UNSUCCESSFUL.
         */
        public int seulex_(int n, S_fp fcn, int ifcn, double x,
            double[] y, double xend, double h, double rtol, double atol,
            J_fp jac, int ijac, M_fp mas, int imas, S_fp solout, int iout)
        {
            int km, km2;
            int ijob;
            double hmax;
            int nmax;
            double thet;
            double wkdec;
            double wkjac;
            double wkfcn;
            int lambda, nsequ;
            double wksol, wkrow;
            bool implct;

            bool autnms;
            double uround;
            
            nstep = 0;
            naccpt = 0;
            nrejct = 0;
            ndec = 0;
            nsol = 0;

            this.rtol = rtol;
            this.atol = atol;

            // NMAX - The maximal number of steps
            nmax = 100000;

            // KM - Maximum number of columns in the extrapolation table.
            km = 12;

            if (km <= 2)
            {
                Console.WriteLine(" CURIOUS INPUT KM=", km);
                return -1;
            }

            // NSEQU - Choice of step size sequence
            //     IF NSEQU = 1 THEN 1,2,3,4,6,8,12,16,24,32,48...
            //     IF NSEQU = 2 THEN 2,3,4,6,8,12,16,24,32,48,64...
            //     IF NSEQU = 3 THEN 1,2,3,4,5,6,7,8,9,10...
            //     IF NSEQU = 4 THEN 2,3,4,5,6,7,8,9,10,11...
            nsequ = 2;

            if (nsequ < 1 || nsequ > 4)
            {
                Console.WriteLine(" CURIOUS INPUT NSEQU=", nsequ);
                return -1;
            }

            // LAMBDA - Parameter for dense output
            lambda = 0;

            if (lambda < 0 || lambda > 1)
            {
                Console.WriteLine(" CURIOUS INPUT LAMBDA=", lambda);
                return -1;
            }

            // UROUND   Smallest number satisfying 1.0+UROUND>1.0
            uround = 1e-16;

            // Maximal step size
            hmax = xend - x;

            // THET - Decides whether the jacobian should be recomputed;
            //     Increase THET, to 0.01 say, when jacobian evaluations are costly.
            //     For small systems THET should be smaller.
            thet = Math.Min(1e-4, rtol);

            // FAC1,FAC2 - Parameters for step size selection
            //     The new step size for the j-th diagonal entry is chosen subject
            //     to the restriction
            //         facmin / FAC2 <= hnew(j) / hold <= 1 / facmin
            //     where
            //         facmin = FAC1 * (1 / (j-1))
            fac1 = 0.1;
            fac2 = 4.0;

            // FAC3, FAC4 - Parameters for the order selection
            //     Order is decreased if    w(k-1) <= w(k) * FAC3
            //     Order is increased if    w(k) <= w(k-1) * FAC4
            fac3 = 0.7;
            fac4 = 0.9;

            // SAFE1, SAFE2 - Safety factors for step size prediction
            //     hnew = h * SAFE1 * (SAFE2 * tol / err)^(1/(j-1))
            safe1 = 0.6;
            safe2 = 0.93;

            // WKFCN,WKJAC,WKDEC,WKSOL  estimated work for  FCN,JAC,DEC,SOL
            wkfcn = 1.0;
            wkjac = 5.0;
            wkdec = 1.0;
            wksol = 1.0;
            wkrow = wkfcn + wksol;

            // Check if tolerances are ok.
            if (atol <= 0.0 || rtol <= uround * 10.0)
            {
                Console.WriteLine(" TOLERANCES ARE TOO SMALL");
                return -1;
            }

            // Autonomous, implicit, banded or not ?
            autnms = ifcn == 0;
            implct = imas != 0;

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
            km2 = km * (km + 1) / 2;
            var yh = new double[n];
            var dy = new double[n];
            var fx = new double[n];
            var yhh = new double[n];
            var dyh = new double[n];
            var del = new double[n];
            var wh = new double[n];
            var scal = new double[n];
            var hh = new double[km];
            var w = new double[km];
            var a = new double[km];
            var _jac = new double[n * n];
            var e = new double[n * n];
            var _mas = new double[n * n];
            var t = new double[n * km];
            var ifac = new double[km];

            // TODO: if dense
            var de = new double[(km + 2) * n];
            var fsafe = new double[km2 * n];

            // Entry points for int workspace
            //var co = new int[n];
            var ip = new int[n];
            var nj = new int[km];
            var iph = new int[km];

            // Call to core integrator
            int idid = seucor_(n, fcn, x, y, xend, hmax, h, km,
                jac, ijac, mas, solout, iout, ijob, nmax, uround,
                nsequ, autnms, implct, yh,
                dy, fx, yhh, dyh, del,
                 wh, scal, hh, w,
                a, _jac, e, _mas, t, ip,
                 nj, thet, wkjac, wkdec, wkrow, km2, ifac,
                 fsafe, lambda, de);

            return idid;
        }

        // ... and here is the core integrator

        int seucor_(int n, S_fp fcn, double x, double[] y,
            double xend, double hmax, double h, int km,
            J_fp jac, int ijac, M_fp mas, S_fp solout, int iout, int ijob,
            int nmax, double uround,
            int nsequ, bool autnms, bool implct, double[] yh,
            double[] dy, double[] fx, double[] yhh, double[] dyh,
            double[] del, double[] wh, double[] scal, double[] hh,
            double[] w, double[] a, double[] fjac, double[] e,
            double[] fmas, double[] t, int[] ip, int[] nj,
            double thet, double wkjac, double wkdec, double wkrow,
            int km2, double[] facul, double[] fsafe,
            int lambda, double[] dens)
        {
            int i1, i2, i3;
            double d1;

            int i, j, k, l, kc = 0, ii, kk, i_n;
            double t1i, err;
            int ipt, klr, krn, lbeg, lend, lrde;
            double delt;
            bool last;
            double xold;
            bool atov;
            double hopt = 0;

            int kopt;
            double xsol, facnj, theta, ysafe, hmaxn;
            int nrsol, irtrn, nsolu;
            bool caljac;
            double dblenj;
            bool reject;
            double factor;
            double errold = 0, posneg;


            // Core integrator for seulex

            if (iout == 2)
            {
                coseu_1.nrd = n;
                // Compute the factorials
                facul[0] = 1.0;
                for (i = 0; i < km - 1; ++i)
                {
                    facul[i + 1] = i * facul[i];
                }
            }
            // Compute mass matrix for implicit case
            if (implct)
            {
                mas(n, fmas, n);
            }

            // Initializations
            lrde = (km + 2) * n;

            // Define the step size sequence
            if (nsequ == 1)
            {
                nj[0] = 1;
                nj[1] = 2;
                nj[2] = 3;
                for (i = 3; i < km; ++i)
                {
                    nj[i] = nj[i - 2] << 1;
                }
            }
            if (nsequ == 2)
            {
                nj[0] = 2;
                nj[1] = 3;
                for (i = 2; i < km; ++i)
                {
                    nj[i] = nj[i - 2] << 1;
                }
            }
            for (i = 0; i < km; ++i)
            {
                if (nsequ == 3)
                {
                    nj[i] = i;
                }
                if (nsequ == 4)
                {
                    nj[i] = i + 1;
                }
            }
            a[0] = wkjac + nj[0] * wkrow + wkdec;
            for (i = 1; i < km; ++i)
            {
                a[i] = a[i - 1] + (nj[i] - 1) * wkrow + wkdec;
            }
            d1 = xend - x;
            posneg = d_sign(1.0, d1);
            k = Math.Max(2, Math.Min(km - 2, (int)(-Math.Log10(rtol + atol) * 0.6 + 1.5)));
            hmaxn = Math.Min(Math.Abs(hmax), Math.Abs(xend - x));
            h = Math.Max(Math.Abs(h), 1e-6);
            h = posneg * Math.Min(h, hmaxn);
            theta = Math.Abs(thet) * 2;
            if (iout != 0)
            {
                irtrn = 1;
                nrsol = 1;
                xsol = x;
                for (i = 0; i < n; ++i)
                {
                    yh[i] = y[i];
                }
                nsolu = n;
                xold = x;
                //solout(nrsol, xold, xsol, yh, dens, lrde, icomp, nrd, nsolu, irtrn);
                if (irtrn < 0)
                {
                    goto L120;
                }
                ipt = 0;
            }
            err = 0.0;
            w[0] = 1e30;
            for (i = 0; i < n; ++i)
            {
                scal[i] = atol + rtol * Math.Abs(y[i]);
            }
            caljac = false;
            reject = false;
            last = false;
            L10:
            if (reject)
            {
                theta = Math.Abs(thet) * 2;
            }
            atov = false;

            // Is xend reached in the next step? 
            if (Math.Abs(xend - x) * 0.1 <= Math.Abs(x) * uround)
            {
                goto L110;
            }
            hopt = h;
            h = posneg * Math.Min(Math.Min(Math.Abs(h), Math.Abs(xend - x)), hmaxn);
            if ((x + h * 1.01 - xend) * posneg > 0.0)
            {
                h = xend - x;
                last = true;
            }
            if (autnms)
            {
                fcn(n, x, y, dy);
            }
            if (theta > thet && !caljac)
            {
                // Computation of the jacobian
                if (ijac == 0)
                {
                    // Compute jacobian matrix numerically 
                    if (!(autnms))
                    {
                        fcn(n, x, y, dy);
                    }

                    // Jacobian is full 
                    for (i = 0; i < n; ++i)
                    {
                        ysafe = y[i];
                        delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                        y[i] = ysafe + delt;
                        fcn(n, x, y, yh);
                        for (j = 0; j < n; ++j)
                        {
                            fjac[j + i * n] = (yh[j] - dy[j]) / delt;
                        }
                        y[i] = ysafe;
                    }
                }
                else
                {
                    // Compute jacobian matrix analytically 
                    jac(n, x, y, fjac, n);
                }
                caljac = true;
            }

            // The first and last step 
            if (nstep == 0 || last)
            {
                ipt = 0;
                ++(nstep);
                for (j = 0; j < k; ++j)
                {
                    kc = j + 1;
                    seul_(j, n, fcn, x, y, dy, fx, fjac,
                        fmas, e, ip, ref h, km, hmaxn, t, scal,
                        nj, hh, w, a, yhh, dyh, del,
                        wh, ref err, ref theta,
                        ref errold,
                        autnms, implct, ref reject,
                        ref atov, fsafe, km2, iout, ipt, ijob);
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

            // Basic integration step 
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
                seul_(j, n, fcn, x, y, dy, fx, fjac,
                    fmas, e, ip,
                    ref h, km, hmaxn, t, scal, nj, hh, w,
                    a, yhh, dyh, del, wh, ref err, ref theta,
                    ref errold, autnms, implct,
                    ref reject, ref atov, fsafe, km2, iout,
                    ipt, ijob);
                if (atov)
                {
                    goto L10;
                }
            }

            // Convergence monitor 
            if (k == 2 || reject)
            {
                goto L50;
            }
            if (err <= 1.0)
            {
                goto L60;
            }
            if (err > (nj[k] * nj[k - 1]) * 4.0)
            {
                goto L100;
            }
            L50:
            seul_(k - 1, n, fcn, x, y, dy, fx, fjac,
                fmas, e, ip, ref h,
                km, hmaxn, t, scal, nj, hh, w, a,
                yhh, dyh, del, wh, ref err, ref theta, ref errold,
                autnms, implct, ref reject, ref atov,
                fsafe, km2, iout, ipt, ijob);
            if (atov)
            {
                goto L10;
            }
            kc = k;
            if (err <= 1.0)
            {
                goto L60;
            }
            // Hope for convergence in line K+1 
            L55:
            if (err > (double)nj[k] * 2.0)
            {
                goto L100;
            }
            kc = k + 1;
            seul_(kc - 1, n, fcn, x, y, dy, fx, fjac,
                fmas, e, ip, ref h,
                km, hmaxn, t, scal, nj, hh, w, a,
                yhh, dyh, del, wh, ref err, ref theta, ref errold,
                autnms, implct, ref reject, ref atov,
                fsafe, km2, iout, ipt, ijob);
            if (atov)
            {
                goto L10;
            }
            if (err > 1.0)
            {
                goto L100;
            }

            // Step is accepted 
            L60:
            xold = x;
            x += h;
            if (iout == 2)
            {
                coseu_1.kright = kc;
                for (i = 0; i < n; ++i)
                {
                    dens[i] = y[i];
                }
            }
            for (i = 0; i < n; ++i)
            {
                t1i = t[i * km];
                scal[i] = atol + rtol * Math.Abs(t1i);
                y[i] = t1i;
            }
            ++(naccpt);
            caljac = false;
            if (iout == 2)
            {
                coseu_1.xold = xold;
                coseu_1.h = h;
                for (i = 0; i < n; ++i)
                {
                    dens[n + i] = y[i];
                }
                i1 = coseu_1.kright - 1;
                for (klr = 0; klr < i1; ++klr)
                {
                    // Compute differences 
                    if (klr >= 2)
                    {
                        i2 = kc;
                        for (kk = klr - 1; kk < i2; ++kk)
                        {
                            lbeg = (kk + 1) * kk / 2;
                            lend = lbeg - kk + 2;
                            for (l = lbeg; l >= lend; --l)
                            {
                                for (i = 0; i < n; ++i)
                                {
                                    fsafe[l + i * km2] -= fsafe[l - 1 + i * km2];
                                }
                            }
                        }
                    }
                    // Compute derivatives at right end
                    for (kk = klr + lambda - 1; kk < kc; ++kk)
                    {
                        facnj = Math.Pow(nj[kk], klr) / facul[klr + 1]; // TODO: pow_di()
                        ipt = (kk + 1) * kk / 2;
                        for (i = 0; i < n; ++i)
                        {
                            krn = (kk - lambda + 1) * n;
                            dens[krn + i] = fsafe[ipt + i * km2] * facnj;
                        }
                    }
                    for (j = klr + lambda; j < kc; ++j)
                    {
                        dblenj = (double)nj[j];
                        i3 = klr + lambda + 1;
                        for (l = j; l >= i3; --l)
                        {
                            factor = dblenj / nj[l - 1] - 1.0;
                            for (i = 0; i < n; ++i)
                            {
                                krn = (l - lambda + 1) * n + i;
                                dens[krn - n] = dens[krn] + (dens[krn] - dens[krn - n]) / factor;
                            }
                        }
                    }
                }
                // Compute the coefficients of the interpolation polynomial 
                for (i_n = 0; i_n < n; ++i_n)
                {
                    for (j = 0; j < coseu_1.kright; ++j)
                    {
                        ii = n * j + i_n;
                        dens[ii] -= dens[ii - n];
                    }
                }
            }
            if (iout != 0)
            {
                irtrn = 1;
                nrsol = naccpt + 1;
                xsol = x;
                for (i = 0; i < n; ++i)
                {
                    yh[i] = y[i];
                }
                nsolu = n;
                //solout(nrsol, xold, xsol, yh, dens, lrde, icomp, nrd, nsolu, irtrn);
                if (irtrn < 0)
                {
                    goto L120;
                }
            }
            // Compute optimal order 
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
            // After a rejected step 
            L80:
            if (reject)
            {
                k = Math.Min(kopt, kc);
                h = posneg * Math.Min(Math.Abs(h), Math.Abs(hh[k - 1]));
                reject = false;
                goto L10;
            }
            // Compute step size for next step 
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

            // Step is rejected 
            L100:
            k = Math.Min(Math.Min(k, kc), km - 1);
            if (k > 2 && w[k - 2] < w[k - 1] * fac3)
            {
                --k;
            }
            ++(nrejct);
            h = posneg * hh[k - 1];
            last = false;
            reject = true;
            if (caljac)
            {
                goto L30;
            }
            goto L10;
            // Solution exit 
            L110:
            h = hopt;
            return 1;
            // Fail exit 
            L120:
            //do_fio("", x, h);
            return -1;
        }

        // This subroutine computes the j-th line of the 
        // extrapolation table and provides an estimate 
        // of the optimal step size 
        int seul_(int jj, int n, S_fp fcn, double x,
            double[] y, double[] dy, double[] fx, double[] fjac,
            double[] fmas, double[] e,
            int[] ip, ref double h, int km, double hmaxn,
            double[] t, double[] scal, int[] nj, double[] hh,
            double[] w, double[] a, double[] yh, double[] dyh,
            double[] del, double[] wh, ref double err,
            ref double theta,
            ref double errold, bool autnms,
            bool implct, ref bool reject, ref bool atov,
            double[] fsafe, int km2, int iout, int ipt,
            int ijob)
        {
            double d1, d2;

            int i, j, l, m;
            double hj;
            int mm;
            double hji;
            int ier = 0;
            double sum, del1, del2;
            double fac, facmin, expo;

            hj = h / nj[jj];
            hji = 1.0 / hj;
            dc_decsol.decomr_(n, fjac, fmas, hji, e, ip, ref ier, ijob);
            if (ier != 0)
            {
                goto L79;
            }
            ++(ndec);

            // Starting procedure 
            if (!(autnms))
            {
                d1 = x + hj;
                fcn(n, d1, y, dy);
            }
            for (i = 0; i < n; ++i)
            {
                yh[i] = y[i];
                del[i] = dy[i];
            }
            dc_decsol.slvseu_(n, e, ip, del, ijob);
            ++(nsol);
            m = nj[jj];
            if (iout == 2 && m == jj + 1)
            {
                ++(ipt);
                for (i = 0; i < n; ++i)
                {
                    fsafe[ipt + i * km2] = del[i];
                }
            }

            // Semi-implicit Euler method 
            if (m > 1)
            {
                for (mm = 0; mm < m - 1; ++mm)
                {
                    for (i = 0; i < n; ++i)
                    {
                        yh[i] += del[i];
                    }
                    if (autnms)
                    {
                        fcn(n, x + hj * (mm + 1), yh, dyh);
                    }
                    else
                    {
                        fcn(n, x + hj * (mm + 2), yh, dyh);
                    }
                    if (mm == 0 && jj < 2)
                    {
                        // Stability check 
                        del1 = 0.0;
                        for (i = 0; i < n; ++i)
                        {
                            /* Computing 2nd power */
                            d1 = del[i] / scal[i];
                            del1 += d1 * d1;
                        }
                        del1 = Math.Sqrt(del1);
                        if (implct)
                        {
                            for (i = 0; i < n; ++i)
                            {
                                wh[i] = del[i];
                            }
                            //if (mlb == nm1)
                            {
                                for (i = 0; i < n; ++i)
                                {
                                    sum = 0.0;
                                    for (j = 0; j < n; ++j)
                                    {
                                        sum += fmas[i + j * n] * wh[j];
                                    }
                                    del[i] = sum;
                                }
                            }
                        }
                        if (!(autnms))
                        {
                            d1 = x + hj;
                            fcn(n, d1, yh, wh);
                            for (i = 0; i < n; ++i)
                            {
                                del[i] = wh[i] - del[i] * hji;
                            }
                        }
                        else
                        {
                            for (i = 0; i < n; ++i)
                            {
                                del[i] = dyh[i] - del[i] * hji;
                            }
                        }
                        dc_decsol.slvseu_(n, e, ip, del, ijob);
                        ++(nsol);
                        del2 = 0.0;
                        for (i = 0; i < n; ++i)
                        {
                            /* Computing 2nd power */
                            d1 = del[i] / scal[i];
                            del2 += d1 * d1;
                        }
                        del2 = Math.Sqrt(del2);
                        theta = del2 / Math.Max(1.0, del1);
                        if (theta > 1.0)
                        {
                            goto L79;
                        }
                    }
                    dc_decsol.slvseu_(n, e, ip, dyh, ijob);
                    ++(nsol);
                    for (i = 0; i < n; ++i)
                    {
                        del[i] = dyh[i];
                    }
                    if (iout == 2 && mm + 1 >= m - jj - 1)
                    {
                        ++(ipt);
                        for (i = 0; i < n; ++i)
                        {
                            fsafe[ipt + i * km2] = del[i];
                        }
                    }
                }
            }
            for (i = 0; i < n; ++i)
            {
                t[jj + i * km] = yh[i] + del[i];
            }

            // Polynomial extrapolation 
            if (jj == 0)
            {
                return 0;
            }
            for (l = jj; l >= 1; --l)
            {
                fac = nj[jj] / (double)nj[l - 1] - 1.0;
                for (i = 0; i < n; ++i)
                {
                    t[l - 1 + i * km] = t[l + i * km] + (t[l + i * km] - t[l - 1 + i * km]) / fac;
                }
            }
            err = 0.0;
            for (i = 0; i < n; ++i)
            {
                /* Computing 2nd power */
                d2 = Math.Min(Math.Abs(t[i * km] - t[i * km + 1]) / scal[i], 1e15);
                err += d2 * d2;
            }
            if (err >= 1e30)
            {
                goto L79;
            }
            err = Math.Sqrt(err / (double)(n));
            if (jj > 1 && err >= errold)
            {
                goto L79;
            }
            d1 = err * 4;
            errold = Math.Max(d1, 1.0);
            // Compute optimal step sizes 
            expo = 1.0 / (jj + 1);
            facmin = Math.Pow(fac1, expo);
            fac = Math.Min(fac2 / facmin, Math.Max(facmin, Math.Pow(err / safe1, expo) / safe2));
            fac = 1.0 / fac;

            hh[jj] = Math.Min(Math.Abs(h) * fac, hmaxn);
            w[jj] = a[jj] / hh[jj];
            return 0;
            L79:
            atov = true;
            h *= 0.5;
            reject = true;
            return 0;
        }

        //  This function can be used for coninuous output in conection 
        //  with the output-subroutine for seulex. it provides an 
        //  approximation to the ii-th component of the solution at x. 
        double contex_(int i, double x, double[] rc)
        {
            // Compute the interpolated value 
            double theta = (x - coseu_1.xold) / coseu_1.h;
            double ret_val = rc[coseu_1.kright * coseu_1.nrd + i];
            for (int j = 1; j < coseu_1.kright; ++j)
            {
                ret_val = rc[(coseu_1.kright + 1 - j) * coseu_1.nrd + i] + ret_val * (theta - 1.0);
            }
            return rc[i] + ret_val * theta;
        }
    }
}
