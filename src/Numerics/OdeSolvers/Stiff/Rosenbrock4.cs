// Based on Fortran code RODAS
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;

    using S_fp = System.Action<int, double, double[], double[]>;
    using J_fp = System.Action<int, double, double[], LinearAlgebra.Double.DenseMatrix>;

    public class Rosenbrock4
    {
        RosenbrockErrorController controller;

        double[] rtol, atol;
        
        double hold;

        struct conros
        {
            public double xold, h;
            public int n;
        };

        conros conros_ = new conros();

        public Rosenbrock4(RosenbrockErrorController controller)
        {
            this.controller = controller;
        }

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c21, c31,
            c32, c41, c42, c43, c51, c52, c53, c54, c61, c62, c63, c64, c65,
            d21, d22, d23, d24, d25, d31, d32, d33, d34, d35;
        double c2, c3, c4, d1, d2, d3, d4, gamma;

        /**
         *     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
         *     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
         *     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4
         *     (WITH STEP SIZE CONTROL).
         *     C.F. SECTIONS IV.7  AND VI.3
         *
         *     AUTHORS: E. HAIRER AND G. WANNER
         *              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
         *              CH-1211 GENEVE 24, SWITZERLAND
         *              E-MAIL:  Ernst.Hairer@math.unige.ch
         *                       Gerhard.Wanner@math.unige.ch
         *
         *     THIS CODE IS PART OF THE BOOK:
         *         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
         *         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
         *         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
         *         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
         *
         *     VERSION OF OCTOBER 28, 1996
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
         *                 RPAR, IPAR (SEE BELOW)
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
         *                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS
         *                 FURNISHED BY THE CALLING PROGRAM.
         *                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
         *                    BE FULL AND THE PARTIAL DERIVATIVES ARE
         *                    STORED IN DFY AS
         *                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
         *
         *     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
         *                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X
         *                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1;
         *                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0).
         *                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM
         *                    SUBROUTINE DFX(N,X,Y,FX,RPAR,IPAR)
         *                    DOUBLE PRECISION X,Y(N),FX(N)
         *                    FX(1)= 0...
         *
         *     MAS         THE MASS-MATRIX M.
         * ----------------------------------------------------------------------
         *
         *     SOPHISTICATED SETTING OF PARAMETERS
         *     -----------------------------------
         *              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
         *              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(4)
         *              AS WELL AS IWORK(1),IWORK(2) DIFFERENT FROM ZERO.
         *              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
         *
         *    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
         *              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
         *
         *    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
         *              IF IWORK(2).EQ.1  METHOD (SEE BOOK, PAGE 452)
         *              IF IWORK(2).EQ.2  SAME METHOD WITH DIFFERENT PARAMETERS
         *              IF IWORK(2).EQ.3  METHOD WITH COEFF. OF GERD STEINEBACH
         *              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1.
         *
         *    IWORK(3)  SWITCH FOR STEP SIZE STRATEGY
         *              IF IWORK(3).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
         *              IF IWORK(3).EQ.2  CLASSICAL APPROACH
         *              THE DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=1.
         *
         *    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
         *
         *    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
         *
         * -----------------------------------------------------------------------
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
         *                   IDID= 1  COMPUTATION SUCCESSFUL,
         *                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
         *                   IDID=-1  INPUT IS NOT CONSISTENT,
         *                   IDID=-2  LARGER NMAX IS NEEDED,
         *                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
         *                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
         *
         *   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
         *                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
         *   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
         *                      OR NUMERICALLY)
         *   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
         *   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
         *   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
         *                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
         *   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
         *   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
         */
        public int rodas_(int n, S_fp fcn, int ifcn, double x, double[] y, double xend, double h, double[] rtol, double[] atol,
            int itol, J_fp jac, S_fp dfx, DenseMatrix mas, double[] work, int[] iwork)
        {
            int ndec, njac, nfcn, nsol, nstep;
            int meth;
            int nmax;
            bool implct;
            bool autnms;
            double uround;

            this.rtol = rtol;
            this.atol = atol;

            // Function Body
            nfcn = 0;
            nstep = 0;
            njac = 0;
            ndec = 0;
            nsol = 0;

            // NMAX , THE MAXIMAL NUMBER OF STEPS
            if (iwork[0] == 0)
            {
                nmax = 100000;
            }
            else
            {
                nmax = iwork[0];
                if (nmax <= 0)
                {
                    Console.WriteLine(" Wrong input IWORK(1)=", iwork[0]);
                    return -1;
                }
            }

            // METH   COEFFICIENTS OF THE METHOD
            if (iwork[1] == 0)
            {
                meth = 1;
            }
            else
            {
                meth = iwork[1];
                if (meth <= 0 || meth >= 4)
                {
                    Console.WriteLine(" Curious input IWORK(2)=", iwork[1]);
                    return -1;
                }
            }

            // UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0
            if (work[0] == 0.0)
            {
                uround = 1e-16;
            }
            else
            {
                uround = work[0];
                if (uround < 1e-16 || uround >= 1.0)
                {
                    Console.WriteLine(" Coefficients have 16 digits, UROUND=", work[0]);
                    return -1;
                }
            }
            
            // Check if tolerances are OK.
            if (itol == 0)
            {
                if (atol[0] <= 0.0 || rtol[0] <= uround * 10.0)
                {
                    Console.WriteLine(" Tolerances are too small");
                    return -1;
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    if (atol[i] <= 0.0 || rtol[i] <= uround * 10.0)
                    {
                        Console.WriteLine(" Tolerances(%d) are too small", i);
                        return -1;
                    }
                }
            }

            // Autonomous, implicit or not ?
            autnms = ifcn == 0;
            implct = mas != null;
            
            // Prepare the entry-points for the arrays in work
            var ynew = new double[n];
            var dy1 = new double[n];
            var dy = new double[n];
            var ak1 = new double[n];
            var ak2 = new double[n];
            var ak3 = new double[n];
            var ak4 = new double[n];
            var ak5 = new double[n];
            var ak6 = new double[n];
            var fx = new double[n];
            var con = new double[4 * n];
            
            // Call to core integrator
            int idid = roscor_(n, fcn, x, y, xend, h, rtol, atol,
                itol, jac, dfx, mas,
                nmax, uround, meth,
                autnms, implct,
                ynew, dy1, dy,
                ak1, ak2, ak3, ak4, ak5, ak6, fx,
                con, ref nfcn,
                ref njac, ref nstep, ref ndec, ref nsol);

            iwork[13] = nfcn;
            iwork[14] = njac;
            iwork[15] = nstep;
            iwork[18] = ndec;
            iwork[19] = nsol;

            return idid;
        }
        
        // ... and here is the core integrator
        int roscor_(int n, S_fp fcn, double x, double[] y,
            double xend, double h, double[] rtol, double[] atol,
            int itol, J_fp jac, S_fp dfx, DenseMatrix mas,
            int nmax, double uround, int meth,
            bool autnms, bool implct,
            double[] ynew, double[] dy1,
            double[] dy, double[] ak1, double[] ak2, double[] ak3,
            double[] ak4, double[] ak5, double[] ak6, double[] fx,
            double[] cont,
            ref int nfcn, ref int njac, ref int nstep, ref int ndec, ref int nsol)
        {
            var lu = new ReusableLU(n);
            
            // Local variables
            int i, j;
            double hd1 = 0, hd2 = 0, hd3 = 0, hd4 = 0;
            double fac, hc21, hc31, hc32, hc41, hc42, hc43, hc51, hc52, hc53, hc54, hc61, hc62;
            int ier = 0;
            double hc63, hc64, hc65, err;
            double delt;
            bool last;
            double hopt = 0;
            
            //int ijob = implct ? 5 : 1;

            double ysafe;
            int nsing;
            double xdelt;
            
            double posneg;

            //
            // Core integrator for RODAS
            //

            var fjac = new DenseMatrix(n);
            var fmodjac = new DenseMatrix(n);

            // Function Body
            conros_.n = n;
            
            // Set the parameters of the method
            rocoe_(meth);

            // Initial preparations
            posneg = d_sign(1.0, xend - x);
            
            if (Math.Abs(h) <= uround * 10.0)
            {
                h = 1e-6;
            }
            h = Math.Min(Math.Abs(h), Math.Abs(xend - x));

            h = d_sign(h, posneg);
            last = false;
            nsing = 0;

            if (autnms)
            {
                hd1 = 0.0;
                hd2 = 0.0;
                hd3 = 0.0;
                hd4 = 0.0;
            }
            
            {
                conros_.xold = x;
                conros_.h = h;
                //solout(naccpt + 1, conros_.xold, x, y, cont, lrc, n, irtrn);
            }

            // Basic integration step
            L1:
            if (nstep > nmax)
            {
                //do_fio("  Exit of RODAS at X= ,e18.4", (x));
                Console.WriteLine(" More than NMAX =" + nmax + "steps are needed");
                return -2;
            }

            if (Math.Abs(h) * 0.1 <= Math.Abs(x) * uround)
            {
                //do_fio("  Exit of RODAS at X= ,e18.4", (x));
                Console.WriteLine(" Step size too small, H=" + h);
                return -3;
            }

            if (last)
            {
                h = hopt;

                return 1;
            }
            hopt = h;
            if ((x + h * 1.0001 - xend) * posneg >= 0.0)
            {

                h = xend - x;
                last = true;
            }

            //
            // Computation of the Jacobian
            //
            fcn(n, x, y, dy1);
            nfcn++;
            njac++;
            if (jac == null)
            {
                // Compute Jacobian matrix numerically
                // Jacobian is full
                for (i = 0; i < n; ++i)
                {
                    ysafe = y[i];
                    delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                    y[i] = ysafe + delt;
                    fcn(n, x, y, ak1);
                    for (j = 0; j < n; ++j)
                    {
                        fjac.At(i, j, (ak1[j] - dy1[j]) / delt);
                    }
                    y[i] = ysafe;
                }
            }
            else
            {
                // Compute jacobian matrix analytically
                jac(n, x, y, fjac);
            }

            if (!(autnms))
            {
                if (dfx == null)
                {
                    // Compute numerically the derivative with respect to x
                    delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(x)));
                    xdelt = x + delt;
                    fcn(n, xdelt, y, ak1);
                    for (j = 0; j < n; ++j)
                    {
                        fx[j] = (ak1[j] - dy1[j]) / delt;
                    }
                }
                else
                {
                    // Compute analytically the derivative with respect to x
                    dfx(n, x, y, fx);
                }
            }
            L2:
            //
            // Compute the stages
            //
            fac = 1.0 / (h * gamma);
            Factorize(n, lu, fmodjac, fjac, mas, fac);
            //dc_decsol.decomr_(n, fjac, ldjac, mas, fac, e, lde, ip, ref ier, ijob, implct, ip);
            
            if (ier != 0) // TODO: check determinant?
            {
                // Singular matrix
                ++nsing;
                if (nsing >= 5)
                {
                    //do_fio("  Exit of RODAS at X= ,e18.4", (x));
                    throw new Exception("Matrix is repeatedly singular.");
                    //return -4;
                }

                h *= 0.5;
                //reject = true; // TODO: controller.Reject();
                last = false;
                goto L2;
            }

            ++(ndec);

            // Prepare for the computation of the 6 stages
            hc21 = c21 / h;
            hc31 = c31 / h;
            hc32 = c32 / h;
            hc41 = c41 / h;
            hc42 = c42 / h;
            hc43 = c43 / h;
            hc51 = c51 / h;
            hc52 = c52 / h;
            hc53 = c53 / h;
            hc54 = c54 / h;
            hc61 = c61 / h;
            hc62 = c62 / h;
            hc63 = c63 / h;
            hc64 = c64 / h;
            hc65 = c65 / h;
            if (!(autnms))
            {
                hd1 = h * d1;
                hd2 = h * d2;
                hd3 = h * d3;
                hd4 = h * d4;
            }

            // THE STAGES
            Solve(n, lu, mas, dy1, ak1, fx, ynew, hd1, false);
            //dc_decsol.slvrod_(n, fjac, ldjac, mas, fac, e, lde, ip, dy1, ak1, fx, ynew, hd1, ijob, false);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a21 * ak1[i];
            }

            fcn(n, x + c2 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc21 * ak1[i];
            }
            Solve(n, lu, mas, dy, ak2, fx, ynew, hd2, true);
            //dc_decsol.slvrod_(n, fjac, ldjac, mas, fac, e, lde, ip, dy, ak2, fx, ynew, hd2, ijob, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a31 * ak1[i] + a32 * ak2[i];
            }

            fcn(n, x + c3 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc31 * ak1[i] + hc32 * ak2[i];
            }
            Solve(n, lu, mas, dy, ak3, fx, ynew, hd3, true);
            //dc_decsol.slvrod_(n, fjac, ldjac, mas, fac, e, lde, ip, dy, ak3, fx, ynew, hd3, ijob, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a41 * ak1[i] + a42 * ak2[i] + a43 * ak3[i];
            }

            fcn(n, x + c4 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc41 * ak1[i] + hc42 * ak2[i] + hc43 * ak3[i];
            }
            Solve(n, lu, mas, dy, ak4, fx, ynew, hd4, true);
            //dc_decsol.slvrod_(n, fjac, ldjac, mas, fac, e, lde, ip, dy, ak4, fx, ynew, hd4, ijob, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a51 * ak1[i] + a52 * ak2[i] + a53 * ak3[i]
                    + a54 * ak4[i];
            }

            fcn(n, x + h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ak6[i] = hc52 * ak2[i] + hc54 * ak4[i] + hc51 * ak1[i] + hc53 * ak3[i];
            }
            Solve(n, lu, mas, dy, ak5, fx, ak6, 0.0, true);
            //dc_decsol.slvrod_(n, fjac, ldjac, mas, fac, e, lde, ip, dy, ak5, fx, ak6, 0.0, ijob, true);

            // ---- Embedded solution
            for (i = 0; i < n; i++)
            {
                ynew[i] += ak5[i];
            }

            fcn(n, x + h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                cont[i] = hc61 * ak1[i] + hc62 * ak2[i] + hc65 * ak5[i] + hc64 * ak4[i] + hc63 * ak3[i];
            }
            Solve(n, lu, mas, dy, ak6, fx, cont, 0.0, true);
            //dc_decsol.slvrod_(n, fjac, ldjac, mas, fac, e, lde, ip, dy, ak6, fx, cont, 0.0, ijob, true);

            // ---- New solution
            for (i = 0; i < n; i++)
            {
                ynew[i] += ak6[i];
            }

            nsol += 6;

            nfcn += 5;

            // ---- Dense output
            //if (dense)
            {
                for (i = 0; i < n; ++i)
                {
                    cont[i] = y[i];
                    cont[i + n * 2] = d21 * ak1[i] + d22 * ak2[i] + d23 * ak3[i] + d24 * ak4[i] + d25 * ak5[i];
                    cont[i + n * 3] = d31 * ak1[i] + d32 * ak2[i] + d33 * ak3[i] + d34 * ak4[i] + d35 * ak5[i];
                }
            }

            ++(nstep);

            hold = h;

            //
            // Error estimation
            //
            err = Error(n, h, y, ynew, ak6);
            
            //
            // Is the error small enough ?
            //
            if (controller.Success(err, posneg, ref h))
            {
                for (i = 0; i < n; ++i)
                {
                    y[i] = ynew[i];
                }
                conros_.xold = x;

                x += hold;

                /*
                if (iout != 0)
                {
                    for (i = 0; i < n; ++i)
                    {
                        cont[conros_.n + i] = y[i];
                    }
                    conros_.h = h;
                    //solout(naccpt + 1, conros_.xold, x, y, cont, lrc, n, irtrn);
                }
                //*/

                goto L1;
            }
            else
            {
                // Step is rejected
                last = false;

                goto L2;
            }
        }
        
        public double Error(int n, double h, double[] y, double[] ynew, double[] yerr)
        {
            double temp, sk, err = 0.0;

            for (int i = 0; i < n; i++)
            {
                //if (itol == 0)
                {
                    sk = atol[0] + rtol[0] * Math.Max(Math.Abs(y[i]), Math.Abs(ynew[i]));
                }
                //else
                //{
                //    sk = atol[i] + rtol[i] * Math.Max(Math.Abs(y[i]), Math.Abs(ynew[i]));
                //}

                temp = yerr[i] / sk;
                err += temp * temp;
            }

            return Math.Sqrt(err / n);
        }

        // This function can be used for continuous output in connection
        // with the output-subroutine for RODAS. It provides an
        // approximation to the i-th component of the solution at x.
        double Interpolate(int i, double x, double[] cont)
        {
            // Local variables
            double s = (x - conros_.xold) / conros_.h;

            return cont[i] * (1 - s) + s * (cont[i + conros_.n] + (1 - s) * (cont[i + (conros_.n << 1)] + s * cont[i + conros_.n * 3]));
        }
        
        int rocoe_(int meth)
        {
            double bet2p, bet3p, bet4p;

            switch (meth)
            {
                case 1:
                    c2 = 0.386;
                    c3 = 0.21;
                    c4 = 0.63;
                    bet2p = 0.0317;
                    bet3p = 0.0635;
                    bet4p = 0.3438;
                    d1 = 0.25;
                    d2 = -0.1043;
                    d3 = 0.1035;
                    d4 = -0.03620000000000023;
                    a21 = 1.544;
                    a31 = 0.9466785280815826;
                    a32 = 0.2557011698983284;
                    a41 = 3.314825187068521;
                    a42 = 2.896124015972201;
                    a43 = 0.9986419139977817;
                    a51 = 1.221224509226641;
                    a52 = 6.019134481288629;
                    a53 = 12.53708332932087;
                    a54 = -0.687886036105895;
                    c21 = -5.6688;
                    c31 = -2.430093356833875;
                    c32 = -0.2063599157091915;
                    c41 = -0.1073529058151375;
                    c42 = -9.594562251023355;
                    c43 = -20.47028614809616;
                    c51 = 7.496443313967647;
                    c52 = -10.24680431464352;
                    c53 = -33.99990352819905;
                    c54 = 11.7089089320616;
                    c61 = 8.083246795921522;
                    c62 = -7.981132988064893;
                    c63 = -31.52159432874371;
                    c64 = 16.31930543123136;
                    c65 = -6.058818238834054;
                    gamma = 0.25;
                    d21 = 10.12623508344586;
                    d22 = -7.487995877610167;
                    d23 = -34.80091861555747;
                    d24 = -7.992771707568823;
                    d25 = 1.025137723295662;
                    d31 = -0.6762803392801253;
                    d32 = 6.087714651680015;
                    d33 = 16.43084320892478;
                    d34 = 24.76722511418386;
                    d35 = -6.594389125716872;
                    break;
                case 2:
                    c2 = 0.3507221;
                    c3 = 0.2557041;
                    c4 = 0.681779;
                    bet2p = 0.0317;
                    bet3p = 0.0047369;
                    bet4p = 0.3438;
                    d1 = 0.25;
                    d2 = -0.06902209999999998;
                    d3 = -9.671999999999459e-4;
                    d4 = -0.08797900000000025;
                    a21 = 1.4028884;
                    a31 = 0.6581212688557198;
                    a32 = -1.320936088384301;
                    a41 = 7.131197445744498;
                    a42 = 16.02964143958207;
                    a43 = -5.561572550509766;
                    a51 = 22.73885722420363;
                    a52 = 67.38147284535289;
                    a53 = -31.2187749303856;
                    a54 = 0.7285641833203814;
                    c21 = -5.1043536;
                    c31 = -2.899967805418783;
                    c32 = 4.040399359702244;
                    c41 = -32.64449927841361;
                    c42 = -99.35311008728094;
                    c43 = 49.99119122405989;
                    c51 = -76.46023087151691;
                    c52 = -278.5942120829058;
                    c53 = 153.9294840910643;
                    c54 = 10.97101866258358;
                    c61 = -76.29701586804983;
                    c62 = -294.2795630511232;
                    c63 = 162.0029695867566;
                    c64 = 23.6516690309527;
                    c65 = -7.652977706771382;
                    gamma = 0.25;
                    d21 = -38.71940424117216;
                    d22 = -135.8025833007622;
                    d23 = 64.51068857505875;
                    d24 = -4.192663174613162;
                    d25 = -2.53193205033506;
                    d31 = -14.99268484949843;
                    d32 = -76.30242396627033;
                    d33 = 58.65928432851416;
                    d34 = 16.61359034616402;
                    d35 = -0.6758691794084156;
                    break;
                case 3:
                    // Coefficients for RODAS with order 4 for linear parabolic problems
                    // Gerd Steinebach (1993)
                    gamma = 0.25;
                    c2 = gamma * 3.0;
                    c3 = 0.21;
                    c4 = 0.63;
                    bet2p = 0.0;
                    bet3p = c3 * c3 * (c3 / 6.0 - gamma / 2.0) / (gamma * gamma);
                    bet4p = 0.3438;
                    d1 = 0.25;
                    d2 = -0.5;
                    d3 = -0.023504;
                    d4 = -0.0362;
                    a21 = 3.0;
                    a31 = 1.831036793486759;
                    a32 = 0.4955183967433795;
                    a41 = 2.304376582692669;
                    a42 = -0.05249275245743001;
                    a43 = -1.176798761832782;
                    a51 = -7.170454962423024;
                    a52 = -4.741636671481785;
                    a53 = -16.31002631330971;
                    a54 = -1.062004044111401;
                    c21 = -12.0;
                    c31 = -8.791795173947035;
                    c32 = -2.207865586973518;
                    c41 = 10.81793056857153;
                    c42 = 6.780270611428266;
                    c43 = 19.5348594464241;
                    c51 = 34.19095006749676;
                    c52 = 15.49671153725963;
                    c53 = 54.7476087596413;
                    c54 = 14.16005392148534;
                    c61 = 34.62605830930532;
                    c62 = 15.30084976114473;
                    c63 = 56.99955578662667;
                    c64 = 18.40807009793095;
                    c65 = -5.714285714285717;

                    d21 = 25.09876703708589;
                    d22 = 11.62013104361867;
                    d23 = 28.49148307714626;
                    d24 = -5.664021568594133;
                    d25 = 0.0;
                    d31 = 1.638054557396973;
                    d32 = -0.7373619806678748;
                    d33 = 8.47791821923899;
                    d34 = 15.9925314877952;
                    d35 = -1.882352941176471;
                    break;
            }
            return 0;
        }

        private void Factorize(int n, ReusableLU lu, DenseMatrix tmp, DenseMatrix jac, DenseMatrix mas, double fac)
        {
            var a = tmp.Values;
            var b = jac.Values;

            bool dae = mas != null; // Not tested.

            if (!dae)
            {
                // M = identity, Jacobian a full matrix
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        a[i * n + j] = -b[i * n + j];
                    }

                    a[i * n + i] += fac;
                }
            }
            else
            {
                // M is a full matrix, Jacobian a full matrix
                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        a[i * n + j] = mas.At(i, j) * fac - b[i * n + j];
                    }
                }
            }

            lu.Compute(tmp);
        }

        void Solve(int n, ReusableLU lu, DenseMatrix mas, double[] dy, double[] ak, double[] fx, double[] ynew, double hd, bool stage1)
        {
            if (hd == 0.0)
            {
                for (int i = 0; i < n; i++)
                {
                    ak[i] = dy[i];
                }
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    ak[i] = dy[i] + hd * fx[i];
                }
            }

            bool dae = mas != null; // Not tested.

            if (!dae)
            {
                // M = identity, Jacobian a full matrix
                if (stage1)
                {
                    for (int i = 0; i < n; i++)
                    {
                        ak[i] += ynew[i];
                    }
                }
            }
            else
            {
                // M is a full matrix, Jacobian a full matrix
                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < n; ++j)
                    {
                        sum += mas.At(i, j) * ynew[j];
                    }
                    ak[i] += sum;
                }
            }

            lu.Solve(ak, ak);
        }
    }
}
