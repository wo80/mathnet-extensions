// Based on Fortran code RODAS
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;

    /// <summary>
    /// Numerical solution of a stiff (or differential algebraic) system of first
    /// order ordinary differential equations. This is an embedded Rosenbrock
    /// method of order (3)4 with step size control and dense output.
    /// 
    /// Authors: E. Hairer and G. Wanner
    ///
    /// This code is part of the book:
    ///         E. Hairer and G. Wanner
    ///         Solving Ordinary Differential Equations II.
    ///         Stiff and differential-algebraic problems. (2nd edition)
    ///         Springer-Verlag (1996)
    /// </summary>
    public class Rosenbrock4
    {
        RosenbrockErrorController controller;

        int n;

        Action<double, double[], double[]> fcn;

        bool autonomous;
        Action<double, double[], DenseMatrix> jac;
        Action<double, double[], double[]> dfx;
        DenseMatrix mas;

        double rtol, atol;

        public int ndec, nsol; // TODO: remove
        
        double told, hold;

        // Workspace

        double[] ynew;
        double[] dy, dy1;
        double[] ak1, ak2, ak3, ak4, ak5, ak6;
        double[] fx, cont;

        ReusableLU lu;
        DenseMatrix fjac;
        DenseMatrix mjac;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="n">dimension of the system</param>
        /// <param name="fcn">subroutine computing the value of f(x,y)</param>
        /// <param name="autonomous">f(x,y) independent of x (autonomous)</param>
        /// <param name="jac">subroutine which computes the partial derivatives of f(x,y) with respect to y</param>
        /// <param name="dfx">subroutine which computes the partial derivatives of f(x,y) with respect to x</param>
        /// <param name="mas">the mass-matrix m.</param>
        /// <param name="rtol">relative error tolerances</param>
        /// <param name="atol">absolute error tolerances</param>
        /// <param name="controller"></param>
        public Rosenbrock4(int n, Action<double, double[], double[]> fcn,
            bool autonomous, Action<double, double[], DenseMatrix> jac,
            Action<double, double[], double[]> dfx, DenseMatrix mas,
            double rtol, double atol,
            RosenbrockErrorController controller)
        {
            this.n = n;
            this.fcn = fcn;
            this.autonomous = autonomous;
            this.jac = jac;
            this.dfx = dfx;
            this.mas = mas;
            this.rtol = rtol;
            this.atol = atol;

            this.controller = controller;

            ynew = new double[n];
            dy1 = new double[n];
            dy = new double[n];
            ak1 = new double[n];
            ak2 = new double[n];
            ak3 = new double[n];
            ak4 = new double[n];
            ak5 = new double[n];
            ak6 = new double[n];
            fx = new double[n];
            cont = new double[4 * n];

            lu = new ReusableLU(n);
            fjac = new DenseMatrix(n);
            mjac = new DenseMatrix(n);
        }

        double sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c21, c31,
            c32, c41, c42, c43, c51, c52, c53, c54, c61, c62, c63, c64, c65,
            d21, d22, d23, d24, d25, d31, d32, d33, d34, d35;
        double c2, c3, c4, d1, d2, d3, d4, gamma;

        /**
         *     output parameters
         *     -----------------
         *     x           x-value where the solution is computed
         *                 (after successful return x=xend)
         *
         *     y(n)        solution at x
         *
         *     h           predicted step size of the last accepted step
         *
         * 
         *     x           initial x-value
         *
         *     y(n)        initial values for y
         *
         *     xend        final x-value (xend-x may be positive or negative)
         *
         *     h           initial step size guess;
         *                 for stiff equations with initial transient,
         *                 h=1.d0/(norm of f'), usually 1.d-2 or 1.d-3, is good.
         *                 this choice is not very important, the code quickly
         *                 adapts its step size (if h=0.d0, the code puts h=1.d-6).
         */
        public int Integrate(double t, double[] y, double tend, double h)
        {
            double eps = Precision.DoublePrecision;

            // NMAX , THE MAXIMAL NUMBER OF STEPS
            int nmax = 100000;

            // METH   COEFFICIENTS OF THE METHOD
            //   1  method (see book, page 452)
            //   2  same method with different parameters
            //   3  method with coeff. of gerd steinebach

            int meth = 1;

            if (meth <= 0 || meth >= 4)
            {
                Console.WriteLine(" Curious input IWORK(2)=", meth);
                return -1;
            }

            // Check if tolerances are OK.
            if (atol <= 0.0 || rtol <= eps * 10.0)
            {
                Console.WriteLine(" Tolerances are too small");
                return -1;
            }
            
            ndec = nsol = 0;

            // Call to core integrator
            int nstep = Integrate(t, y, tend, h, nmax, meth, true);

            if (nstep > nmax)
            {
                Console.WriteLine(" More than NMAX =" + nmax + "steps are needed");
            }

            return nstep;
        }

        // ... and here is the core integrator
        private int Integrate(double t, double[] y, double tend, double h, int nmax, int meth, bool dense)
        {
            int n = this.n;

            var fcn = this.fcn;
            var jac = this.jac;
            var dfx = this.dfx;
            var mas = this.mas;

            bool last;
            bool autonomous = this.autonomous;
            
            int nstep = 0;

            double posneg, eps = Precision.DoublePrecision;
            
            // Set the parameters of the method
            LoadMethod(meth);

            // Initial preparations
            posneg = sign(1.0, tend - t);

            if (Math.Abs(h) <= eps * 10.0)
            {
                h = 1e-6;
            }
            h = Math.Min(Math.Abs(h), Math.Abs(tend - t));

            h = sign(h, posneg);
            last = false;
            
            // Basic integration step
            while (nstep < nmax)
            {
                if (last)
                {
                    return nstep;
                }

                if ((t + h * 1.0001 - tend) * posneg >= 0.0)
                {
                    h = tend - t;
                    last = true;
                }

                Step(ref h, ref t, y, posneg, dense, ref nstep);
            }

            return nstep;
        }

        private void Step(ref double h, ref double t, double[] y, double posneg, bool dense, ref int nstep)
        {
            int n = this.n;

            var fcn = this.fcn;
            var jac = this.jac;
            var dfx = this.dfx;
            var mas = this.mas;

            bool autonomous = this.autonomous;
            
            int nsing = 0;

            double err, eps = Precision.DoublePrecision;

            fcn(t, y, dy1);

            // Computation of the Jacobian
            if (jac == null)
            {
                // Compute Jacobian matrix numerically
                for (int i = 0; i < n; ++i)
                {
                    double ysafe = y[i];
                    double delta = Math.Sqrt(eps * Math.Max(1e-5, Math.Abs(ysafe)));
                    y[i] = ysafe + delta;
                    fcn(t, y, ak1);
                    for (int j = 0; j < n; ++j)
                    {
                        fjac.At(i, j, (ak1[j] - dy1[j]) / delta);
                    }
                    y[i] = ysafe;
                }
            }
            else
            {
                jac(t, y, fjac);
            }

            // Computation of derivative with respect to x
            if (!autonomous)
            {
                if (dfx == null)
                {
                    // Compute the derivative with respect to x numerically
                    double delt = Math.Sqrt(eps * Math.Max(1e-5, Math.Abs(t)));
                    fcn(t + delt, y, ak1);
                    for (int i = 0; i < n; ++i)
                    {
                        fx[i] = (ak1[i] - dy1[i]) / delt;
                    }
                }
                else
                {
                    dfx(t, y, fx);
                }
            }

            // Loop until success or stepsize too small.
            while (true)
            {
                if (Math.Abs(h) * 0.1 <= Math.Abs(t) * eps)
                {
                    throw new NumericalBreakdownException("Step size too small, h=" + h);
                }

                Factorize(n, lu, mjac, fjac, mas, 1.0 / (h * gamma));

                if (lu.Determinant < eps) // TODO: remove?
                {
                    // Singular matrix
                    if (++nsing >= 5)
                    {
                        throw new Exception("Matrix is repeatedly singular.");
                        //return -4;
                    }

                    h *= 0.5;
                    //reject = true; // TODO: controller.Reject();
                    //last = false;

                    continue;
                }

                Step(h, t, y);

                nstep++;

                hold = h;
                
                // Error estimation
                err = Error(n, h, y, ynew, ak6);
                
                // Is the error small enough ?
                if (controller.Success(err, posneg, ref h))
                {
                    if (dense)
                    {
                        PrepareInterpolation(y);
                    }

                    for (int i = 0; i < n; ++i)
                    {
                        y[i] = ynew[i];
                    }

                    told = t; // Used for interpolation.

                    t += hold;
                    
                    hold = h; // Used for interpolation.

                    return;
                }

                // Step is rejected
                //last = false;
            }
        }

        private void PrepareInterpolation(double[] y)
        {
            for (int i = 0; i < n; ++i)
            {
                cont[i] = y[i];
                cont[i + n] = ynew[i];
                cont[i + n * 2] = d21 * ak1[i] + d22 * ak2[i] + d23 * ak3[i] + d24 * ak4[i] + d25 * ak5[i];
                cont[i + n * 3] = d31 * ak1[i] + d32 * ak2[i] + d33 * ak3[i] + d34 * ak4[i] + d35 * ak5[i];
            }
        }

        private void Step(double h, double x, double[] y)
        {
            int i;

            double hd1 = 0, hd2 = 0, hd3 = 0, hd4 = 0;
            double hc21, hc31, hc32, hc41, hc42, hc43, hc51, hc52, hc53, hc54, hc61, hc62, hc63, hc64, hc65;

            var fcn = this.fcn;
            var mas = this.mas;

            // Compute the stages

            ndec++;

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

            if (!autonomous)
            {
                hd1 = h * d1;
                hd2 = h * d2;
                hd3 = h * d3;
                hd4 = h * d4;
            }

            Solve(n, lu, mas, dy1, ak1, fx, ynew, hd1, false);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a21 * ak1[i];
            }

            fcn(x + c2 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc21 * ak1[i];
            }
            Solve(n, lu, mas, dy, ak2, fx, ynew, hd2, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a31 * ak1[i] + a32 * ak2[i];
            }

            fcn(x + c3 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc31 * ak1[i] + hc32 * ak2[i];
            }
            Solve(n, lu, mas, dy, ak3, fx, ynew, hd3, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a41 * ak1[i] + a42 * ak2[i] + a43 * ak3[i];
            }

            fcn(x + c4 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc41 * ak1[i] + hc42 * ak2[i] + hc43 * ak3[i];
            }
            Solve(n, lu, mas, dy, ak4, fx, ynew, hd4, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a51 * ak1[i] + a52 * ak2[i] + a53 * ak3[i]
                    + a54 * ak4[i];
            }

            fcn(x + h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ak6[i] = hc52 * ak2[i] + hc54 * ak4[i] + hc51 * ak1[i] + hc53 * ak3[i];
            }
            Solve(n, lu, mas, dy, ak5, fx, ak6, 0.0, true);

            // Embedded solution
            for (i = 0; i < n; i++)
            {
                ynew[i] += ak5[i];
            }

            fcn(x + h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                cont[i] = hc61 * ak1[i] + hc62 * ak2[i] + hc65 * ak5[i] + hc64 * ak4[i] + hc63 * ak3[i];
            }
            Solve(n, lu, mas, dy, ak6, fx, cont, 0.0, true);

            // New solution
            for (i = 0; i < n; i++)
            {
                ynew[i] += ak6[i];
            }

            nsol += 6;
        }

        public double Error(int n, double h, double[] y, double[] ynew, double[] yerr)
        {
            double temp, sk, err = 0.0;

            for (int i = 0; i < n; i++)
            {
                sk = atol + rtol * Math.Max(Math.Abs(y[i]), Math.Abs(ynew[i]));

                temp = yerr[i] / sk;
                err += temp * temp;
            }

            return Math.Sqrt(err / n);
        }

        // This function can be used for continuous output in connection
        // with the output-subroutine for RODAS. It provides an
        // approximation to the i-th component of the solution at x.
        public double Interpolate(int i, double t)
        {
            // Local variables
            double s = (t - told) / hold;

            return cont[i] * (1 - s) + s * (cont[i + n] + (1 - s) * (cont[i + 2 * n] + s * cont[i + 3 * n]));
        }

        private int LoadMethod(int meth)
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
