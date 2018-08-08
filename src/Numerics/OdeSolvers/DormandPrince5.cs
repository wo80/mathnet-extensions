
namespace MathNet.Numerics.OdeSolvers
{
    using System;

    /// <summary>
    /// Numerical solution of a system of first order ordinary differential equations y'=f(x,y).
    /// This is an explicit runge-kutta method of order (4)5 due to Dormand & Prince (with stepsize
    /// control and dense output)
    ///
    /// Authors: E. Hairer and G. Wanner
    ///
    /// This code is described in:
    ///         E. Hairer, S.P. Norsett and G. Wanner
    ///         Solving Ordinary Differential Equations I. Nonstiff Problems (2nd edition)
    ///         Springer-Verlag (1993)
    /// </summary>
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
        
        IErrorController controller;
        StiffnessChecker stiff;

        int n;

        Action<double, double[], double[]> fcn;

        double told, dtold;
        double[] dxdt, dxdtnew, xout, xtemp, xerr;
        double[] k2, k3, k4, k5, k6;
        double[] r;

        double rtol, atol;

        /// <summary>
        /// Initializes a new instance of the <see cref="DormandPrince5"/> class.
        /// </summary>
        /// <param name="n">Dimension of the system.</param>
        /// <param name="controller">The error controller.</param>
        /// <param name="stiff">The stiffness detector.</param>
        public DormandPrince5(int n, IErrorController controller, StiffnessChecker stiff)
        {
            this.n = n;

            this.controller = controller;
            this.stiff = stiff;

            xout = new double[n];
            xtemp = new double[n];
            xerr = new double[n];
            dxdt = new double[n];
            dxdtnew = new double[n];

            k2 = new double[n];
            k3 = new double[n];
            k4 = new double[n];
            k5 = new double[n];
            k6 = new double[n];

            r = new double[5 * n];
        }

        double sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="fcn">Function computing the value of f(t,x).</param>
        /// <param name="t">Initial t-value.</param>
        /// <param name="x">Initial values for x.</param>
        /// <param name="tend">Final t-value (<c>tend-t</c> may be positive or negative).</param>
        /// <param name="rtol">Relative error tolerance.</param>
        /// <param name="atol">Absolute error tolerance.</param>
        /// <returns>Number of computed steps.</returns>
        public int Integrate(Action<double, double[], double[]> fcn, double t, double[] x, double tend,
            double rtol, double atol)
        {
            this.fcn = fcn;
            
            // NMAX the maximal number of steps
            int nmax = 100000;

            // Initial step size
            double h = 0.0;
            
            this.rtol = rtol;
            this.atol = atol;

            fcn(t, x, dxdt);

            if (h == 0.0)
            {
                double hmax = tend - t;
                double posneg = sign(1.0, hmax);
                h = AdaptiveIntegrator.Initialize(fcn, 5, t, x, tend, posneg, k2, k3, dxdt, Math.Abs(hmax), rtol, atol);
            }

            // Call to core integrator
            int nstep = Integrate(t, x, tend, h, nmax);

            if (nstep >= nmax)
            {
                Console.WriteLine(" More than NMAX =" + nmax + " steps are needed");
            }

            return nstep;
        }
        
        // Core integrator for DOPRI5
        int Integrate(double t, double[] x, double tend, double dt, int nmax)
        {
            int nstep = 0;
            int irtrn = 0;
            double posneg = sign(1.0, tend - t);

            // UROUND smallest number satisfying 1.0 + UROUND > 1.0
            double uround = Precision.DoublePrecision;

            // Initial preparations
            bool last = false;

            // Basic integration step
            while (nstep < nmax)
            {
                if (Math.Abs(dt) * 0.1 <= Math.Abs(t) * uround)
                {
                    throw new NumericalBreakdownException("Step size too small, h=" + dt);
                }

                if ((t + dt * 1.01 - tend) * posneg > 0.0)
                {
                    dt = tend - t;
                    last = true;
                }

                nstep++;

                if (irtrn > 1)
                {
                    fcn(t, x, dxdt);
                }

                Step(t, dt, x);

                // Error estimation
                double err = Error(dt, x);

                dtold = dt;
                told = t;

                // Computation of HNEW
                if (controller.Success(err, posneg, ref dt))
                {
                    // Stiffness detection
                    if (!stiff.Check(controller.Accepted, dtold, dxdtnew, k6, xout, xtemp))
                    {
                        throw new Exception(" The problem seems to become stiff at t = " + t);
                    }

                    //if (dense)
                    {
                        PrepareInterpolation(dtold, x);
                    }

                    for (int i = 0; i < n; ++i)
                    {
                        dxdt[i] = dxdtnew[i];
                        x[i] = xout[i];
                    }

                    t = t + dtold;

                    // Normal exit
                    if (last)
                    {
                        return nstep;
                    }
                }
                else
                {
                    last = false;
                }
            }

            return nstep;
        }

        private void PrepareInterpolation(double dt, double[] x)
        {
            int n = this.n;

            for (int i = 0; i < n; ++i)
            {
                double dx = xout[i] - x[i];
                double bspl = dt * dxdt[i] - dx;

                r[i] = x[i];
                r[n + i] = dx;
                r[n * 2 + i] = bspl;
                r[n * 3 + i] = -(dt) * dxdtnew[i] + dx - bspl;
                r[n * 4 + i] = dt * (d1 * dxdt[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * dxdtnew[i]);
            }
        }

        private void Step(double t, double dt, double[] x)
        {
            int i;

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
                xtemp[i] = x[i] + dt * (a61 * dxdt[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            }

            double th = t + dt;
            fcn(th, xtemp, k6);
            for (i = 0; i < n; ++i)
            {
                xout[i] = x[i] + dt * (a71 * dxdt[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            }

            fcn(th, xout, dxdtnew);
            for (i = 0; i < n; ++i)
            {
                xerr[i] = dt * (e1 * dxdt[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * dxdtnew[i]);
            }
        }

        double Error(double dt, double[] x)
        {
            double err = 0.0, sk, temp;

            for (int i = 0; i < n; ++i)
            {
                sk = atol + rtol * Math.Max(Math.Abs(x[i]), Math.Abs(xtemp[i]));

                temp = xerr[i] / sk;
                err += temp * temp;
            }

            return Math.Sqrt(err / n);
        }

        /// <summary>
        /// Provides an approximation to the i-th component of the solution at t.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="t"></param>
        /// <param name="con"></param>
        /// <returns></returns>
        public double Interpolate(int i, double t, double[] con)
        {
            int n = this.n;

            double s = (t - told) / dtold;
            double s1 = 1.0 - s;

            return con[i] + s * (con[n + i] + s1 * (con[2 * n + i] + s * (con[3 * n + i] + s1 * con[4 * n + i])));
        }
    }
}
