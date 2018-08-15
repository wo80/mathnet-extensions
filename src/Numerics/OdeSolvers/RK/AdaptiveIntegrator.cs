
namespace MathNet.Numerics.OdeSolvers.RK
{
    using System;

    public class AdaptiveIntegrator
    {
        private static double sign(double a, double b)
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
        public static int Integrate(int n, IRungeKuttaStepper rk, Action<double, double[], double[]> fcn, double t, double[] x, double tend, double rtol, double atol)
        {
            // NMAX the maximal number of steps
            int nmax = 100000;

            // Initial step size
            double h = 0.0;

            var dxdt = new double[n];

            fcn(t, x, dxdt);

            double hmax = tend - t;
            double posneg = sign(1.0, hmax);

            if (h == 0.0)
            {
                h = Initialize(n, fcn, rk.Order, t, x, dxdt, posneg, Math.Abs(hmax), rtol, atol);
            }

            int step = 0;

            // EPS smallest number satisfying 1.0 + EPS > 1.0
            double eps = Precision.DoublePrecision;

            // Initial preparations
            bool last = false;

            // Basic integration step
            while (step < nmax)
            {
                if (Math.Abs(h) * 0.1 <= Math.Abs(t) * eps)
                {
                    throw new NumericalBreakdownException("Step size too small, h=" + h);
                }

                if ((t + h * 1.01 - tend) * posneg > 0.0)
                {
                    h = tend - t;
                    last = true;
                }
                
                // Call to core integrator
                rk.Integrate(ref t, ref h, x, ref step, nmax, posneg, false);

                Console.WriteLine("Y = [{0}, {1}]", x[0], x[1]);

                // Normal exit
                if (last)
                {
                    break;
                }
            }

            if (step >= nmax)
            {
                Console.WriteLine(" More than NMAX =" + nmax + " steps are needed");
            }

            return step;
        }

        /// <summary>
        /// Computation of an initial step size guess
        /// </summary>
        /// <param name="n"></param>
        /// <param name="fcn"></param>
        /// <param name="order"></param>
        /// <param name="t">The initial time.</param>
        /// <param name="x">The initial value.</param>
        /// <param name="dxdt">The initial function value.</param>
        /// <param name="posneg"></param>
        /// <param name="hmax"></param>
        /// <param name="rtol"></param>
        /// <param name="atol"></param>
        /// <returns></returns>
        public static double Initialize(int n, Action<double, double[], double[]> fcn, int order, double t, double[] x, double[] dxdt,
            double posneg, double hmax, double rtol, double atol)
        {
            double temp;
            double sk, h, h1, der2, der12;

            var f0 = new double[n];
            var f1 = new double[n];

            // Compute a first guess for explicit euler as
            //   H = 0.01 * NORM (Y0) / NORM (F0)
            // The increment for explicit euler is small compared to the solution.
            
            double dnf = 0.0;
            double dny = 0.0;

            for (int i = 0; i < n; ++i)
            {
                sk = atol + rtol * Math.Abs(x[i]);

                temp = dxdt[i] / sk;
                dnf += temp * temp;

                temp = x[i] / sk;
                dny += temp * temp;
            }

            if (dnf <= 1.0e-10 || dny <= 1.0e-10)
            {
                h = 1.0e-6;
            }
            else
            {
                h = Math.Sqrt(dny / dnf) * 0.01;
            }

            h = Math.Min(h, hmax);
            h = sign(h, posneg);

            // Perform an explicit euler step
            for (int i = 0; i < n; ++i)
            {
                f0[i] = x[i] + h * dxdt[i];
            }

            fcn(t + h, f0, f1);

            // Estimate the second derivative of the solution
            der2 = 0.0;

            for (int i = 0; i < n; ++i)
            {
                sk = atol + rtol * Math.Abs(x[i]);

                temp = (f1[i] - f0[i]) / sk;
                der2 += temp * temp;
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
                h1 = Math.Pow(0.01 / der12, 1.0 / order);
            }

            h = Math.Min(Math.Min(Math.Abs(h) * 100, h1), hmax);

            return sign(h, posneg);
        }
    }
}
