
namespace MathNet.Numerics.Integration
{
    using System;

    /// <summary>
    /// Adaptive integrator using Gauss-Lobatto method with a Kronrod extension.
    /// </summary>
    /// <remarks>
    /// Reference:
    /// W. Gander and W. Gautschi: Adaptive Quadrature - Revisited
    /// 
    /// Matlab code available on https://www.inf.ethz.ch/personal/gander/
    /// </remarks>
    public class AdaptiveGaussLobatto
    {
        const double Epsilon = 2.22044604925031308e-16;

        private static readonly double alpha = Math.Sqrt(2.0 / 3.0);
        private static readonly double beta = 1.0 / Math.Sqrt(5.0);

        private const double x1 = 0.942882415695479719056351758432;
        private const double x2 = 0.641853342345781305781235541329;
        private const double x3 = 0.236383199662149880282223773493;

        //Abscissas for Gauss-Lobatto-Kronrod quadrature.
        private static readonly double[] x = { 0.0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1 };

        private const double a1 = 0.0158271919734801830871699867333;
        private const double a2 = 0.0942738402188500455312825050771;
        private const double a3 = 0.1550719873365853962536359798021;
        private const double a4 = 0.1888215739601824544200053393730;
        private const double a5 = 0.1997734052268585267920680220665;
        private const double a6 = 0.2249264653333395270160176879964;
        private const double a7 = 0.2426110719014077337996409579033;

        double tol;
        
        // Maximum number of function evaluations.
        int count, limit = 100000;

        /// <summary>
        /// Indicates wether integration stopped because of loss of precision.
        /// </summary>
        public bool OutOfTolerance
        {
            get;
            private set;
        }

        /// <summary>
        /// Adaptive Lobatto integration.
        /// </summary>
        /// <param name="tol">Desired tolerance.</param>
        /// <param name="limit">The maximum number of function evaluations.</param>
        public AdaptiveGaussLobatto(double tolerance, int limit = 0)
        {
            this.tol = tolerance < Epsilon ? Epsilon : tolerance;
            this.limit = limit > 0 ? limit : 100000;
        }

        /// <summary>
        /// Adaptive approximation of the definite integral by the Gauss-Lobatto rule.
        /// </summary>
        /// <param name="func">The analytic smooth function to integrate.</param>
        /// <param name="a">Left endpoint of the integration interval (inclusive and finite).</param>
        /// <param name="b">Right endpoint of the integration interval (inclusive and finite).</param>
        /// <param name="tolerance">The expected accuracy of the approximation.</param>
        /// <param name="limit">The maximum number of function evaluations.</param>
        /// <returns>Approximation of the finite integral in the given interval.</returns>
        public static double Integrate(Func<double, double> func, double a, double b, double tolerance, int limit = 0)
        {
            var lobatto = new AdaptiveGaussLobatto(tolerance, limit);

            return lobatto.Integrate(func, a, b);
        }

        /// <summary>
        /// Adaptive approximation of the definite integral by the Gauss-Lobatto rule.
        /// </summary>
        public double Integrate(Func<double, double> func, double a, double b)
        {
            OutOfTolerance = false;

            double m, h, fa, fb, i1, i2, est, erri1, erri2, r;
            double[] y = new double[13];

            m = 0.5 * (a + b);
            h = 0.5 * (b - a);

            fa = y[0] = func(a);
            for (int i = 1; i < 12; i++)
            {
                y[i] = func(m + x[i] * h);
            }
            fb = y[12] = func(b);

            // 4-point Gauss-Lobatto formula.
            i2 = (h / 6) * (y[0] + y[12] + 5 * (y[4] + y[8]));

            // 7-point Kronrod extension.
            i1 = (h / 1470) * (77 * (y[0] + y[12]) + 432 * (y[2] + y[10]) + 625 * (y[4] + y[8]) + 672 * y[6]);
            
            // 13-point Kronrod extension.
            est = h * (a1 * (y[0] + y[12]) + a2 * (y[1] + y[11]) + a3 * (y[2] + y[10]) +
                a4 * (y[3] + y[9]) + a5 * (y[4] + y[8]) + a6 * (y[5] + y[7]) + a7 * y[6]);

            erri1 = Math.Abs(i1 - est);
            erri2 = Math.Abs(i2 - est);

            r = (erri2 != 0.0) ? erri1 / erri2 : 1.0;

            // Error of i1 will be sufficiently small, so we can increase tolerance.
            if (r > 0.0 && r < 1.0)
            {
                tol = tol / r;
            }

            est = Math.Abs(est) * tol / Epsilon;

            if (est == 0.0)
            {
                est = b - a;
            }

            count = 13;

            return AdaptLobatto(func, a, b, fa, fb, est);
        }

        double AdaptLobatto(Func<double, double> func, double a, double b, double fa, double fb, double est)
        {
            double m, h, mll, ml, mr, mrr, fmll, fml, fm, fmrr, fmr, i1, i2;

            m = 0.5 * (a + b);
            h = 0.5 * (b - a);

            mll = m - alpha * h;
            ml = m - beta * h;
            mr = m + beta * h;
            mrr = m + alpha * h;

            fmll = func(mll);
            fml = func(ml);
            fm = func(m);
            fmr = func(mr);
            fmrr = func(mrr);

            count += 5;

            // 4-point Gauss-Lobatto formula.
            i2 = h / 6 * (fa + fb + 5 * (fml + fmr));

            if (double.IsNaN(i2) || double.IsInfinity(i2))
            {
                return i2;
            }

            // 7-point Kronrod extension.
            i1 = h / 1470 * (77 * (fa + fb) + 432 * (fmll + fmrr) + 625 * (fml + fmr) + 672 * fm);

            if (double.IsNaN(i1) || double.IsInfinity(i1))
            {
                return i1;
            }

            // TODO: might get optimized by compiler. Should we explicitly check for machine epsilon?

            // Should we stop the recursion?
            if ((est + (i1 - i2)) == est || mll <= a || b <= mrr || count > limit)
            {
                if ((mll <= a || b <= mrr))
                {
                    // Interval contains no more machine numbers.
                    OutOfTolerance = true;
                }

                return i1;
            }
            else
            {
                // Subdivide interval.
                return AdaptLobatto(func, a, mll, fa, fmll, est) +
                    AdaptLobatto(func, mll, ml, fmll, fml, est) +
                    AdaptLobatto(func, ml, m, fml, fm, est) +
                    AdaptLobatto(func, m, mr, fm, fmr, est) +
                    AdaptLobatto(func, mr, mrr, fmr, fmrr, est) +
                    AdaptLobatto(func, mrr, b, fmrr, fb, est);
            }
        }
    }
}