
namespace MathNet.Numerics.Integration
{
    using System;

    /// <summary>
    /// Adaptive quadrature with mixed quadrature rules.
    /// </summary>
    /// <remarks>
    /// Application of Mixed Quadrature Rules in the Adaptive Quadrature Routine 
    /// Debasish Das, Rajani B. Dash 
    /// 
    /// Gen. Math. Notes, Vol. 18, No. 1, September, 2013, pp. 46-63 
    /// ISSN 2219-7184; Copyright © ICSRS Publication, 2013
    /// 
    /// Available online at http://www.geman.in
    /// </remarks>
    public class AdaptiveMixedDasDash
    {
        const double Epsilon = 2.22044604925031308e-16;

        const double o = 0.81649658092772615; // sqrt(2) / sqrt(3)
        const double p = 0.70710678118654746; // 1 / sqrt(2)
        const double q = 0.65465367070797709; // sqrt(3) / sqrt(7)
        const double r = 0.44721359549995793; // 1 / sqrt(5)

        double tol;

        int limit;
        
        /// <summary>
        /// Adaptive Lobatto integration.
        /// </summary>
        /// <param name="tolerance">Desired tolerance.</param>
        /// <param name="limit">The maximum depth of recursion.</param>
        public AdaptiveMixedDasDash(double tolerance, int limit = 0)
        {
            this.tol = tolerance < Epsilon ? Epsilon : tolerance;
            this.limit = limit > 0 ? limit : 10000;
        }

        /// <summary>
        /// Adaptive approximation of the definite integral by mixed Quadrature Rules.
        /// </summary>
        /// <param name="func">The analytic smooth function to integrate.</param>
        /// <param name="a">Left endpoint of the integration interval (inclusive and finite).</param>
        /// <param name="b">Right endpoint of the integration interval (inclusive and finite).</param>
        /// <param name="tolerance">The expected accuracy of the approximation.</param>
        /// <param name="limit">The maximum number of function evaluations.</param>
        /// <returns>Approximation of the finite integral in the given interval.</returns>
        public static double Integrate(Func<double, double> func, double a, double b, double tolerance, int limit = 0)
        {
            var mixed = new AdaptiveMixedDasDash(tolerance, limit);

            return mixed.Integrate(func, a, b);
        }

        public double Integrate(Func<double, double> func, double a, double b)
        {
            return IntegrateAdapt(func, a, b, this.tol);
        }

        /// <summary>
        /// Approximate integral by adaptive subdivision.
        /// </summary>
        private double IntegrateAdapt(Func<double, double> f, double a, double b, double tol)
        {
            double m, int1, int2;

            m = (a + b) / 2;

            // Low order approximation.
            int1 = REC_Integrand_R1(f, a, b);

            if (double.IsNaN(int1) || double.IsInfinity(int1))
            {
                return int1;
            }

            // High order approximation.
            int2 = REC_Integrand_R2(f, a, m) + REC_Integrand_R3(f, m, b);

            if (double.IsNaN(int2) || double.IsInfinity(int2))
            {
                return int2;
            }

            if (Math.Abs(int2 - int1) < tol / 2)
            {
                // Error criterion satisfied.
                return int2;
            }

            if (--this.limit < 0)
            {
                return int2;
            }

            // Error criterion not satisfied.
            return (IntegrateAdapt(f, a, m, tol) + IntegrateAdapt(f, m, b, tol));
        }

        /// <summary>
        /// Mixed quadrature rule RL4CC5 of precision 7 using convex combination of the Lobatto
        /// 4-point the and Clenshaw-Curtis 5-point rule (each of precision 5).
        /// </summary>
        private double REC_Integrand_R1(Func<double, double> f, double a, double b)
        {
            double h, m, RL4CC5;
            double x1, x2, x5, x6, x7, x10, x11;

            h = (b - a) / 2;
            m = (a + b) / 2;

            x1 = f(a);
            x2 = f(b);
            x5 = f(m - p * h);
            x6 = f(m + p * h);
            x7 = f(m);
            x10 = f(m - r * h);
            x11 = f(m + r * h);

            RL4CC5 = h / 630 * (57 * (x1 + x2) + 256 * (x5 + x6) + 125 * (x10 + x11) + 384 * x7);

            return RL4CC5;
        }

        /// <summary>
        /// Mixed quadrature rule RL4CC5L5 of precision 9 using linear combination of the Lobatto
        /// 5-point rule RL5 and the mixed quadrature rule RL4CC5 (each of precision 7).
        /// </summary>
        private double REC_Integrand_R2(Func<double, double> f, double a, double b)
        {
            double h, m;
            double x1, x2, x5, x6, x7, x8, x9, x10, x11, RL4CC5L5;

            h = (b - a) / 2;
            m = (a + b) / 2;

            x1 = f(a);
            x2 = f(b);
            x5 = f(m - p * h);
            x6 = f(m + p * h);
            x7 = f(m);
            x8 = f(m - q * h);
            x9 = f(m + q * h);
            x10 = f(m - r * h);
            x11 = f(m + r * h);

            RL4CC5L5 = h / 1890 * (129 * (x1 + x2) + 2560 * (x5 + x6)
                - 2401 * (x8 + x9) + 1250 * (x10 + x11) + 704 * x7);

            return RL4CC5L5;
        }

        /// <summary>
        /// Mixed quadrature rule RL4CC5L5KEL4 of precision 11 using linear combination of the mixed
        /// quadrature rule RL4CC5L5 and the Kronrod extension of the Lobatto 4-point rule (each of
        /// precision 9).
        /// </summary>
        private double REC_Integrand_R3(Func<double, double> f, double a, double b)
        {
            double h, m;
            double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, RL4CC5L5KEL4;

            h = (b - a) / 2;
            m = (a + b) / 2;

            x1 = f(a);
            x2 = f(b);
            x3 = f(m - o * h);
            x4 = f(m + o * h);
            x5 = f(m - p * h);
            x6 = f(m + p * h);
            x7 = f(m);
            x8 = f(m - q * h);
            x9 = f(m + q * h);
            x10 = f(m - r * h);
            x11 = f(m + r * h);

            RL4CC5L5KEL4 = h / 727650 * (35175 * (x1 + x2) + 268272 * (x3 + x4)
                - 250880 * (x5 + x6) + 235298 * (x8 + x9)
                + 265625 * (x10 + x11) + 348320 * x7);

            return RL4CC5L5KEL4;
        }
    }
}
