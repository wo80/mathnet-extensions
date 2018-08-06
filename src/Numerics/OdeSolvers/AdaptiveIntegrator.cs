
namespace MathNet.Numerics.OdeSolvers
{
    using System;

    public class AdaptiveIntegrator
    {
        public static double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }
        
        /* ----  Computation of an initial step size guess */

        public static double Initialize(Action<double, double[], double[]> fcn, int iord, double t, double[] x, double tend, double posneg,
            double[] f0, double[] f1, double[] y1, double hmax, double[] rtol, double[] atol, int itol)
        {
            double temp;

            int n = x.Length;

            /* Local variables */
            double h;
            double h1, sk, dnf, dny, der2, der12, atoli, rtoli;

            // Compute a first guess for explicit euler as
            //   H = 0.01 * NORM (Y0) / NORM (F0)
            // The increment for explicit euler is small
            // compared to the solution

            /* Function Body */
            dnf = 0.0;
            dny = 0.0;
            atoli = atol[0];
            rtoli = rtol[0];
            if (itol == 0)
            {
                for (int i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Abs(x[i]);

                    temp = f0[i] / sk;
                    dnf += temp * temp;

                    temp = x[i] / sk;
                    dny += temp * temp;
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Abs(x[i]);

                    temp = f0[i] / sk;
                    dnf += temp * temp;

                    temp = x[i] / sk;
                    dny += temp * temp;
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

            // Perform an explicit euler step
            for (int i = 0; i < n; ++i)
            {
                y1[i] = x[i] + h * f0[i];
            }

            fcn(t + h, y1, f1);

            // Estimate the second derivative of the solution
            der2 = 0.0;
            if (itol == 0)
            {
                for (int i = 0; i < n; ++i)
                {
                    sk = atoli + rtoli * Math.Abs(x[i]);

                    temp = (f1[i] - f0[i]) / sk;
                    der2 += temp * temp;
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    sk = atol[i] + rtol[i] * Math.Abs(x[i]);

                    temp = (f1[i] - f0[i]) / sk;
                    der2 += temp * temp;
                }
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
                h1 = Math.Pow(0.01 / der12, 1.0 / iord);
            }

            h = Math.Min(Math.Min(Math.Abs(h) * 100, h1), hmax);

            return d_sign(h, posneg);
        }
    }
}
