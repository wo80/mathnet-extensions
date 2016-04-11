
namespace MathNet.Numerics.LinearAlgebra.Double
{
    using System;

    /// <summary>
    /// DenseVector extensions methods.
    /// </summary>
    public static class DenseVectorExtensions
    {
        /// <summary>
        /// Clear this vector.
        /// </summary>
        public static void Clear(this DenseVector vector)
        {
            var x = (double[])vector;

            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Scale this vector.
        /// </summary>
        public static void Scale(this DenseVector vector, double alpha)
        {
            var vx = (double[])vector;

            int length = vx.Length;

            for (int i = 0; i < length; i++)
            {
                vx[i] *= alpha;
            }
        }

        /// <summary>
        /// Computes the L2 norm of a vector avoiding overflow.
        /// </summary>
        public static double L2NormRobust(this DenseVector x)
        {
            var vx = (double[])x;

            int length = vx.Length;

            double scale = 0.0, ssq = 1.0;

            for (int i = 0; i < length; ++i)
            {
                if (vx[i] != 0.0)
                {
                    double absxi = Math.Abs(vx[i]);
                    if (scale < absxi)
                    {
                        ssq = 1.0 + ssq * (scale / absxi) * (scale / absxi);
                        scale = absxi;
                    }
                    else
                    {
                        ssq += (absxi / scale) * (absxi / scale);
                    }
                }
            }

            return scale * Math.Sqrt(ssq);
        }
    }
}
