
namespace MathNet.Numerics.LinearAlgebra.Double
{
    using MathNet.Numerics.Threading;
    using System;

    /// <summary>
    /// Vector helper methods.
    /// </summary>
    public static class VectorHelper
    {
        /// <summary>
        /// Clone the given vector.
        /// </summary>
        public static DenseVector Clone(DenseVector vector)
        {
            return DenseVector.OfArray((double[])vector);
        }

        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(DenseVector src, DenseVector dst)
        {
            var x = (double[])src;
            var y = (double[])dst;

            Array.Copy(x, 0, y, 0, x.Length);
        }

        /// <summary>
        /// Set all vector elements to zero.
        /// </summary>
        public static void Clear(DenseVector vector)
        {
            var x = (double[])vector;

            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Scale given vector, x = alpha * x.
        /// </summary>
        public static void Scale(double alpha, DenseVector x)
        {
            var vx = (double[])x;

            int length = vx.Length;

            for (int i = 0; i < length; i++)
            {
                vx[i] *= alpha;
            }
        }

        /// <summary>
        /// Copy scaled vector, target = alpha * x.
        /// </summary>
        public static void Scale(double alpha, DenseVector x, DenseVector target)
        {
            var vx = (double[])x;
            var vt = (double[])target;

            int length = vx.Length;

            for (int i = 0; i < length; i++)
            {
                vt[i] = alpha * vx[i];
            }
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static double DotProduct(DenseVector x, DenseVector y)
        {
            var vx = (double[])x;
            var vy = (double[])y;

            int length = vx.Length;

            double result = 0.0;

            for (int i = 0; i < length; i++)
            {
                result += vx[i] * vy[i];
            }

            return result;
        }

        /// <summary>
        /// Compute sum of scaled vectors, target = alpha * x + y.
        /// </summary>
        public static void Add(double alpha, DenseVector x, DenseVector y, DenseVector target)
        {
            var vx = (double[])x;
            var vy = (double[])y;
            var vt = (double[])target;

            int length = vx.Length;

            CommonParallel.For(0, length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    vt[i] = alpha * vx[i] + vy[i];
                }
            });
        }

        /// <summary>
        /// Compute sum of scaled vectors, target = alpha * x + beta * y.
        /// </summary>
        public static void Add(double alpha, DenseVector x, double beta, DenseVector y, DenseVector target)
        {
            var vx = (double[])x;
            var vy = (double[])y;
            var vt = (double[])target;

            int length = vx.Length;

            CommonParallel.For(0, length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    vt[i] = alpha * vx[i] + beta * vy[i];
                }
            });
        }

        /// <summary>
        /// Compute sum of scaled vectors, target = alpha * x + beta * y + z.
        /// </summary>
        public static void Add(double alpha, DenseVector x, double beta, DenseVector y, DenseVector z,
            DenseVector target)
        {
            var vx = (double[])x;
            var vy = (double[])y;
            var vz = (double[])z;
            var vt = (double[])target;

            int length = vx.Length;

            CommonParallel.For(0, length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    vt[i] = alpha * vx[i] + beta * vy[i] + vz[i];
                }
            });
        }

        /// <summary>
        /// Compute sum of scaled vectors, target = alpha * x + beta * y + gamma * z.
        /// </summary>
        public static void Add(double alpha, DenseVector x, double beta, DenseVector y, 
            double gamma, DenseVector z, DenseVector target)
        {
            var vx = (double[])x;
            var vy = (double[])y;
            var vz = (double[])z;
            var vt = (double[])target;

            int length = vx.Length;

            CommonParallel.For(0, length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    vt[i] = alpha * vx[i] + beta * vy[i] + gamma * vz[i];
                }
            });
        }
    }
}
