
namespace MathNet.Numerics.LinearAlgebra.Complex
{
    using MathNet.Numerics.Threading;
    using System;
    using System.Numerics;

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
            return DenseVector.OfArray((Complex[])vector);
        }

        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(DenseVector src, DenseVector dst)
        {
            var x = (Complex[])src;
            var y = (Complex[])dst;

            Array.Copy(x, 0, y, 0, x.Length);
        }

        /// <summary>
        /// Set all vector elements to zero.
        /// </summary>
        public static void Clear(DenseVector vector)
        {
            var x = (Complex[])vector;

            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Scale given vector, x = alpha * x.
        /// </summary>
        public static void Scale(Complex alpha, DenseVector x)
        {
            var vx = (Complex[])x;

            int length = vx.Length;

            for (int i = 0; i < length; i++)
            {
                vx[i] *= alpha;
            }
        }

        /// <summary>
        /// Copy scaled vector, target = alpha * x.
        /// </summary>
        public static void Scale(Complex alpha, DenseVector x, DenseVector target)
        {
            var vx = (Complex[])x;
            var vt = (Complex[])target;

            int length = vx.Length;

            for (int i = 0; i < length; i++)
            {
                vt[i] = alpha * vx[i];
            }
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static Complex DotProduct(DenseVector x, DenseVector y)
        {
            var vx = (Complex[])x;
            var vy = (Complex[])y;

            int length = vx.Length;

            Complex result = 0.0;

            for (int i = 0; i < length; i++)
            {
                result += Complex.Conjugate(vx[i]) * vy[i];
            }

            return result;
        }

        /// <summary>
        /// Compute sum of scaled vectors, target = alpha * x + y.
        /// </summary>
        public static void Add(Complex alpha, DenseVector x, DenseVector y, DenseVector target)
        {
            var vx = (Complex[])x;
            var vy = (Complex[])y;
            var vt = (Complex[])target;

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
        public static void Add(Complex alpha, DenseVector x, Complex beta, DenseVector y, DenseVector target)
        {
            var vx = (Complex[])x;
            var vy = (Complex[])y;
            var vt = (Complex[])target;

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
        public static void Add(Complex alpha, DenseVector x, Complex beta, DenseVector y, DenseVector z,
            DenseVector target)
        {
            var vx = (Complex[])x;
            var vy = (Complex[])y;
            var vz = (Complex[])z;
            var vt = (Complex[])target;

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
        public static void Add(Complex alpha, DenseVector x, Complex beta, DenseVector y, 
            Complex gamma, DenseVector z, DenseVector target)
        {
            var vx = (Complex[])x;
            var vy = (Complex[])y;
            var vz = (Complex[])z;
            var vt = (Complex[])target;

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
