
namespace MathNet.Numerics.LinearAlgebra.Complex
{
    using System;
    using System.Numerics;

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
            var x = (Complex[])vector;

            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Scale this vector.
        /// </summary>
        public static void Scale(this DenseVector vector, Complex alpha)
        {
            var vx = (Complex[])vector;

            int length = vx.Length;

            for (int i = 0; i < length; i++)
            {
                vx[i] *= alpha;
            }
        }
    }
}
