
namespace MathNet.Numerics.LinearAlgebra
{
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;

    internal static class VectorExtensions
    {
        public static T[] Data<T>(this Vector<T> vector)
            where T : struct, IEquatable<T>, IFormattable
        {
            return ((DenseVectorStorage<T>)vector.Storage).Data;
        }

        /// <summary>
        /// Copy a vector to an array.
        /// </summary>
        /// <param name="src">The source vector.</param>
        /// <param name="dst">The destination array.</param>
        public static void CopyTo<T>(this Vector<T> src, T[] dst)
            where T : struct, IEquatable<T>, IFormattable
        {
            int length = dst.Length;

            var dense = src.Storage as DenseVectorStorage<T>;

            if (dense != null)
            {
                Array.Copy(dense.Data, 0, dst, 0, length);
                return;
            }

            for (int i = 0; i < length; i++)
            {
                dst[i] = src.At(i);
            }
        }
    }
}
