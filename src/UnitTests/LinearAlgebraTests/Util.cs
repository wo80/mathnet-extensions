
namespace MathNet.Numerics.Extensions.UnitTests.LinearAlgebraTests
{
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;

    static class Util
    {
        public static T[] GetData<T>(Vector<T> vector)
            where T : struct, IEquatable<T>, IFormattable
        {
            return ((DenseVectorStorage<T>)vector.Storage).Data;
        }

        public static bool Equals(Vector<double> v, double[] y, double eps)
        {
            var x = GetData(v);

            for (int i = 0; i < v.Count; i++)
            {
                if (Math.Abs(x[i] - y[i]) > eps)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Returns a permutation vector of length n.
        /// </summary>
        /// <param name="n">Length of the permutation.</param>
        /// <param name="seed">0: identity, -1: reverse, seed > 0: random permutation</param>
        /// <returns>Permutation vector.</returns>
        public static int[] Permutation(int n, int seed = 0)
        {
            int i, j, tmp;
            int[] p = new int[n];

            if (seed == 0)
            {
                // Identity
                for (i = 0; i < n; i++)
                {
                    p[i] = i;
                }
            }
            else
            {
                // Inverse
                for (i = 0; i < n; i++)
                {
                    p[i] = n - i - 1;
                }
            }

            if (seed > 0)
            {
                // Randomize permutation.
                var rand = new System.Random(seed);

                for (i = 0; i < n; i++)
                {
                    // j = rand integer in range k to n-1
                    j = i + (rand.Next() % (n - i));

                    // swap p[k] and p[j]
                    tmp = p[j];
                    p[j] = p[i];
                    p[i] = tmp;
                }
            }

            return p;
        }
    }
}
