
namespace ExtensionsBenchmark
{
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;
    using System.Diagnostics;

    static class Util
    {
        #region Stopwatch

        private static Stopwatch timer = new Stopwatch();

        public static void Tic()
        {
            timer.Restart();
        }

        public static double Toc()
        {
            timer.Stop();

            return TimeSpan.FromTicks(timer.ElapsedTicks).TotalMilliseconds;
        }

        #endregion

        public static T[] GetData<T>(Vector<T> vector)
            where T : struct, IEquatable<T>, IFormattable
        {
            return ((DenseVectorStorage<T>)vector.Storage).Data;
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
                var rand = new Random(seed);

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

        /// <summary>
        /// Returns a permutation vector of length n.
        /// </summary>
        /// <param name="n">Length of the permutation.</param>
        /// <param name="rand">random source</param>
        /// <returns>Permutation vector.</returns>
        public static int[] Permutation(int n, Random rand)
        {
            int i, j, tmp;
            int[] p = new int[n];

            for (i = 0; i < n; i++)
            {
                p[i] = n - i - 1;
            }

            for (i = 0; i < n; i++)
            {
                // j = rand integer in range k to n-1
                j = i + (rand.Next() % (n - i));

                // swap p[k] and p[j]
                tmp = p[j];
                p[j] = p[i];
                p[i] = tmp;
            }

            return p;
        }
    }
}
