
namespace MathNet.Numerics.LinearAlgebra
{
    using System;

    internal static class Helper
    {
        const int InsertionSortThreshold = 10;

        public static void SortIndices<T>(int n, T[] ax, int[] ap, int[] ai)
            where T : struct, IEquatable<T>, IFormattable
        {
            int from, to, size, p, q, idx;
            T val;

            for (int i = 0; i < n; i++)
            {
                from = ap[i];
                to = ap[i + 1] - 1;

                size = to - from + 1;

                if (size > InsertionSortThreshold)
                {
                    // Quicksort
                    Array.Sort(ai, ax, from, size);
                }
                else
                {
                    // Insertion sort
                    for (p = from + 1; p <= to; p++)
                    {
                        idx = ai[p];
                        val = ax[p];
                        q = p - 1;
                        while (q >= from && ai[q] > idx)
                        {
                            ai[q + 1] = ai[q];
                            ax[q + 1] = ax[q];
                            q--;
                        }
                        ai[q + 1] = idx;
                        ax[q + 1] = val;
                    }
                }
            }
        }

        /// <summary>
        /// Cumulative sum of given array.
        /// </summary>
        /// <param name="sum">Output: cumulative sum of counts</param>
        /// <param name="counts">input array, overwritten with sum</param>
        /// <param name="size">length of counts</param>
        /// <returns>sum[size] (non-zeros)</returns>
        public static int CumulativeSum(int[] sum, int[] counts, int size)
        {
            int i, nz = 0;

            for (i = 0; i < size; i++)
            {
                sum[i] = nz;
                nz += counts[i];
                counts[i] = sum[i]; // also copy p[0..n-1] back into c[0..n-1]
            }

            sum[size] = nz;

            return nz;
        }

        /// <summary>
        /// Returns the positions of the diagonal elements of a sparse matrix.
        /// </summary>
        /// <param name="n">Row count</param>
        /// <param name="ap">Row pointers</param>
        /// <param name="ai">Column indices</param>
        /// <param name="throwOnMissingDiag"></param>
        /// <returns></returns>
        public static int[] FindDiagonalIndices(int n, int[] ap, int[] ai, bool throwOnMissingDiag = false)
        {
            int[] diag = new int[n];

            for (int i = 0; i < n; i++)
            {
                diag[i] = Array.BinarySearch(ai, ap[i], ap[i + 1] - ap[i], i);

                if (diag[i] < 0 && throwOnMissingDiag)
                {
                    throw new Exception("Missing diagonal entry on row " + (i + 1));
                }
            }

            return diag;
        }
    }
}
