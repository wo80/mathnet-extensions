
namespace MathNet.Numerics.LinearAlgebra.Double
{
    using System;
    using MathNet.Numerics.LinearAlgebra.Storage;

    internal static class StorageConverter
    {
        /// <summary>
        /// Convert a coordinate storage to compressed sparse row (CSR) format.
        /// </summary>
        /// <param name="storage">Coordinate storage.</param>
        /// <param name="cleanup">Remove and sum duplicate entries.</param>
        /// <returns>Compressed sparse row storage.</returns>
        public static SparseMatrix ToSparseMatrix(CoordinateStorage<double> coo, bool cleanup = true)
        {
            int nrows = coo.RowCount;
            int ncols = coo.ColumnCount;

            var rowind = coo.RowIndices;
            var colind = coo.ColumnIndices;
            var values = coo.Values;

            int p, k, nz = coo.NonZerosCount;

            var counts = new int[nrows];

            for (k = 0; k < nz; k++)
            {
                // Determine row-lengths.
                counts[rowind[k]]++;
            }

            var result = new SparseMatrix(nrows, ncols);
            var storage = result.Storage as SparseCompressedRowMatrixStorage<double>;

            var ap = storage.RowPointers;

            // Get row pointers (starting position of each row).
            int valueCount = Helper.CumulativeSum(ap, counts, nrows);

            var ai = new int[valueCount];
            var ax = new double[valueCount];

            // Fill in output matrix.
            for (k = 0; k < nz; k++)
            {
                p = counts[rowind[k]]++;
                ai[p] = colind[k];
                ax[p] = values[k];
            }

            Helper.SortIndices(nrows, ax, ap, ai);

            if (cleanup)
            {
                Cleanup(ncols, nrows, ap, ref ai, ref ax);
            }

            storage.ColumnIndices = ai;
            storage.Values = ax;

            return result;
        }

        public static void Cleanup(int ncols, int nrows, int[] ap, ref int[] ai, ref double[] ax)
        {
            int i, j, p, q, nz = 0;
            int[] marker = new int[ncols];

            for (j = 0; j < ncols; j++)
            {
                marker[j] = -1; // Column j not yet seen
            }

            for (i = 0; i < nrows; i++)
            {
                q = nz; // Row i will start at q
                for (p = ap[i]; p < ap[i + 1]; p++)
                {
                    j = ai[p]; // A(i,j) is nonzero
                    if (marker[j] >= q)
                    {
                        ax[marker[j]] += ax[p]; // A(i,j) is a duplicate
                    }
                    else
                    {
                        marker[j] = nz; // Record where column j occurs
                        ai[nz] = j; // keep A(i,j)
                        ax[nz] = ax[p];

                        nz += 1;
                    }
                }
                ap[i] = q; // Record start of row i
            }

            ap[nrows] = nz;

            // Remove extra space from arrays.
            Array.Resize<double>(ref ax, nz);
            Array.Resize<int>(ref ai, nz);
        }
    }
}
