
namespace MathNet.Numerics.LinearAlgebra.Storage
{
    using System;

    /// <summary>
    /// Storage extension methods.
    /// </summary>
    public static class StorageExtensions
    {
        #region SparseCompressedRowMatrixStorage extensions

        /// <summary>
        /// Returns the indices of the diagonal entries in the sparse matrix storage.
        /// </summary>
        /// <param name="storage"></param>
        /// <param name="throwOnMissingDiag"></param>
        /// <returns>Indices of the diagonal entries.</returns>
        public static int[] FindDiagonalIndices<T>(this SparseCompressedRowMatrixStorage<T> storage,
            bool throwOnMissingDiag = false)
            where T : struct, IEquatable<T>, IFormattable
        {
            return Helper.FindDiagonalIndices(storage.RowCount, storage.RowPointers,
                storage.ColumnIndices, throwOnMissingDiag);
        }

        /// <summary>
        /// Permute the rows of a matrix in CSR format, B = P * A, where P represents
        /// a permutation matrix. 
        /// </summary>
        /// <param name="storage">Input storage.</param>
        /// <param name="target">Target storage.</param>
        /// <param name="perm">Permutation array of length RowCount.</param>
        /// <param name="sorted">Ensure that column indices will be sorted.</param>
        /// <remarks>
        /// The permutation P is defined through the array perm: for each j,
        /// perm(j) represents the destination row number of row number j:
        /// 
        /// a(i,j) in the original matrix becomes a(perm(i),j) in the output matrix.
        /// </remarks>
        public static void PermuteRows<T>(this SparseCompressedRowMatrixStorage<T> storage,
            SparseCompressedRowMatrixStorage<T> target, int[] perm, bool sorted = true)
            where T : struct, IEquatable<T>, IFormattable
        {
            int k, n = storage.RowCount;

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            var bx = target.Values;
            var bp = target.RowPointers;
            var bi = target.ColumnIndices;

            // Determine pointers for output matix. 
            for (int i = 0; i < n; i++)
            {
                k = perm[i];
                bp[k + 1] = ap[i + 1] - ap[i];
            }

            // Get pointers from lengths
            bp[0] = 0;
            for (int i = 0; i < n; i++)
            {
                bp[i + 1] += bp[i];
            }

            // Copying
            for (int i = 0; i < n; i++)
            {
                // Old row = i, new row = perm(i), k = new pointer
                k = bp[perm[i]];
                for (int j = ap[i]; j < ap[i + 1]; j++)
                {
                    bi[k] = ai[j];
                    bx[k] = ax[j];
                    k = k + 1;
                }
            }

            if (sorted)
            {
                Helper.SortIndices(storage.RowCount, bx, bp, bi);
            }
        }

        /// <summary>
        /// Permutes the columns of a matrix in CSR format, B = A * P, where P represents
        /// a permutation matrix.
        /// </summary>
        /// <param name="storage">Input storage.</param>
        /// <param name="target">Target storage.</param>
        /// <param name="perm">Permutation array of length ColumnCount.</param>
        /// <param name="copy">Copy matrix values (not needed if used 'in place').</param>
        /// <param name="sorted">Ensure that column indices will be sorted.</param>
        /// <remarks>
        /// The permutation matrix P maps column j into column perm(j), i.e., 
        /// on return a(i,j) in the original matrix becomes a(i,perm(j)) in the
        /// output matrix.
        /// 
        /// Notes:
        /// 
        /// 1. This routine is in place: aj, bj can be the same.
        /// 2. If the matrix is initially sorted (by increasing column number) 
        ///    then bx, bi, bj may not be on return.
        /// </remarks>
        public static void PermuteColumns<T>(this SparseCompressedRowMatrixStorage<T> storage,
            SparseCompressedRowMatrixStorage<T> target, int[] perm, bool copy = false,
            bool sorted = true)
            where T : struct, IEquatable<T>, IFormattable
        {
            int n = storage.RowCount;

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            var bx = target.Values;
            var bp = target.RowPointers;
            var bi = target.ColumnIndices;

            int i, nnz = ap[n];

            for (i = 0; i < nnz; i++)
            {
                bi[i] = perm[ai[i]];
            }

            if (copy)
            {
                Array.Copy(ax, bx, nnz);
                Array.Copy(ap, bp, n);
            }

            if (sorted)
            {
                Helper.SortIndices(storage.RowCount, bx, bp, bi);
            }
        }

        /// <summary>
        /// Extract row from storage.
        /// </summary>
        /// <param name="storage"></param>
        /// <param name="rowIndex">The row index to extract.</param>
        /// <param name="target">Dense array.</param>
        /// <returns>The requested row.</returns>
        public static void GetRow<T>(this SparseCompressedRowMatrixStorage<T> storage,
            int rowIndex, T[] target, ExistingData existingData = ExistingData.Clear)
            where T : struct, IEquatable<T>, IFormattable
        {
            if (target.Length != storage.ColumnCount)
            {
                throw new Exception();
            }

            if (existingData == ExistingData.Clear)
            {
                Array.Clear(target, 0, target.Length);
            }

            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int rowEnd = ap[rowIndex + 1];

            for (int k = ap[rowIndex]; k < rowEnd; k++)
            {
                target[ai[k]] = ax[k];
            }
        }

        /// <summary>
        /// Extract rows from given storage.
        /// </summary>
        /// <param name="storage"></param>
        /// <param name="rowIndices">The rows to extract.</param>
        /// <param name="target">The target storage (will be resized, if needed).</param>
        /// <returns>The requested sub-matrix.</returns>
        public static void GetRows<T>(this SparseCompressedRowMatrixStorage<T> storage,
            int[] rowIndices, SparseCompressedRowMatrixStorage<T> target)
            where T : struct, IEquatable<T>, IFormattable
        {
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int count = rowIndices.Length;

            if (target.RowCount != count || storage.ColumnCount != target.ColumnCount)
            {
                throw new Exception();
            }

            int i, needed = 0;

            // Compute space needed.
            for (int k = 0; k < count; k++)
            {
                i = rowIndices[k];

                needed += ap[i + 1] - ap[i];
            }

            var cp = target.RowPointers;
            var ci = target.ColumnIndices;
            var cx = target.Values;

            // We assume that cx and ci have the same length.
            int nz = cx.Length;

            if (nz < needed)
            {
                // Resize workspace.
                Array.Resize(ref ci, needed);
                Array.Resize(ref cx, needed);
            }

            nz = 0;

            for (int k = 0; k < count; k++)
            {
                i = rowIndices[k];

                int rowEnd = ap[i + 1];

                cp[k] = nz;

                // Copy row.
                for (int j = ap[i]; j < rowEnd; j++)
                {
                    cx[nz] = ax[j];
                    ci[nz] = ai[j];

                    nz = nz + 1;
                }
            }

            cp[count] = nz;

            target.ColumnIndices = ci;
            target.Values = cx;
        }

        /// <summary>
        /// Extract columns from given storage.
        /// </summary>
        /// <param name="storage"></param>
        /// <param name="columnIndices">The columns to extract.</param>
        /// <param name="target">The target storage (will be resized, if needed).</param>
        /// <returns>The requested sub-matrix.</returns>
        public static void GetColumns<T>(this SparseCompressedRowMatrixStorage<T> storage,
            int[] columnIndices, SparseCompressedRowMatrixStorage<T> target)
            where T : struct, IEquatable<T>, IFormattable
        {
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int count = columnIndices.Length;
            int rowCount = storage.RowCount;

            if (target.ColumnCount != count || rowCount != target.RowCount)
            {
                throw new Exception();
            }

            int i, j, nz, needed = 0;

            var columnCounts = new int[storage.ColumnCount];

            nz = storage.ValueCount;

            // Count columns.
            for (i = 0; i < nz; i++)
            {
                columnCounts[ai[i]]++;
            }

            // Compute space needed.
            for (int k = 0; k < count; k++)
            {
                needed += columnCounts[columnIndices[k]];
            }

            var cp = target.RowPointers;
            var ci = target.ColumnIndices;
            var cx = target.Values;

            // We assume that cx and ci have the same length.
            nz = cx.Length;

            if (nz < needed)
            {
                // Resize workspace.
                Array.Resize(ref ci, needed);
                Array.Resize(ref cx, needed);
            }

            int column, rowStart, rowEnd;

            nz = 0;

            for (int k = 0; k < rowCount; k++)
            {
                rowStart = ap[k];
                rowEnd = ap[k + 1];

                cp[k] = nz;

                // Look for columns in current row.
                for (j = 0; j < count; j++)
                {
                    column = columnIndices[j];

                    i = Array.BinarySearch(ai, rowStart, rowEnd - rowStart, column);

                    if (i >= 0)
                    {
                        ci[nz] = j;
                        cx[nz] = ax[i];

                        nz = nz + 1;
                    }
                }
            }

            cp[rowCount] = nz;

            target.ColumnIndices = ci;
            target.Values = cx;
        }

        /// <summary>
        /// Replaces the rows of the storage with the given rows.
        /// </summary>
        /// <param name="storage"></param>
        /// <param name="rowIndices">The indices of the rows to replace.</param>
        /// <param name="source">The storage containing the rows to copy.</param>
        public static void SetRows<T>(this SparseCompressedRowMatrixStorage<T> storage,
            int[] rowIndices, SparseCompressedRowMatrixStorage<T> source)
            where T : struct, IEquatable<T>, IFormattable
        {
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int i, n = ap.Length; 
            int rowCount = rowIndices.Length;

            int valueCount = 0;

            // Count the number of nonzeros in given rows of input matrix.
            for (int k = 0; k < rowCount; k++)
            {
                i = rowIndices[k];

                valueCount += ap[i + 1] - ap[i];
            }

            int nz = storage.ValueCount;
            int size = nz + (source.ValueCount - valueCount);

            // Check if the matrix storage has to be resized.
            if (ai.Length < size)
            {
                Array.Resize(ref ai, size);
                Array.Resize(ref ax, size);

                storage.ColumnIndices = ai;
                storage.Values = ax;
            }

            size = ai.Length;

            var cp = source.RowPointers;
            var ci = source.ColumnIndices;
            var cx = source.Values;

            int count, length, rowStart, rowEnd;

            // Replace the rows.
            for (int k = 0; k < rowCount; k++)
            {
                i = rowIndices[k];

                rowStart = ap[i];
                rowEnd = ap[i + 1];

                // Length of the row to insert.
                length = cp[k + 1] - cp[k];

                // Number of elements to shift.
                count = nz - rowEnd;

                // Copy all data behind the current block.
                Array.Copy(ai, rowEnd, ai, rowStart + length, count);
                Array.Copy(ax, rowEnd, ax, rowStart + length, count);

                // Replace row.
                Array.Copy(ci, cp[k], ai, rowStart, length);
                Array.Copy(cx, cp[k], ax, rowStart, length);

                // The difference of the row lengths.
                count = length - (rowEnd - rowStart);

                // Fix row pointers.
                for (int j = i + 1; j < n; j++)
                {
                    ap[j] += count;
                }

                nz += count;
            }
        }

        /// <summary>
        /// Drops entries from a sparse matrix.
        /// </summary>
        /// <param name="storage">The sparse storage.</param>
        /// <param name="func">Drop element a_{i,j} if func(i, j, aij) is false.</param>
        /// <param name="sorted">Ensure that column indices will be sorted.</param>
        /// <returns>New number of entries in A.</returns>
        public static int Keep<T>(this SparseCompressedRowMatrixStorage<T> storage,
            Func<int, int, T, bool> func, bool sorted = true)
            where T : struct, IEquatable<T>, IFormattable
        {
            int rowCount = storage.RowCount;

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            int i, j, nz = 0;

            for (i = 0; i < rowCount; i++)
            {
                j = ap[i];

                // Record new location of row i.
                ap[i] = nz;

                for (; j < ap[i + 1]; j++)
                {
                    if (func(i, ai[j], ax[j]))
                    {
                        // Keep A(i,j).
                        ax[nz] = ax[j];
                        ai[nz] = ai[j];
                        nz++;
                    }
                }
            }

            // Record new nonzero count.
            ap[rowCount] = nz;

            if (sorted)
            {
                Helper.SortIndices(storage.RowCount, ax, ap, ai);
            }
            
            // Remove extra space.
            Array.Resize(ref storage.Values, nz);
            Array.Resize(ref storage.ColumnIndices, nz);

            return nz;
        }

        #endregion
    }
}
