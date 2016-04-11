
namespace MathNet.Numerics.LinearAlgebra.Storage
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// Modified compressed sparse row storage (MSR).
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <remarks>
    /// </remarks>
    public class ModifiedSparseRowStorage<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// The number of rows.
        /// </summary>
        protected int rowCount;

        /// <summary>
        /// The number of columns.
        /// </summary>
        protected int columnCount;

        /// <summary>
        /// Combined row pointers and column indices.
        /// </summary>
        /// <remarks>
        /// The row pointers are stored at positions (1:n) and column indices to off-diagonal
        /// elements at positions (n+1:nnz).
        /// </remarks>
        public int[] Indices;

        /// <summary>
        /// Storage values.
        /// </summary>
        /// <remarks>
        /// The diagonal is stored at positions (1:n). Position (n+1) is unused.
        /// </remarks>
        public T[] Values;

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int RowCount
        {
            get { return rowCount; }
        }

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int ColumnCount
        {
            get { return columnCount; }
        }

        /// <summary>
        /// Gets the number of non-zero entries.
        /// </summary>
        public int ValueCount
        {
            get { return Indices[rowCount]; }
        }

        /// <summary>
        /// Initializes a new instance of the ModifiedSparseRowStorage class.
        /// </summary>
        public ModifiedSparseRowStorage(int rowCount, int columnCount)
            : this(rowCount, columnCount, 0)
        {
        }

        /// <summary>
        /// Initializes a new instance of the ModifiedSparseRowStorage class.
        /// </summary>
        public ModifiedSparseRowStorage(int rowCount, int columnCount, int valueCount)
        {
            if (rowCount <= 0)
            {
                throw new ArgumentOutOfRangeException("rowCount");
            }

            if (columnCount <= 0)
            {
                throw new ArgumentOutOfRangeException("columnCount");
            }

            this.rowCount = rowCount;
            this.columnCount = columnCount;

            this.Indices = new int[rowCount + 1 + valueCount];
            this.Values = new T[rowCount + 1 + valueCount];
        }

        /// <summary>
        /// Enumerates the sparse matrix storage.
        /// </summary>
        /// <returns>Enumeration of all storage entries (i, j, a_{ij}).</returns>
        public IEnumerable<Tuple<int, int, T>> Enumerate()
        {
            int n = this.RowCount;

            var ax = this.Values;
            var ai = this.Indices;

            // The diagonal entries.
            for (int i = 0; i < n; i++)
            {
                yield return new Tuple<int, int, T>(i, i, ax[i]);
            }

            // The off-diagonal entries (row-wise).
            for (int i = 0; i < n; i++)
            {
                int end = ai[i + 1];

                for (var j = ai[i]; j < end; j++)
                {
                    yield return new Tuple<int, int, T>(i, ai[j], ax[j]);
                }
            }
        }
    }
}
