
namespace MathNet.Numerics.LinearAlgebra.Storage
{
    using System;
    using System.Numerics;

    /// <summary>
    /// Coordinate storage sparse matrix format.
    /// </summary>
    public class CoordinateStorage<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        private static readonly T Zero = default(T);

        private int nrows;
        private int ncols;
        private int nz; // Number of entries in triplet matrix
        private int nzmax; // Maximum number of entries

        private int[] rowind; // Row indices (size nzmax)
        private int[] colind; // Column indices (size nzmax)
        private T[] values; // Numerical values (size nzmax)

        /// <summary>
        /// Row indices (size = NonZerosCount)
        /// </summary>
        public int[] RowIndices => rowind;

        /// <summary>
        /// Column indices (size = NonZerosCount)
        /// </summary>
        public int[] ColumnIndices => colind;

        /// <summary>
        /// Numerical values (size = NonZerosCount)
        /// </summary>
        public T[] Values => values;

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int RowCount => nrows;

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int ColumnCount => ncols;

        /// <summary>
        /// Gets the number of non-zero entries.
        /// </summary>
        public int NonZerosCount => nz;

        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateStorage{T}"/> class.
        /// </summary>
        public CoordinateStorage(int rowCount, int columnCount)
            : this(rowCount, columnCount, 4)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateStorage{T}"/> class.
        /// </summary>
        public CoordinateStorage(int rowCount, int columnCount, int nzmax)
        {
            if (rowCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(rowCount), Resources.MatrixDimensionNonNegative);
            }

            if (columnCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(columnCount), Resources.MatrixDimensionNonNegative);
            }

            if (nzmax < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(nzmax), Resources.ValueNonNegative);
            }

            nrows = rowCount;
            ncols = columnCount;

            this.nzmax = (nzmax = Math.Max(nzmax, 1));
            this.nz = 0;

            rowind = new int[nzmax];
            colind = new int[nzmax];
            values = new T[nzmax];
        }

        /// <summary>
        /// Adds an entry to the storage.
        /// </summary>
        /// <param name="i">Row index of new entry</param>
        /// <param name="j">Column index of new entry</param>
        /// <param name="value">Numerical value of new entry</param>
        public void At(int i, int j, T value)
        {
            if (value.Equals(Zero))
            {
                return;
            }

            if (i < 0 || i >= nrows)
            {
                throw new ArgumentOutOfRangeException(nameof(i));
            }

            if (j < 0 || j >= ncols)
            {
                throw new ArgumentOutOfRangeException(nameof(j));
            }

            if (nz >= nzmax)
            {
                Resize(2 * nzmax);
            }

            rowind[nz] = i;
            colind[nz] = j;
            values[nz] = value;

            nz += 1;
        }

        /// <summary>
        /// Filter matrix values.
        /// </summary>
        /// <param name="func">Filter function returning true if value should be kept,
        /// false if value should be discarded.</param>
        /// <returns>New number of non-zeros.</returns>
        /// <remarks>
        /// Filter function arguments:
        /// 
        /// 1 = Row index i
        /// 2 = Column index j
        /// 3 = Value of entry (i,j)
        /// 
        /// Element a_{i,j} is dropped, if func(i, j, aij) returns false.
        /// </remarks>
        public int Keep(Func<int, int, T, bool> func)
        {
            int k = 0;

            for (int i = 0; i < nz; i++)
            {
                int ai = rowind[i];
                int aj = colind[i];
                var ax = values[i];

                if (func(ai, aj, ax))
                {
                    // Keep A(i,j).
                    rowind[k] = ai;
                    colind[k] = aj;
                    values[k] = ax;
                    k++;
                }
            }

            return nz = k;
        }

        /// <summary>
        /// Remove all values from the storage (without freeing the memory).
        /// </summary>
        public void Clear()
        {
            Array.Clear(rowind, 0, nzmax);
            Array.Clear(colind, 0, nzmax);
            Array.Clear(values, 0, nzmax);

            nz = 0;
        }

        /// <summary>
        /// Transpose the coordinate storage inline.
        /// </summary>
        public void TransposeInline()
        {
            // Transposing is just a matter of switching row and column indices.
            var temp = colind;
            colind = rowind;
            rowind = temp;
        }

        /// <summary>
        /// Convert the coordinate storage to a sparse matrix using <see cref="SparseCompressedRowMatrixStorage{T}"/>.
        /// </summary>
        /// <param name="cleanup">Remove and sum duplicate entries.</param>
        /// <returns></returns>
        public Matrix<T> ToSparseMatrix(bool cleanup = true)
        {
            if (typeof(T) == typeof(double))
            {
                return LinearAlgebra.Double.StorageConverter
                    .ToSparseMatrix(this as CoordinateStorage<double>, cleanup) as Matrix<T>;
            }

            if (typeof(T) == typeof(Complex))
            {
                return LinearAlgebra.Complex.StorageConverter
                    .ToSparseMatrix(this as CoordinateStorage<Complex>, cleanup) as Matrix<T>;
            }

            throw new NotImplementedException();
        }

        /// <summary>
        /// Resize the storage arrays of the sparse matrix.
        /// </summary>
        /// <param name="size">The new size of the storage arrays.</param>
        /// <remarks>
        /// Use size = 0 to automatically resize to non-zeros count.
        /// </remarks>
        protected void Resize(int size)
        {
            if (size <= 0)
            {
                size = nz;
            }

            Array.Resize(ref rowind, size);
            Array.Resize(ref colind, size);
            Array.Resize(ref values, size);

            this.nzmax = size;
        }
    }
}
