
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
        public int[] RowIndices
        {
            get { return rowind; }
        }

        /// <summary>
        /// Column indices (size = NonZerosCount)
        /// </summary>
        public int[] ColumnIndices
        {
            get { return colind; }
        }

        /// <summary>
        /// Numerical values (size = NonZerosCount)
        /// </summary>
        public T[] Values
        {
            get { return values; }
        }

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int RowCount
        {
            get { return nrows; }
        }

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int ColumnCount
        {
            get { return ncols; }
        }

        /// <summary>
        /// Gets the number of non-zero entries.
        /// </summary>
        public int NonZerosCount
        {
            get { return nz; }
        }

        public CoordinateStorage(int rowCount, int columnCount, bool alloc = true)
            : this(rowCount, columnCount, 4, true)
        {
        }

        /// <summary>
        /// Initializes a new instance of the CoordinateStorage class.
        /// </summary>
        public CoordinateStorage(int rowCount, int columnCount, int nzmax, bool alloc = true)
        {
            this.nrows = rowCount;
            this.ncols = columnCount;

            if (alloc)
            {
                this.nzmax = (nzmax = Math.Max(nzmax, 1));
                this.nz = 0;

                this.rowind = new int[nzmax];
                this.colind = new int[nzmax];
                this.values = new T[nzmax];
            }
        }

        /// <summary>
        /// Adds an entry. Memory and dimension of the matrix are increased if necessary.
        /// </summary>
        /// <param name="i">Row index of new entry</param>
        /// <param name="j">Column index of new entry</param>
        /// <param name="value">Numerical value of new entry</param>
        /// <returns>True if successful, false otherwise</returns>
        public void At(int i, int j, T value)
        {
            if (i < 0 || j < 0)
            {
                return;
            }

            if (value.Equals(Zero))
            {
                return;
            }

            if (nz >= nzmax)
            {
                this.Resize(2 * nzmax);
            }

            if (i < 0 || i >= nrows)
            {
                throw new ArgumentOutOfRangeException("i");
            }

            if (j < 0 || j >= ncols)
            {
                throw new ArgumentOutOfRangeException("j");
            }

            rowind[nz] = i;
            colind[nz] = j;
            values[nz] = value;

            nz += 1;
        }

        public Matrix<T> ToSparseMatrix(bool cleanup = true)
        {
            if (typeof(T) == typeof(double))
            {
                return MathNet.Numerics.LinearAlgebra.Double.StorageConverter
                    .ToSparseMatrix(this as CoordinateStorage<double>, cleanup) as Matrix<T>;
            }

            if (typeof(T) == typeof(Complex))
            {
                return MathNet.Numerics.LinearAlgebra.Complex.StorageConverter
                    .ToSparseMatrix(this as CoordinateStorage<Complex>, cleanup) as Matrix<T>;
            }

            throw new NotImplementedException();
        }

        /// <summary>
        /// Resize the storage arrays of the sparse matrix.
        /// </summary>
        /// <param name="size">The new size of Values and ColumnIndices arrays.</param>
        /// <remarks>
        /// Use size = 0 to automatically resize to non-zeros count.
        /// </remarks>
        protected void Resize(int size)
        {
            if (size <= 0)
            {
                size = this.nz;
            }

            Array.Resize<int>(ref this.rowind, size);
            Array.Resize<int>(ref this.colind, size);
            Array.Resize<T>(ref this.values, size);

            this.nzmax = size;
        }
    }
}
