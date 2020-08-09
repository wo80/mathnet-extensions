
namespace MathNet.Numerics.LinearAlgebra.Double
{
    using System;
    using MathNet.Numerics.LinearAlgebra.Storage;

    /// <summary>
    /// SparseMatrix extension methods.
    /// </summary>
    public static class SparseMatrixExtensions
    {
        #region Linear Algebra (Vector)

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = A*x + y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void Multiply(this SparseMatrix matrix, Vector<double> x, Vector<double> y)
        {
            Multiply(matrix, x.Data(), y.Data());
        }

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = A*x + y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void Multiply(this SparseMatrix matrix, double[] x, double[] y)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            int nrows = matrix.RowCount;

            double yi;

            for (int i = 0; i < nrows; i++)
            {
                yi = 0.0;

                // Compute the inner product of row i with vector x
                for (int k = ap[i]; k < ap[i + 1]; k++)
                {
                    yi += ax[k] * x[ai[k]];
                }

                // Store result in y(i) 
                y[i] += yi;
            }
        }

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = alpha*A*x + beta*y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="alpha">Scalar to multiply with matrix.</param>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="beta">Scalar to multiply with vector y.</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void Multiply(this SparseMatrix matrix, double alpha, Vector<double> x,
            double beta, Vector<double> y)
        {
            Multiply(matrix, alpha, x.Data(), beta, y.Data());
        }

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = alpha*A*x + beta*y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="alpha">Scalar to multiply with matrix.</param>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="beta">Scalar to multiply with vector y.</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void Multiply(this SparseMatrix matrix, double alpha, double[] x, double beta, double[] y)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            int nrows = matrix.RowCount;

            double yi;

            int end, start = ap[0];

            for (int i = 0; i < nrows; i++)
            {
                end = ap[i + 1];

                yi = beta * y[i];
                for (int k = start; k < end; k++)
                {
                    yi += alpha * ax[k] * x[ai[k]];
                }
                y[i] = yi;

                start = end;
            }
        }

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = A'*x + y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void TransposeMultiply(this SparseMatrix matrix, Vector<double> x, Vector<double> y)
        {
            TransposeMultiply(matrix, x.Data(), y.Data());
        }

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = A'*x + y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void TransposeMultiply(this SparseMatrix matrix, double[] x, double[] y)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            int nrows = matrix.RowCount;

            int end;

            for (int j = 0; j < nrows; j++)
            {
                end = ap[j + 1];

                // Loop over the rows
                for (int k = ap[j]; k < end; k++)
                {
                    y[ai[k]] += x[j] * ax[k];
                }
            }
        }

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = alpha*A'*x + beta*y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="alpha">Scalar to multiply with matrix.</param>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="beta">Scalar to multiply with vector y.</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void TransposeMultiply(this SparseMatrix matrix, double alpha, Vector<double> x,
            double beta, Vector<double> y)
        {
            Multiply(matrix, alpha, x.Data(), beta, y.Data());
        }

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = alpha*A'*x + beta*y. 
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="alpha">Scalar to multiply with matrix.</param>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="beta">Scalar to multiply with vector y.</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        /// <remarks>
        /// Input values of vector y will be accumulated.
        /// </remarks>
        public static void TransposeMultiply(this SparseMatrix matrix, double alpha, double[] x, double beta, double[] y)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            int nrows = matrix.RowCount;
            int ncols = matrix.ColumnCount;

            // Scale y by beta
            for (int j = 0; j < ncols; j++)
            {
                y[j] = beta * y[j];
            }

            int end;
            double xi;

            // For each column j of A^T, do a daxpy / scatter
            for (int i = 0; i < nrows; i++)
            {
                xi = alpha * x[i];

                end = ap[i + 1];

                for (int k = ap[i]; k < end; k++)
                {
                    y[ai[k]] += ax[k] * xi;
                }
            }
        }

        #endregion

        #region Matrix Operations

        /// <summary>
        /// Compute the transpose of the matix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="target">The transposed matrix.</param>
        public static void FastTranspose(this SparseMatrix matrix, SparseMatrix target)
        {
            if (target == null)
            {
                throw new ArgumentNullException("target");
            }

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            int rowCount = matrix.RowCount;
            int columnCount = matrix.ColumnCount;

            if (rowCount != target.ColumnCount || columnCount != target.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixDimensions);
            }

            int i, j, p;

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            storage = target.Storage as SparseCompressedRowMatrixStorage<double>;

            var cp = storage.RowPointers;
            var ci = storage.ColumnIndices;
            var cx = storage.Values;

            if (cx.Length < ax.Length)
            {
                // The only purpose of this method, is NOT to allocate new memory
                // for the transpose, but ...
                storage.Values = cx = new double[ax.Length];
                storage.ColumnIndices = ci = new int[ai.Length];
            }

            int[] w = new int[columnCount];

            for (p = 0; p < ap[rowCount]; p++)
            {
                // Column counts
                w[ai[p]]++;
            }

            // Column pointers
            Helper.CumulativeSum(cp, w, columnCount);

            for (i = 0; i < rowCount; i++)
            {
                for (p = ap[i]; p < ap[i + 1]; p++)
                {
                    j = w[ai[p]]++;

                    // Place A(i,j) as entry C(j,i)
                    ci[j] = i;
                    cx[j] = ax[p];
                }
            }
        }

        /// <summary>
        /// Permute the rows of a sparse matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="perm">Permutation vector.</param>
        public static void FastPermuteRows(this SparseMatrix matrix, int[] perm)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            // TODO: clone needed?
            var temp = matrix.Clone();

            var target = temp.Storage as SparseCompressedRowMatrixStorage<double>;

            if (target == null)
            {
                throw new Exception();
            }

            storage.PermuteRows(target, perm);

            target.CopyTo(storage);
        }

        /// <summary>
        /// Permute the columns of a sparse matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="perm">Permutation vector.</param>
        public static void FastPermuteColumns(this SparseMatrix matrix, int[] perm)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            storage.PermuteColumns(storage, perm);
        }

        /// <summary>
        /// Compute the sum of this matrix and another matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="other">The other matrix.</param>
        /// <returns></returns>
        public static Matrix<double> FastAdd(this SparseMatrix matrix, SparseMatrix other)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            storage = other.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var bx = storage.Values;
            var bp = storage.RowPointers;
            var bi = storage.ColumnIndices;

            int nrows_A = matrix.RowCount;
            int ncols_A = matrix.ColumnCount;

            int nrows_B = other.RowCount;
            int ncols_B = other.ColumnCount;

            if (nrows_A != nrows_B || ncols_A != ncols_B)
            {
                throw new Exception(); // TODO: Exception text
            }

            int ia, ib, ic, jcol, count = 0;
            int pos;
            int[] marker;

            marker = new int[ncols_A];

            var result = new SparseMatrix(nrows_A, ncols_B);

            storage = result.Storage as SparseCompressedRowMatrixStorage<double>;

            var cp = storage.RowPointers;

            for (ia = 0; ia < ncols_A; ia++)
            {
                marker[ia] = -1;
            }

            count = 0;
            cp[0] = 0;

            // Count non-zeros of each row.
            for (ic = 0; ic < nrows_A; ic++)
            {
                for (ia = ap[ic]; ia < ap[ic + 1]; ia++)
                {
                    jcol = ai[ia];
                    marker[jcol] = ic;
                    count++;
                }

                for (ib = bp[ic]; ib < bp[ic + 1]; ib++)
                {
                    jcol = bi[ib];
                    if (marker[jcol] != ic)
                    {
                        marker[jcol] = ic;
                        count++;
                    }
                }

                cp[ic + 1] = count;
            }

            var ci = new int[count];
            var cx = new double[count];

            for (ia = 0; ia < ncols_A; ia++)
            {
                marker[ia] = -1;
            }

            pos = 0;

            // Add matrices.
            for (ic = 0; ic < nrows_A; ic++)
            {
                for (ia = ap[ic]; ia < ap[ic + 1]; ia++)
                {
                    jcol = ai[ia];
                    ci[pos] = jcol;
                    cx[pos] = ax[ia];
                    marker[jcol] = pos;
                    pos++;
                }

                for (ib = bp[ic]; ib < bp[ic + 1]; ib++)
                {
                    jcol = bi[ib];
                    if (marker[jcol] < cp[ic])
                    {
                        ci[pos] = jcol;
                        cx[pos] = bx[ib];
                        marker[jcol] = pos;
                        pos++;
                    }
                    else
                    {
                        cx[marker[jcol]] += bx[ib];
                    }
                }
            }

            storage.Values = cx;
            storage.ColumnIndices = ci;

            Helper.SortIndices(nrows_A, cx, cp, ci);

            return result;
        }

        /// <summary>
        /// Computes the Kronecker product of this matrix with given matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="other">The other matrix.</param>
        /// <returns>Kronecker product.</returns>
        public static SparseMatrix FastKroneckerProduct(this SparseMatrix matrix, SparseMatrix other)
        {
            var storageA = matrix.Storage as SparseCompressedRowMatrixStorage<double>;
            var storageB = other.Storage as SparseCompressedRowMatrixStorage<double>;

            var ap = storageA.RowPointers;
            var aj = storageA.ColumnIndices;
            var ax = storageA.Values;

            var bp = storageB.RowPointers;
            var bj = storageB.ColumnIndices;
            var bx = storageB.Values;

            int rowsA = matrix.RowCount;
            int rowsB = other.RowCount;

            var counts = new int[rowsA * rowsB];

            int k = 0;

            // Count non-zeros in each row of kron(A, B).
            for (int i = 0; i < rowsA; i++)
            {
                for (int j = 0; j < rowsB; j++)
                {
                    counts[k++] = (ap[i + 1] - ap[i]) * (bp[j + 1] - bp[j]);
                }
            }

            int colsA = matrix.ColumnCount;
            int colsB = other.ColumnCount;

            var C = new SparseMatrix(rowsA * rowsB, colsA * colsB);
            var storageC = C.Storage as SparseCompressedRowMatrixStorage<double>;

            var cp = storageC.RowPointers;

            int nnz = Helper.CumulativeSum(cp, counts, counts.Length);

            var cj = storageC.ColumnIndices = new int[nnz];
            var cx = storageC.Values = new double[nnz];

            k = 0;

            // For each row in A ...
            for (int ia = 0; ia < rowsA; ia++)
            {
                // ... and each row in B ...
                for (int ib = 0; ib < rowsB; ib++)
                {
                    // ... get element a_{ij}
                    for (int j = ap[ia]; j < ap[ia + 1]; j++)
                    {
                        var idx = aj[j];
                        var aij = ax[j];

                        // ... and multiply it with current row of B
                        for (int s = bp[ib]; s < bp[ib + 1]; s++)
                        {
                            cj[k] = (idx * colsB) + bj[s];
                            cx[k] = aij * bx[s];
                            k++;
                        }
                    }
                }
            }

            return C;
        }

        /// <summary>
        /// Adds a scalar to the diagonal entries of a sparse matrix A = A + s*I.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="value">Scalar to add to the diagonal entries.</param>
        /// <remarks>
        /// The matrix may be expanded slightly to allow for additions of nonzero elements
        /// to previously non-existing diagonals.
        /// </remarks>
        public static void FastAddDiagonalScalar(this SparseMatrix matrix, double value)
        {
            int nrows = matrix.RowCount;

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            int i, j, k;

            var diagind = Helper.FindDiagonalIndices(nrows, ap, ai);

            int icount = 0;

            for (j = 0; j < nrows; j++)
            {
                if (diagind[j] < 0)
                {
                    ++icount;
                }
                else
                {
                    ax[diagind[j]] += value;
                }
            }

            // If no diagonal elements to insert in data structure, return.
            if (icount == 0)
            {
                return;
            }

            // Shift the nonzero elements if needed, to allow for created
            // diagonal elements.
            int count = ap[nrows] + icount;

            // Resize storage accordingly.
            Array.Resize(ref ai, count);
            Array.Resize(ref ax, count);

            int rowStart, rowEnd;

            bool test;

            // Copy rows backward.
            for (i = nrows - 1; i >= 0; i--)
            {
                // Go through row i.
                rowStart = ap[i];
                rowEnd = ap[i + 1];

                ap[i + 1] = count;
                test = (diagind[i] < 0);

                for (k = rowEnd - 1; k >= rowStart; k--)
                {
                    j = ai[k];
                    if (test && j < i)
                    {
                        test = false;
                        count--;

                        ax[count] = value;
                        ai[count] = i;
                        diagind[i] = count;
                    }
                    count--;
                    ax[count] = ax[k];
                    ai[count] = j;
                }

                // Diagonal element has not been added yet.
                if (test)
                {
                    count--;
                    ax[count] = value;
                    ai[count] = i;
                    diagind[i] = count;
                }
            }

            ap[0] = count;

            // Update storage references.
            storage.Values = ax;
            storage.ColumnIndices = ai;
        }

        /// <summary>
        /// Adds a diagonal matrix to a general sparse matrix B = A + Diag.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="diag">Array containing the matrix diagonal.</param>
        /// <param name="result">The resulting sparse matrix.</param>
        /// <remarks>
        /// The matrix may be expanded slightly to allow for additions of nonzero elements
        /// to previously non-existing diagonals.
        /// </remarks>
        public static void FastAddDiagonalMatrix(this SparseMatrix matrix, Vector diag, SparseMatrix result)
        {
            FastAddDiagonalMatrix(matrix, diag.Data(), result);
        }

        /// <summary>
        /// Adds a diagonal matrix to a general sparse matrix B = A + Diag.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="diag">Array containing the matrix diagonal.</param>
        /// <param name="result">The resulting sparse matrix.</param>
        /// <remarks>
        /// The matrix may be expanded slightly to allow for additions of nonzero elements
        /// to previously non-existing diagonals.
        /// </remarks>
        public static void FastAddDiagonalMatrix(this SparseMatrix matrix, double[] diag, SparseMatrix result)
        {
            int nrows = matrix.RowCount;

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            storage = result.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var bx = storage.Values;
            var bi = storage.ColumnIndices;
            var bp = storage.RowPointers;

            int i, j, k;

            // Copy int arrays into result data structure if required.
            if (!ReferenceEquals(matrix, result))
            {
                Array.Copy(ap, bp, nrows + 1);
                Array.Copy(ai, bi, ap[nrows]);
                Array.Copy(ax, bx, ap[nrows]);
            }

            // Get positions of diagonal elements in data structure.
            var diagind = Helper.FindDiagonalIndices(nrows, ap, ai);

            // Count number of holes in diagonal and add diag(*) elements to
            // valid diagonal entries.
            int icount = 0;

            for (j = 0; j < nrows; j++)
            {
                if (diagind[j] < 0)
                {
                    icount++;
                }
                else
                {
                    bx[diagind[j]] = ax[diagind[j]] + diag[j];
                }
            }

            // If no diagonal elements to insert, return.
            if (icount == 0)
            {
                return;
            }

            // Shift the nonzero elements if needed, to allow for created
            // diagonal elements.
            int k0 = bp[nrows] + icount;

            // Resize storage accordingly.
            Array.Resize(ref bi, k0);
            Array.Resize(ref bx, k0);

            int rowStart, rowEnd;
            bool test;

            // Copy rows backward.
            for (i = nrows - 1; i >= 0; i--)
            {
                // Go through row i.
                rowStart = bp[i];
                rowEnd = bp[i + 1];

                bp[i + 1] = k0;
                test = diagind[i] < 0;

                for (k = rowEnd - 1; k >= rowStart; k--)
                {
                    j = bi[k];
                    if (test && j < i)
                    {
                        test = false;
                        k0--;
                        bx[k0] = diag[i];
                        bi[k0] = i;
                        diagind[i] = k0;
                    }
                    k0--;
                    bx[k0] = ax[k];
                    bi[k0] = j;
                }

                // Diagonal element has not been added yet.
                if (test)
                {
                    k0--;
                    bx[k0] = diag[i];
                    bi[k0] = i;
                    diagind[i] = k0;
                }
            }

            bp[0] = k0;

            // Update storage references.
            storage.Values = bx;
            storage.ColumnIndices = bi;
        }

        #endregion

        #region Special Matrix Operations

        /// <summary>
        /// Make the matrix lower triangular.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="strict">If true, the diagonal will be excluded.</param>
        public static void FastLowerTriangle(this SparseMatrix matrix, bool strict = false)
        {
            if (strict)
            {
                Keep(matrix, (i, j, a) => { return i > j; });
            }
            else
            {
                Keep(matrix, (i, j, a) => { return i >= j; });
            }
        }

        /// <summary>
        /// Make the matrix upper triangular.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="strict">If true, the diagonal will be excluded.</param>
        public static void FastUpperTriangle(this SparseMatrix matrix, bool strict = false)
        {
            if (strict)
            {
                Keep(matrix, (i, j, a) => { return i < j; });
            }
            else
            {
                Keep(matrix, (i, j, a) => { return i <= j; });
            }
        }

        /// <summary>
        /// Gets the norms of each row of A.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="norm">Norm indicator (1 = 1-norm, 2 = 2-norm, 0 = max norm).</param>
        /// <returns>Vector containing the norms (output).</returns>
        public static double[] FastRowNorms(this SparseMatrix matrix, int norm)
        {
            var result = new double[matrix.RowCount];

            FastRowNorms(matrix, norm, result);

            return result;
        }

        /// <summary>
        /// Gets the norms of each row of A.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="norm">Norm indicator (1 = 1-norm, 2 = 2-norm, 0 = max norm).</param>
        /// <param name="result">Vector containing the norms (output).</param>
        public static void FastRowNorms(this SparseMatrix matrix, int norm, double[] result)
        {
            int i, k, k1, k2;
            double scal;

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;

            int nrows = matrix.RowCount;

            if (norm == 0)
            {
                for (i = 0; i < nrows; i++)
                {
                    // Compute the norm of each element.
                    scal = 0.0;
                    k1 = ap[i];
                    k2 = ap[i + 1];

                    for (k = k1; k < k2; k++)
                    {
                        scal = Math.Max(scal, Math.Abs(ax[k]));
                    }

                    result[i] = scal;
                }
            }
            else if (norm == 1)
            {

                for (i = 0; i < nrows; i++)
                {
                    // Compute the norm of each element.
                    scal = 0.0;
                    k1 = ap[i];
                    k2 = ap[i + 1];

                    for (k = k1; k < k2; k++)
                    {
                        scal += Math.Abs(ax[k]);
                    }

                    result[i] = scal;
                }
            }
            else if (norm == 2)
            {

                for (i = 0; i < nrows; i++)
                {
                    // Compute the norm of each element.
                    scal = 0.0;
                    k1 = ap[i];
                    k2 = ap[i + 1];

                    for (k = k1; k < k2; k++)
                    {
                        scal += ax[k] * ax[k];
                    }

                    result[i] = Math.Sqrt(scal);
                }
            }
            else
            {
                throw new ArgumentException("Norm not supported.");
            }
        }

        /// <summary>
        /// Gets the norms of each column of A.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="norm">Norm indicator (1 = 1-norm, 2 = 2-norm, 0 = max norm).</param>
        /// <returns>Vector containing the norms (output).</returns>
        public static double[] FastColumnNorms(this SparseMatrix matrix, int norm)
        {
            var result = new double[matrix.ColumnCount];

            FastColumnNorms(matrix, norm, result);

            return result;
        }

        /// <summary>
        /// Gets the norms of each column of A.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="norm">Norm indicator (1 = 1-norm, 2 = 2-norm, 0 = max norm).</param>
        /// <param name="result">Vector containing the norms (output).</param>
        public static void FastColumnNorms(this SparseMatrix matrix, int norm, double[] result)
        {
            int i, j;

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ai = storage.ColumnIndices;

            int ncols = matrix.ColumnCount;
            int nnz = storage.ValueCount;

            Array.Clear(result, 0, result.Length);

            if (norm == 0)
            {
                for (i = 0; i < nnz; i++)
                {
                    j = ai[i];

                    result[j] = Math.Max(result[j], Math.Abs(ax[i]));
                }
            }
            else if (norm == 1)
            {
                for (i = 0; i < nnz; i++)
                {
                    result[ai[i]] += Math.Abs(ax[i]);
                }
            }
            else if (norm == 2)
            {
                for (i = 0; i < nnz; i++)
                {
                    result[ai[i]] += ax[i] * ax[i];
                }

                for (i = 0; i < ncols; i++)
                {
                    result[i] = Math.Sqrt(result[i]);
                }
            }
            else
            {
                throw new ArgumentException("Norm not supported.");
            }
        }

        /// <summary>
        /// Scales the matrix rows, i.e. performs the matrix by matrix product B = Diag * A  (in place).
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="diag">Array representing the diagonal matrix.</param>
        /// <param name="result">Resulting matrix (can be the same as this instance).</param>
        public static void FastScaleRows(this SparseMatrix matrix, Vector diag, SparseMatrix result)
        {
            FastScaleRows(matrix, diag.Data(), result);
        }

        /// <summary>
        /// Scales the matrix rows, i.e. performs the matrix by matrix product B = Diag * A  (in place).
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="diag">Array representing the diagonal matrix.</param>
        /// <param name="result">Resulting matrix (can be the same as this instance).</param>
        public static void FastScaleRows(this SparseMatrix matrix, double[] diag, SparseMatrix result)
        {
            int n = matrix.RowCount;

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            storage = result.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var bx = storage.Values;
            var bp = storage.RowPointers;
            var bi = storage.ColumnIndices;

            int i, k, k1, k2;
            double scal;

            for (i = 0; i < n; i++)
            {
                k1 = ap[i];
                k2 = ap[i + 1];

                scal = diag[i];

                // Normalize each row.
                for (k = k1; k < k2; k++)
                {
                    bx[k] = ax[k] * scal;
                }
            }

            if (!ReferenceEquals(matrix, result))
            {
                Array.Copy(ap, bp, n + 1);
                Array.Copy(ai, bi, ap[n + 1]);
            }
        }

        /// <summary>
        /// Scales the matrix columns, i.e. performs the matrix by matrix product B = A * Diag  (in place).
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="diag">Array representing the diagonal matrix.</param>
        /// <param name="result">Resulting matrix (can be the same as this instance).</param>
        public static void FastScaleColumns(this SparseMatrix matrix, Vector diag, SparseMatrix result)
        {
            FastScaleColumns(matrix, diag.Data(), result);
        }

        /// <summary>
        /// Scales the matrix columns, i.e. performs the matrix by matrix product B = A * Diag  (in place).
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="diag">Array representing the diagonal matrix.</param>
        /// <param name="result">Resulting matrix (can be the same as this instance).</param>
        public static void FastScaleColumns(this SparseMatrix matrix, double[] diag, SparseMatrix result)
        {
            int n = matrix.RowCount;

            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            storage = result.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var bx = storage.Values;
            var bp = storage.RowPointers;
            var bi = storage.ColumnIndices;

            int i, k, k1, k2;

            for (i = 0; i < n; i++)
            {
                k1 = ap[i];
                k2 = ap[i + 1];

                // Scale each element.
                for (k = k1; k < k2; k++)
                {
                    bx[k] = ax[k] * diag[ai[k]];
                }
            }

            if (!ReferenceEquals(matrix, result))
            {
                Array.Copy(ap, bp, n + 1);
                Array.Copy(ai, bi, ap[n + 1]);
            }
        }

        /// <summary>
        /// Scales the rows of A such that their norms are 1.0 on return.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="norm">Norm indicator (1 = 1-norm, 2 = 2-norm, 0 = max norm)</param>
        public static void FastNormalizeRows(this SparseMatrix matrix, int norm)
        {
            int n = matrix.RowCount;

            var diag = new double[n];

            FastRowNorms(matrix, norm, diag);

            for (int i = 0; i < n; i++)
            {
                if (diag[i] != 0.0)
                {
                    diag[i] = 1.0 / diag[i];
                }
            }

            FastScaleRows(matrix, diag, matrix);
        }

        /// <summary>
        /// Scales the columns of A such that their norms are 1.0 on return.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="norm">Norm indicator (1 = 1-norm, 2 = 2-norm, 0 = max norm)</param>
        public static void FastNormalizeColumns(this SparseMatrix matrix, int norm)
        {
            int n = matrix.ColumnCount;

            var diag = new double[n];

            FastColumnNorms(matrix, norm, diag);

            for (int i = 0; i < n; i++)
            {
                if (diag[i] != 0.0)
                {
                    diag[i] = 1.0 / diag[i];
                }
            }

            FastScaleColumns(matrix, diag, matrix);
        }

        #endregion

        /// <summary>
        /// Creates a dense vector that contains the values from the requested matrix row.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="rowIndex">The row to copy.</param>
        /// <returns>The requested matrix row.</returns>
        public static Vector<double> FastGetRow(this SparseMatrix matrix, int rowIndex)
        {
            var result = new DenseVector(matrix.ColumnCount);

            FastGetRow(matrix, rowIndex, result);

            return result;
        }

        /// <summary>
        /// Copies a row into the given vector.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="rowIndex">The row to copy.</param>
        /// <param name="target"></param>
        /// <param name="existingData"></param>
        /// <returns>The requested matrix row.</returns>
        public static void FastGetRow(this SparseMatrix matrix, int rowIndex, DenseVector target,
            ExistingData existingData = ExistingData.Clear)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            storage.GetRow(rowIndex, (double[])target);
        }

        /// <summary>
        /// Returns a sparse sub-matrix containing the requested rows.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="rowIndices">The rows to extract.</param>
        /// <returns>The requested sub-matrix.</returns>
        public static SparseMatrix FastGetRows(this SparseMatrix matrix, int[] rowIndices)
        {
            var target = new SparseMatrix(rowIndices.Length, matrix.ColumnCount);

            matrix.FastGetRows(rowIndices, target);

            return target;
        }

        /// <summary>
        /// Extracts a sparse sub-matrix containing the requested rows.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="rowIndices">The rows to extract.</param>
        /// <param name="target">The target matrix.</param>
        /// <returns>The requested sub-matrix.</returns>
        public static void FastGetRows(this SparseMatrix matrix, int[] rowIndices,
            SparseMatrix target)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var targetStorage = target.Storage as SparseCompressedRowMatrixStorage<double>;

            if (targetStorage == null)
            {
                throw new Exception();
            }

            storage.GetRows(rowIndices, targetStorage);
        }

        /// <summary>
        /// Returns a sparse sub-matrix containing the requested columns.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="columnIndices">The columns to extract.</param>
        /// <returns>The requested sub-matrix.</returns>
        public static SparseMatrix FastGetColumns(this SparseMatrix matrix, int[] columnIndices)
        {
            var target = new SparseMatrix(matrix.RowCount, columnIndices.Length);

            matrix.FastGetColumns(columnIndices, target);

            return target;
        }

        /// <summary>
        /// Extracts a sparse sub-matrix containing the requested columns.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="columnIndices">The columns to extract.</param>
        /// <param name="target">The target matrix.</param>
        /// <returns>The requested sub-matrix.</returns>
        public static void FastGetColumns(this SparseMatrix matrix, int[] columnIndices,
            SparseMatrix target)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var targetStorage = target.Storage as SparseCompressedRowMatrixStorage<double>;

            if (targetStorage == null)
            {
                throw new Exception();
            }

            storage.GetColumns(columnIndices, targetStorage);
        }

        /// <summary>
        /// Replaces the rows of the matrix with the given rows.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="rowIndices">The indices of the rows to replace.</param>
        /// <param name="source">The sparse matrix containing the rows to copy.</param>
        public static void FastSetRows(this SparseMatrix matrix,
            int[] rowIndices, SparseMatrix source)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            var sourceStorage = source.Storage as SparseCompressedRowMatrixStorage<double>;

            if (sourceStorage == null)
            {
                throw new Exception();
            }

            storage.SetRows(rowIndices, sourceStorage);
        }

        #region Utility Functions

        /// <summary>
        /// Returns the diagonal entries of the sparse matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <returns></returns>
        public static double[] FastDiagonal(this SparseMatrix matrix)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            int i, n = storage.RowCount;

            var ax = storage.Values;
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;

            var diag = new double[n];

            // Find the indices to the diagonal entries
            for (int k = 0; k < n; k++)
            {
                i = Array.BinarySearch(ai, ap[k], ap[k + 1] - ap[k], k);

                if (i >= 0)
                {
                    diag[k] = ax[i];
                }
            }

            return diag;
        }

        /// <summary>
        /// Removes numerically zero entries from a matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="tolerance">Drop tolerance.</param>
        /// <returns>New number of non-zero entries in A.</returns>
        public static int DropZeros(this SparseMatrix matrix, double tolerance = 0.0)
        {
            Func<int, int, double, bool> func;

            if (tolerance <= 0.0)
            {
                func = (i, j, aij) =>
                {
                    return (aij != 0.0);
                };
            }
            else
            {
                func = (i, j, aij) =>
                {
                    return (Math.Abs(aij) > tolerance);
                };
            }

            return Keep(matrix, func);
        }

        /// <summary>
        /// Drops entries from a sparse matrix.
        /// </summary>
        /// <param name="matrix">The sparse matrix.</param>
        /// <param name="func">Drop element a_{i,j} if func(i, j, aij) is false.</param>
        /// <returns>New number of entries in A.</returns>
        public static int Keep(this SparseMatrix matrix, Func<int, int, double, bool> func)
        {
            var storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            if (storage == null)
            {
                throw new Exception();
            }

            return storage.Keep(func);
        }

        #endregion
    }
}
