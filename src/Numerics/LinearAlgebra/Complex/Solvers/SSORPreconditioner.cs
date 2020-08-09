// -----------------------------------------------------------------------
// <copyright file="SSORPreconditioner.cs">
// Copyright (c) 2003-2006 Bjørn-Ove Heimsund
// Copyright (c) 2006-2014 Samuel Halliday
// Copyright (c) 2014 Christian Woltering, C# version
// </copyright>
// -----------------------------------------------------------------------

#if !MIT_ONLY

namespace MathNet.Numerics.LinearAlgebra.Complex.Solvers
{
    using MathNet.Numerics.LinearAlgebra.Solvers;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;
    using System.Numerics;

    /// <summary>
    /// Symmetric successive overrelaxation (SSOR) preconditioner. 
    /// </summary>
    /// <remarks>
    /// Preconditioner for iterative methods that solve symmetric, positive definite problems.
    /// 
    /// From matrix-toolkits-java (LGPL): https://github.com/fommil/matrix-toolkits-java
    /// </remarks>
    public sealed class SSORPreconditioner : IPreconditioner<Complex>
    {
        private SparseCompressedRowMatrixStorage<Complex> storage;
        private int[] diag; // Position of diagonal entries.
        private Complex[] work;

        private double omega;
        private bool reverse;

        /// <summary>
        /// Gets or sets the overrelaxation parameter (between 0 and 2).
        /// </summary>
        public double Omega
        {
            get { return omega; }
            set
            {
                if (value < 0 || value > 2)
                {
                    throw new ArgumentException("Omega must be between 0 and 2");
                }

                omega = value;
            }
        }

        /// <summary>
        /// Indicates, whether the reverse (backward) sweep is to be done.
        /// Without this, the method is SOR instead of SSOR.
        /// </summary>
        public bool Reverse
        {
            get { return reverse; }
            set { reverse = value; }
        }

        public bool IsInitialized
        {
            get;
            private set;
        }

        public void Initialize(Matrix<Complex> matrix)
        {
            this.storage = matrix.Storage as SparseCompressedRowMatrixStorage<Complex>;

            if (storage == null)
            {
                throw new ArgumentException(Resources.MatrixMustBeSparse, "matrix");
            }

            int rowCount = storage.RowCount;

            if (rowCount != storage.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare, "matrix");
            }

            this.reverse = true;
            this.omega = 1.0;
            this.work = new Complex[rowCount];
            this.diag = storage.FindDiagonalIndices(true);

            IsInitialized = true;
        }

        /// <summary>
        /// Approximates the solution to the matrix equation <b>Ax = b</b>.
        /// </summary>
        /// <param name="input">The right hand side vector b.</param>
        /// <param name="result">The left hand side vector x.</param>
        public void Approximate(Vector<Complex> input, Vector<Complex> result)
        {
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int end, rowCount = storage.RowCount;
            Complex sum;

            result.CopyTo(work);

            // Forward sweep.
            for (int i = 0; i < rowCount; i++)
            {
                sum = 0.0;

                end = diag[i];

                for (int j = ap[i]; j < end; j++)
                {
                    sum += ax[j] * work[ai[j]];
                }

                end = ap[i + 1];

                for (int j = diag[i] + 1; j < end; j++)
                {
                    sum += ax[j] * result[ai[j]];
                }

                sum = (input[i] - sum) / ax[diag[i]];

                work[i] = result[i] + omega * (sum - result[i]);
            }

            // Stop here if the reverse sweep was not requested
            if (!reverse)
            {
                result.SetValues(work);
                return;
            }

            // The backward sweep.
            for (int i = rowCount - 1; i >= 0; i--)
            {
                sum = 0.0;

                end = diag[i];

                for (int j = ap[i]; j < end; j++)
                {
                    sum += ax[j] * work[ai[j]];
                }

                end = ap[i + 1];

                for (int j = diag[i] + 1; j < end; j++)
                {
                    sum += ax[j] * result[ai[j]];
                }

                sum = (input[i] - sum) / ax[diag[i]];

                result[i] = work[i] + omega * (sum - work[i]);
            }
        }
    }
}

#endif