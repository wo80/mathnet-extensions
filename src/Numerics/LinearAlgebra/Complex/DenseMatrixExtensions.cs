
namespace MathNet.Numerics.LinearAlgebra.Complex
{
    using System;
    using System.Numerics;
    using MathNet.Numerics.LinearAlgebra.Storage;

    public static class DenseMatrixExtensions
    {
        /// <summary>
        /// Compute the eigenvalues for a symmetric, generalized eigenvalue problem A*x = lambda*B*x.
        /// </summary>
        /// <param name="A">Symmetric matrix.</param>
        /// <param name="B">Symmetric, positive definite matrix.</param>
        /// <returns>A vector of eigenvalues.</returns>
        /// <remarks>
        /// See http://www.cmth.ph.ic.ac.uk/people/a.mackinnon/Lectures/compphys/node72.html
        /// and http://www.netlib.org/lapack/lug/node54.html
        /// </remarks>
        public static Vector<Complex> GeneralizedEigenvalues(this DenseMatrix A, DenseMatrix B)
        {
            // Cholesky factor of B.
            var L = B.Cholesky().Factor;

            // Compute L^-1.
            InvertLowerTriangle((DenseMatrix)L);

            // Compute L^-t = (L^-1)^t.
            var Lt = L.Transpose();

            // Build L^-1 * A * L^-t
            A.Multiply(Lt, Lt);
            L.Multiply(Lt, Lt);

            var evd = Lt.Evd(Symmetricity.Hermitian);

            return evd.EigenValues;
        }

        /// <summary>
        /// Compute the eigenvalues for a symmetric, generalized eigenvalue problem A*x = lambda*B*x.
        /// </summary>
        /// <param name="A">Symmetric matrix.</param>
        /// <param name="B">Symmetric, positive definite matrix.</param>
        /// <param name="E">Matrix containing the eigenvectors on return.</param>
        /// <returns>A vector of eigenvalues.</returns>
        /// <remarks>
        /// See http://www.cmth.ph.ic.ac.uk/people/a.mackinnon/Lectures/compphys/node72.html
        /// and http://www.netlib.org/lapack/lug/node54.html
        /// </remarks>
        public static Vector<Complex> GeneralizedEigenvalues(this DenseMatrix A, DenseMatrix B, DenseMatrix E)
        {
            // Cholesky factor of B.
            var L = B.Cholesky().Factor;

            // Compute L^-1.
            InvertLowerTriangle((DenseMatrix)L);

            // Compute L^-t = (L^-1)^t.
            var Lt = L.Transpose();

            // Save L^-t for recovery of eigenvectors.
            var copy = (DenseMatrix)Lt.Clone();

            // Build L^-1 * A * L^-t
            A.Multiply(Lt, Lt);
            L.Multiply(Lt, Lt);

            var evd = Lt.Evd(Symmetricity.Hermitian);

            // Recover eigenvectors.
            copy.Multiply(evd.EigenVectors, E);

            return evd.EigenValues;
        }

        /// <summary>
        /// Tranpose the matrix (in place).
        /// </summary>
        public static void TransposeInline(this DenseMatrix matrix)
        {
            int n = ValidateSquareMatrix(matrix);

            var data = GetData(matrix);

            Complex a, at;

            for (int j = 0; j < n; j++)
            {
                for (int i = j + 1; i < n; i++)
                {
                    a = data[i + j * n];
                    at = data[j + i * n];

                    data[i + j * n] = at;
                    data[j + i * n] = a;
                }
            }
        }

        /// <summary>
        /// Symmetrize the matrix (in place).
        /// </summary>
        /// <remarks>
        /// Computes (A + A') / 2.
        /// </remarks>
        public static void Symmetrize(this DenseMatrix matrix)
        {
            int n = ValidateSquareMatrix(matrix);

            var data = GetData(matrix);

            Complex value;

            for (int j = 0; j < n; j++)
            {
                for (int i = j + 1; i < n; i++)
                {
                    value = (data[i + j * n] + data[j + i * n]) * 0.5;

                    data[i + j * n] = value;
                    data[j + i * n] = value;
                }
            }
        }

        /// <summary>
        /// Clear the lower triangle of the matrix (in place).
        /// </summary>
        public static void ClearLowerTriangle(this DenseMatrix matrix)
        {
            int n = ValidateSquareMatrix(matrix);

            var data = GetData(matrix);

            for (int j = 0; j < n; j++)
            {
                for (int i = j + 1; i < n; i++)
                {
                    data[i + j * n] = 0.0;
                }
            }
        }

        /// <summary>
        /// Clear the upper triangle of the matrix (in place).
        /// </summary>
        public static void ClearUpperTriangle(this DenseMatrix matrix)
        {
            int n = ValidateSquareMatrix(matrix);

            var data = GetData(matrix);

            for (int j = 0; j < n; j++)
            {
                for (int i = j + 1; i < n; i++)
                {
                    data[j + i * n] = 0.0;
                }
            }
        }

        /// <summary>
        /// Compute the inverse of the lower triangle (in place).
        /// </summary>
        public static void InvertLowerTriangle(this DenseMatrix matrix)
        {
            int n = matrix.RowCount;
            Complex sum;

            for (int j = 0; j < n; j++)
            {
                matrix.At(j, j, 1.0 / matrix.At(j, j));

                for (int i = j + 1; i < n; i++)
                {
                    sum = 0.0;
                    for (int k = j; k < i; k++)
                    {

                        sum -= matrix.At(i, k) * matrix.At(k, j);
                    }

                    matrix.At(i, j, sum / matrix.At(i, i));
                }
            }
        }

        /// <summary>
        /// Compute the inverse of the upper triangle (in place).
        /// </summary>
        public static void InvertUpperTriangle(this DenseMatrix matrix)
        {
            int n = matrix.RowCount;
            Complex sum;

            for (int j = n - 1; j > -1; j--)
            {
                matrix.At(j, j, 1.0 / matrix.At(j, j));

                for (int i = j - 1; i > -1; i--)
                {
                    sum = 0.0;
                    for (int k = j; k > i; k--)
                    {

                        sum -= matrix.At(i, k) * matrix.At(k, j);
                    }

                    matrix.At(i, j, sum / matrix.At(i, i));
                }
            }
        }

        private static int ValidateSquareMatrix(DenseMatrix matrix)
        {
            int n = matrix.RowCount;

            if (matrix.ColumnCount != n)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare);
            }

            return n;
        }

        private static Complex[] GetData(Matrix<Complex> matrix)
        {
            return ((DenseColumnMajorMatrixStorage<Complex>)matrix.Storage).Data;
        }
    }
}
