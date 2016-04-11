
namespace MathNet.Numerics.LinearAlgebra.Double
{
    using System;
    using MathNet.Numerics.LinearAlgebra.Storage;

    public static class MatrixComparer
    {
        public static bool Equals(Matrix<double> A, Matrix<double> B, double eps = 0.0)
        {
            if (A == null || B == null)
            {
                return false;
            }

            if (A is DenseMatrix && B is DenseMatrix)
            {
                return Equals((DenseMatrix)A, (DenseMatrix)B, eps);
            }

            if (A is SparseMatrix && B is SparseMatrix)
            {
                return Equals((SparseMatrix)A, (SparseMatrix)B, eps);
            }

            return false;
        }

        #region DenseMatrix comparison

        private static bool Equals(DenseMatrix A, DenseMatrix B, double eps = 0.0)
        {
            if (A == null || B == null)
            {
                return false;
            }

            var sa = A.Storage as DenseColumnMajorMatrixStorage<double>;
            var sb = B.Storage as DenseColumnMajorMatrixStorage<double>;

            if (sa == null || sb == null)
            {
                return A.Equals(B);
            }

            int nnz = sa.Data.Length;

            return CompareValues(nnz, sa.Data, sb.Data, eps);
        }

        #endregion

        #region SparseMatrix comparison

        private static bool Equals(SparseMatrix A, SparseMatrix B, double eps = 0.0)
        {
            if (A == null || B == null)
            {
                return false;
            }

            var sa = A.Storage as SparseCompressedRowMatrixStorage<double>;
            var sb = B.Storage as SparseCompressedRowMatrixStorage<double>;

            if (sa == null || sb == null)
            {
                return A.Equals(B);
            }

            int nnz = sa.ValueCount;

            if (!CompareStructure(nnz, sa.RowPointers, sa.ColumnIndices, sb.RowPointers, sb.ColumnIndices))
            {
                return false;
            }

            return CompareValues(nnz, sa.Values, sb.Values, eps);
        }

        private static bool CompareStructure(int nnz, int[] ap, int[] ai, int[] bp, int[] bi, bool structure = false)
        {
            int n = ap.Length;

            if (n != bp.Length)
            {
                return false;
            }

            for (int i = 0; i < n; i++)
            {
                if (ap[i] != bp[i])
                {
                    return false;
                }
            }

            n = structure ? ai.Length : nnz;

            if (structure && n != bi.Length)
            {
                return false;
            }

            for (int i = 0; i < n; i++)
            {
                if (ai[i] != bi[i])
                {
                    return false;
                }
            }

            return true;
        }

        private static bool CompareValues(int nnz, double[] ax, double[] bx, double eps, bool structure = false)
        {
            int n = structure ? ax.Length : nnz;

            if (structure && n != bx.Length)
            {
                return false;
            }

            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(ax[i] - bx[i]) > eps)
                {
                    return false;
                }
            }

            return true;
        }

        #endregion
    }
}
