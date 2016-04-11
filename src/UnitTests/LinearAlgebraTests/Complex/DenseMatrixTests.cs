
namespace MathNet.Numerics.Extensions.UnitTests.LinearAlgebraTests.Complex
{
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Complex;
    using NUnit.Framework;
    using System;
    using System.Linq;
    using System.Numerics;

    public class DenseMatrixTests
    {
        private const int N = 10;
        private const int RAND_SEED = 51793;

        private const double EPS = 1e-12;

        [Test]
        public void TestTransposeInline()
        {
            var A = (DenseMatrix)DenseMatrix.Build.Random(N, N, RAND_SEED);

            var S = A.Transpose();

            A.TransposeInline();

            Assert.IsTrue(S.Equals(A));
        }

        [Test]
        public void TestSymmetrize()
        {
            var A = (DenseMatrix)DenseMatrix.Build.Random(N, N, RAND_SEED);

            var S = (A + A.Transpose()) * 0.5;

            A.Symmetrize();

            Assert.IsTrue(S.Equals(A));
        }

        [Test]
        public void TestClearUpperTriangle()
        {
            var A = (DenseMatrix)DenseMatrix.Build.Random(N, N, RAND_SEED);

            A.ClearUpperTriangle();

            Complex sum = 0.0;

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = i + 1; j < A.ColumnCount; j++)
                {
                    sum += A[i, j];
                }
            }

            Assert.IsTrue(Complex.Abs(sum) < 1e-20);
        }

        [Test]
        public void TestClearLowerTriangle()
        {
            var A = (DenseMatrix)DenseMatrix.Build.Random(N, N, RAND_SEED);

            A.ClearLowerTriangle();

            Complex sum = 0.0;

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = i + 1; j < A.ColumnCount; j++)
                {
                    sum += A[j, i];
                }
            }

            Assert.IsTrue(Complex.Abs(sum) < 1e-20);
        }

        [Test]
        public void TestInvertUpperTriangle()
        {
            var A = DenseMatrix.Build.Random(N, N, RAND_SEED);

            var U = (DenseMatrix)A.UpperTriangle();
            var S = (DenseMatrix)U.Inverse();

            U.InvertUpperTriangle();

            // Shouldn't use random matrix: can be ill-conditioned.
            Assert.IsTrue(MatrixComparer.Equals(S, U, 1000 * EPS));
        }

        [Test]
        public void TestInvertLowerTriangle()
        {
            var A = DenseMatrix.Build.Random(N, N, RAND_SEED);

            var L = (DenseMatrix)A.LowerTriangle();
            var S = (DenseMatrix)L.Inverse();

            L.InvertLowerTriangle();

            // Shouldn't use random matrix: can be ill-conditioned.
            Assert.IsTrue(MatrixComparer.Equals(S, L, 1000 * EPS));
        }

        [Test]
        public void TestSymEigs()
        {
            var A = DenseMatrix.OfColumnMajor(3, 3, new Complex[]
            {
                4, 4, 8, 4, 0, 4, 8, 4, 4
            });

            var B = DenseMatrix.OfColumnMajor(3, 3, new Complex[]
            {
                2, -1, 0, -1, 2, -1, 0, -1, 2
            });

            var BA = B.Inverse() * A;

            var evd = BA.Evd(Symmetricity.Asymmetric);

            var e1 = evd.EigenValues.OrderBy(a => a.Magnitude).ToArray();
            var e2 = A.GeneralizedEigenvalues(B).OrderBy(a => a.Magnitude).ToArray();

            for (int i = 0; i < 3; i++)
            {
                var v1 = e1[i];
                var v2 = e2[i];

                Assert.IsTrue(Math.Abs(v1.Imaginary) < EPS);
                Assert.IsTrue(Math.Abs(v2.Imaginary) < EPS);

                Assert.IsTrue(Math.Abs(v1.Real - v2.Real) < EPS);
            }
        }
    }
}
