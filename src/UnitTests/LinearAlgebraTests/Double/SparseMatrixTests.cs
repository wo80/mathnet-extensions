using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using NUnit.Framework;

namespace MathNet.Numerics.Extensions.UnitTests.LinearAlgebraTests.Double
{
    public class SparseMatrixTests
    {
        private const double EPS = 1e-12;
        private const int RAND_SEED = 230183;

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestMultiply1(int size)
        {
            var A = MatrixLoader.A(size);

            int rows = A.RowCount;
            int cols = A.ColumnCount;

            var x = Vector.Build.Dense(cols, 1.0);
            var y = new double[rows];

            A.Multiply(Util.GetData(x), y);

            var z = A * x;

            Assert.IsTrue(Helper.VectorEquals(z, y, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestMultiply2(int size)
        {
            var A = MatrixLoader.A(size);

            int rows = A.RowCount;
            int cols = A.ColumnCount;

            var x = Vector.Build.Dense(cols, 1.0);
            var y = new double[rows];

            A.Multiply(2.0, Util.GetData(x), 0.0, y);

            var z = 2.0 * A * x;

            Assert.IsTrue(Helper.VectorEquals(z, y, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestTransposeMultiply1(int size)
        {
            var A = MatrixLoader.A(size);

            int rows = A.RowCount;
            int cols = A.ColumnCount;

            var x = Vector.Build.Dense(rows, 1.0);
            var y = new double[cols];

            A.TransposeMultiply(Util.GetData(x), y);

            var z = A.Transpose() * x;

            Assert.IsTrue(Helper.VectorEquals(z, y, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestTransposeMultiply2(int size)
        {
            var A = MatrixLoader.A(size);

            int rows = A.RowCount;
            int cols = A.ColumnCount;

            var x = Vector.Build.Dense(rows, 1.0);
            var y = new double[cols];

            A.TransposeMultiply(2.0, Util.GetData(x), 0.0, y);

            var z = 2.0 * A.Transpose() * x;

            Assert.IsTrue(Helper.VectorEquals(z, y, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestPermuteRows(int size)
        {
            var A = MatrixLoader.A(size);

            var p = Util.Permutation(A.RowCount, RAND_SEED);

            var B = A.Clone() as SparseMatrix;
            var C = A.Clone() as SparseMatrix;

            B.PermuteRows(new Permutation(p));
            C.FastPermuteRows(p);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestPermuteColumns(int size)
        {
            var A = MatrixLoader.A(size);

            var p = Util.Permutation(A.ColumnCount, RAND_SEED);

            var B = A.Clone() as SparseMatrix;
            var C = A.Clone() as SparseMatrix;

            B.PermuteColumns(new Permutation(p));
            C.FastPermuteColumns(p);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestAdd(int size)
        {
            var A = MatrixLoader.A(size);
            var B = MatrixLoader.B(size, false);

            var C = A.Add(B);
            var D = A.FastAdd(B);

            Assert.IsTrue(MatrixComparer.Equals(C, D, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestTranspose(int size)
        {
            var A = MatrixLoader.A(size);

            var B = A.Transpose();
            var C = CreateSparse.Zeros(A.ColumnCount, A.RowCount, A.NonZerosCount);
            A.FastTranspose(C);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        public void TestLowerTriangle(int size)
        {
            var A = MatrixLoader.A(size);

            var B = A.LowerTriangle() as SparseMatrix;
            var C = A.Clone() as SparseMatrix;

            C.FastLowerTriangle();

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        public void TestUpperTriangle(int size)
        {
            var A = MatrixLoader.A(size);

            var B = A.UpperTriangle() as SparseMatrix;
            var C = A.Clone() as SparseMatrix;

            C.FastUpperTriangle();

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestRowNorms(int size)
        {
            var A = MatrixLoader.A(size);

            Vector<double> b;
            var a = new double[A.RowCount];

            b = A.RowNorms(1);
            A.FastRowNorms(1, a);
            Assert.IsTrue(Util.Equals(b, a, 0.0));

            b = A.RowNorms(2);
            A.FastRowNorms(2, a);
            Assert.IsTrue(Util.Equals(b, a, 0.0));

            b = A.RowNorms(double.PositiveInfinity);
            A.FastRowNorms(0, a);
            Assert.IsTrue(Util.Equals(b, a, 0.0));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestColumnNorms(int size)
        {
            var A = MatrixLoader.A(size);

            Vector<double> b;
            var a = new double[A.ColumnCount];

            b = A.ColumnNorms(1);
            A.FastColumnNorms(1, a);
            Assert.IsTrue(Util.Equals(b, a, 0.0));

            b = A.ColumnNorms(2);
            A.FastColumnNorms(2, a);
            Assert.IsTrue(Util.Equals(b, a, 0.0));

            b = A.ColumnNorms(double.PositiveInfinity);
            A.FastColumnNorms(0, a);
            Assert.IsTrue(Util.Equals(b, a, 0.0));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestScaleRows(int size)
        {
            var A = MatrixLoader.A(size);

            var a = A.RowNorms(1);
            var b = Util.GetData(a);
            var D = SparseMatrix.OfDiagonalVector(a);

            var B = D.Multiply(A);
            A.FastScaleRows(b, A);

            Assert.IsTrue(MatrixComparer.Equals(A, B, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestScaleColumns(int size)
        {
            var A = MatrixLoader.A(size);

            var a = A.ColumnNorms(1);
            var b = Util.GetData(a);
            var D = SparseMatrix.OfDiagonalVector(a);

            var B = A.Multiply(D);
            A.FastScaleColumns(b, A);

            Assert.IsTrue(MatrixComparer.Equals(A, B, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestNormalizeRows(int size)
        {
            var A = MatrixLoader.A(size);

            var B = A.NormalizeRows(1);
            A.FastNormalizeRows(1);

            Assert.IsTrue(MatrixComparer.Equals(A, B, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestNormalizeColumns(int size)
        {
            var A = MatrixLoader.A(size);

            var B = A.NormalizeColumns(1);
            A.FastNormalizeColumns(1);

            Assert.IsTrue(MatrixComparer.Equals(A, B, EPS));
        }

        [TestCase(55)]
        public void TestAddDiagonalScalar(int size)
        {
            var A = MatrixLoader.A(size);

            int rows = A.RowCount;
            int cols = A.ColumnCount;

            var a = DenseVector.Create(rows, 1.0);
            var D = DiagonalMatrix.OfDiagonal(rows, cols, a);

            var B = A.Add(D);
            var C = A.Clone() as SparseMatrix;

            C.FastAddDiagonalScalar(1.0);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        public void TestAddDiagonalMatrix(int size)
        {
            var A = MatrixLoader.A(size);

            int rows = A.RowCount;
            int cols = A.ColumnCount;

            var a = DenseVector.Create(rows, (i) => 1.0 + i);
            var D = DiagonalMatrix.OfDiagonal(rows, cols, a);

            var B = A.Add(D);
            var C = A.Clone() as SparseMatrix;

            C.FastAddDiagonalMatrix(Util.GetData(a), C);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [Test]
        public void TestAddDiagonalMatrixNew()
        {
            var diag = DenseVector.Create(10, 1.0);

            var A = CreateSparse.Random(10, 10, 0.8);

            A.FastLowerTriangle();
            A.SetDiagonal(diag);

            var empty = CreateSparse.Zeros(A.RowCount, A.ColumnCount, A.NonZerosCount);

            // Add diagonal and save to new matrix.
            A.FastAddDiagonalMatrix(diag.Negate().ToArray(), empty);

            var d = empty.FastDiagonal();

            for (int i = 0; i < 10; i++)
            {
                Assert.AreEqual(d[i], 0.0);
            }
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestKroneckerProduct(int size)
        {
            var A = MatrixLoader.A(size);

            var E = DenseMatrix.OfColumnMajor(2, 2, new[] { 1.0, 1.0, 1.0, 1.0 });
            var S = SparseMatrix.OfMatrix(E);

            var B = S.KroneckerProduct(A);
            var C = S.FastKroneckerProduct(A);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestDropZeros(int size)
        {
            var A = MatrixLoader.A(size);

            double threshold = 1.0e-5;

            var B = A.Clone() as SparseMatrix;

            foreach (var item in A.EnumerateIndexed(Zeros.AllowSkip))
            {
                int i = item.Item1;
                int j = item.Item2;

                if (i != j)
                {
                    // Make all off-diagonal entries small.
                    B[i, j] = threshold / 10;
                }
            }

            var C = B.Clone() as SparseMatrix;

            B.CoerceZero(threshold);
            C.DropZeros(threshold);

            Assert.IsTrue(MatrixComparer.Equals(B, C, EPS));
        }

        [TestCase(55)]
        [TestCase(57)]
        [TestCase(75)]
        public void TestKeep(int size)
        {
            var A = MatrixLoader.A(size);

            var B = A.Clone() as SparseMatrix;

            A.Keep((i, j, a) => i != j);

            foreach (var item in A.EnumerateIndexed(Zeros.AllowSkip))
            {
                // No more diagonal entries.
                Assert.IsTrue(item.Item1 != item.Item2);
            }
        }
    }
}
