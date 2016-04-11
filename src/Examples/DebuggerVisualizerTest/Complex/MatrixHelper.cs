
namespace DebuggerVisualizerTest.Complex
{
    using System.Numerics;

    using MathNet.MatrixDebuggerVisualizer.Complex;
    using MathNet.MatrixDebuggerVisualizer.Services;
    using MathNet.Numerics.LinearAlgebra.Complex;
    using MathNet.Numerics.Data.Text;

    static class MatrixHelper
    {
        public static IStorageAdapter LoadMatrix(string file)
        {
            var A = (SparseMatrix)MatrixMarketReader.ReadMatrix<Complex>(file);

            return new SparseStorageAdapter(A);
        }

        public static IStorageAdapter CreateSpecial(bool laplace, int nx, int ny)
        {
            SparseMatrix A;

            if (laplace)
            {
                A = CreateSparse.Laplacian(nx, ny);
            }
            else
            {
                A = CreateSparse.Wathen(nx, ny);
            }

            return new SparseStorageAdapter(A);
        }

        public static IStorageAdapter CreateRandom(int rowCount, int columnCount,
            double density, bool symmetric)
        {
            var A = symmetric ? CreateSparse.RandomSymmetric(rowCount, density) :
                CreateSparse.Random(rowCount, columnCount, density);

            var norms = ToComplexArray((2 * A.RowNorms(1)).ToArray());

            A.FastAddDiagonalMatrix(norms, A);

            return new SparseStorageAdapter(A);
        }

        public static IStorageAdapter CreateRandomDense(int rowCount, int columnCount)
        {
            var A = (DenseMatrix)DenseMatrix.Build.Random(rowCount, columnCount);

            return new DenseStorageAdapter(A);
        }

        private static Complex[] ToComplexArray(double[] arr)
        {
            int n = arr.Length;

            var temp = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                temp[i] = arr[i];
            }

            return temp;
        }
    }
}
