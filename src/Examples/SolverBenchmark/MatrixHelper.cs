
namespace SolverBenchmark
{
    using System.Numerics;

    using CreateDouble = MathNet.Numerics.LinearAlgebra.Double.CreateSparse;
    using CreateComplex = MathNet.Numerics.LinearAlgebra.Complex.CreateSparse;

    using DoubleExtensions = MathNet.Numerics.LinearAlgebra.Double.SparseMatrixExtensions;
    using ComplexExtensions = MathNet.Numerics.LinearAlgebra.Complex.SparseMatrixExtensions;

    static class MatrixHelper
    {
        public static object CreateLaplacian(string type, int nx, int ny)
        {
            var doubleType = !type.StartsWith("Complex");

            if (doubleType)
            {
                return MathNet.Numerics.LinearAlgebra.Double.CreateSparse.Laplacian(nx, ny);
            }
            else
            {
                return MathNet.Numerics.LinearAlgebra.Complex.CreateSparse.Laplacian(nx, ny);
            }
        }

        public static object CreateRandom(string type, int rowCount, int columnCount,
            double density, bool symmetric)
        {
            var doubleType = !type.StartsWith("Complex");

            if (doubleType)
            {
                var A = symmetric ? CreateDouble.RandomSymmetric(rowCount, density) :
                    CreateDouble.Random(rowCount, columnCount, density);

                var norms = (2 * A.RowNorms(1)).ToArray();

                DoubleExtensions.FastAddDiagonalMatrix(A, norms, A);

                return A;
            }
            else
            {
                var A = symmetric ? CreateComplex.RandomSymmetric(rowCount, density) :
                    CreateComplex.Random(rowCount, columnCount, density);

                var norms = ToComplexArray((2 * A.RowNorms(1)).ToArray());

                ComplexExtensions.FastAddDiagonalMatrix(A, norms, A);

                return A;
            }
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
