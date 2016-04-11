
namespace ExtensionsBenchmark.Benchmark.Double
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using System;

    internal static class Create
    {
        public static SparseMatrix SparseMatrix(int rows, int columns, double density, bool symmetric)
        {
            if (symmetric)
            {
                return CreateSparse.RandomSymmetric(rows, density);
            }

            return CreateSparse.Random(rows, columns, density);
        }

        public static DenseVector Vector(int size, double value)
        {
            return DenseVector.Create(size, value);
        }

        public static double[] Array(int size, double value)
        {
            var a = new double[size];

            for (int i = 0; i < size; i++)
            {
                a[i] = value;
            }
            return a;
        }
    }
}
