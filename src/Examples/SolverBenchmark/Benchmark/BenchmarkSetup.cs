
namespace SolverBenchmark.Benchmark
{
    using MathNet.Numerics.Data.Text;
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;
    using System.Numerics;
    using System.Threading;

    public class BenchmarkSetup
    {
        public string Name { get; set; }

        public object Matrix { get; set; }

        public int RowCount { get; set; }

        public int ColumnCount { get; set; }

        public int NonZeros { get; set; }

        public bool Symmetric { get; set; }

        public int IterationsLimit { get; set; }

        public double Tolerance { get; set; }

        public object Solver { get; set; }

        public object Preconditioner { get; set; }

        public CancellationToken CancellationToken { get; set; }

        // Controlled by context.
        public BenchmarkResult Result { get; internal set; }

        public void SetMatrix<T>(object matrix)
            where T : struct, IEquatable<T>, IFormattable
        {
            var A = matrix as Matrix<T>;

            this.Matrix = matrix;
            this.RowCount = A.RowCount;
            this.ColumnCount = A.ColumnCount;
            this.Symmetric = A.IsHermitian();

            var storage = A.Storage as SparseCompressedRowMatrixStorage<T>;

            if (storage == null)
            {
                throw new Exception("Expected sparse matrix storage.");
            }

            this.NonZeros = storage.ValueCount;
        }

        public bool CheckForNull()
        {
            return Matrix != null && Solver != null && Preconditioner != null;
        }
    }
}
