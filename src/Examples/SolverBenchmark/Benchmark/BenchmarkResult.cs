
namespace SolverBenchmark.Benchmark
{
    using MathNet.Numerics.LinearAlgebra.Solvers;

    public class BenchmarkResult
    {
        /// <summary>
        /// Size of the problem.
        /// </summary>
        public int Size;

        /// <summary>
        /// Execution time of the iterative solver.
        /// </summary>
        public long Time;

        /// <summary>
        /// Total number of the iterations.
        /// </summary>
        public int IterationCount;

        /// <summary>
        /// Residual history.
        /// </summary>
        public double[] ResidualHistory;

        /// <summary>
        /// Relative error of the result.
        /// </summary>
        public double Error;

        /// <summary>
        /// Status of the iterative solver.
        /// </summary>
        public IterationStatus Status;
    }
}
