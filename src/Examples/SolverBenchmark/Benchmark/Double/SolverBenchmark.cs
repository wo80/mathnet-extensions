
namespace SolverBenchmark.Benchmark.Double
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Solvers;

    class SolverBenchmark : BenchmarkBase
    {
        public override string Name
        {
            get { return string.Empty; }
        }

        public override BenchmarkResult Run(BenchmarkSetup config)
        {
            var result = new BenchmarkResult();

            if (!(config.Solver is IIterativeSolver<double>))
            {
                return result;
            }

            if (!(config.Preconditioner is IPreconditioner<double>))
            {
                return result;
            }

            var matrix = config.Matrix as SparseMatrix;

            if (matrix == null)
            {
                return result;
            }

            int n = matrix.ColumnCount;

            var s = DenseVector.Create(n, (i) => 1.0 + ((double)i / (n - 1)));
            var b = matrix.Multiply(s) as DenseVector;
            var x = DenseVector.Create(n, 0.0);

            var monitor = new IterationMonitor<double>();

            var iterator = new Iterator<double>(
                new IterationCountStopCriterion<double>(config.IterationsLimit),
                //new ResidualStopCriterion<double>(config.Tolerance, b.L2Norm()),
                new ResidualStopCriterion<double>(config.Tolerance),
                new DivergenceStopCriterion<double>(),
                new CancellationStopCriterion<double>(config.CancellationToken),
                monitor
            );

            var solver = (IIterativeSolver<double>)config.Solver;
            var preconditioner = (IPreconditioner<double>)config.Preconditioner;

            preconditioner.Initialize(matrix);

            Util.Tic();
            solver.Solve(matrix, b, x, iterator, preconditioner);

            var residual = b.Clone();

            matrix.Multiply(-1.0, x, 1.0, residual);

            result.Time = Util.Toc();
            result.Status = iterator.Status;
            result.IterationCount = monitor.LastIteration;
            result.ResidualHistory = monitor.ResidualNormHistory.ToArray();
            result.Error = residual.L2Norm() / b.L2Norm();

            return result;
        }
    }
}
