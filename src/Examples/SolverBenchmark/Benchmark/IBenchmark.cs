
namespace SolverBenchmark.Benchmark
{
    public interface IBenchmark
    {
        string Name { get; }

        BenchmarkResult Run(BenchmarkSetup config);
    }
}
