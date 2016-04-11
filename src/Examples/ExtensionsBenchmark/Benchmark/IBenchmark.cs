using System.Collections.Generic;

namespace ExtensionsBenchmark.Benchmark
{
    public interface IBenchmark
    {
        string Name { get; }

        BenchmarkResult Run(BenchmarkSetup config);

        List<BenchmarkSetup> GetConfigurations();
    }
}
