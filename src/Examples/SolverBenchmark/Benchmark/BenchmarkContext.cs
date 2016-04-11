
namespace SolverBenchmark.Benchmark
{
    using System;
    using System.Threading.Tasks;

    public class BenchmarkContext
    {
        IBenchmark benchmark;

        public string Name { get { return benchmark.Name; } }

        public BenchmarkSetup Configuration { get; private set; }

        public BenchmarkContext(BenchmarkSetup setup, string type)
        {
            if (type.Equals("double", StringComparison.OrdinalIgnoreCase))
            {
                this.benchmark = new Double.SolverBenchmark();
            }
            else
            {
                this.benchmark = new Complex.SolverBenchmark();
            }

            Configuration = setup;
        }

        public async Task<bool> Run(IProgress<int> progress)
        {
            return await Task<bool>.Run(() =>
            {
                var config = Configuration;

                config.Result = benchmark.Run(config);
                //progress.Report(i++);

                return true;
            });
        }

        public long TotalTime()
        {
            return Configuration.Result.Time;
        }
    }
}
