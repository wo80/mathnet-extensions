using System.Collections.Generic;

namespace ExtensionsBenchmark.Benchmark
{
    abstract class BenchmarkBase : IBenchmark
    {
        public abstract string Name
        {
            get;
        }

        public abstract BenchmarkResult Run(BenchmarkSetup config);

        public abstract List<BenchmarkSetup> GetConfigurations();

        public void Initialize()
        {
        }

        public override string ToString()
        {
            return Name;
        }

        protected List<BenchmarkSetup> SetupTests(int start, int step, int count, int repeat = 1)
        {
            var tests = new List<BenchmarkSetup>();

            for (int i = 0; i < count; i++)
            {
                int size = start + i * step;

                tests.Add(new BenchmarkSetup(size, size, repeat));
            }

            return tests;
        }
    }
}
