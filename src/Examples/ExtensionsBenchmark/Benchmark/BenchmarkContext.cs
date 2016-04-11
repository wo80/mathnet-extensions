using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace ExtensionsBenchmark.Benchmark
{
    public class BenchmarkContext
    {
        IBenchmark test;

        public string Name { get { return test.Name; } }

        public List<BenchmarkSetup> Configurations { get; private set; }

        /// <summary>
        /// Gets or sets a factor for the problem size (matrix dimensions).
        /// </summary>
        public double Factor { get; set; }

        public BenchmarkContext(IBenchmark test)
        {
            this.test = test;

            Factor = 1;
            Configurations = test.GetConfigurations();
        }

        public async Task<bool> RunTests(double density, bool symmetric, IProgress<int> progress)
        {
            return await Task<bool>.Run(() =>
            {
                int count = Configurations.Count;
                int i = 1;

                foreach (var config in Configurations)
                {
                    config.Density = density;
                    config.Symmetric = symmetric;

                    config.RowCount = (int)(Factor * config.RowCount);
                    config.ColumnCount = (int)(Factor * config.ColumnCount);

                    config.Result = test.Run(config);
                    progress.Report(i++);
                }

                return true;
            });
        }

        public Tuple<long, long> TotalTime()
        {
            return new Tuple<long, long>(
                Configurations.Aggregate(0L, (i, s) => i + (s.Result == null ? 0 : s.Result.Time1)),
                Configurations.Aggregate(0L, (i, s) => i + (s.Result == null ? 0 : s.Result.Time2))
            );
        }

        public bool Success()
        {
            return Configurations.All(s => s.Result == null ? false : s.Result.Success);
        }
    }
}
