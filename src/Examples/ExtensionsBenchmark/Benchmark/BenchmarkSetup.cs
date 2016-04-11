
namespace ExtensionsBenchmark.Benchmark
{
    public class BenchmarkSetup
    {
        // Controlled by test and context.
        public int RowCount { get; internal set; }

        // Controlled by test and context.
        public int ColumnCount { get; internal set; }

        // Controlled by test.
        public int Repeat { get; private set; }

        // Controlled by user.
        public double Density { get; set; }

        // Controlled by user.
        public bool Symmetric { get; set; }

        // Controlled by context.
        public BenchmarkResult Result { get; internal set; }

        public BenchmarkSetup(int rows, int columns, int repeat)
        {
            RowCount = rows;
            ColumnCount = columns;
            Repeat = repeat;

            Density = 0.1;
            Symmetric = false;
        }
    }
}
