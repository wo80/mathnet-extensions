
namespace ExtensionsBenchmark.Benchmark
{
    public class BenchmarkResult
    {
        /// <summary>
        /// Size of the problem.
        /// </summary>
        public int Size;

        /// <summary>
        /// Execution time of the original MathNet code.
        /// </summary>
        public long Time1;

        /// <summary>
        /// Execution time of the extension code.
        /// </summary>
        public long Time2;

        /// <summary>
        /// MathNet and extension produced the same result.
        /// </summary>
        public bool Success;
    }
}
