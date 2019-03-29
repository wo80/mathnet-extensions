
namespace SolverBenchmark
{
    using System;
    using System.Diagnostics;

    static class Util
    {
        #region Stopwatch

        private static Stopwatch timer = new Stopwatch();

        public static void Tic()
        {
            timer.Restart();
        }

        public static double Toc()
        {
            timer.Stop();

            return TimeSpan.FromTicks(timer.ElapsedTicks).TotalMilliseconds;
        }

        #endregion

        public static string GetPreconditionerName(Type type)
        {
            var name = type.Name.Replace("Preconditioner", string.Empty);

            if (name.StartsWith("Unit"))
            {
                // Special name for UnitPreconditioner.
                name = "(None)";
            }

            return name;
        }
    }
}
