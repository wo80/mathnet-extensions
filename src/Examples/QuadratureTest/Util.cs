
namespace QuadratureTest
{
    using System;
    using System.Diagnostics;

    static class Util
    {
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
    }
}
