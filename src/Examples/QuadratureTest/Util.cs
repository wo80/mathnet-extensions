
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

        public static long Toc()
        {
            timer.Stop();

            return timer.ElapsedMilliseconds;
        }
    }
}
