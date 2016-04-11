
namespace ExtensionsBenchmark.Benchmark
{
    using MathNet.Numerics.LinearAlgebra;
    using System;
    using Z = System.Numerics.Complex;

    static class Helper
    {
        public static bool Equals(Vector<double> v, double[] y, double eps)
        {
            var x = Util.GetData(v);

            for (int i = 0; i < v.Count; i++)
            {
                if (Math.Abs(x[i] - y[i]) > eps)
                {
                    return false;
                }
            }
            return true;
        }

        public static bool Equals(Vector<Z> v, Z[] y, double eps)
        {
            var x = Util.GetData(v);

            for (int i = 0; i < v.Count; i++)
            {
                if (Math.Abs(x[i].Real - y[i].Real) > eps ||
                    Math.Abs(x[i].Imaginary - y[i].Imaginary) > eps)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
