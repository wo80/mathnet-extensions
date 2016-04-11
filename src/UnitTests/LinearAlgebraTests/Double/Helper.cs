
namespace MathNet.Numerics.Extensions.UnitTests.LinearAlgebraTests.Double
{
    using MathNet.Numerics.LinearAlgebra;
    using System;

    static class Helper
    {
        public static bool VectorEquals(Vector<double> v, double[] y, double eps)
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
    }
}
