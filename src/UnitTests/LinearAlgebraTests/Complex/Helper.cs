

namespace MathNet.Numerics.Extensions.UnitTests.LinearAlgebraTests.Complex
{
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Complex;
    using System.Numerics;

    static class Helper
    {
        public static bool VectorEquals(Vector<Complex> v, Complex[] y, double eps)
        {
            var x = Util.GetData(v);

            for (int i = 0; i < v.Count; i++)
            {
                if (Complex.Abs(x[i] - y[i]) > eps)
                {
                    return false;
                }
            }
            return true;
        }

        public static SparseMatrix CreateDiagonalMatrix(Vector<double> a)
        {
            var c = DenseVector.Create(a.Count, (i) => a[i]);

            return SparseMatrix.OfDiagonalVector(c);
        }

        internal static Vector<Complex> ToComplexVector(double[] a)
        {
            return Vector<Complex>.Build.Dense(a.Length, i => a[i]);
        }

        internal static Vector<Complex> ToComplexVector(Vector<double> a)
        {
            return Vector<Complex>.Build.Dense(a.Count, i => a[i]);
        }
    }
}
