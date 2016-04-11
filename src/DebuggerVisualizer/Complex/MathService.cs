
namespace MathNet.MatrixDebuggerVisualizer.Complex
{
    using MathNet.MatrixDebuggerVisualizer.Services;
    using System;
    using System.Numerics;

    public class MathService : IMathService<Complex>
    {
        public double Abs(Complex value)
        {
            return Complex.Abs(value);
        }

        public double Sum2(Complex a, Complex b)
        {
            var z = a + b;

            return (z * z).Real;
        }

        public double Difference2(Complex a, Complex b)
        {
            var z = a - b;

            return (z * z).Real;
        }

        public double Square(Complex a)
        {
            return (a * Complex.Conjugate(a)).Real;
        }
    }
}
