
namespace MathNet.MatrixDebuggerVisualizer.Complex32
{
    using MathNet.Numerics;
    using MathNet.MatrixDebuggerVisualizer.Services;
    using System;

    public class MathService : IMathService<Complex32>
    {
        public double Abs(Complex32 value)
        {
            return Complex32.Abs(value);
        }

        public double Sum2(Complex32 a, Complex32 b)
        {
            var z = a + b;

            return (z * z).Real;
        }

        public double Difference2(Complex32 a, Complex32 b)
        {
            var z = a - b;

            return (z * z).Real;
        }

        public double Square(Complex32 a)
        {
            return (a * Complex32.Conjugate(a)).Real;
        }
    }
}
