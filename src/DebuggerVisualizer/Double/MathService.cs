
namespace MathNet.MatrixDebuggerVisualizer.Double
{
    using MathNet.MatrixDebuggerVisualizer.Services;
    using System;

    public class MathService : IMathService<double>
    {
        public double Abs(double value)
        {
            return Math.Abs(value);
        }

        public double Sum2(double a, double b)
        {
            return (a + b) * (a + b);
        }

        public double Difference2(double a, double b)
        {
            return (a - b) * (a - b);
        }

        public double Square(double a)
        {
            return a * a;
        }
    }
}
