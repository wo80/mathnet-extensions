
namespace MathNet.MatrixDebuggerVisualizer.Single
{
    using MathNet.MatrixDebuggerVisualizer.Services;
    using System;

    public class MathService : IMathService<float>
    {
        public double Abs(float value)
        {
            return Math.Abs(value);
        }

        public double Sum2(float a, float b)
        {
            return (a + b) * (a + b);
        }

        public double Difference2(float a, float b)
        {
            return (a - b) * (a - b);
        }

        public double Square(float a)
        {
            return a * a;
        }
    }
}
