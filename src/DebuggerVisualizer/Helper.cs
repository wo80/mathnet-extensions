
namespace MathNet.MatrixDebuggerVisualizer
{
    using MathNet.Numerics;
    using System;

    public static class Helper
    {
        public static int SizeOf<T>()
            where T : struct, IEquatable<T>, IFormattable
        {
            if (typeof(T) == typeof(double))
            {
                return Constants.SizeOfDouble;
            }

            if (typeof(T) == typeof(System.Numerics.Complex))
            {
                return Constants.SizeOfComplex;
            }

            if (typeof(T) == typeof(float))
            {
                return Constants.SizeOfFloat;
            }

            if (typeof(T) == typeof(MathNet.Numerics.Complex32))
            {
                return Constants.SizeOfComplex32;
            }

            return 0;
        }

        public static string GetTotalBytesString(long bytes)
        {
            if (bytes < 1024)
            {
                return string.Format("{0} bytes", bytes);
            }

            if (bytes < 1024 * 1024)
            {
                return string.Format("{0:0.0} KB", bytes / 1024.0);
            }

            if (bytes < 1024 * 1024 * 1024)
            {
                return string.Format("{0:0.0} MB", bytes / (1024.0 * 1024.0));
            }

            return string.Format("{0:0.0} GB", bytes / (1024.0 * 1024.0 * 1024.0));
        }
    }
}
