
namespace DebuggerVisualizerTest
{
    using MathNet.MatrixDebuggerVisualizer.Services;

    static class MatrixHelper
    {
        public static IStorageAdapter LoadMatrix(bool doubleType, string file)
        {
            if (doubleType)
            {
                return Double.MatrixHelper.LoadMatrix(file);
            }
            else
            {
                return Complex.MatrixHelper.LoadMatrix(file);
            }
        }

        public static IStorageAdapter CreateSpecial(bool doubleType, bool laplace, int nx, int ny)
        {
            if (doubleType)
            {
                return Double.MatrixHelper.CreateSpecial(laplace, nx, ny);
            }
            else
            {
                return Complex.MatrixHelper.CreateSpecial(laplace, nx, ny);
            }
        }

        public static IStorageAdapter CreateRandom(bool doubleType, int rowCount, int columnCount,
            double density, bool symmetric)
        {
            if (doubleType)
            {
                return Double.MatrixHelper.CreateRandom(rowCount, columnCount, density, symmetric);
            }
            else
            {
                return Complex.MatrixHelper.CreateRandom(rowCount, columnCount, density, symmetric);
            }
        }

        public static IStorageAdapter CreateRandomDense(bool doubleType, int rowCount, int columnCount)
        {
            if (doubleType)
            {
                return Double.MatrixHelper.CreateRandomDense(rowCount, columnCount);
            }
            else
            {
                return Complex.MatrixHelper.CreateRandomDense(rowCount, columnCount);
            }
        }
    }
}
