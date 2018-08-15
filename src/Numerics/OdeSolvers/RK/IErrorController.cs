
namespace MathNet.Numerics.OdeSolvers.RK
{
    public interface IErrorController
    {
        double NextStepSize { get; }

        int Rejected { get; }

        int Accepted { get; }

        void Reset();

        bool Success(double err, double posneg, ref double h);
    }
}
