
namespace MathNet.Numerics.OdeSolvers.RK
{
    public interface IRungeKuttaStepper
    {
        int Order { get; }

        double Integrate(ref double t, ref double dt, double[] x, ref int step, int nmax, double posneg, bool dense);

        double Error(double dt, double[] x);

        double Interpolate(int i, double t, double[] con);
    }
}
