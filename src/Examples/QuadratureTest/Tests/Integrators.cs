using MathNet.Numerics.Integration;

namespace QuadratureTest.Tests
{
    public class TrapeziumIntegrator : IIntegrator
    {
        public string Name
        {
            get { return "Adaptive Trapezium Rule"; }
        }

        public double Integrate(ITestFunction testFunction, double tolerance)
        {
            return NewtonCotesTrapeziumRule.IntegrateAdaptive(testFunction.Eval,
                testFunction.A, testFunction.B, tolerance);
        }

        public override string ToString()
        {
            return Name;
        }
    }

    public class GaussLobattoIntegrator : IIntegrator
    {
        public string Name
        {
            get { return "Adaptive Gauss-Lobatto Rule"; }
        }

        public double Integrate(ITestFunction testFunction, double tolerance)
        {
            return AdaptiveGaussLobatto.Integrate(testFunction.Eval,
                testFunction.A, testFunction.B, tolerance);
        }

        public override string ToString()
        {
            return Name;
        }
    }

    public class MixedGaussIntegrator : IIntegrator
    {
        public string Name
        {
            get { return "Adaptive Mixed Gauss Rule"; }
        }

        public double Integrate(ITestFunction testFunction, double tolerance)
        {
            return AdaptiveMixedDasDash.Integrate(testFunction.Eval,
                testFunction.A, testFunction.B, tolerance, 1000);
        }

        public override string ToString()
        {
            return Name;
        }
    }
}
