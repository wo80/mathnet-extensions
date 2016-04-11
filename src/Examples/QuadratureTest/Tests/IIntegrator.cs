
namespace QuadratureTest.Tests
{
    public interface IIntegrator
    {
        string Name { get; }

        double Integrate(ITestFunction testFunction, double tolerance);
    }
}
