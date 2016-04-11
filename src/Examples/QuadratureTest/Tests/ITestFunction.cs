
namespace QuadratureTest.Tests
{
    public interface ITestFunction
    {
        /// <summary>
        /// Gets the number of calls to the eval function.
        /// </summary>
        int Count { get; }

        /// <summary>
        /// Gets the left endpoint of the interval.
        /// </summary>
        double A { get; }

        /// <summary>
        /// Gets the right endpoint of the interval.
        /// </summary>
        double B { get; }

        /// <summary>
        /// Gets the exact value of the integral.
        /// </summary>
        double ExactValue { get; }

        /// <summary>
        /// Evaluates the function at given point.
        /// </summary>
        double Eval(double x);
    }
}
