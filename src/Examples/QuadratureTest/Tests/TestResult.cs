
namespace QuadratureTest.Tests
{
    public class TestResult
    {
        public string Name { get; set; }

        public double Value { get; set; }

        public double Error { get; set; }

        public double ExactValue { get; set; }

        public int Count { get; set; }

        public bool HasInvalidValue()
        {
            return double.IsInfinity(Value) || double.IsNaN(Value);
        }
    }
}
