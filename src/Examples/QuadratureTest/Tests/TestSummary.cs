using System.Collections.Generic;

namespace QuadratureTest.Tests
{
    public class TestSummary
    {
        public string TestSet { get; private set; }

        public string Integrator { get; private set; }

        public double TargetError { get; private set; }

        public int FailedCount { get; set; }

        public int TotalCount { get; set; }

        public long FunctionCalls { get; set; }

        public long Time { get; set; }

        public List<TestResult> Results { get; private set; }

        public TestSummary(ITestSet tests, IIntegrator integrator, double targetError)
        {
            this.Results = new List<TestResult>();

            this.TestSet = tests.Name;
            this.Integrator = integrator.Name;
            this.TargetError = targetError;
        }

        internal void Add(TestResult result)
        {
            this.Results.Add(result);

            this.FunctionCalls += result.Count;
        }
    }
}
