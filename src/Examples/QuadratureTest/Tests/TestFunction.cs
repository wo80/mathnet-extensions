
namespace QuadratureTest.Tests
{
    using System;

    public abstract class TestFunction : ITestFunction
    {
        protected int count;

        protected double a;
        protected double b;
        protected double exact;

        public int Count
        {
            get { return count; }
        }

        public double A
        {
            get { return a; }
        }

        public double B
        {
            get { return b; }
        }

        public double ExactValue
        {
            get { return exact; }
        }

        public abstract double Eval(double x);

        public TestFunction()
            : this(0.0, 1.0, 0.0)
        {
        }

        public TestFunction(double exact)
            : this(0.0, 1.0, exact)
        {
        }

        public TestFunction(double a, double b, double exact)
        {
            this.count = 0;

            this.a = a;
            this.b = b;
            this.exact = exact;
        }
    }
}
