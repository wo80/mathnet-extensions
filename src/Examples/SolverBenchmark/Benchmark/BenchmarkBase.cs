
namespace SolverBenchmark.Benchmark
{
    using System;
    using System.Globalization;
    using System.Linq;
    using System.Reflection;

    abstract class BenchmarkBase : IBenchmark
    {
        public abstract string Name
        {
            get;
        }

        public abstract BenchmarkResult Run(BenchmarkSetup config);

        public virtual void Initialize()
        {
        }

        protected object CreateInstance(Type type)
        {
            return ReflectionHelper.CreateInstance(type);
        }

        public override string ToString()
        {
            return Name;
        }
    }
}
