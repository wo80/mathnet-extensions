using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ExtensionsBenchmark.Benchmark
{
    public interface IBenchmarkCollection
    {
        string Name { get; }

        List<IBenchmark> Items { get; }
    }
}
