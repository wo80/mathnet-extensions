
namespace QuadratureTest.Tests
{
    using System.Collections.Generic;

    public interface ITestSet
    {
        string Name { get; }

        IEnumerable<ITestFunction> GetTestFunctions();
    }
}
