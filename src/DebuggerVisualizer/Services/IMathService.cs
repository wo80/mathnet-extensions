using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.MatrixDebuggerVisualizer.Services
{
    public interface IMathService<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        double Abs(T value);
        double Sum2(T a, T b);
        double Difference2(T a, T b);
        double Square(T a);
    }
}
