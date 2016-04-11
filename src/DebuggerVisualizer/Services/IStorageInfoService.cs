using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.MatrixDebuggerVisualizer.Services
{
    public interface IStorageInfoService
    {
        StorageInfo StorageInfo { get; }

        void Update(int level);
    }
}
