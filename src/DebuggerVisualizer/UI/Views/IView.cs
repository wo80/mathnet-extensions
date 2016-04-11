using MathNet.MatrixDebuggerVisualizer.Services;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public interface IView
    {
        IStorageAdapter StorageAdapter { get; set; }

        void ActivateView();
    }
}
