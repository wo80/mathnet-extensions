
namespace MathNet.MatrixDebuggerVisualizer.UI
{
    using MathNet.MatrixDebuggerVisualizer.Services;
    using MathNet.MatrixDebuggerVisualizer.UI.Views;
    using System;
    using System.Windows.Forms;

    public partial class SparseMatrixVisualizerForm : Form
    {
        public SparseMatrixVisualizerForm()
        {
            InitializeComponent();
        }

        public void SetStorageAdapter(IStorageAdapter adapter)
        {
            storageView.StorageAdapter = adapter;
            infoView.StorageAdapter = adapter;
            spyView.StorageAdapter = adapter;
            denseView.StorageAdapter = adapter;

            // Since TabIndexChanged doesn't fire on load, initialize view manually.
            storageView.ActivateView();
        }

        private void TabControl_TabIndexChanged(object sender, EventArgs e)
        {
            var view = tabControl.SelectedTab.Controls[0] as IView;

            if (view != null)
            {
                view.ActivateView();
            }
        }
    }
}
