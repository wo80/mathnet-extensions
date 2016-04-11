using MathNet.MatrixDebuggerVisualizer.Services;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public partial class SparseSpyView : UserControl, IView
    {
        bool initialized = false;

        public IStorageAdapter StorageAdapter { get; set; }

        public SparseSpyView()
        {
            InitializeComponent();

            toolStrip1.Renderer = CustomToolStripRenderer.Instance;
        }

        public void ActivateView()
        {
            if (!initialized)
            {
                matrixSpyControl.Initialize(StorageAdapter);
                initialized = true;
            }
        }

        private void btnOverview_Click(object sender, EventArgs e)
        {
            matrixSpyControl.Reset();
        }

        private void btnZoomIn_Click(object sender, EventArgs e)
        {
            matrixSpyControl.Zoom(1, 0.5f, 0.5f);
        }

        private void btnZoomOut_Click(object sender, EventArgs e)
        {
            matrixSpyControl.Zoom(-1, 0.5f, 0.5f);
        }
    }
}
