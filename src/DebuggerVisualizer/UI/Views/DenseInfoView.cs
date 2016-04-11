using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra.Double;
using System.IO;
using System.Globalization;
using MathNet.MatrixDebuggerVisualizer.Services;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public partial class DenseInfoView : UserControl, IView
    {
        private static readonly NumberFormatInfo NumberFormat = CultureInfo.InvariantCulture.NumberFormat;

        bool initialized = false;

        public IStorageAdapter StorageAdapter { get; set; }

        public DenseInfoView()
        {
            InitializeComponent();

            toolStrip1.Renderer = CustomToolStripRenderer.Instance;
        }

        public void ActivateView()
        {
            if (!initialized)
            {
                InitializeView();
                initialized = true;
            }
        }

        private void InitializeView()
        {
            if (StorageAdapter == null)
            {
                return;
            }

            var info = StorageAdapter.GetStorageInfo(0);

            var matrix = StorageAdapter.GetMatrix();
            var storage = StorageAdapter.GetStorage();

            lbObject.Text = matrix.GetType().Name;
            lbStorageType.Text = storage.GetType().Name.Replace("`1", string.Empty);
            lbDataType.Text = StorageAdapter.GetDataTypeName();

            int rows = info.RowCount;
            int cols = info.ColumnCount;
            int values = info.ValueCount;

            lbSize.Text = string.Format("{0} x {1}", rows, cols);
            lbValueCount.Text = string.Format("{0} (100 %)", info.ValueCount);
            lbBytes.Text = Helper.GetTotalBytesString(info.TotalBytes);
        }

        private void btnSave_Click(object sender, EventArgs e)
        {
            var dlg = new SaveFileDialog();

            dlg.Filter = "MATLAB file (*.mat)|*.mat|MatrixMarket file (*.mtx)|*.mtx|CSV file (*.csv)|*.csv";

            if (dlg.ShowDialog() == DialogResult.OK)
            {
                StorageAdapter.Save(dlg.FileName);
            }
        }
    }
}
