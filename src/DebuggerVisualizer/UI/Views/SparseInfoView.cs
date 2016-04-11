using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.MatrixDebuggerVisualizer.Services;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public partial class SparseInfoView : UserControl, IView
    {
        bool initialized = false;

        public IStorageAdapter StorageAdapter { get; set; }

        public SparseInfoView()
        {
            InitializeComponent();

            toolStrip.Renderer = CustomToolStripRenderer.Instance;
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
            lbValueCount.Text = string.Format("{0} ({1:0.00} %)", values, (values * 100f) / (rows * cols));
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

        private void lbUpdate0_Click(object sender, EventArgs e)
        {
            UpdateLevel1();
        }

        private void lbUpdate1_Click(object sender, EventArgs e)
        {
            UpdateLevel2();
        }

        private void lbUpdate2_Click(object sender, EventArgs e)
        {
            UpdateLevel3();
        }

        private void lbUpdate3_Click(object sender, EventArgs e)
        {
            UpdateLevel4();
        }

        private void btnRefresh_Click(object sender, EventArgs e)
        {
            UpdateLevel1();
            UpdateLevel2();
            UpdateLevel3();
            UpdateLevel4();
        }

        private bool IsSquareMatrix()
        {
            return StorageAdapter.RowCount == StorageAdapter.ColumnCount;
        }

        private void UpdateLevel1()
        {
            var info = StorageAdapter.GetStorageInfo(1);

            lbNzUpper.Text = string.Format("{0} / {1}", info.NonZerosUpper, info.NonZerosLower);
            lbNzDiag.Text = info.NonZerosDiagonal.ToString();

            lbRowWeight.Text = string.Format("{0} / {1}", info.MinRowLength, info.MaxRowLength);
            lbNzRow.Text = string.Format("\u00D8 {0:0.0}   (\u03C3 = {1:0.0})",
                info.NonZerosPerRow, info.NonZerosPerRowDev);

            lbColumnWeight.Text = string.Format("{0} / {1}", info.MinColumnLength, info.MaxColumnLength);
            lbNzColumn.Text = string.Format("\u00D8 {0:0.0}   (\u03C3 = {1:0.0})",
                info.NonZerosPerColumn, info.NonZerosPerColumnDev);

            int n = Math.Min(info.RowCount, info.ColumnCount);

            if (info.NonZerosDiagonal < n)
            {
                this.lbInfo1.Visible = true;
                this.toolTip.SetToolTip(this.lbInfo1,
                    string.Format("Matrix has {0} zeros on diagonal.", n - info.NonZerosDiagonal));
            }

            if (info.ZeroRowsCount != 0)
            {
                this.lbInfo2.Visible = true;
                this.toolTip.SetToolTip(this.lbInfo2,
                    string.Format("Matrix has {0} empty rows.", info.ZeroRowsCount));
            }

            if (info.ZeroColumnsCount != 0)
            {
                this.lbInfo3.Visible = true;
                this.toolTip.SetToolTip(this.lbInfo3,
                    string.Format("Matrix has {0} empty columns.", info.ZeroRowsCount));
            }
        }

        private void UpdateLevel2()
        {
            if (!IsSquareMatrix()) return;

            var info = StorageAdapter.GetStorageInfo(2);

            // Diagonals
            var diag = info.DiagonalsCount;
            int ndiag = 0;

            for (int i = 0; i < info.RowCount + info.ColumnCount; i++)
            {
                if (diag[i] != 0) ndiag++;
            }

            lbNumDiag.Text = ndiag.ToString();
            lbDistDiag.Text = string.Format("\u00D8 {0:0.0}   (\u03C3 = {1:0.0})",
                info.DiagonalDistanceAverage, info.DiagonalDistanceDeviation);

            lbBandUpper.Text = string.Format("{0} / {1}", info.UpperBandwidth, info.LowerBandwidth);
            lbBandAvg.Text = string.Format("\u00D8 {0:0.0}   (max = {1})",
                info.AverageBandwidth, info.MaximumBandwidth);

            lbBand90.Text = info.Band90.ToString();
            lbBand80.Text = info.Band80.ToString();
        }

        private void UpdateLevel3()
        {
            if (!IsSquareMatrix()) return;

            var info = StorageAdapter.GetStorageInfo(3);

            lbNormMax.Text = info.MaxAbsoluteValue.ToString("0.0");
            lbNormFrob.Text = info.FrobeniusNorm.ToString("0.0");
            lbNorm1.Text = info.OneNorm.ToString("0.0");
            lbNormInf.Text = info.InfNorm.ToString("0.0");

            int dd = info.DiagonallyDominantRows;

            lbDomRows.Text = dd.ToString();
            lbDomRowsPcnt.Text = string.Format("{0:0.0} %", 100.0 * dd / info.RowCount);

            dd = info.DiagonallyDominantColumns;

            lbDomColumns.Text = dd.ToString();
            lbDomColumnsPcnt.Text = string.Format("{0:0.0} %", 100.0 * dd / info.ColumnCount);
        }

        private void UpdateLevel4()
        {
            if (!IsSquareMatrix()) return;

            var info = StorageAdapter.GetStorageInfo(4);

            int imatch = info.ElementsInSymmetry;

            lbNumSym.Text = imatch.ToString();
            lbNumSymPcnt.Text = string.Format("{0:0.0} %", 100 * imatch / (double)info.ValueCount);
            lbFrobSym.Text = info.FrobeniusNormSymPart.ToString("0.0");
            lbFrobNonSym.Text = info.FrobeniusNormNonSymPart.ToString("0.0");
        }
    }
}
