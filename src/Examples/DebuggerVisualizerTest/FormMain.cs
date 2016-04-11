using MathNet.MatrixDebuggerVisualizer.UI;
using System;
using System.Windows.Forms;

namespace DebuggerVisualizerTest
{
    public partial class FormMain : Form
    {
        public FormMain()
        {
            InitializeComponent();
        }

        private void FormMain_Load(object sender, EventArgs e)
        {
            comboSparseSize.SelectedIndex = 0;
            comboWathen.SelectedIndex = 0;
            comboDenseSize.SelectedIndex = 0;
        }

        private void btnFileOpen_Click(object sender, EventArgs e)
        {
            var ofd = new OpenFileDialog();

            ofd.Filter = "MatrixMarket file (*.mtx)|*.mtx";

            if (ofd.ShowDialog() == DialogResult.OK)
            {
                comboFile.Items.Insert(0, ofd.FileName);
                comboFile.SelectedIndex = 0;
            }
        }

        private void btnSparseFile_Click(object sender, EventArgs e)
        {
            if (string.IsNullOrEmpty(comboFile.Text))
            {
                return;
            }

            var form = new SparseMatrixVisualizerForm();

            var A = MatrixHelper.LoadMatrix(radioDouble.Checked, comboFile.Text);

            form.SetStorageAdapter(A);
            form.ShowDialog();
        }

        private void btnSparseRandom_Click(object sender, EventArgs e)
        {
            var form = new SparseMatrixVisualizerForm();

            int rows, cols;
            bool symmetric;

            if (!Util.GetSparseRandomSize(comboSparseSize.Text, out rows, out cols, out symmetric))
            {
                return;
            }

            var A = MatrixHelper.CreateRandom(radioDouble.Checked, rows, cols, 0.05, symmetric);

            form.SetStorageAdapter(A);
            form.ShowDialog();
        }

        private void btnSpecial_Click(object sender, EventArgs e)
        {
            var form = new SparseMatrixVisualizerForm();

            int nx, ny;

            bool laplacian;

            if (!Util.GetSpecialSize(comboWathen.Text, out nx, out ny, out laplacian))
            {
                return;
            }

            var A = MatrixHelper.CreateSpecial(radioDouble.Checked, laplacian, nx, ny);

            form.SetStorageAdapter(A);
            form.ShowDialog();
        }

        private void btnDense_Click(object sender, EventArgs e)
        {
            var form = new DenseMatrixVisualizerForm();

            int rows, cols;

            if (!Util.GetDenseRandomSize(comboDenseSize.Text, 10, out rows, out cols))
            {
                return;
            }

            var A = MatrixHelper.CreateRandomDense(radioDouble.Checked, rows, cols);

            form.SetStorageAdapter(A);
            form.ShowDialog();
        }
    }
}
