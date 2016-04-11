using SolverBenchmark.Benchmark;
using System;
using System.IO;
using System.Numerics;
using System.Text.RegularExpressions;
using System.Windows.Forms;
using MatrixMarketReader = MathNet.Numerics.Data.Text.MatrixMarketReader;

namespace SolverBenchmark
{
    public partial class FormMatrix : Form
    {
        string typeName;

        public FormMatrix()
        {
            InitializeComponent();
        }

        public FormMatrix(string typeName)
            : this()
        {
            lbType.Text = typeName;
            this.typeName = typeName;
        }

        private void FormMatrix_Load(object sender, EventArgs e)
        {
            cbLaplace.SelectedIndex = 0;
            cbRandom.SelectedIndex = 0;
        }

        private void btnFile_Click(object sender, EventArgs e)
        {
            var ofd = new OpenFileDialog();

            ofd.Filter = "MatrixMarket file (*.mtx)|*.mtx";

            if (ofd.ShowDialog() == DialogResult.OK)
            {
                tbFile.Text = ofd.FileName;
            }
        }

        public object GetMatrix()
        {
            if (rbLaplace.Checked)
            {
                return LoadLaplace();
            }

            if (rbRandom.Checked)
            {
                return LoadRandom();
            }

            throw new NotImplementedException("Loading matix from file not implemented.");
        }

        private object LoadLaplace()
        {
            int nx, ny;

            GetLaplacianSize(cbLaplace.Text, out nx, out ny);

            return MatrixHelper.CreateLaplacian(typeName, nx, ny);
        }

        private object LoadRandom()
        {
            int rows, cols;
            bool symmetric;

            GetSparseRandomSize(cbRandom.Text, out rows, out cols, out symmetric);

            if (typeName == "Double")
            {
                return MatrixHelper.CreateRandom(typeName, rows, cols, 0.05, symmetric);
            }
            else
            {
                return MatrixHelper.CreateRandom(typeName, rows, cols, 0.05, symmetric);
            }
        }

        private void LoadFile()
        {
            if (!File.Exists(tbFile.Text))
            {
                MessageBox.Show("File doesn't exist.");
                return;
            }

            if (typeName == "Double")
            {
                var A = MatrixMarketReader.ReadMatrix<double>(tbFile.Text);

                // ...
            }
            else
            {
                var A = MatrixMarketReader.ReadMatrix<Complex>(tbFile.Text);

                // ...
            }
        }

        private void GetDenseRandomSize(string text, int value, out int rows, out int cols)
        {
            rows = value;
            cols = value;

            var match = Regex.Match(text, "(?<rows>\\d+)\\s?x\\s?(?<cols>\\d+)");

            if (match.Success)
            {
                int.TryParse(match.Groups["rows"].Value.Trim(), out rows);
                int.TryParse(match.Groups["cols"].Value.Trim(), out cols);
            }
        }

        private void GetSparseRandomSize(string text, out int rows, out int cols, out bool symmetric)
        {
            GetDenseRandomSize(text, 100, out rows, out cols);

            symmetric = text.ToLower().Contains("symmetric");
        }

        private void GetLaplacianSize(string text, out int nx, out int ny)
        {
            nx = 20;
            ny = 20;

            var match = Regex.Match(text, "nx = (?<nx>\\d+), ny = (?<ny>\\d+)");

            if (match.Success)
            {
                int.TryParse(match.Groups["nx"].Value.Trim(), out nx);
                int.TryParse(match.Groups["ny"].Value.Trim(), out ny);
            }
        }
    }
}
