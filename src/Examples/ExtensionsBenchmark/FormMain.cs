using ExtensionsBenchmark.Benchmark;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ExtensionsBenchmark
{
    public partial class FormMain : Form
    {
        IBenchmarkCollection tests;

        bool busy;

        public FormMain()
        {
            InitializeComponent();
        }

        private void FormMain_Load(object sender, EventArgs e)
        {
            comboType.SelectedIndex = 0;
            comboSymmetricity.SelectedIndex = 0;
            comboSparsity.SelectedIndex = 1;

            LoadTests();
        }

        private async void listBox_DoubleClick(object sender, EventArgs e)
        {
            var selection = listBox.SelectedItem;

            if (selection == null) return;

            var test = selection as IBenchmark;

            if (test == null) return;

            double factor = 1.0;
            double density = GetSparsity();
            bool symmetric = comboSymmetricity.SelectedIndex == 1;

            if (busy) return;

            busy = true;
            await resultsView1.RunTest(test, factor, density, symmetric);
            busy = false;
        }

        private void comboType_SelectedIndexChanged(object sender, EventArgs e)
        {
            LoadTests();
        }

        private async void btnCopy_Click(object sender, EventArgs e)
        {
            btnCopy.Enabled = false;

            var report = await resultsView1.GetReport();

            Clipboard.SetText(report);

            btnCopy.Enabled = true;
        }

        private void LoadTests()
        {
            if (comboType.SelectedIndex == 0)
            {
                tests = new Benchmark.Double.SparseMatrixTest();
            }
            else
            {
                tests = new Benchmark.Complex.SparseMatrixTest();
            }

            listBox.Items.Clear();
            listBox.Items.AddRange(tests.Items.ToArray());
        }

        private double GetSparsity()
        {
            switch (comboSparsity.SelectedIndex)
            {
                case 0:
                    return 0.02;
                case 1:
                    return 0.05;
                case 2:
                    return 0.10;
                case 3:
                    return 0.15;
                case 4:
                default:
                    break;
            }

            return 0.05;
        }
    }
}
