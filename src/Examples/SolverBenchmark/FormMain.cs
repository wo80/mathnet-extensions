using SolverBenchmark.Benchmark;
using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace SolverBenchmark
{
    public partial class FormMain : Form
    {
        string typeName;
        bool busy;

        ObjectCache cache;
        object matrix;

        CancellationTokenSource cts;

        public FormMain()
        {
            InitializeComponent();
        }

        private async void FormMain_Load(object sender, EventArgs e)
        {
            cache = new ObjectCache();

            comboType.SelectedIndex = 0;
            comboTolerance.SelectedIndex = 2;

            lbInfo.Text = "Initializing ...";
            await LoadTypes();
            lbInfo.Text = "Ready";
        }

        private async void comboType_SelectedIndexChanged(object sender, EventArgs e)
        {
            typeName = comboType.SelectedItem.ToString();

            await LoadTypes();
        }

        private async Task LoadTypes()
        {
            List<Type> solvers;
            List<Type> preconditioners;

            if (comboType.SelectedIndex == 0)
            {
                solvers = await ReflectionHelper.GetSolvers<double>();
                preconditioners = await ReflectionHelper.GetPreconditioners<double>();
            }
            else
            {
                solvers = await ReflectionHelper.GetSolvers<Complex>();
                preconditioners = await ReflectionHelper.GetPreconditioners<Complex>();
            }

            comboSolver.Enabled = true;
            comboSolver.Items.Clear();

            foreach (var item in solvers)
            {
                comboSolver.Items.Add(new TypeNameProvider(item));
            }

            comboSolver.SelectedIndex = 0;

            comboPrecon.Enabled = true;
            comboPrecon.Items.Clear();

            foreach (var item in preconditioners)
            {
                comboPrecon.Items.Add(new TypeNameProvider(item));
            }

            comboPrecon.SelectedIndex = 0;
        }

        private void btnReset_Click(object sender, EventArgs e)
        {
            cache.Refresh();
            comboType.SelectedIndex = 0;

            benchmarkResultsView1.ClearResults();
        }

        private void btnLoad_Click(object sender, EventArgs e)
        {
            try
            {
                var form = new FormMatrix(typeName);

                if (form.ShowDialog() == DialogResult.OK)
                {
                    matrix = form.GetMatrix();
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }
        }

        private async void btnRun_Click(object sender, EventArgs e)
        {
            if (busy)
            {
                if (cts != null) cts.Cancel();

                return;
            }

            cts = new CancellationTokenSource();

            var solver = comboSolver.SelectedItem as TypeNameProvider;
            var preconditioner = comboPrecon.SelectedItem as TypeNameProvider;

            var setup = new BenchmarkSetup();

            setup.CancellationToken = cts.Token;

            if (matrix == null)
            {
                matrix = MatrixHelper.CreateLaplacian(typeName, 50, 50);
            }

            if (typeName == "Double")
            {
                setup.SetMatrix<double>(matrix);
            }
            else
            {
                setup.SetMatrix<Complex>(matrix);
            }

            setup.Solver = cache.GetSolver(typeName, solver.Type);
            setup.Preconditioner = cache.GetPreconditioner(typeName, preconditioner.Type);
            setup.Tolerance = GetTolerance();
            setup.IterationsLimit = GetLimit();

            btnRun.Image = Properties.Resources.stop;

            try
            {
                busy = true;
                await benchmarkResultsView1.Run(setup, typeName);
            }
            catch (Exception ex)
            {
                lbInfo.Text = ex.Message;
            }
            finally
            {
                busy = false;
            }

            btnRun.Image = Properties.Resources.play;
        }

        private double GetTolerance()
        {
            double tol;

            if (double.TryParse(comboTolerance.Text, out tol))
            {
                return tol;
            }

            return 1e-8;
        }

        private int GetLimit()
        {
            int limit;

            if (int.TryParse(tbLimit.Text, out limit))
            {
                return limit;
            }

            return 200;
        }

        private void btnPropSolver_Click(object sender, EventArgs e)
        {
            var t = comboSolver.SelectedItem as TypeNameProvider;
            var solver = cache.GetSolver(typeName, t.Type);

            ShowProperties(solver);
        }

        private void btnPropPrecon_Click(object sender, EventArgs e)
        {
            var t = comboPrecon.SelectedItem as TypeNameProvider;
            var precon = cache.GetPreconditioner(typeName, t.Type);

            ShowProperties(precon);
        }

        private void ShowProperties(object obj)
        {
            if (!ReflectionHelper.HasProperties(obj))
            {
                MessageBox.Show("The selected object has no properties to configure.");
            }
            else
            {
                var form = new FormProperties();

                form.SetObject("", obj);
                form.ShowDialog();
            }
        }

        class TypeNameProvider
        {
            public Type Type;
            public string Name;

            public TypeNameProvider(Type type)
            {
                this.Type = type;
                this.Name = Util.GetPreconditionerName(type);
            }

            public override string ToString()
            {
                return this.Name;
            }
        }
    }
}
