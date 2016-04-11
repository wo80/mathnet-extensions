using QuadratureTest.Tests;
using System;
using System.Globalization;
using System.Windows.Forms;

namespace QuadratureTest
{
    public partial class FormMain : Form
    {
        public FormMain()
        {
            InitializeComponent();
        }

        private void FormMain_Load(object sender, EventArgs e)
        {
            comboTestSet.Items.Add(new DasDashTestSet());
            comboTestSet.Items.Add(new EngelsTestSet());

            comboIntegrator.Items.Add(new TrapeziumIntegrator());
            comboIntegrator.Items.Add(new GaussLobattoIntegrator());
            comboIntegrator.Items.Add(new MixedGaussIntegrator());

            comboTestSet.SelectedIndex = 0;
            comboIntegrator.SelectedIndex = 0;
            comboError.SelectedIndex = 0;
        }

        private void btnRun_Click(object sender, EventArgs e)
        {
            var tests = comboTestSet.SelectedItem as ITestSet;
            var integrator = comboIntegrator.SelectedItem as IIntegrator;

            double error = double.Parse(comboError.SelectedItem.ToString(), NumberFormatInfo.InvariantInfo);

            testResultsView1.Run(tests, integrator, error);
        }
    }
}
