using QuadratureTest.Tests;
using System;
using System.Drawing;
using System.Windows.Forms;

namespace QuadratureTest
{
    public partial class TestResultsView : UserControl
    {
        public TestResultsView()
        {
            InitializeComponent();
        }

        public void Run(ITestSet tests, IIntegrator integrator, double targetError)
        {
            listView1.Items.Clear();

            var summary = new TestSummary(tests, integrator, targetError);

            var functions = tests.GetTestFunctions();

            int total = 0, failed = 0;

            Util.Tic();

            foreach (var item in functions)
            {
                total++;

                try
                {
                    var value = integrator.Integrate(item, targetError);
                    Console.WriteLine(value);

                    var result = GetTestResult(item, value);

                    ShowTestResult(result, targetError);

                    summary.Add(result);

                    if (result.HasInvalidValue())
                    {
                        failed++;
                    }
                }
                catch (Exception e)
                {
                    ProcessError(item, e);

                    failed++;
                }
            }

            summary.Time = Util.Toc();

            summary.TotalCount = total;
            summary.FailedCount = failed;

            ShowTestSummary(summary);
        }

        private void ShowTestSummary(TestSummary summary)
        {
            var lvi = new ListViewItem(new string[]
            {
                summary.Integrator,
                summary.TestSet,
                summary.TargetError.ToString("0.0e00"),
                summary.TotalCount.ToString(),
                summary.FailedCount.ToString(),
                summary.Time.ToString(),
                summary.FunctionCalls.ToString()
            });

            listView2.Items.Add(lvi);
        }

        private TestResult GetTestResult(ITestFunction item, double value)
        {
            var result = new TestResult();

            double error = Math.Abs(value - item.ExactValue);

            result.Name = item.GetType().Name;
            result.Value = value;
            result.Error = error;
            result.ExactValue = item.ExactValue;
            result.Count = item.Count;

            return result;
        }

        private void ShowTestResult(TestResult result, double targetError)
        {
            double relativeError =  result.Error / Math.Abs(result.ExactValue);

            if (result.ExactValue == 0)
            {
                relativeError = result.Error;
            }

            double errorRatio = result.Error / targetError;

            var lvi = new ListViewItem(new string[]
            {
                result.Name,
                result.Value.ToString(),
                result.ExactValue.ToString(),
                result.Error.ToString("0.0e00"),
                relativeError.ToString("0.0e00"),
                result.Count.ToString()
            });

            if (result.HasInvalidValue())
            {
                lvi.ForeColor = Color.Red;
            }
            else if (errorRatio > 2.5)
            {
                lvi.ForeColor = Color.Peru;
            }

            listView1.Items.Add(lvi);
        }

        private void ProcessError(ITestFunction item, Exception e)
        {
            var lvi = new ListViewItem(new string[]
            {
                item.GetType().Name,
                "FAILED",
                String.Empty,
                String.Empty,
                String.Empty,
                String.Empty,
                String.Empty
            });

            lvi.ForeColor = Color.Red;

            listView1.Items.Add(lvi);
        }
    }
}
