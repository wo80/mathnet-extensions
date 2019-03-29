using ExtensionsBenchmark.Benchmark;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ExtensionsBenchmark
{
    public partial class BenchmarkView : UserControl
    {
        private OxyPlot.WindowsForms.PlotView plot;

        EnvironmentHelper env;

        public BenchmarkView()
        {
            InitializeComponent();
            InitializeOxyPlot();
        }

        private void listView1_SelectedIndexChanged(object sender, EventArgs e)
        {
            var items = listView1.SelectedItems;

            if (items.Count > 0)
            {
                var context = items[0].Tag as BenchmarkContext;

                if (context != null)
                {
                    VisualizeTestResult(context);
                }
            }
        }

        public async Task RunTest(IBenchmark test, double factor, double density, bool symmetric)
        {
            var context = new BenchmarkContext(test);

            var progress = new Progress<int>(i => { });

            context.Factor = factor;
            await context.RunTests(density, symmetric, progress);

            UpdateResultsList(context, factor, density, symmetric);
            VisualizeTestResult(context);
        }

        public async Task<string> GetReport()
        {
            return await Task<string>.Run(() =>
            {
                if (env == null)
                {
                    env = EnvironmentHelper.GetCurrentInfo();
                }

                var report = new Report();

                report.Append(env);

                var results = new List<BenchmarkContext>();

                foreach (ListViewItem item in listView1.Items)
                {
                    var context = item.Tag as BenchmarkContext;

                    if (context != null)
                    {
                        results.Add(context);
                    }
                }

                report.Append(results);

                return report.ToString();
            });
        }

        private void UpdateResultsList(BenchmarkContext context, double factor, double density,
            bool symmetric)
        {
            var time = context.TotalTime();

            bool success = context.Success();

            var lvi = new ListViewItem(new string[]
            {
                context.Name,
                symmetric ? "Yes" : "No",
                success ? "Yes" : "No",
                time.Item1.ToString("0") + "ms",
                time.Item2.ToString("0") + "ms",
                Report.CompareResult((double)time.Item1 / time.Item2)
            });

            if (!success)
            {
                lvi.ForeColor = Color.Red;
            }

            lvi.Tag = context;

            listView1.Items.Add(lvi);
        }

        private void VisualizeTestResult(BenchmarkContext context)
        {
            var pm = new PlotModel
            {
                Title = context.Name,
                PlotType = PlotType.XY,
                Background = OxyColors.White
            };

            var ls1 = new LineSeries();
            var ls2 = new LineSeries();

            ls1.Title = "MathNet";
            ls2.Title = "Extensions";

            int dmin = int.MaxValue;
            int dmax = 0;

            double tmax = 0;

            foreach (var item in context.Configurations)
            {
                var result = item.Result;

                dmin = Math.Min(item.RowCount, dmin);
                dmax = Math.Max(item.RowCount, dmax);

                tmax = Math.Max(result.Time1, Math.Max(result.Time2, tmax));

                ls1.Points.Add(new DataPoint(item.RowCount, result.Time1));
                ls2.Points.Add(new DataPoint(item.RowCount, result.Time2));
            }

            pm.Axes.Add(new LinearAxis()
            {
                Position = AxisPosition.Bottom,
                IsPanEnabled = false,
                IsZoomEnabled = false,
                Minimum = dmin,
                Maximum = dmax
            });

            pm.Axes.Add(new LinearAxis()
            {
                Position = AxisPosition.Left,
                IsPanEnabled = false,
                IsZoomEnabled = false,
                Minimum = 0,
                Maximum = tmax
            });

            pm.Series.Add(ls1);
            pm.Series.Add(ls2);

            plot.Model = pm;
        }

        private void InitializeOxyPlot()
        {
            this.plot = new OxyPlot.WindowsForms.PlotView();

            this.plot.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot.Location = new System.Drawing.Point(0, 0);
            this.plot.Margin = new System.Windows.Forms.Padding(0);
            this.plot.Name = "plot";
            this.plot.BackColor = Color.White;
            this.plot.Size = new System.Drawing.Size(500, 250);
            this.plot.TabIndex = 0;

            this.splitContainer.Panel2.Controls.Add(this.plot);
        }
    }
}
