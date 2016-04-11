using SolverBenchmark.Benchmark;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using System;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Collections.Generic;

namespace SolverBenchmark
{
    public partial class BenchmarkView : UserControl
    {
        bool ignoreCheckedEvent = false;

        private OxyPlot.WindowsForms.PlotView plot;

        public BenchmarkView()
        {
            InitializeComponent();
            InitializeOxyPlot();
        }

        public async Task Run(BenchmarkSetup setup, string type)
        {
            if (!setup.CheckForNull())
            {
                return;
            }

            var context = new BenchmarkContext(setup, type);

            var progress = new Progress<int>(i => { });

            await context.Run(progress);

            AddResultToList(context);

            VisualizeTestResults();
        }

        public void ClearResults()
        {
            plot.Model = null;
            listView1.Items.Clear();
        }

        private void VisualizeTestResults()
        {
            var list = new List<BenchmarkContext>();

            foreach (ListViewItem item in listView1.CheckedItems)
            {
                list.Add((BenchmarkContext)item.Tag);
            }

            VisualizeTestResult(list);
        }

        private void AddResultToList(BenchmarkContext context)
        {
            ignoreCheckedEvent = true;

            var time = context.TotalTime();

            var config = context.Configuration;
            var result = config.Result;

            int m = config.RowCount;
            int n = config.ColumnCount;

            var item = new ListViewItem(new string[]
            {
                config.Solver.GetType().Name,
                Util.GetPreconditionerName(config.Preconditioner.GetType()),
                m.ToString(),
                n.ToString(),
                String.Format(NumberFormatInfo.InvariantInfo,
                    "{0:0.00}%", (double)config.NonZeros / (m * n)),
                config.Symmetric ? "Yes" : "No",
                result.IterationCount.ToString(),
                result.Time.ToString() + "ms",
                result.Status.ToString(),
                result.Error.ToString("0.0E00")
            });

            item.Checked = true;
            item.Tag = context;

            listView1.Items.Add(item);

            ignoreCheckedEvent = false;
        }

        private LineSeries GetLineSeries(BenchmarkContext context,
            ref double min, ref double max, ref int iterations)
        {
            var config = context.Configuration;

            iterations = Math.Max(iterations, config.Result.IterationCount);

            var ls = new LineSeries();

            ls.Title = config.Solver.GetType().Name + " "
                + Util.GetPreconditionerName(config.Preconditioner.GetType());

            int i = 0;

            var result = config.Result;

            foreach (var value in result.ResidualHistory)
            {
                min = Math.Min(value, min);
                max = Math.Max(value, max);

                ls.Points.Add(new DataPoint(i++, value));
            }

            return ls;
        }

        private void VisualizeTestResult(IEnumerable<BenchmarkContext> list)
        {
            var pm = new PlotModel
            {
                Title = string.Empty,
                PlotType = PlotType.XY,
                Background = OxyColors.White
            };

            double min = double.MaxValue;
            double max = 0;
            int iterations = 0;

            foreach (var context in list)
            {
                var ls = GetLineSeries(context, ref min, ref max, ref iterations);

                pm.Series.Add(ls);
            }

            pm.Axes.Add(new LinearAxis()
            {
                Position = AxisPosition.Bottom,
                IsPanEnabled = false,
                IsZoomEnabled = false,
                Minimum = 0,
                Maximum = iterations
            });

            pm.Axes.Add(new LogarithmicAxis()
            {
                UseSuperExponentialFormat = true,
                Position = AxisPosition.Left,
                IsPanEnabled = false,
                IsZoomEnabled = false,
                Minimum = min,
                Maximum = max
            });

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

        private void listView1_ItemChecked(object sender, ItemCheckedEventArgs e)
        {
            if (!ignoreCheckedEvent)
            {
                VisualizeTestResults();
            }
        }
    }
}
