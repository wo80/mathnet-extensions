using MathNet.MatrixDebuggerVisualizer.Services;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Windows.Forms;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public partial class DenseStorageView : UserControl, IView
    {
        private const int MAX_ROWS = 1000;
        private const int MAX_COLUMNS = 10;
        private const string VALUE_FORMAT = "0.00000000";

        private static readonly NumberFormatInfo NumberFormat = CultureInfo.InvariantCulture.NumberFormat;

        bool initialized = false;
        int currentColumn = 0;

        public IStorageAdapter StorageAdapter { get; set; }

        public DenseStorageView()
        {
            InitializeComponent();

            toolStrip1.Renderer = CustomToolStripRenderer.Instance;

            ListViewHelper.EnableDoubleBuffer(listView);
        }

        public void ActivateView()
        {
            if (StorageAdapter == null)
            {
                return;
            }

            if (!initialized)
            {
                CreateListViewColumns();
                UpdateListView();
            }
        }

        private void CreateListViewColumns()
        {
            var info = StorageAdapter.GetStorageInfo(0);

            int rows = Math.Min(MAX_ROWS, info.RowCount);
            int columns = Math.Min(MAX_COLUMNS, info.ColumnCount);

            listView.Columns.Add("index", string.Empty, 40);

            for (int j = 0; j < columns; j++)
            {
                listView.Columns.Add("col" + j, j.ToString(), 80, HorizontalAlignment.Right, -1);
            }

            initialized = true;
        }

        private void UpdateListView()
        {
            if (!initialized) return;

            UpdateButtons();

            tbCurrentColumn.Text = currentColumn.ToString();

            var info = StorageAdapter.GetStorageInfo(0);

            int rows = Math.Min(MAX_ROWS, info.RowCount);
            int columns = Math.Min(MAX_COLUMNS, info.ColumnCount);

            var list = new List<ListViewItem>(rows);

            // First column displays the row index, so we have (columns + 1) entries.
            var array = new string[columns + 1];

            int start = currentColumn;

            if (start + columns > info.ColumnCount)
            {
                start = info.ColumnCount - columns;
            }

            for (int j = 0; j < columns; j++)
            {
                listView.Columns[j + 1].Text = (start + j).ToString();
            }

            int end = start + columns - 1;

            for (int i = 0; i < rows; i++)
            {
                array[0] = i.ToString();

                foreach (var item in StorageAdapter.EnumerateRow(i))
                {
                    int j = item.Column;
                    double value = item.Value;

                    if (j > end)
                    {
                        break;
                    }

                    if (j < start)
                    {
                        continue;
                    }

                    array[j + 1 - start] = ((Math.Abs(value) < 2e-16) ? "0"
                        : value.ToString(VALUE_FORMAT, NumberFormat));
                }

                list.Add(new ListViewItem(array));
            }

            listView.Items.Clear();
            listView.Items.AddRange(list.ToArray());
        }

        private void UpdateButtons()
        {
            var info = StorageAdapter.GetStorageInfo(0);

            btnFrist.Enabled = currentColumn > 0;
            btnPrevious.Enabled = currentColumn > 0;
            btnNext.Enabled = currentColumn < info.ColumnCount - MAX_COLUMNS;
            btnLast.Enabled = currentColumn < info.ColumnCount - MAX_COLUMNS;
        }

        private void btnFrist_Click(object sender, EventArgs e)
        {
            currentColumn = 0;

            UpdateListView();
        }

        private void btnPrevious_Click(object sender, EventArgs e)
        {
            if (currentColumn > 0)
            {
                currentColumn--;
                UpdateListView();
            }
        }

        private void btnNext_Click(object sender, EventArgs e)
        {
            var info = StorageAdapter.GetStorageInfo(0);

            if (currentColumn < info.ColumnCount - MAX_COLUMNS)
            {
                currentColumn++;
                UpdateListView();
            }
        }

        private void btnLast_Click(object sender, EventArgs e)
        {
            var info = StorageAdapter.GetStorageInfo(0);

            if (currentColumn < info.ColumnCount - MAX_COLUMNS)
            {
                currentColumn = Math.Max(0, info.ColumnCount - MAX_COLUMNS);
                UpdateListView();
            }
        }
    }
}
