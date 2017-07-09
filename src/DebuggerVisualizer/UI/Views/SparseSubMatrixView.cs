using MathNet.MatrixDebuggerVisualizer.Services;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Windows.Forms;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public partial class SparseSubMatrixView : UserControl, IView
    {
        private const int MAX_ROWS = 1000;
        private const int MAX_COLUMNS = 10;
        private const string VALUE_FORMAT = "0.00000000";

        private static readonly NumberFormatInfo NumberFormat = CultureInfo.InvariantCulture.NumberFormat;

        bool initialized = false;

        int currentRow = 0;
        int currentColumn = 0;
        int currentWidth = 10;
        int currentHeight = 10;

        public IStorageAdapter StorageAdapter { get; set; }

        public SparseSubMatrixView()
        {
            InitializeComponent();

            toolStrip1.Renderer = CustomToolStripRenderer.Instance;
            comboPlaces.SelectedIndex = 0;

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
                GetSelectedRange();
                CreateListViewColumns();
                UpdateListView();
            }
        }

        private void btnUpdate_Click(object sender, EventArgs e)
        {
            ValidateRange();

            if (currentWidth + 1 != listView.Columns.Count)
            {
                CreateListViewColumns();
            }

            UpdateListView();
        }

        private void CreateListViewColumns()
        {
            var info = StorageAdapter.GetStorageInfo(0);

            int rows = Math.Min(currentHeight, info.RowCount);
            int columns = Math.Min(currentWidth, info.ColumnCount);

            listView.Items.Clear();
            listView.Columns.Clear();

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

            if (!ValidateRange()) return;

            var info = StorageAdapter.GetStorageInfo(0);

            int rows = Math.Min(currentHeight, info.RowCount);
            int columns = Math.Min(currentWidth, info.ColumnCount);

            var list = new List<ListViewItem>(rows);

            // Store the complete sub-matrix as string array.
            var array = new string[rows][];

            int rowStart = currentRow;
            int columnStart = currentColumn;

            if (rowStart + rows > info.RowCount)
            {
                rowStart = info.RowCount - rows;
            }

            if (columnStart + columns > info.ColumnCount)
            {
                columnStart = info.ColumnCount - columns;
            }

            for (int j = 0; j < columns; j++)
            {
                listView.Columns[j + 1].Text = (columnStart + j).ToString();
            }

            for (int i = 0; i < rows; i++)
            {
                // First column displays the row index, so we have (columns + 1) entries.
                var s = array[i] = new string[columns + 1];

                s[0] = i.ToString();

                for (int j = 1; j < columns + 1; j++)
                {
                    s[j] = "0"; // Display default value.
                }
            }

            string valueFormat = GetNumberFormat();

            foreach (var item in StorageAdapter.EnumerateSubmatrix(
                rowStart, rowStart + rows,
                columnStart, columnStart + columns - 1))
            {
                int i = item.Row;
                int j = item.Column;
                double value = item.Value;

                var s = array[i - rowStart];

                s[j + 1 - columnStart] = ((Math.Abs(value) < 2e-16) ? "0"
                    : value.ToString(valueFormat, NumberFormat));
            }

            for (int i = 0; i < rows; i++)
            {
                list.Add(new ListViewItem(array[i]));
            }

            listView.Items.Clear();
            listView.Items.AddRange(list.ToArray());
        }

        private bool ValidateRange()
        {
            var info = StorageAdapter.GetStorageInfo(0);

            GetSelectedRange();

            // Rows:

            if (currentRow < 0)
            {
                currentRow = 0;
            }

            if (currentRow >= info.RowCount)
            {
                currentRow = info.RowCount;
            }

            // Columns:

            if (currentColumn < 0)
            {
                currentColumn = 0;
            }

            if (currentColumn >= info.ColumnCount)
            {
                currentColumn = info.ColumnCount;
            }

            // Width:

            if (currentWidth < 1)
            {
                currentWidth = 1;
            }

            if (currentWidth > MAX_COLUMNS)
            {
                currentWidth = MAX_COLUMNS;
            }

            // Height:

            if (currentHeight < 1)
            {
                currentHeight = 1;
            }

            if (currentHeight > MAX_ROWS)
            {
                currentHeight = MAX_ROWS;
            }

            // Update text boxes:
            tbRow.Text = currentRow.ToString();
            tbColumn.Text = currentColumn.ToString();
            tbWidth.Text = currentWidth.ToString();
            tbHeight.Text = currentHeight.ToString();

            return true;
        }

        private string GetNumberFormat()
        {
            switch (comboPlaces.SelectedIndex)
            {
                case 0:
                    return VALUE_FORMAT;
                case 1:
                    return "0.0";
                case 2:
                    return "0.00";
                case 3:
                    return "0.000";
                case 4:
                    return "0.0000";
                case 5:
                    return "0.00000";
                default:
                    break;
            }

            return VALUE_FORMAT;
        }

        private void GetSelectedRange()
        {
            if (!int.TryParse(tbRow.Text, out currentRow))
            {
                currentRow = 0;
            }

            if (!int.TryParse(tbColumn.Text, out currentColumn))
            {
                currentColumn = 0;
            }

            if (!int.TryParse(tbWidth.Text, out currentWidth))
            {
                currentWidth = 10;
            }

            if (!int.TryParse(tbHeight.Text, out currentHeight))
            {
                currentHeight = 10;
            }
        }
    }
}
