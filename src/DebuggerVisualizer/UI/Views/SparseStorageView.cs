using MathNet.MatrixDebuggerVisualizer.Services;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    public partial class SparseStorageView : UserControl, IView
    {
        bool initialized = false;
        int current;

        public IStorageAdapter StorageAdapter { get; set; }

        public SparseStorageView()
        {
            InitializeComponent();

            toolStrip1.Renderer = CustomToolStripRenderer.Instance;

            ListViewHelper.EnableDoubleBuffer(listViewValues);
        }

        public void ActivateView()
        {
            if (!initialized)
            {
                DisplayRow(0);
                initialized = true;
            }
        }

        private void CheckRange()
        {
            bool enable = (current > 0);

            btnFirst.Enabled = enable;
            btnPrevious.Enabled = enable;

            enable = (current < StorageAdapter.RowCount - 1);

            btnNext.Enabled = enable;
            btnLast.Enabled = enable;
        }

        private void DisplayRow(int i)
        {
            if (StorageAdapter == null)
            {
                return;
            }

            if (i < 0 || i >= StorageAdapter.RowCount)
            {
                return;
            }

            current = i;

            tbCurrentRow.Text = i.ToString();

            listViewValues.Items.Clear();

            foreach (var item in StorageAdapter.EnumerateRow(i))
            {
                listViewValues.Items.Add(new ListViewItem(new string[]
                {
                    string.Empty,
                    item.Row.ToString(),
                    item.Column.ToString(),
                    item.Value.ToString(),
                    item.Z != 0.0 ? item.Z.ToString() : string.Empty
                }));
            }

            CheckRange();
        }

        private void First_Click(object sender, EventArgs e)
        {
            DisplayRow(0);
        }

        private void Previous_Click(object sender, EventArgs e)
        {
            DisplayRow(current - 1);
        }

        private void Next_Click(object sender, EventArgs e)
        {
            DisplayRow(current + 1);
        }

        private void Last_Click(object sender, EventArgs e)
        {
            DisplayRow(StorageAdapter.RowCount - 1);
        }

        private void tbCurrentRow_KeyUp(object sender, KeyEventArgs e)
        {
            if (e.KeyCode == Keys.Enter)
            {
                int i = 0;

                if (int.TryParse(tbCurrentRow.Text, out i))
                {
                    DisplayRow(i);
                }
            }
        }
    }
}
