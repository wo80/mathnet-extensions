namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    partial class SparseStorageView
    {
        /// <summary> 
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary> 
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Component Designer generated code

        /// <summary> 
        /// Required method for Designer support - do not modify 
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.toolStrip1 = new System.Windows.Forms.ToolStrip();
            this.btnFirst = new System.Windows.Forms.ToolStripButton();
            this.btnPrevious = new System.Windows.Forms.ToolStripButton();
            this.tbCurrentRow = new System.Windows.Forms.ToolStripTextBox();
            this.btnNext = new System.Windows.Forms.ToolStripButton();
            this.btnLast = new System.Windows.Forms.ToolStripButton();
            this.listViewValues = new System.Windows.Forms.ListView();
            this.column0 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.column1 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.column2 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.column3 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.column4 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.toolStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // toolStrip1
            // 
            this.toolStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.btnFirst,
            this.btnPrevious,
            this.tbCurrentRow,
            this.btnNext,
            this.btnLast});
            this.toolStrip1.Location = new System.Drawing.Point(0, 0);
            this.toolStrip1.Name = "toolStrip1";
            this.toolStrip1.Size = new System.Drawing.Size(500, 25);
            this.toolStrip1.TabIndex = 0;
            this.toolStrip1.Text = "toolStrip1";
            // 
            // btnFirst
            // 
            this.btnFirst.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnFirst.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.first;
            this.btnFirst.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnFirst.Name = "btnFirst";
            this.btnFirst.Size = new System.Drawing.Size(23, 22);
            this.btnFirst.Text = "First";
            this.btnFirst.Click += new System.EventHandler(this.First_Click);
            // 
            // btnPrevious
            // 
            this.btnPrevious.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnPrevious.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.previous;
            this.btnPrevious.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnPrevious.Name = "btnPrevious";
            this.btnPrevious.Size = new System.Drawing.Size(23, 22);
            this.btnPrevious.Text = "Previous";
            this.btnPrevious.Click += new System.EventHandler(this.Previous_Click);
            // 
            // tbCurrentRow
            // 
            this.tbCurrentRow.MaxLength = 16;
            this.tbCurrentRow.Name = "tbCurrentRow";
            this.tbCurrentRow.Size = new System.Drawing.Size(60, 25);
            this.tbCurrentRow.TextBoxTextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            this.tbCurrentRow.KeyUp += new System.Windows.Forms.KeyEventHandler(this.tbCurrentRow_KeyUp);
            // 
            // btnNext
            // 
            this.btnNext.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnNext.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.next;
            this.btnNext.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnNext.Name = "btnNext";
            this.btnNext.Size = new System.Drawing.Size(23, 22);
            this.btnNext.Text = "Next";
            this.btnNext.Click += new System.EventHandler(this.Next_Click);
            // 
            // btnLast
            // 
            this.btnLast.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnLast.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.last;
            this.btnLast.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnLast.Name = "btnLast";
            this.btnLast.Size = new System.Drawing.Size(23, 22);
            this.btnLast.Text = "Last";
            this.btnLast.Click += new System.EventHandler(this.Last_Click);
            // 
            // listViewValues
            // 
            this.listViewValues.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.listViewValues.Columns.AddRange(new System.Windows.Forms.ColumnHeader[] {
            this.column0,
            this.column1,
            this.column2,
            this.column3,
            this.column4});
            this.listViewValues.Dock = System.Windows.Forms.DockStyle.Fill;
            this.listViewValues.FullRowSelect = true;
            this.listViewValues.GridLines = true;
            this.listViewValues.HeaderStyle = System.Windows.Forms.ColumnHeaderStyle.Nonclickable;
            this.listViewValues.Location = new System.Drawing.Point(0, 25);
            this.listViewValues.Name = "listViewValues";
            this.listViewValues.Size = new System.Drawing.Size(500, 375);
            this.listViewValues.TabIndex = 1;
            this.listViewValues.UseCompatibleStateImageBehavior = false;
            this.listViewValues.View = System.Windows.Forms.View.Details;
            // 
            // column0
            // 
            this.column0.Text = "";
            this.column0.Width = 20;
            // 
            // column1
            // 
            this.column1.Text = "Row";
            this.column1.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // column2
            // 
            this.column2.Text = "Column";
            this.column2.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // column3
            // 
            this.column3.Text = "Value";
            this.column3.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.column3.Width = 160;
            // 
            // column4
            // 
            this.column4.Text = "";
            this.column4.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.column4.Width = 160;
            // 
            // StorageView
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.listViewValues);
            this.Controls.Add(this.toolStrip1);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.Name = "StorageView";
            this.Size = new System.Drawing.Size(500, 400);
            this.toolStrip1.ResumeLayout(false);
            this.toolStrip1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ToolStrip toolStrip1;
        private System.Windows.Forms.ToolStripButton btnFirst;
        private System.Windows.Forms.ListView listViewValues;
        private System.Windows.Forms.ColumnHeader column0;
        private System.Windows.Forms.ColumnHeader column1;
        private System.Windows.Forms.ColumnHeader column2;
        private System.Windows.Forms.ColumnHeader column3;
        private System.Windows.Forms.ColumnHeader column4;
        private System.Windows.Forms.ToolStripButton btnPrevious;
        private System.Windows.Forms.ToolStripTextBox tbCurrentRow;
        private System.Windows.Forms.ToolStripButton btnNext;
        private System.Windows.Forms.ToolStripButton btnLast;
    }
}
