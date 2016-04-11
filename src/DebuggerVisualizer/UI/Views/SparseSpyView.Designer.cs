namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    partial class SparseSpyView
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
            this.btnOverview = new System.Windows.Forms.ToolStripButton();
            this.btnZoomIn = new System.Windows.Forms.ToolStripButton();
            this.btnZoomOut = new System.Windows.Forms.ToolStripButton();
            this.matrixSpyControl = new MathNet.MatrixDebuggerVisualizer.UI.Controls.MatrixSpyControl();
            this.toolStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // toolStrip1
            // 
            this.toolStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.btnOverview,
            this.btnZoomIn,
            this.btnZoomOut});
            this.toolStrip1.Location = new System.Drawing.Point(0, 0);
            this.toolStrip1.Name = "toolStrip1";
            this.toolStrip1.Size = new System.Drawing.Size(500, 25);
            this.toolStrip1.TabIndex = 0;
            this.toolStrip1.Text = "toolStrip1";
            // 
            // btnOverview
            // 
            this.btnOverview.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnOverview.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.arrow_in;
            this.btnOverview.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnOverview.Name = "btnOverview";
            this.btnOverview.Size = new System.Drawing.Size(23, 22);
            this.btnOverview.Text = "Show complete matrix";
            this.btnOverview.Click += new System.EventHandler(this.btnOverview_Click);
            // 
            // btnZoomIn
            // 
            this.btnZoomIn.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnZoomIn.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.zoom_in;
            this.btnZoomIn.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnZoomIn.Name = "btnZoomIn";
            this.btnZoomIn.Size = new System.Drawing.Size(23, 22);
            this.btnZoomIn.Text = "Zoom In";
            this.btnZoomIn.Click += new System.EventHandler(this.btnZoomIn_Click);
            // 
            // btnZoomOut
            // 
            this.btnZoomOut.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnZoomOut.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.zoom_out;
            this.btnZoomOut.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnZoomOut.Name = "btnZoomOut";
            this.btnZoomOut.Size = new System.Drawing.Size(23, 22);
            this.btnZoomOut.Text = "Zoom Out";
            this.btnZoomOut.Click += new System.EventHandler(this.btnZoomOut_Click);
            // 
            // matrixSpyControl
            // 
            this.matrixSpyControl.BackColor = System.Drawing.Color.White;
            this.matrixSpyControl.Dock = System.Windows.Forms.DockStyle.Fill;
            this.matrixSpyControl.ForeColor = System.Drawing.Color.Black;
            this.matrixSpyControl.Location = new System.Drawing.Point(0, 25);
            this.matrixSpyControl.Name = "matrixSpyControl";
            this.matrixSpyControl.Size = new System.Drawing.Size(500, 375);
            this.matrixSpyControl.TabIndex = 1;
            // 
            // SpyView
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.matrixSpyControl);
            this.Controls.Add(this.toolStrip1);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.Name = "SpyView";
            this.Size = new System.Drawing.Size(500, 400);
            this.toolStrip1.ResumeLayout(false);
            this.toolStrip1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ToolStrip toolStrip1;
        private System.Windows.Forms.ToolStripButton btnOverview;
        private Controls.MatrixSpyControl matrixSpyControl;
        private System.Windows.Forms.ToolStripButton btnZoomIn;
        private System.Windows.Forms.ToolStripButton btnZoomOut;
    }
}
