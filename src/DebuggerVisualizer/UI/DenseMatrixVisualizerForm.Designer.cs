namespace MathNet.MatrixDebuggerVisualizer.UI
{
    partial class DenseMatrixVisualizerForm
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

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.tabControl = new MathNet.MatrixDebuggerVisualizer.UI.CustomTabControl();
            this.tabPage1 = new System.Windows.Forms.TabPage();
            this.storageView = new MathNet.MatrixDebuggerVisualizer.UI.Views.DenseStorageView();
            this.tabPage2 = new System.Windows.Forms.TabPage();
            this.infoView = new MathNet.MatrixDebuggerVisualizer.UI.Views.DenseInfoView();
            this.tabControl.SuspendLayout();
            this.tabPage1.SuspendLayout();
            this.tabPage2.SuspendLayout();
            this.SuspendLayout();
            // 
            // tabControl
            // 
            this.tabControl.Controls.Add(this.tabPage1);
            this.tabControl.Controls.Add(this.tabPage2);
            this.tabControl.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tabControl.Location = new System.Drawing.Point(0, 0);
            this.tabControl.Name = "tabControl";
            this.tabControl.Padding = new System.Drawing.Point(10, 5);
            this.tabControl.SelectedIndex = 0;
            this.tabControl.Size = new System.Drawing.Size(692, 523);
            this.tabControl.TabIndex = 0;
            this.tabControl.SelectedIndexChanged += new System.EventHandler(this.TabControl_TabIndexChanged);
            // 
            // tabPage1
            // 
            this.tabPage1.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(41)))), ((int)(((byte)(57)))), ((int)(((byte)(85)))));
            this.tabPage1.Controls.Add(this.storageView);
            this.tabPage1.Location = new System.Drawing.Point(4, 29);
            this.tabPage1.Name = "tabPage1";
            this.tabPage1.Size = new System.Drawing.Size(684, 490);
            this.tabPage1.TabIndex = 0;
            this.tabPage1.Text = "Values";
            // 
            // storageView
            // 
            this.storageView.BackColor = System.Drawing.Color.WhiteSmoke;
            this.storageView.Dock = System.Windows.Forms.DockStyle.Fill;
            this.storageView.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.storageView.Location = new System.Drawing.Point(0, 0);
            this.storageView.Name = "storageView";
            this.storageView.Size = new System.Drawing.Size(684, 490);
            this.storageView.StorageAdapter = null;
            this.storageView.TabIndex = 0;
            // 
            // tabPage2
            // 
            this.tabPage2.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(41)))), ((int)(((byte)(57)))), ((int)(((byte)(85)))));
            this.tabPage2.Controls.Add(this.infoView);
            this.tabPage2.Location = new System.Drawing.Point(4, 29);
            this.tabPage2.Name = "tabPage2";
            this.tabPage2.Size = new System.Drawing.Size(684, 490);
            this.tabPage2.TabIndex = 1;
            this.tabPage2.Text = "Info";
            // 
            // infoView
            // 
            this.infoView.BackColor = System.Drawing.Color.WhiteSmoke;
            this.infoView.Dock = System.Windows.Forms.DockStyle.Fill;
            this.infoView.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.infoView.ForeColor = System.Drawing.Color.White;
            this.infoView.Location = new System.Drawing.Point(0, 0);
            this.infoView.Name = "infoView";
            this.infoView.Size = new System.Drawing.Size(684, 490);
            this.infoView.StorageAdapter = null;
            this.infoView.TabIndex = 0;
            // 
            // DenseMatrixVisualizerForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(692, 523);
            this.Controls.Add(this.tabControl);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.MinimumSize = new System.Drawing.Size(600, 400);
            this.Name = "DenseMatrixVisualizerForm";
            this.Text = "Dense Matrix Visualizer";
            this.tabControl.ResumeLayout(false);
            this.tabPage1.ResumeLayout(false);
            this.tabPage2.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private CustomTabControl tabControl;
        private System.Windows.Forms.TabPage tabPage1;
        private System.Windows.Forms.TabPage tabPage2;
        private Views.DenseInfoView infoView;
        private Views.DenseStorageView storageView;

    }
}