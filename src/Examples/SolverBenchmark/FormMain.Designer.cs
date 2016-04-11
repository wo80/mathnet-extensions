namespace SolverBenchmark
{
    partial class FormMain
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
            this.toolStrip1 = new System.Windows.Forms.ToolStrip();
            this.toolStripLabel1 = new System.Windows.Forms.ToolStripLabel();
            this.comboType = new System.Windows.Forms.ToolStripComboBox();
            this.toolStripSeparator1 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripLabel2 = new System.Windows.Forms.ToolStripLabel();
            this.comboSolver = new System.Windows.Forms.ToolStripComboBox();
            this.toolStripSeparator2 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripLabel3 = new System.Windows.Forms.ToolStripLabel();
            this.comboPrecon = new System.Windows.Forms.ToolStripComboBox();
            this.toolStripSeparator3 = new System.Windows.Forms.ToolStripSeparator();
            this.statusStrip1 = new System.Windows.Forms.StatusStrip();
            this.labelInfo = new System.Windows.Forms.ToolStripStatusLabel();
            this.btnPropSolver = new System.Windows.Forms.ToolStripButton();
            this.btnPropPrecon = new System.Windows.Forms.ToolStripButton();
            this.btnLoad = new System.Windows.Forms.ToolStripButton();
            this.btnRun = new System.Windows.Forms.ToolStripButton();
            this.btnReset = new System.Windows.Forms.ToolStripButton();
            this.toolStripSeparator4 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripLabel4 = new System.Windows.Forms.ToolStripLabel();
            this.comboTolerance = new System.Windows.Forms.ToolStripComboBox();
            this.benchmarkResultsView1 = new SolverBenchmark.BenchmarkView();
            this.toolStripLabel5 = new System.Windows.Forms.ToolStripLabel();
            this.tbLimit = new System.Windows.Forms.ToolStripTextBox();
            this.toolStripSeparator5 = new System.Windows.Forms.ToolStripSeparator();
            this.lbInfo = new System.Windows.Forms.ToolStripStatusLabel();
            this.toolStrip1.SuspendLayout();
            this.statusStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // toolStrip1
            // 
            this.toolStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.btnReset,
            this.toolStripSeparator4,
            this.toolStripLabel1,
            this.comboType,
            this.toolStripSeparator1,
            this.toolStripLabel2,
            this.comboSolver,
            this.btnPropSolver,
            this.toolStripSeparator2,
            this.toolStripLabel3,
            this.comboPrecon,
            this.btnPropPrecon,
            this.toolStripSeparator3,
            this.toolStripLabel4,
            this.comboTolerance,
            this.toolStripLabel5,
            this.tbLimit,
            this.toolStripSeparator5,
            this.btnLoad,
            this.btnRun});
            this.toolStrip1.Location = new System.Drawing.Point(0, 0);
            this.toolStrip1.Name = "toolStrip1";
            this.toolStrip1.Size = new System.Drawing.Size(842, 25);
            this.toolStrip1.TabIndex = 0;
            this.toolStrip1.Text = "toolStrip1";
            // 
            // toolStripLabel1
            // 
            this.toolStripLabel1.Name = "toolStripLabel1";
            this.toolStripLabel1.Size = new System.Drawing.Size(35, 22);
            this.toolStripLabel1.Text = "Type:";
            // 
            // comboType
            // 
            this.comboType.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboType.Items.AddRange(new object[] {
            "Double",
            "Complex"});
            this.comboType.Name = "comboType";
            this.comboType.Size = new System.Drawing.Size(75, 25);
            this.comboType.SelectedIndexChanged += new System.EventHandler(this.comboType_SelectedIndexChanged);
            // 
            // toolStripSeparator1
            // 
            this.toolStripSeparator1.Name = "toolStripSeparator1";
            this.toolStripSeparator1.Size = new System.Drawing.Size(6, 25);
            // 
            // toolStripLabel2
            // 
            this.toolStripLabel2.Name = "toolStripLabel2";
            this.toolStripLabel2.Size = new System.Drawing.Size(41, 22);
            this.toolStripLabel2.Text = "Solver:";
            // 
            // comboSolver
            // 
            this.comboSolver.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboSolver.Enabled = false;
            this.comboSolver.Name = "comboSolver";
            this.comboSolver.Size = new System.Drawing.Size(100, 25);
            // 
            // toolStripSeparator2
            // 
            this.toolStripSeparator2.Name = "toolStripSeparator2";
            this.toolStripSeparator2.Size = new System.Drawing.Size(6, 25);
            // 
            // toolStripLabel3
            // 
            this.toolStripLabel3.Name = "toolStripLabel3";
            this.toolStripLabel3.Size = new System.Drawing.Size(80, 22);
            this.toolStripLabel3.Text = "Preconditioner:";
            // 
            // comboPrecon
            // 
            this.comboPrecon.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboPrecon.Enabled = false;
            this.comboPrecon.Name = "comboPrecon";
            this.comboPrecon.Size = new System.Drawing.Size(100, 25);
            // 
            // toolStripSeparator3
            // 
            this.toolStripSeparator3.Name = "toolStripSeparator3";
            this.toolStripSeparator3.Size = new System.Drawing.Size(6, 25);
            // 
            // statusStrip1
            // 
            this.statusStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.lbInfo});
            this.statusStrip1.Location = new System.Drawing.Point(0, 551);
            this.statusStrip1.Name = "statusStrip1";
            this.statusStrip1.Size = new System.Drawing.Size(842, 22);
            this.statusStrip1.TabIndex = 1;
            this.statusStrip1.Text = "statusStrip1";
            // 
            // labelInfo
            // 
            this.labelInfo.Name = "labelInfo";
            this.labelInfo.Size = new System.Drawing.Size(11, 17);
            this.labelInfo.Text = "-";
            // 
            // btnPropSolver
            // 
            this.btnPropSolver.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnPropSolver.Image = global::SolverBenchmark.Properties.Resources.properties;
            this.btnPropSolver.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnPropSolver.Name = "btnPropSolver";
            this.btnPropSolver.Size = new System.Drawing.Size(23, 22);
            this.btnPropSolver.Text = "Solver Properties";
            this.btnPropSolver.Click += new System.EventHandler(this.btnPropSolver_Click);
            // 
            // btnPropPrecon
            // 
            this.btnPropPrecon.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnPropPrecon.Image = global::SolverBenchmark.Properties.Resources.properties;
            this.btnPropPrecon.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnPropPrecon.Name = "btnPropPrecon";
            this.btnPropPrecon.Size = new System.Drawing.Size(23, 22);
            this.btnPropPrecon.Text = "Preconditioner Properties";
            this.btnPropPrecon.Click += new System.EventHandler(this.btnPropPrecon_Click);
            // 
            // btnLoad
            // 
            this.btnLoad.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnLoad.Image = global::SolverBenchmark.Properties.Resources.open;
            this.btnLoad.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnLoad.Name = "btnLoad";
            this.btnLoad.Size = new System.Drawing.Size(23, 22);
            this.btnLoad.Text = "Load Matrix File";
            this.btnLoad.Click += new System.EventHandler(this.btnLoad_Click);
            // 
            // btnRun
            // 
            this.btnRun.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnRun.Image = global::SolverBenchmark.Properties.Resources.play;
            this.btnRun.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnRun.Name = "btnRun";
            this.btnRun.Size = new System.Drawing.Size(23, 22);
            this.btnRun.Text = "Run Benchmark";
            this.btnRun.Click += new System.EventHandler(this.btnRun_Click);
            // 
            // btnReset
            // 
            this.btnReset.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnReset.Image = global::SolverBenchmark.Properties.Resources.refresh;
            this.btnReset.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnReset.Name = "btnReset";
            this.btnReset.Size = new System.Drawing.Size(23, 22);
            this.btnReset.Text = "Reset to defaults";
            this.btnReset.Click += new System.EventHandler(this.btnReset_Click);
            // 
            // toolStripSeparator4
            // 
            this.toolStripSeparator4.Name = "toolStripSeparator4";
            this.toolStripSeparator4.Size = new System.Drawing.Size(6, 25);
            // 
            // toolStripLabel4
            // 
            this.toolStripLabel4.Name = "toolStripLabel4";
            this.toolStripLabel4.Size = new System.Drawing.Size(58, 22);
            this.toolStripLabel4.Text = "Tolerance:";
            // 
            // comboTolerance
            // 
            this.comboTolerance.AutoSize = false;
            this.comboTolerance.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboTolerance.DropDownWidth = 60;
            this.comboTolerance.Items.AddRange(new object[] {
            "1e-4",
            "1e-6",
            "1e-8",
            "1e-10",
            "1e-12"});
            this.comboTolerance.Name = "comboTolerance";
            this.comboTolerance.Size = new System.Drawing.Size(60, 25);
            // 
            // benchmarkResultsView1
            // 
            this.benchmarkResultsView1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.benchmarkResultsView1.Location = new System.Drawing.Point(0, 25);
            this.benchmarkResultsView1.Name = "benchmarkResultsView1";
            this.benchmarkResultsView1.Size = new System.Drawing.Size(842, 526);
            this.benchmarkResultsView1.TabIndex = 2;
            // 
            // toolStripLabel5
            // 
            this.toolStripLabel5.Name = "toolStripLabel5";
            this.toolStripLabel5.Size = new System.Drawing.Size(32, 22);
            this.toolStripLabel5.Text = "Limit:";
            // 
            // tbLimit
            // 
            this.tbLimit.Name = "tbLimit";
            this.tbLimit.Size = new System.Drawing.Size(40, 25);
            this.tbLimit.Text = "200";
            this.tbLimit.TextBoxTextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // toolStripSeparator5
            // 
            this.toolStripSeparator5.Name = "toolStripSeparator5";
            this.toolStripSeparator5.Size = new System.Drawing.Size(6, 25);
            // 
            // lbInfo
            // 
            this.lbInfo.Name = "lbInfo";
            this.lbInfo.Size = new System.Drawing.Size(11, 17);
            this.lbInfo.Text = "-";
            // 
            // FormMain
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(842, 573);
            this.Controls.Add(this.benchmarkResultsView1);
            this.Controls.Add(this.statusStrip1);
            this.Controls.Add(this.toolStrip1);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.Name = "FormMain";
            this.Text = "MathNet.Numerics Solver Benchmark";
            this.Load += new System.EventHandler(this.FormMain_Load);
            this.toolStrip1.ResumeLayout(false);
            this.toolStrip1.PerformLayout();
            this.statusStrip1.ResumeLayout(false);
            this.statusStrip1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ToolStrip toolStrip1;
        private System.Windows.Forms.ToolStripLabel toolStripLabel1;
        private System.Windows.Forms.ToolStripComboBox comboType;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator1;
        private System.Windows.Forms.ToolStripLabel toolStripLabel2;
        private System.Windows.Forms.ToolStripComboBox comboSolver;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator2;
        private System.Windows.Forms.ToolStripLabel toolStripLabel3;
        private System.Windows.Forms.ToolStripComboBox comboPrecon;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator3;
        private System.Windows.Forms.ToolStripButton btnLoad;
        private System.Windows.Forms.ToolStripButton btnRun;
        private System.Windows.Forms.StatusStrip statusStrip1;
        private System.Windows.Forms.ToolStripStatusLabel labelInfo;
        private BenchmarkView benchmarkResultsView1;
        private System.Windows.Forms.ToolStripButton btnPropPrecon;
        private System.Windows.Forms.ToolStripButton btnPropSolver;
        private System.Windows.Forms.ToolStripButton btnReset;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator4;
        private System.Windows.Forms.ToolStripLabel toolStripLabel4;
        private System.Windows.Forms.ToolStripComboBox comboTolerance;
        private System.Windows.Forms.ToolStripLabel toolStripLabel5;
        private System.Windows.Forms.ToolStripTextBox tbLimit;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator5;
        private System.Windows.Forms.ToolStripStatusLabel lbInfo;
    }
}