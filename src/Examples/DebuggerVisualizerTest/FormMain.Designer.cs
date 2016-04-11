namespace DebuggerVisualizerTest
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
            this.btnSparseFile = new System.Windows.Forms.Button();
            this.comboFile = new System.Windows.Forms.ComboBox();
            this.btnDense = new System.Windows.Forms.Button();
            this.comboDenseSize = new System.Windows.Forms.ComboBox();
            this.label1 = new System.Windows.Forms.Label();
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.groupBox2 = new System.Windows.Forms.GroupBox();
            this.btnWathen = new System.Windows.Forms.Button();
            this.btnSparseRandom = new System.Windows.Forms.Button();
            this.label7 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.comboWathen = new System.Windows.Forms.ComboBox();
            this.comboSparseSize = new System.Windows.Forms.ComboBox();
            this.label3 = new System.Windows.Forms.Label();
            this.btnFileOpen = new System.Windows.Forms.Button();
            this.label8 = new System.Windows.Forms.Label();
            this.radioDouble = new System.Windows.Forms.RadioButton();
            this.radioComplex = new System.Windows.Forms.RadioButton();
            this.groupBox1.SuspendLayout();
            this.groupBox2.SuspendLayout();
            this.SuspendLayout();
            // 
            // btnSparseFile
            // 
            this.btnSparseFile.Location = new System.Drawing.Point(367, 15);
            this.btnSparseFile.Name = "btnSparseFile";
            this.btnSparseFile.Size = new System.Drawing.Size(75, 23);
            this.btnSparseFile.TabIndex = 0;
            this.btnSparseFile.Text = "Show";
            this.btnSparseFile.UseVisualStyleBackColor = true;
            this.btnSparseFile.Click += new System.EventHandler(this.btnSparseFile_Click);
            // 
            // comboFile
            // 
            this.comboFile.FormattingEnabled = true;
            this.comboFile.Location = new System.Drawing.Point(119, 17);
            this.comboFile.Name = "comboFile";
            this.comboFile.Size = new System.Drawing.Size(208, 21);
            this.comboFile.TabIndex = 1;
            // 
            // btnDense
            // 
            this.btnDense.Location = new System.Drawing.Point(367, 21);
            this.btnDense.Name = "btnDense";
            this.btnDense.Size = new System.Drawing.Size(75, 23);
            this.btnDense.TabIndex = 0;
            this.btnDense.Text = "Show";
            this.btnDense.UseVisualStyleBackColor = true;
            this.btnDense.Click += new System.EventHandler(this.btnDense_Click);
            // 
            // comboDenseSize
            // 
            this.comboDenseSize.FormattingEnabled = true;
            this.comboDenseSize.Items.AddRange(new object[] {
            "5 x 5",
            "10 x 10",
            "10 x 20",
            "20 x 10",
            "50 x 50"});
            this.comboDenseSize.Location = new System.Drawing.Point(119, 23);
            this.comboDenseSize.Name = "comboDenseSize";
            this.comboDenseSize.Size = new System.Drawing.Size(242, 21);
            this.comboDenseSize.TabIndex = 1;
            // 
            // label1
            // 
            this.label1.Location = new System.Drawing.Point(12, 26);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(101, 13);
            this.label1.TabIndex = 2;
            this.label1.Text = "Random:";
            this.label1.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.btnDense);
            this.groupBox1.Controls.Add(this.label1);
            this.groupBox1.Controls.Add(this.comboDenseSize);
            this.groupBox1.Location = new System.Drawing.Point(12, 171);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Size = new System.Drawing.Size(448, 62);
            this.groupBox1.TabIndex = 3;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "Dense Matrix";
            // 
            // groupBox2
            // 
            this.groupBox2.Controls.Add(this.comboFile);
            this.groupBox2.Controls.Add(this.btnWathen);
            this.groupBox2.Controls.Add(this.btnSparseRandom);
            this.groupBox2.Controls.Add(this.label7);
            this.groupBox2.Controls.Add(this.label5);
            this.groupBox2.Controls.Add(this.comboWathen);
            this.groupBox2.Controls.Add(this.comboSparseSize);
            this.groupBox2.Controls.Add(this.label3);
            this.groupBox2.Controls.Add(this.btnFileOpen);
            this.groupBox2.Controls.Add(this.btnSparseFile);
            this.groupBox2.Location = new System.Drawing.Point(12, 46);
            this.groupBox2.Name = "groupBox2";
            this.groupBox2.Size = new System.Drawing.Size(448, 119);
            this.groupBox2.TabIndex = 4;
            this.groupBox2.TabStop = false;
            this.groupBox2.Text = "Sparse Matrix";
            // 
            // btnWathen
            // 
            this.btnWathen.Location = new System.Drawing.Point(367, 83);
            this.btnWathen.Name = "btnWathen";
            this.btnWathen.Size = new System.Drawing.Size(75, 23);
            this.btnWathen.TabIndex = 0;
            this.btnWathen.Text = "Show";
            this.btnWathen.UseVisualStyleBackColor = true;
            this.btnWathen.Click += new System.EventHandler(this.btnSpecial_Click);
            // 
            // btnSparseRandom
            // 
            this.btnSparseRandom.Location = new System.Drawing.Point(367, 49);
            this.btnSparseRandom.Name = "btnSparseRandom";
            this.btnSparseRandom.Size = new System.Drawing.Size(75, 23);
            this.btnSparseRandom.TabIndex = 0;
            this.btnSparseRandom.Text = "Show";
            this.btnSparseRandom.UseVisualStyleBackColor = true;
            this.btnSparseRandom.Click += new System.EventHandler(this.btnSparseRandom_Click);
            // 
            // label7
            // 
            this.label7.Location = new System.Drawing.Point(9, 88);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(104, 13);
            this.label7.TabIndex = 2;
            this.label7.Text = "Special:";
            this.label7.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label5
            // 
            this.label5.Location = new System.Drawing.Point(6, 54);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(107, 13);
            this.label5.TabIndex = 2;
            this.label5.Text = "Random:";
            this.label5.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // comboWathen
            // 
            this.comboWathen.FormattingEnabled = true;
            this.comboWathen.Items.AddRange(new object[] {
            "Laplacian: nx = 20, ny = 20",
            "Laplacian: nx = 20, ny = 40",
            "Laplacian: nx = 40, ny = 20",
            "Laplacian: nx = 50, ny = 50",
            "Wathen: nx = 5, ny = 5",
            "Wathen: nx = 10, ny = 5",
            "Wathen: nx = 10, ny = 10",
            "Wathen: nx = 20, ny = 20",
            "Wathen: nx = 20, ny = 40",
            "Wathen: nx = 50, ny = 50"});
            this.comboWathen.Location = new System.Drawing.Point(119, 85);
            this.comboWathen.Name = "comboWathen";
            this.comboWathen.Size = new System.Drawing.Size(242, 21);
            this.comboWathen.TabIndex = 1;
            // 
            // comboSparseSize
            // 
            this.comboSparseSize.FormattingEnabled = true;
            this.comboSparseSize.Items.AddRange(new object[] {
            "500 x 500",
            "1000 x 1000",
            "1000 x 2000",
            "2000 x 1000",
            "5000 x 5000",
            "2500 x 2500 (Symmetric)",
            "5000 x 5000 (Symmetric)"});
            this.comboSparseSize.Location = new System.Drawing.Point(119, 51);
            this.comboSparseSize.Name = "comboSparseSize";
            this.comboSparseSize.Size = new System.Drawing.Size(242, 21);
            this.comboSparseSize.TabIndex = 1;
            // 
            // label3
            // 
            this.label3.Location = new System.Drawing.Point(9, 20);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(104, 13);
            this.label3.TabIndex = 2;
            this.label3.Text = "File:";
            this.label3.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // btnFileOpen
            // 
            this.btnFileOpen.Location = new System.Drawing.Point(333, 15);
            this.btnFileOpen.Name = "btnFileOpen";
            this.btnFileOpen.Size = new System.Drawing.Size(28, 23);
            this.btnFileOpen.TabIndex = 0;
            this.btnFileOpen.Text = "...";
            this.btnFileOpen.UseVisualStyleBackColor = true;
            this.btnFileOpen.Click += new System.EventHandler(this.btnFileOpen_Click);
            // 
            // label8
            // 
            this.label8.Location = new System.Drawing.Point(15, 9);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(110, 13);
            this.label8.TabIndex = 2;
            this.label8.Text = "Select Type:";
            this.label8.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // radioDouble
            // 
            this.radioDouble.AutoSize = true;
            this.radioDouble.Checked = true;
            this.radioDouble.Location = new System.Drawing.Point(131, 7);
            this.radioDouble.Name = "radioDouble";
            this.radioDouble.Size = new System.Drawing.Size(63, 17);
            this.radioDouble.TabIndex = 5;
            this.radioDouble.TabStop = true;
            this.radioDouble.Text = "Double";
            this.radioDouble.UseVisualStyleBackColor = true;
            // 
            // radioComplex
            // 
            this.radioComplex.AutoSize = true;
            this.radioComplex.Location = new System.Drawing.Point(200, 7);
            this.radioComplex.Name = "radioComplex";
            this.radioComplex.Size = new System.Drawing.Size(69, 17);
            this.radioComplex.TabIndex = 5;
            this.radioComplex.Text = "Complex";
            this.radioComplex.UseVisualStyleBackColor = true;
            // 
            // FormMain
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(472, 249);
            this.Controls.Add(this.radioComplex);
            this.Controls.Add(this.radioDouble);
            this.Controls.Add(this.groupBox2);
            this.Controls.Add(this.groupBox1);
            this.Controls.Add(this.label8);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Name = "FormMain";
            this.Text = "Math.NET Matrix DebuggerVisualizer Test";
            this.Load += new System.EventHandler(this.FormMain_Load);
            this.groupBox1.ResumeLayout(false);
            this.groupBox2.ResumeLayout(false);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button btnSparseFile;
        private System.Windows.Forms.ComboBox comboFile;
        private System.Windows.Forms.Button btnDense;
        private System.Windows.Forms.ComboBox comboDenseSize;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.GroupBox groupBox2;
        private System.Windows.Forms.Button btnWathen;
        private System.Windows.Forms.Button btnSparseRandom;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.ComboBox comboWathen;
        private System.Windows.Forms.ComboBox comboSparseSize;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Button btnFileOpen;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.RadioButton radioDouble;
        private System.Windows.Forms.RadioButton radioComplex;
    }
}

