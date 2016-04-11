namespace SolverBenchmark
{
    partial class FormMatrix
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
            this.label1 = new System.Windows.Forms.Label();
            this.lbType = new System.Windows.Forms.Label();
            this.rbLaplace = new System.Windows.Forms.RadioButton();
            this.rbRandom = new System.Windows.Forms.RadioButton();
            this.btnCancel = new System.Windows.Forms.Button();
            this.btnOK = new System.Windows.Forms.Button();
            this.cbLaplace = new System.Windows.Forms.ComboBox();
            this.cbRandom = new System.Windows.Forms.ComboBox();
            this.rbFile = new System.Windows.Forms.RadioButton();
            this.btnFile = new System.Windows.Forms.Button();
            this.tbFile = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(12, 9);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(74, 13);
            this.label1.TabIndex = 0;
            this.label1.Text = "Current type:";
            // 
            // lbType
            // 
            this.lbType.AutoSize = true;
            this.lbType.Location = new System.Drawing.Point(92, 9);
            this.lbType.Name = "lbType";
            this.lbType.Size = new System.Drawing.Size(45, 13);
            this.lbType.TabIndex = 1;
            this.lbType.Text = "Double";
            // 
            // rbLaplace
            // 
            this.rbLaplace.AutoSize = true;
            this.rbLaplace.Checked = true;
            this.rbLaplace.Location = new System.Drawing.Point(15, 35);
            this.rbLaplace.Name = "rbLaplace";
            this.rbLaplace.Size = new System.Drawing.Size(73, 17);
            this.rbLaplace.TabIndex = 2;
            this.rbLaplace.TabStop = true;
            this.rbLaplace.Text = "Laplacian";
            this.rbLaplace.UseVisualStyleBackColor = true;
            // 
            // rbRandom
            // 
            this.rbRandom.AutoSize = true;
            this.rbRandom.Location = new System.Drawing.Point(15, 62);
            this.rbRandom.Name = "rbRandom";
            this.rbRandom.Size = new System.Drawing.Size(68, 17);
            this.rbRandom.TabIndex = 2;
            this.rbRandom.Text = "Random";
            this.rbRandom.UseVisualStyleBackColor = true;
            // 
            // btnCancel
            // 
            this.btnCancel.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.btnCancel.DialogResult = System.Windows.Forms.DialogResult.Cancel;
            this.btnCancel.Location = new System.Drawing.Point(257, 140);
            this.btnCancel.Name = "btnCancel";
            this.btnCancel.Size = new System.Drawing.Size(75, 23);
            this.btnCancel.TabIndex = 3;
            this.btnCancel.Text = "Cancel";
            this.btnCancel.UseVisualStyleBackColor = true;
            // 
            // btnOK
            // 
            this.btnOK.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.btnOK.DialogResult = System.Windows.Forms.DialogResult.OK;
            this.btnOK.Location = new System.Drawing.Point(176, 140);
            this.btnOK.Name = "btnOK";
            this.btnOK.Size = new System.Drawing.Size(75, 23);
            this.btnOK.TabIndex = 3;
            this.btnOK.Text = "OK";
            this.btnOK.UseVisualStyleBackColor = true;
            // 
            // cbLaplace
            // 
            this.cbLaplace.FormattingEnabled = true;
            this.cbLaplace.Items.AddRange(new object[] {
            "nx = 20, ny = 20",
            "nx = 20, ny = 40",
            "nx = 40, ny = 20",
            "nx = 50, ny = 50"});
            this.cbLaplace.Location = new System.Drawing.Point(95, 34);
            this.cbLaplace.Name = "cbLaplace";
            this.cbLaplace.Size = new System.Drawing.Size(237, 21);
            this.cbLaplace.TabIndex = 4;
            // 
            // cbRandom
            // 
            this.cbRandom.FormattingEnabled = true;
            this.cbRandom.Items.AddRange(new object[] {
            "1000 x 1000 (Symmetric)",
            "2500 x 2500 (Symmetric)",
            "5000 x 5000 (Symmetric)",
            "1000 x 1000",
            "2500 x 2500",
            "5000 x 5000"});
            this.cbRandom.Location = new System.Drawing.Point(95, 61);
            this.cbRandom.Name = "cbRandom";
            this.cbRandom.Size = new System.Drawing.Size(237, 21);
            this.cbRandom.TabIndex = 4;
            // 
            // rbFile
            // 
            this.rbFile.AutoSize = true;
            this.rbFile.Location = new System.Drawing.Point(15, 89);
            this.rbFile.Name = "rbFile";
            this.rbFile.Size = new System.Drawing.Size(43, 17);
            this.rbFile.TabIndex = 2;
            this.rbFile.Text = "File";
            this.rbFile.UseVisualStyleBackColor = true;
            // 
            // btnFile
            // 
            this.btnFile.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.btnFile.DialogResult = System.Windows.Forms.DialogResult.Cancel;
            this.btnFile.Location = new System.Drawing.Point(257, 88);
            this.btnFile.Name = "btnFile";
            this.btnFile.Size = new System.Drawing.Size(75, 23);
            this.btnFile.TabIndex = 3;
            this.btnFile.Text = "Choose";
            this.btnFile.UseVisualStyleBackColor = true;
            this.btnFile.Click += new System.EventHandler(this.btnFile_Click);
            // 
            // tbFile
            // 
            this.tbFile.Location = new System.Drawing.Point(95, 88);
            this.tbFile.Name = "tbFile";
            this.tbFile.Size = new System.Drawing.Size(154, 22);
            this.tbFile.TabIndex = 5;
            // 
            // FormMatrix
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(344, 175);
            this.ControlBox = false;
            this.Controls.Add(this.tbFile);
            this.Controls.Add(this.cbRandom);
            this.Controls.Add(this.cbLaplace);
            this.Controls.Add(this.btnOK);
            this.Controls.Add(this.btnFile);
            this.Controls.Add(this.btnCancel);
            this.Controls.Add(this.rbFile);
            this.Controls.Add(this.rbRandom);
            this.Controls.Add(this.rbLaplace);
            this.Controls.Add(this.lbType);
            this.Controls.Add(this.label1);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Name = "FormMatrix";
            this.Text = "Choose Matrix";
            this.Load += new System.EventHandler(this.FormMatrix_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label lbType;
        private System.Windows.Forms.RadioButton rbLaplace;
        private System.Windows.Forms.RadioButton rbRandom;
        private System.Windows.Forms.Button btnCancel;
        private System.Windows.Forms.Button btnOK;
        private System.Windows.Forms.ComboBox cbLaplace;
        private System.Windows.Forms.ComboBox cbRandom;
        private System.Windows.Forms.RadioButton rbFile;
        private System.Windows.Forms.Button btnFile;
        private System.Windows.Forms.TextBox tbFile;
    }
}