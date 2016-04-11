namespace MathNet.MatrixDebuggerVisualizer.UI.Views
{
    partial class SparseInfoView
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
            this.components = new System.ComponentModel.Container();
            this.toolStrip = new System.Windows.Forms.ToolStrip();
            this.btnSave = new System.Windows.Forms.ToolStripButton();
            this.btnRefresh = new System.Windows.Forms.ToolStripButton();
            this.label1 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.lbObject = new System.Windows.Forms.Label();
            this.lbStorageType = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.lbSize = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.lbValueCount = new System.Windows.Forms.Label();
            this.panel1 = new System.Windows.Forms.Panel();
            this.lbInfo3 = new System.Windows.Forms.Label();
            this.lbInfo2 = new System.Windows.Forms.Label();
            this.lbInfo1 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label20 = new System.Windows.Forms.Label();
            this.lbNzRowDev = new System.Windows.Forms.Label();
            this.lbNzDiag = new System.Windows.Forms.Label();
            this.label16 = new System.Windows.Forms.Label();
            this.label22 = new System.Windows.Forms.Label();
            this.label24 = new System.Windows.Forms.Label();
            this.label14 = new System.Windows.Forms.Label();
            this.lbNzUpper = new System.Windows.Forms.Label();
            this.lbColumnWeight = new System.Windows.Forms.Label();
            this.lbRowWeight = new System.Windows.Forms.Label();
            this.lbNzColumn = new System.Windows.Forms.Label();
            this.lbNzRow = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.label37 = new System.Windows.Forms.Label();
            this.lbDistDiag = new System.Windows.Forms.Label();
            this.label35 = new System.Windows.Forms.Label();
            this.lbNumDiag = new System.Windows.Forms.Label();
            this.panel2 = new System.Windows.Forms.Panel();
            this.lbUpdate1 = new System.Windows.Forms.Label();
            this.label23 = new System.Windows.Forms.Label();
            this.label27 = new System.Windows.Forms.Label();
            this.label13 = new System.Windows.Forms.Label();
            this.label25 = new System.Windows.Forms.Label();
            this.label29 = new System.Windows.Forms.Label();
            this.lbBand80 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.lbBand90 = new System.Windows.Forms.Label();
            this.lbBandAvg = new System.Windows.Forms.Label();
            this.lbBandUpper = new System.Windows.Forms.Label();
            this.panel3 = new System.Windows.Forms.Panel();
            this.lbUpdate3 = new System.Windows.Forms.Label();
            this.label32 = new System.Windows.Forms.Label();
            this.lbNumSymPcnt = new System.Windows.Forms.Label();
            this.lbNumSym = new System.Windows.Forms.Label();
            this.label42 = new System.Windows.Forms.Label();
            this.lbFrobNonSym = new System.Windows.Forms.Label();
            this.label39 = new System.Windows.Forms.Label();
            this.lbFrobSym = new System.Windows.Forms.Label();
            this.panel4 = new System.Windows.Forms.Panel();
            this.panel5 = new System.Windows.Forms.Panel();
            this.lbUpdate2 = new System.Windows.Forms.Label();
            this.label21 = new System.Windows.Forms.Label();
            this.lbDomRowsPcnt = new System.Windows.Forms.Label();
            this.lbDomColumnsPcnt = new System.Windows.Forms.Label();
            this.label17 = new System.Windows.Forms.Label();
            this.label34 = new System.Windows.Forms.Label();
            this.lbDomRows = new System.Windows.Forms.Label();
            this.lbNormInf = new System.Windows.Forms.Label();
            this.lbNorm1 = new System.Windows.Forms.Label();
            this.label31 = new System.Windows.Forms.Label();
            this.label19 = new System.Windows.Forms.Label();
            this.lbNormMax = new System.Windows.Forms.Label();
            this.lbDomColumns = new System.Windows.Forms.Label();
            this.label26 = new System.Windows.Forms.Label();
            this.lbNormFrob = new System.Windows.Forms.Label();
            this.toolTip = new System.Windows.Forms.ToolTip(this.components);
            this.label11 = new System.Windows.Forms.Label();
            this.lbBytes = new System.Windows.Forms.Label();
            this.label15 = new System.Windows.Forms.Label();
            this.lbDataType = new System.Windows.Forms.Label();
            this.toolStrip.SuspendLayout();
            this.panel1.SuspendLayout();
            this.panel2.SuspendLayout();
            this.panel3.SuspendLayout();
            this.panel4.SuspendLayout();
            this.panel5.SuspendLayout();
            this.SuspendLayout();
            // 
            // toolStrip
            // 
            this.toolStrip.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.btnSave,
            this.btnRefresh});
            this.toolStrip.Location = new System.Drawing.Point(0, 0);
            this.toolStrip.Name = "toolStrip";
            this.toolStrip.Size = new System.Drawing.Size(600, 25);
            this.toolStrip.TabIndex = 0;
            // 
            // btnSave
            // 
            this.btnSave.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnSave.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.save;
            this.btnSave.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnSave.Name = "btnSave";
            this.btnSave.Size = new System.Drawing.Size(23, 22);
            this.btnSave.Text = "Save";
            this.btnSave.Click += new System.EventHandler(this.btnSave_Click);
            // 
            // btnRefresh
            // 
            this.btnRefresh.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.btnRefresh.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.refresh;
            this.btnRefresh.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.btnRefresh.Name = "btnRefresh";
            this.btnRefresh.Size = new System.Drawing.Size(23, 22);
            this.btnRefresh.Click += new System.EventHandler(this.btnRefresh_Click);
            // 
            // label1
            // 
            this.label1.BackColor = System.Drawing.Color.WhiteSmoke;
            this.label1.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label1.ForeColor = System.Drawing.Color.DimGray;
            this.label1.Location = new System.Drawing.Point(5, 39);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(69, 13);
            this.label1.TabIndex = 2;
            this.label1.Text = "Object:";
            this.label1.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // label3
            // 
            this.label3.BackColor = System.Drawing.Color.WhiteSmoke;
            this.label3.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label3.ForeColor = System.Drawing.Color.DimGray;
            this.label3.Location = new System.Drawing.Point(5, 57);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(69, 13);
            this.label3.TabIndex = 2;
            this.label3.Text = "Storage:";
            this.label3.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // lbObject
            // 
            this.lbObject.ForeColor = System.Drawing.Color.Black;
            this.lbObject.Location = new System.Drawing.Point(80, 39);
            this.lbObject.Name = "lbObject";
            this.lbObject.Size = new System.Drawing.Size(285, 13);
            this.lbObject.TabIndex = 2;
            this.lbObject.Text = "-";
            // 
            // lbStorageType
            // 
            this.lbStorageType.AutoEllipsis = true;
            this.lbStorageType.ForeColor = System.Drawing.Color.Black;
            this.lbStorageType.Location = new System.Drawing.Point(80, 57);
            this.lbStorageType.Name = "lbStorageType";
            this.lbStorageType.Size = new System.Drawing.Size(285, 13);
            this.lbStorageType.TabIndex = 2;
            this.lbStorageType.Text = "-";
            // 
            // label4
            // 
            this.label4.BackColor = System.Drawing.Color.WhiteSmoke;
            this.label4.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label4.ForeColor = System.Drawing.Color.DimGray;
            this.label4.Location = new System.Drawing.Point(377, 39);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(82, 13);
            this.label4.TabIndex = 2;
            this.label4.Text = "Size:";
            this.label4.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // lbSize
            // 
            this.lbSize.ForeColor = System.Drawing.Color.Black;
            this.lbSize.Location = new System.Drawing.Point(465, 39);
            this.lbSize.Name = "lbSize";
            this.lbSize.Size = new System.Drawing.Size(128, 13);
            this.lbSize.TabIndex = 2;
            this.lbSize.Text = "-";
            // 
            // label6
            // 
            this.label6.BackColor = System.Drawing.Color.WhiteSmoke;
            this.label6.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label6.ForeColor = System.Drawing.Color.DimGray;
            this.label6.Location = new System.Drawing.Point(371, 57);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(88, 13);
            this.label6.TabIndex = 2;
            this.label6.Text = "Value Count:";
            this.label6.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // lbValueCount
            // 
            this.lbValueCount.ForeColor = System.Drawing.Color.Black;
            this.lbValueCount.Location = new System.Drawing.Point(465, 57);
            this.lbValueCount.Name = "lbValueCount";
            this.lbValueCount.Size = new System.Drawing.Size(128, 13);
            this.lbValueCount.TabIndex = 2;
            this.lbValueCount.Text = "-";
            // 
            // panel1
            // 
            this.panel1.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.panel1.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(242)))), ((int)(((byte)(242)))), ((int)(((byte)(242)))));
            this.panel1.Controls.Add(this.lbInfo3);
            this.panel1.Controls.Add(this.lbInfo2);
            this.panel1.Controls.Add(this.lbInfo1);
            this.panel1.Controls.Add(this.label5);
            this.panel1.Controls.Add(this.label20);
            this.panel1.Controls.Add(this.lbNzRowDev);
            this.panel1.Controls.Add(this.lbNzDiag);
            this.panel1.Controls.Add(this.label16);
            this.panel1.Controls.Add(this.label22);
            this.panel1.Controls.Add(this.label24);
            this.panel1.Controls.Add(this.label14);
            this.panel1.Controls.Add(this.lbNzUpper);
            this.panel1.Controls.Add(this.lbColumnWeight);
            this.panel1.Controls.Add(this.lbRowWeight);
            this.panel1.Controls.Add(this.lbNzColumn);
            this.panel1.Controls.Add(this.lbNzRow);
            this.panel1.Controls.Add(this.label8);
            this.panel1.Controls.Add(this.label2);
            this.panel1.Controls.Add(this.label9);
            this.panel1.Location = new System.Drawing.Point(0, 0);
            this.panel1.Name = "panel1";
            this.panel1.Size = new System.Drawing.Size(597, 76);
            this.panel1.TabIndex = 3;
            // 
            // lbInfo3
            // 
            this.lbInfo3.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.lbInfo3.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.warning;
            this.lbInfo3.Location = new System.Drawing.Point(577, 44);
            this.lbInfo3.Name = "lbInfo3";
            this.lbInfo3.Size = new System.Drawing.Size(16, 16);
            this.lbInfo3.TabIndex = 1;
            this.lbInfo3.Visible = false;
            // 
            // lbInfo2
            // 
            this.lbInfo2.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.lbInfo2.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.warning;
            this.lbInfo2.Location = new System.Drawing.Point(577, 28);
            this.lbInfo2.Name = "lbInfo2";
            this.lbInfo2.Size = new System.Drawing.Size(16, 16);
            this.lbInfo2.TabIndex = 1;
            this.lbInfo2.Visible = false;
            // 
            // lbInfo1
            // 
            this.lbInfo1.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.lbInfo1.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.warning;
            this.lbInfo1.Location = new System.Drawing.Point(577, 12);
            this.lbInfo1.Name = "lbInfo1";
            this.lbInfo1.Size = new System.Drawing.Size(16, 16);
            this.lbInfo1.TabIndex = 1;
            this.lbInfo1.Visible = false;
            // 
            // label5
            // 
            this.label5.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.refresh_small;
            this.label5.Location = new System.Drawing.Point(3, 8);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(16, 16);
            this.label5.TabIndex = 1;
            this.label5.Click += new System.EventHandler(this.lbUpdate0_Click);
            // 
            // label20
            // 
            this.label20.ForeColor = System.Drawing.Color.Black;
            this.label20.Location = new System.Drawing.Point(381, 46);
            this.label20.Name = "label20";
            this.label20.Size = new System.Drawing.Size(80, 13);
            this.label20.TabIndex = 0;
            this.label20.Text = "Nonzeros:";
            this.label20.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNzRowDev
            // 
            this.lbNzRowDev.ForeColor = System.Drawing.Color.Black;
            this.lbNzRowDev.Location = new System.Drawing.Point(381, 28);
            this.lbNzRowDev.Name = "lbNzRowDev";
            this.lbNzRowDev.Size = new System.Drawing.Size(80, 13);
            this.lbNzRowDev.TabIndex = 0;
            this.lbNzRowDev.Text = "Nonzeros:";
            this.lbNzRowDev.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNzDiag
            // 
            this.lbNzDiag.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lbNzDiag.ForeColor = System.Drawing.Color.Black;
            this.lbNzDiag.Location = new System.Drawing.Point(467, 11);
            this.lbNzDiag.Name = "lbNzDiag";
            this.lbNzDiag.Size = new System.Drawing.Size(104, 13);
            this.lbNzDiag.TabIndex = 0;
            this.lbNzDiag.Text = "-";
            // 
            // label16
            // 
            this.label16.ForeColor = System.Drawing.Color.Black;
            this.label16.Location = new System.Drawing.Point(381, 11);
            this.label16.Name = "label16";
            this.label16.Size = new System.Drawing.Size(80, 13);
            this.label16.TabIndex = 0;
            this.label16.Text = "Diagonal:";
            this.label16.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label22
            // 
            this.label22.AutoSize = true;
            this.label22.ForeColor = System.Drawing.Color.DimGray;
            this.label22.Location = new System.Drawing.Point(34, 46);
            this.label22.Name = "label22";
            this.label22.Size = new System.Drawing.Size(60, 13);
            this.label22.TabIndex = 0;
            this.label22.Text = "COLUMNS";
            this.label22.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label24
            // 
            this.label24.AutoSize = true;
            this.label24.ForeColor = System.Drawing.Color.DimGray;
            this.label24.Location = new System.Drawing.Point(34, 11);
            this.label24.Name = "label24";
            this.label24.Size = new System.Drawing.Size(66, 13);
            this.label24.TabIndex = 0;
            this.label24.Text = "NONZEROS";
            this.label24.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.ForeColor = System.Drawing.Color.DimGray;
            this.label14.Location = new System.Drawing.Point(34, 28);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(40, 13);
            this.label14.TabIndex = 0;
            this.label14.Text = "ROWS";
            this.label14.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNzUpper
            // 
            this.lbNzUpper.ForeColor = System.Drawing.Color.Black;
            this.lbNzUpper.Location = new System.Drawing.Point(255, 11);
            this.lbNzUpper.Name = "lbNzUpper";
            this.lbNzUpper.Size = new System.Drawing.Size(120, 13);
            this.lbNzUpper.TabIndex = 0;
            this.lbNzUpper.Text = "-";
            this.lbNzUpper.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbColumnWeight
            // 
            this.lbColumnWeight.ForeColor = System.Drawing.Color.Black;
            this.lbColumnWeight.Location = new System.Drawing.Point(255, 46);
            this.lbColumnWeight.Name = "lbColumnWeight";
            this.lbColumnWeight.Size = new System.Drawing.Size(120, 13);
            this.lbColumnWeight.TabIndex = 0;
            this.lbColumnWeight.Text = "-";
            this.lbColumnWeight.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbRowWeight
            // 
            this.lbRowWeight.ForeColor = System.Drawing.Color.Black;
            this.lbRowWeight.Location = new System.Drawing.Point(255, 28);
            this.lbRowWeight.Name = "lbRowWeight";
            this.lbRowWeight.Size = new System.Drawing.Size(120, 13);
            this.lbRowWeight.TabIndex = 0;
            this.lbRowWeight.Text = "-";
            this.lbRowWeight.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNzColumn
            // 
            this.lbNzColumn.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lbNzColumn.ForeColor = System.Drawing.Color.Black;
            this.lbNzColumn.Location = new System.Drawing.Point(467, 46);
            this.lbNzColumn.Name = "lbNzColumn";
            this.lbNzColumn.Size = new System.Drawing.Size(104, 13);
            this.lbNzColumn.TabIndex = 0;
            this.lbNzColumn.Text = "-";
            // 
            // lbNzRow
            // 
            this.lbNzRow.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lbNzRow.ForeColor = System.Drawing.Color.Black;
            this.lbNzRow.Location = new System.Drawing.Point(467, 28);
            this.lbNzRow.Name = "lbNzRow";
            this.lbNzRow.Size = new System.Drawing.Size(104, 13);
            this.lbNzRow.TabIndex = 0;
            this.lbNzRow.Text = "-";
            // 
            // label8
            // 
            this.label8.ForeColor = System.Drawing.Color.Black;
            this.label8.Location = new System.Drawing.Point(129, 11);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(120, 13);
            this.label8.TabIndex = 0;
            this.label8.Text = "Upper/Lower:";
            this.label8.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label2
            // 
            this.label2.ForeColor = System.Drawing.Color.Black;
            this.label2.Location = new System.Drawing.Point(132, 46);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(117, 13);
            this.label2.TabIndex = 0;
            this.label2.Text = "Shortest/Longest:";
            this.label2.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label9
            // 
            this.label9.ForeColor = System.Drawing.Color.Black;
            this.label9.Location = new System.Drawing.Point(129, 28);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(120, 13);
            this.label9.TabIndex = 0;
            this.label9.Text = "Shortest/Longest:";
            this.label9.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label37
            // 
            this.label37.ForeColor = System.Drawing.Color.Black;
            this.label37.Location = new System.Drawing.Point(381, 7);
            this.label37.Name = "label37";
            this.label37.Size = new System.Drawing.Size(80, 13);
            this.label37.TabIndex = 0;
            this.label37.Text = "Distance:";
            this.label37.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbDistDiag
            // 
            this.lbDistDiag.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lbDistDiag.ForeColor = System.Drawing.Color.Black;
            this.lbDistDiag.Location = new System.Drawing.Point(467, 7);
            this.lbDistDiag.Name = "lbDistDiag";
            this.lbDistDiag.Size = new System.Drawing.Size(126, 13);
            this.lbDistDiag.TabIndex = 0;
            this.lbDistDiag.Text = "-";
            // 
            // label35
            // 
            this.label35.ForeColor = System.Drawing.Color.Black;
            this.label35.Location = new System.Drawing.Point(135, 7);
            this.label35.Name = "label35";
            this.label35.Size = new System.Drawing.Size(114, 13);
            this.label35.TabIndex = 0;
            this.label35.Text = "Total:";
            this.label35.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNumDiag
            // 
            this.lbNumDiag.ForeColor = System.Drawing.Color.Black;
            this.lbNumDiag.Location = new System.Drawing.Point(255, 7);
            this.lbNumDiag.Name = "lbNumDiag";
            this.lbNumDiag.Size = new System.Drawing.Size(120, 13);
            this.lbNumDiag.TabIndex = 0;
            this.lbNumDiag.Text = "-";
            this.lbNumDiag.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // panel2
            // 
            this.panel2.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.panel2.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(242)))), ((int)(((byte)(242)))), ((int)(((byte)(242)))));
            this.panel2.Controls.Add(this.lbUpdate1);
            this.panel2.Controls.Add(this.label23);
            this.panel2.Controls.Add(this.label27);
            this.panel2.Controls.Add(this.label13);
            this.panel2.Controls.Add(this.label25);
            this.panel2.Controls.Add(this.label29);
            this.panel2.Controls.Add(this.lbBand80);
            this.panel2.Controls.Add(this.label12);
            this.panel2.Controls.Add(this.lbBand90);
            this.panel2.Controls.Add(this.label37);
            this.panel2.Controls.Add(this.lbBandAvg);
            this.panel2.Controls.Add(this.lbBandUpper);
            this.panel2.Controls.Add(this.lbDistDiag);
            this.panel2.Controls.Add(this.lbNumDiag);
            this.panel2.Controls.Add(this.label35);
            this.panel2.Location = new System.Drawing.Point(0, 82);
            this.panel2.Name = "panel2";
            this.panel2.Size = new System.Drawing.Size(597, 96);
            this.panel2.TabIndex = 3;
            // 
            // lbUpdate1
            // 
            this.lbUpdate1.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.refresh_small;
            this.lbUpdate1.Location = new System.Drawing.Point(3, 6);
            this.lbUpdate1.Name = "lbUpdate1";
            this.lbUpdate1.Size = new System.Drawing.Size(16, 16);
            this.lbUpdate1.TabIndex = 1;
            this.lbUpdate1.Click += new System.EventHandler(this.lbUpdate1_Click);
            // 
            // label23
            // 
            this.label23.ForeColor = System.Drawing.Color.Black;
            this.label23.Location = new System.Drawing.Point(384, 25);
            this.label23.Name = "label23";
            this.label23.Size = new System.Drawing.Size(77, 13);
            this.label23.TabIndex = 0;
            this.label23.Text = "Average:";
            this.label23.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label27
            // 
            this.label27.ForeColor = System.Drawing.Color.Black;
            this.label27.Location = new System.Drawing.Point(3, 66);
            this.label27.Name = "label27";
            this.label27.Size = new System.Drawing.Size(246, 13);
            this.label27.TabIndex = 0;
            this.label27.Text = "80% of matrix is in the band of width:";
            this.label27.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.ForeColor = System.Drawing.Color.DimGray;
            this.label13.Location = new System.Drawing.Point(34, 25);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(72, 13);
            this.label13.TabIndex = 0;
            this.label13.Text = "BANDWIDTH";
            this.label13.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label25
            // 
            this.label25.ForeColor = System.Drawing.Color.Black;
            this.label25.Location = new System.Drawing.Point(3, 48);
            this.label25.Name = "label25";
            this.label25.Size = new System.Drawing.Size(246, 13);
            this.label25.TabIndex = 0;
            this.label25.Text = "90% of matrix is in the band of width:";
            this.label25.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label29
            // 
            this.label29.AutoSize = true;
            this.label29.ForeColor = System.Drawing.Color.DimGray;
            this.label29.Location = new System.Drawing.Point(34, 7);
            this.label29.Name = "label29";
            this.label29.Size = new System.Drawing.Size(68, 13);
            this.label29.TabIndex = 0;
            this.label29.Text = "DIAGONALS";
            this.label29.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbBand80
            // 
            this.lbBand80.ForeColor = System.Drawing.Color.Black;
            this.lbBand80.Location = new System.Drawing.Point(255, 66);
            this.lbBand80.Name = "lbBand80";
            this.lbBand80.Size = new System.Drawing.Size(120, 13);
            this.lbBand80.TabIndex = 0;
            this.lbBand80.Text = "-";
            this.lbBand80.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label12
            // 
            this.label12.ForeColor = System.Drawing.Color.Black;
            this.label12.Location = new System.Drawing.Point(3, 25);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(246, 13);
            this.label12.TabIndex = 0;
            this.label12.Text = "Upper/Lower:";
            this.label12.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbBand90
            // 
            this.lbBand90.ForeColor = System.Drawing.Color.Black;
            this.lbBand90.Location = new System.Drawing.Point(255, 48);
            this.lbBand90.Name = "lbBand90";
            this.lbBand90.Size = new System.Drawing.Size(120, 13);
            this.lbBand90.TabIndex = 0;
            this.lbBand90.Text = "-";
            this.lbBand90.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbBandAvg
            // 
            this.lbBandAvg.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lbBandAvg.ForeColor = System.Drawing.Color.Black;
            this.lbBandAvg.Location = new System.Drawing.Point(467, 25);
            this.lbBandAvg.Name = "lbBandAvg";
            this.lbBandAvg.Size = new System.Drawing.Size(126, 13);
            this.lbBandAvg.TabIndex = 0;
            this.lbBandAvg.Text = "-";
            // 
            // lbBandUpper
            // 
            this.lbBandUpper.ForeColor = System.Drawing.Color.Black;
            this.lbBandUpper.Location = new System.Drawing.Point(255, 25);
            this.lbBandUpper.Name = "lbBandUpper";
            this.lbBandUpper.Size = new System.Drawing.Size(120, 13);
            this.lbBandUpper.TabIndex = 0;
            this.lbBandUpper.Text = "-";
            this.lbBandUpper.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // panel3
            // 
            this.panel3.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.panel3.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(242)))), ((int)(((byte)(242)))), ((int)(((byte)(242)))));
            this.panel3.Controls.Add(this.lbUpdate3);
            this.panel3.Controls.Add(this.label32);
            this.panel3.Controls.Add(this.lbNumSymPcnt);
            this.panel3.Controls.Add(this.lbNumSym);
            this.panel3.Controls.Add(this.label42);
            this.panel3.Controls.Add(this.lbFrobNonSym);
            this.panel3.Controls.Add(this.label39);
            this.panel3.Controls.Add(this.lbFrobSym);
            this.panel3.Location = new System.Drawing.Point(0, 316);
            this.panel3.Name = "panel3";
            this.panel3.Size = new System.Drawing.Size(600, 70);
            this.panel3.TabIndex = 4;
            // 
            // lbUpdate3
            // 
            this.lbUpdate3.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.refresh_small;
            this.lbUpdate3.Location = new System.Drawing.Point(3, 4);
            this.lbUpdate3.Name = "lbUpdate3";
            this.lbUpdate3.Size = new System.Drawing.Size(16, 16);
            this.lbUpdate3.TabIndex = 1;
            this.lbUpdate3.Click += new System.EventHandler(this.lbUpdate3_Click);
            // 
            // label32
            // 
            this.label32.ForeColor = System.Drawing.Color.Black;
            this.label32.Location = new System.Drawing.Point(25, 6);
            this.label32.Name = "label32";
            this.label32.Size = new System.Drawing.Size(224, 13);
            this.label32.TabIndex = 0;
            this.label32.Text = "Matching elements in symmetry:";
            this.label32.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNumSymPcnt
            // 
            this.lbNumSymPcnt.ForeColor = System.Drawing.Color.DimGray;
            this.lbNumSymPcnt.Location = new System.Drawing.Point(384, 6);
            this.lbNumSymPcnt.Name = "lbNumSymPcnt";
            this.lbNumSymPcnt.Size = new System.Drawing.Size(50, 13);
            this.lbNumSymPcnt.TabIndex = 0;
            this.lbNumSymPcnt.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNumSym
            // 
            this.lbNumSym.ForeColor = System.Drawing.Color.Black;
            this.lbNumSym.Location = new System.Drawing.Point(255, 6);
            this.lbNumSym.Name = "lbNumSym";
            this.lbNumSym.Size = new System.Drawing.Size(120, 13);
            this.lbNumSym.TabIndex = 0;
            this.lbNumSym.Text = "-";
            this.lbNumSym.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label42
            // 
            this.label42.ForeColor = System.Drawing.Color.Black;
            this.label42.Location = new System.Drawing.Point(3, 25);
            this.label42.Name = "label42";
            this.label42.Size = new System.Drawing.Size(246, 13);
            this.label42.TabIndex = 0;
            this.label42.Text = "Frobenius norm of symmetric part:";
            this.label42.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbFrobNonSym
            // 
            this.lbFrobNonSym.ForeColor = System.Drawing.Color.Black;
            this.lbFrobNonSym.Location = new System.Drawing.Point(255, 43);
            this.lbFrobNonSym.Name = "lbFrobNonSym";
            this.lbFrobNonSym.Size = new System.Drawing.Size(120, 13);
            this.lbFrobNonSym.TabIndex = 0;
            this.lbFrobNonSym.Text = "-";
            this.lbFrobNonSym.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label39
            // 
            this.label39.ForeColor = System.Drawing.Color.Black;
            this.label39.Location = new System.Drawing.Point(3, 43);
            this.label39.Name = "label39";
            this.label39.Size = new System.Drawing.Size(246, 13);
            this.label39.TabIndex = 0;
            this.label39.Text = "Frobenius norm of non-symmetric part:";
            this.label39.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbFrobSym
            // 
            this.lbFrobSym.ForeColor = System.Drawing.Color.Black;
            this.lbFrobSym.Location = new System.Drawing.Point(255, 25);
            this.lbFrobSym.Name = "lbFrobSym";
            this.lbFrobSym.Size = new System.Drawing.Size(120, 13);
            this.lbFrobSym.TabIndex = 0;
            this.lbFrobSym.Text = "-";
            this.lbFrobSym.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // panel4
            // 
            this.panel4.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.panel4.AutoScroll = true;
            this.panel4.Controls.Add(this.panel1);
            this.panel4.Controls.Add(this.panel2);
            this.panel4.Controls.Add(this.panel5);
            this.panel4.Controls.Add(this.panel3);
            this.panel4.Location = new System.Drawing.Point(0, 103);
            this.panel4.Name = "panel4";
            this.panel4.Size = new System.Drawing.Size(600, 397);
            this.panel4.TabIndex = 5;
            // 
            // panel5
            // 
            this.panel5.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.panel5.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(242)))), ((int)(((byte)(242)))), ((int)(((byte)(242)))));
            this.panel5.Controls.Add(this.lbUpdate2);
            this.panel5.Controls.Add(this.label21);
            this.panel5.Controls.Add(this.lbDomRowsPcnt);
            this.panel5.Controls.Add(this.lbDomColumnsPcnt);
            this.panel5.Controls.Add(this.label17);
            this.panel5.Controls.Add(this.label34);
            this.panel5.Controls.Add(this.lbDomRows);
            this.panel5.Controls.Add(this.lbNormInf);
            this.panel5.Controls.Add(this.lbNorm1);
            this.panel5.Controls.Add(this.label31);
            this.panel5.Controls.Add(this.label19);
            this.panel5.Controls.Add(this.lbNormMax);
            this.panel5.Controls.Add(this.lbDomColumns);
            this.panel5.Controls.Add(this.label26);
            this.panel5.Controls.Add(this.lbNormFrob);
            this.panel5.Location = new System.Drawing.Point(0, 184);
            this.panel5.Name = "panel5";
            this.panel5.Size = new System.Drawing.Size(600, 126);
            this.panel5.TabIndex = 4;
            // 
            // lbUpdate2
            // 
            this.lbUpdate2.Image = global::MathNet.MatrixDebuggerVisualizer.Properties.Resources.refresh_small;
            this.lbUpdate2.Location = new System.Drawing.Point(3, 5);
            this.lbUpdate2.Name = "lbUpdate2";
            this.lbUpdate2.Size = new System.Drawing.Size(16, 16);
            this.lbUpdate2.TabIndex = 1;
            this.lbUpdate2.Click += new System.EventHandler(this.lbUpdate2_Click);
            // 
            // label21
            // 
            this.label21.ForeColor = System.Drawing.Color.Black;
            this.label21.Location = new System.Drawing.Point(25, 61);
            this.label21.Name = "label21";
            this.label21.Size = new System.Drawing.Size(224, 13);
            this.label21.TabIndex = 0;
            this.label21.Text = "Infinity-norm of A:";
            this.label21.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbDomRowsPcnt
            // 
            this.lbDomRowsPcnt.ForeColor = System.Drawing.Color.DimGray;
            this.lbDomRowsPcnt.Location = new System.Drawing.Point(384, 84);
            this.lbDomRowsPcnt.Name = "lbDomRowsPcnt";
            this.lbDomRowsPcnt.Size = new System.Drawing.Size(50, 13);
            this.lbDomRowsPcnt.TabIndex = 0;
            this.lbDomRowsPcnt.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbDomColumnsPcnt
            // 
            this.lbDomColumnsPcnt.ForeColor = System.Drawing.Color.DimGray;
            this.lbDomColumnsPcnt.Location = new System.Drawing.Point(384, 102);
            this.lbDomColumnsPcnt.Name = "lbDomColumnsPcnt";
            this.lbDomColumnsPcnt.Size = new System.Drawing.Size(50, 13);
            this.lbDomColumnsPcnt.TabIndex = 0;
            this.lbDomColumnsPcnt.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label17
            // 
            this.label17.ForeColor = System.Drawing.Color.Black;
            this.label17.Location = new System.Drawing.Point(25, 43);
            this.label17.Name = "label17";
            this.label17.Size = new System.Drawing.Size(224, 13);
            this.label17.TabIndex = 0;
            this.label17.Text = "1-norm of A:";
            this.label17.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label34
            // 
            this.label34.ForeColor = System.Drawing.Color.Black;
            this.label34.Location = new System.Drawing.Point(3, 84);
            this.label34.Name = "label34";
            this.label34.Size = new System.Drawing.Size(246, 13);
            this.label34.TabIndex = 0;
            this.label34.Text = "Weakly diagonally dominant rows:";
            this.label34.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbDomRows
            // 
            this.lbDomRows.ForeColor = System.Drawing.Color.Black;
            this.lbDomRows.Location = new System.Drawing.Point(255, 84);
            this.lbDomRows.Name = "lbDomRows";
            this.lbDomRows.Size = new System.Drawing.Size(120, 13);
            this.lbDomRows.TabIndex = 0;
            this.lbDomRows.Text = "-";
            this.lbDomRows.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNormInf
            // 
            this.lbNormInf.ForeColor = System.Drawing.Color.Black;
            this.lbNormInf.Location = new System.Drawing.Point(255, 61);
            this.lbNormInf.Name = "lbNormInf";
            this.lbNormInf.Size = new System.Drawing.Size(120, 13);
            this.lbNormInf.TabIndex = 0;
            this.lbNormInf.Text = "-";
            this.lbNormInf.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNorm1
            // 
            this.lbNorm1.ForeColor = System.Drawing.Color.Black;
            this.lbNorm1.Location = new System.Drawing.Point(255, 43);
            this.lbNorm1.Name = "lbNorm1";
            this.lbNorm1.Size = new System.Drawing.Size(120, 13);
            this.lbNorm1.TabIndex = 0;
            this.lbNorm1.Text = "-";
            this.lbNorm1.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label31
            // 
            this.label31.ForeColor = System.Drawing.Color.Black;
            this.label31.Location = new System.Drawing.Point(3, 102);
            this.label31.Name = "label31";
            this.label31.Size = new System.Drawing.Size(246, 13);
            this.label31.TabIndex = 0;
            this.label31.Text = "Weakly diagonally dominant columns:";
            this.label31.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label19
            // 
            this.label19.ForeColor = System.Drawing.Color.Black;
            this.label19.Location = new System.Drawing.Point(25, 7);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(224, 13);
            this.label19.TabIndex = 0;
            this.label19.Text = "Maximum absolute value of A:";
            this.label19.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNormMax
            // 
            this.lbNormMax.ForeColor = System.Drawing.Color.Black;
            this.lbNormMax.Location = new System.Drawing.Point(255, 7);
            this.lbNormMax.Name = "lbNormMax";
            this.lbNormMax.Size = new System.Drawing.Size(120, 13);
            this.lbNormMax.TabIndex = 0;
            this.lbNormMax.Text = "-";
            this.lbNormMax.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbDomColumns
            // 
            this.lbDomColumns.ForeColor = System.Drawing.Color.Black;
            this.lbDomColumns.Location = new System.Drawing.Point(255, 102);
            this.lbDomColumns.Name = "lbDomColumns";
            this.lbDomColumns.Size = new System.Drawing.Size(120, 13);
            this.lbDomColumns.TabIndex = 0;
            this.lbDomColumns.Text = "-";
            this.lbDomColumns.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label26
            // 
            this.label26.ForeColor = System.Drawing.Color.Black;
            this.label26.Location = new System.Drawing.Point(28, 25);
            this.label26.Name = "label26";
            this.label26.Size = new System.Drawing.Size(221, 13);
            this.label26.TabIndex = 0;
            this.label26.Text = "Frobenius norm of A:";
            this.label26.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // lbNormFrob
            // 
            this.lbNormFrob.ForeColor = System.Drawing.Color.Black;
            this.lbNormFrob.Location = new System.Drawing.Point(255, 25);
            this.lbNormFrob.Name = "lbNormFrob";
            this.lbNormFrob.Size = new System.Drawing.Size(120, 13);
            this.lbNormFrob.TabIndex = 0;
            this.lbNormFrob.Text = "-";
            this.lbNormFrob.TextAlign = System.Drawing.ContentAlignment.TopRight;
            // 
            // label11
            // 
            this.label11.BackColor = System.Drawing.Color.WhiteSmoke;
            this.label11.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label11.ForeColor = System.Drawing.Color.DimGray;
            this.label11.Location = new System.Drawing.Point(374, 75);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(85, 13);
            this.label11.TabIndex = 2;
            this.label11.Text = "Total Bytes:";
            this.label11.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // lbBytes
            // 
            this.lbBytes.ForeColor = System.Drawing.Color.Black;
            this.lbBytes.Location = new System.Drawing.Point(465, 75);
            this.lbBytes.Name = "lbBytes";
            this.lbBytes.Size = new System.Drawing.Size(128, 13);
            this.lbBytes.TabIndex = 2;
            this.lbBytes.Text = "-";
            // 
            // label15
            // 
            this.label15.BackColor = System.Drawing.Color.WhiteSmoke;
            this.label15.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label15.ForeColor = System.Drawing.Color.DimGray;
            this.label15.Location = new System.Drawing.Point(5, 75);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(69, 13);
            this.label15.TabIndex = 2;
            this.label15.Text = "Data Type:";
            this.label15.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // lbDataType
            // 
            this.lbDataType.AutoEllipsis = true;
            this.lbDataType.ForeColor = System.Drawing.Color.Black;
            this.lbDataType.Location = new System.Drawing.Point(80, 75);
            this.lbDataType.Name = "lbDataType";
            this.lbDataType.Size = new System.Drawing.Size(285, 13);
            this.lbDataType.TabIndex = 2;
            this.lbDataType.Text = "-";
            // 
            // SparseInfoView
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.BackColor = System.Drawing.Color.WhiteSmoke;
            this.Controls.Add(this.panel4);
            this.Controls.Add(this.lbDataType);
            this.Controls.Add(this.label15);
            this.Controls.Add(this.lbStorageType);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.lbBytes);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.lbValueCount);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.lbSize);
            this.Controls.Add(this.lbObject);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.toolStrip);
            this.Font = new System.Drawing.Font("Segoe UI", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.ForeColor = System.Drawing.Color.White;
            this.Name = "SparseInfoView";
            this.Size = new System.Drawing.Size(600, 500);
            this.toolStrip.ResumeLayout(false);
            this.toolStrip.PerformLayout();
            this.panel1.ResumeLayout(false);
            this.panel1.PerformLayout();
            this.panel2.ResumeLayout(false);
            this.panel2.PerformLayout();
            this.panel3.ResumeLayout(false);
            this.panel4.ResumeLayout(false);
            this.panel5.ResumeLayout(false);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ToolStrip toolStrip;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label lbObject;
        private System.Windows.Forms.Label lbStorageType;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label lbSize;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label lbValueCount;
        private System.Windows.Forms.ToolStripButton btnSave;
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.Panel panel2;
        private System.Windows.Forms.Label lbNzRowDev;
        private System.Windows.Forms.Label lbNzUpper;
        private System.Windows.Forms.Label lbNzRow;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label lbNzDiag;
        private System.Windows.Forms.Label label16;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.Label lbRowWeight;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.Label label37;
        private System.Windows.Forms.Label lbDistDiag;
        private System.Windows.Forms.Label label35;
        private System.Windows.Forms.Label lbNumDiag;
        private System.Windows.Forms.Label label23;
        private System.Windows.Forms.Label label27;
        private System.Windows.Forms.Label label25;
        private System.Windows.Forms.Label lbBand80;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.Label lbBand90;
        private System.Windows.Forms.Label lbBandAvg;
        private System.Windows.Forms.Label lbBandUpper;
        private System.Windows.Forms.Panel panel3;
        private System.Windows.Forms.Label label32;
        private System.Windows.Forms.Label lbNumSym;
        private System.Windows.Forms.Label label42;
        private System.Windows.Forms.Label lbFrobNonSym;
        private System.Windows.Forms.Label label39;
        private System.Windows.Forms.Label lbFrobSym;
        private System.Windows.Forms.Label lbUpdate1;
        private System.Windows.Forms.Label lbUpdate3;
        private System.Windows.Forms.Panel panel4;
        private System.Windows.Forms.ToolTip toolTip;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label lbBytes;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.Label lbDataType;
        private System.Windows.Forms.Label label20;
        private System.Windows.Forms.Label lbColumnWeight;
        private System.Windows.Forms.Label lbNzColumn;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label22;
        private System.Windows.Forms.Label label24;
        private System.Windows.Forms.Label label29;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.Panel panel5;
        private System.Windows.Forms.Label lbUpdate2;
        private System.Windows.Forms.Label label17;
        private System.Windows.Forms.Label lbNorm1;
        private System.Windows.Forms.Label label19;
        private System.Windows.Forms.Label lbNormMax;
        private System.Windows.Forms.Label label26;
        private System.Windows.Forms.Label lbNormFrob;
        private System.Windows.Forms.Label label21;
        private System.Windows.Forms.Label lbNormInf;
        private System.Windows.Forms.Label label34;
        private System.Windows.Forms.Label lbDomRows;
        private System.Windows.Forms.Label label31;
        private System.Windows.Forms.Label lbDomColumns;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label lbInfo3;
        private System.Windows.Forms.Label lbInfo2;
        private System.Windows.Forms.Label lbInfo1;
        private System.Windows.Forms.Label lbNumSymPcnt;
        private System.Windows.Forms.Label lbDomRowsPcnt;
        private System.Windows.Forms.Label lbDomColumnsPcnt;
        private System.Windows.Forms.ToolStripButton btnRefresh;
    }
}
