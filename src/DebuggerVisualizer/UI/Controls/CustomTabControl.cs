
namespace MathNet.MatrixDebuggerVisualizer.UI
{
    using System;
    using System.ComponentModel;
    using System.Drawing;
    using System.Windows.Forms;
    using System.Drawing.Text;
    using MathNet.MatrixDebuggerVisualizer.UI.Themes;

    // Original code on CodeProject: A .NET Flat TabControl (CustomDraw), Oscar Londono

    /// <summary>
    /// Summary description for FlatTabControl.
    /// </summary>
    public class CustomTabControl : System.Windows.Forms.TabControl
    {
        #region Designer

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
            components = new System.ComponentModel.Container();
        }

        #endregion

        #endregion
        
        private const int margin = 5;

        private Color foreColor = BlueTheme.TabControl.ForegroundNormal;
        private Color foreColorActive = BlueTheme.TabControl.ForegroundActive;
        private Color backColor = BlueTheme.TabControl.BackgroundNormal;
        private Color backColorActive = BlueTheme.TabControl.BackgroundActive;

        /// <summary>
        /// Initializes a new instance of the <see cref="CustomTabControl" /> control.
        /// </summary>
		public CustomTabControl()
		{
			// This call is required by the Windows.Forms Form Designer.
			InitializeComponent();

            base.Multiline = false;

			// double buffering
			this.SetStyle(ControlStyles.UserPaint, true);
			this.SetStyle(ControlStyles.AllPaintingInWmPaint, true);
			this.SetStyle(ControlStyles.DoubleBuffer, true);
			this.SetStyle(ControlStyles.ResizeRedraw, true);
			this.SetStyle(ControlStyles.SupportsTransparentBackColor, true);

            this.SelectedIndexChanged += (obj, evt) => { Invalidate(); };
		}

        #region Properties

        new public TabAlignment Alignment
        {
            get { return base.Alignment; }
            set
            {
                TabAlignment ta = value;
                if ((ta != TabAlignment.Top) && (ta != TabAlignment.Bottom))
                {
                    ta = TabAlignment.Top;
                }

                base.Alignment = ta;
            }
        }

        #endregion

		protected override void OnPaint(PaintEventArgs e)
		{
			base.OnPaint(e); 
			
			DrawControl(e.Graphics);
		}

        private void DrawControl(Graphics g)
		{
            if (!Visible)
            {
                return;
            }

			var client = this.ClientRectangle;
            var tab = this.DisplayRectangle;

			// Fill client area
            using (var b = new SolidBrush(backColor))
            {
                g.FillRectangle(b, client);
            }

			int width = tab.Width + margin;

			// Clip region for drawing tabs
			Region clip = g.Clip;
			Rectangle region = new Rectangle(tab.Left, client.Top, width - margin, client.Height);

			g.SetClip(region);

			// Draw tabs
            for (int i = 0; i < this.TabCount; i++)
            {
                DrawTab(g, this.TabPages[i], i);
            }

			g.Clip = clip;
		}

        private void DrawTab(Graphics g, TabPage tabPage, int index)
		{
			var bounds = this.GetTabRect(index);

			bool selected = (this.SelectedIndex == index);

            using (var b = new SolidBrush(selected ? backColorActive : backColor))
            {
                // Fill this tab with background color
                g.FillRectangle(b, bounds);
            }

            if (selected)
            {
                // Clear bottom lines
                Pen pen = new Pen(tabPage.BackColor);

                switch (this.Alignment)
                {
                    case TabAlignment.Top:
                        g.DrawLine(pen, bounds.Left, bounds.Bottom, bounds.Right - 1, bounds.Bottom);
                        g.DrawLine(pen, bounds.Left, bounds.Bottom + 1, bounds.Right - 1, bounds.Bottom + 1);
                        break;

                    case TabAlignment.Bottom:
                        g.DrawLine(pen, bounds.Left, bounds.Top, bounds.Right - 1, bounds.Top);
                        g.DrawLine(pen, bounds.Left, bounds.Top - 1, bounds.Right - 1, bounds.Top - 1);
                        g.DrawLine(pen, bounds.Left, bounds.Top - 2, bounds.Right - 1, bounds.Top - 2);
                        break;
                }

                pen.Dispose();
            }

			// Draw string
			var format = new StringFormat();

			format.Alignment = StringAlignment.Center;  
			format.LineAlignment = StringAlignment.Center;

            g.TextRenderingHint = TextRenderingHint.ClearTypeGridFit;

            using (var b = new SolidBrush(selected ? foreColorActive : foreColor))
            {
                g.DrawString(tabPage.Text, Font, b, bounds, format);
            }
		}
	}
}
