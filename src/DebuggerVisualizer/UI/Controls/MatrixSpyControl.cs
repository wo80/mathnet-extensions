using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Drawing;
using System.Drawing.Drawing2D;
using MathNet.MatrixDebuggerVisualizer.Services;

namespace MathNet.MatrixDebuggerVisualizer.UI.Controls
{
    class MatrixSpyControl : System.Windows.Forms.Control
    {
        // Rendering stuff
        private BufferedGraphics buffer;
        private BufferedGraphicsContext context;
        private Zoom zoom;
        private Timer timer;

        string text = String.Empty;

        bool initialized = false;

        IStorageAdapter storage;

        int rows, columns;

        public MatrixSpyControl()
        {
            //this.SetStyle(ControlStyles.UserPaint, true);
            //this.SetStyle(ControlStyles.OptimizedDoubleBuffer, false);
            this.SetStyle(ControlStyles.ResizeRedraw, true);

            this.BackColor = Color.White;

            context = BufferedGraphicsManager.Current;

            timer = new Timer();
            timer.Interval = 5 * 1000;
            timer.Tick += (sender, e) =>
            {
                timer.Stop();
                text = String.Empty;
                this.Invalidate();
            };
        }

        /// <summary>
        /// Initialize the graphics buffer (should be called in the forms load event).
        /// </summary>
        public void Initialize(IStorageAdapter storage)
        {
            this.storage = storage;

            var info = storage.GetStorageInfo(0);

            rows = info.RowCount;
            columns = info.ColumnCount;


            zoom = new Zoom();
            zoom.Initialize(this.ClientRectangle, rows, columns);

            InitializeBuffer();

            this.Invalidate();
        }

        public void Zoom(int delta, float x, float y)
        {
            if (!initialized) return;

            if (zoom.Update(delta, x, y))
            {
                // Redraw
                this.Render();
            }
        }

        public void Reset()
        {
            zoom.Reset();
            this.Render();
        }

        public override void Refresh()
        {
            this.Render();
        }

        private void InitializeBuffer()
        {
            if (this.Width > 0 && this.Height > 0)
            {
                if (buffer != null)
                {
                    if (this.ClientRectangle == buffer.Graphics.VisibleClipBounds)
                    {
                        this.Invalidate();

                        // Bounds didn't change. Probably we just restored the
                        // window from minimized state.
                        return;
                    }

                    buffer.Dispose();
                }

                buffer = context.Allocate(this.CreateGraphics(), this.ClientRectangle);
                initialized = true;

                this.Render();
            }
        }

        #region Rendering

        private void Render()
        {
            text = String.Empty;

            var g = buffer.Graphics;

            g.Clear(this.BackColor);

            if (!initialized || storage == null)
            {
                return;
            }

            int x = 0, y = 0;

            var view = zoom.Viewport;

            int colStart = Math.Max(0, view.Left);
            int rowStart = Math.Max(0, view.Top);
            int colEnd = Math.Min(columns, view.Right);
            int rowEnd = Math.Min(rows, view.Bottom);

            int size = zoom.BlockSize;

            if (rowStart < rowEnd && colStart < colEnd && zoom.Level > 2)
            {
                foreach (var item in storage.EnumerateSubmatrix(rowStart, rowEnd, colStart, colEnd))
                {
                    x = item.Column;
                    y = item.Row;

                    zoom.WorldToScreen(ref x, ref y);

                    g.FillRectangle(Brushes.Blue, x, y, size, size);
                }
            }
            else
            {
                foreach (var item in storage.Enumerate())
                {
                    x = item.Column;
                    y = item.Row;

                    zoom.WorldToScreen(ref x, ref y);

                    g.FillRectangle(Brushes.Blue, x, y, size, size);
                }
            }

            this.Invalidate();
        }

        #endregion

        #region Protected overrides

        protected override void OnKeyUp(KeyEventArgs e)
        {
            bool update = false;

            if (e.KeyCode == Keys.Left)
            {
                update = zoom.Translate(-1, 0);
            }
            else if (e.KeyCode == Keys.Right)
            {
                update = zoom.Translate(1, 0);
            }
            else if (e.KeyCode == Keys.Down)
            {
                update = zoom.Translate(0, 1);
            }
            else if (e.KeyCode == Keys.Up)
            {
                update = zoom.Translate(0, -1);
            }

            if (update)
            {
                this.Render();
            }
        }

        protected override void OnMouseWheel(MouseEventArgs e)
        {
            //base.OnMouseWheel(e);

            if (!initialized) return;

            float x = e.X / ((float)this.ClientRectangle.Width);
            float y = e.Y / ((float)this.ClientRectangle.Height);

            this.Zoom(e.Delta, x, y);
        }

        protected override void OnMouseClick(MouseEventArgs e)
        {
            base.OnMouseClick(e);

            if (!initialized) return;

            // Enables mouse wheel.
            this.Focus();

            if (e.Button == MouseButtons.Middle)
            {
                this.Reset();
            }
            /*
            else if (button == MouseButtons.Left)
            {
                timer.Stop();

                var nfi = System.Globalization.CultureInfo.InvariantCulture.NumberFormat;

                PointF c = new PointF(x / this.Width, y / this.Height);
                zoom.ScreenToWorld(ref c);
                coordinate = String.Format(nfi, "X:{0} Y:{1}", c.X, c.Y);

                this.Invalidate();

                timer.Start();
            }
            //*/
        }

        protected override void OnResize(EventArgs e)
        {
            base.OnResize(e);

            if (!initialized) return;

            zoom.Initialize(this.ClientRectangle, rows, columns);
            InitializeBuffer();
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            if (!initialized)
            {
                base.OnPaint(e);
                return;
            }

            buffer.Render();

            /*
            if (!String.IsNullOrEmpty(coordinate) && Renderer.Context.HasData)
            {
                Graphics g = e.Graphics;
                g.TextRenderingHint = TextRenderingHint.ClearTypeGridFit;
                g.DrawString(coordinate, this.Font, Brushes.White, 10, 10);
            }
            //*/
        }

        protected override void OnPaintBackground(PaintEventArgs pevent)
        {
            // Do nothing
            if (!initialized)
            {
                base.OnPaintBackground(pevent);
            }
        }

        #endregion
    }
}
