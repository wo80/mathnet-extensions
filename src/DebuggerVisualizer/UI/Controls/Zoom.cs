using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.MatrixDebuggerVisualizer.UI.Controls
{
    /// <summary>
    /// Helper class for martix spy renderer.
    /// </summary>
    class Zoom
    {
        // We have a user control with bounds
        //
        //    (0, 0) (ClientWidth, ClientHeight)
        //
        // and a matrix with dimension
        //
        //    (0, 0) (RowCount - 1, ColumnCount - 1)
        //
        // Ininitially, we want to render the complete non-zero pattern of
        // the matrix centered on the control.
        //
        // The user should be able to zoom into the matrix and select matrix
        // entries with the mouse.
        //
        // The zoom should stop, if one matrix entry takes an area of 5x5 pixels.
        
        // The screen.
        Rectangle screen;

        // The matrix.
        Rectangle matrix;

        /// <summary>
        /// Gets or sets the current viewport (visible submatrix).
        /// </summary>
        public Rectangle Viewport { get; set; }

        /// <summary>
        /// Gets the number of pixels that one matrix entry should take on the screen.
        /// </summary>
        public int BlockSize
        {
            get { return Math.Max(1, (int)(screen.Width / (1.3f * Viewport.Width))); }
        }

        /// <summary>
        /// Gets the zoom level.
        /// </summary>
        public int Level { get; private set; }

        /// <summary>
        /// Gets or sets a clip margin (default is 5% of viewport width on each side).
        /// </summary>
        public float ClipMargin { get; set; }

        int maxZoomLevel = 50;

        public void Initialize(Rectangle screen, int rows, int columns)
        {
            this.screen = screen;

            this.Level = 1;

            // Add a margin so there's some space around the border
            float worldMargin = (columns < rows) ? rows * 0.05f : columns * 0.05f;

            // Get the initial viewport (complete mesh centered on the screen)
            float screenRatio = screen.Width / (float)screen.Height;
            float worldRatio = columns / (float)rows;

            float scale = (columns + worldMargin) / screen.Width;

            if (screenRatio > worldRatio)
            {
                scale = (rows + worldMargin) / screen.Height;
            }

            float centerX = columns / 2f;
            float centerY = rows / 2f;

            // TODO: Add initial margin
            this.Viewport = new Rectangle(
                (int)(centerX - screen.Width * scale / 2),
                (int)(centerY - screen.Height * scale / 2),
                (int)(screen.Width * scale),
                (int)(screen.Height * scale));

            this.ClipMargin = this.Viewport.Width * 0.05f;

            this.matrix = this.Viewport;
        }

        /// <summary>
        /// Zoom in or out of the viewport.
        /// </summary>
        /// <param name="dx">Relative x point position</param>
        /// <param name="dy">Relative y point position</param>
        public bool Translate(int dx, int dy)
        {
            if (Level < 2)
            {
                return false;
            }

            // Current viewport
            float x = Viewport.X;
            float y = Viewport.Y;

            float size = 1f / BlockSize;

            if (dx > 0)
            {
                x += 50 * size;
            }
            else if (dx < 0)
            {
                x -= 50 * size;
            }

            if (dy > 0)
            {
                y += 50 * size;
            }
            else if (dy < 0)
            {
                y -= 50 * size;
            }

            // Set new viewport
            Viewport = new Rectangle((int)x, (int)y, Viewport.Width, Viewport.Height);

            return true;
        }

        /// <summary>
        /// Zoom in or out of the viewport.
        /// </summary>
        /// <param name="amount">Zoom amount</param>
        /// <param name="focusX">Relative x point position</param>
        /// <param name="focusY">Relative y point position</param>
        public bool Update(int amount, float focusX, float focusY)
        {
            float width, height;

            if (amount > 0) // Zoom in
            {
                if (BlockSize >= 10)
                {
                    return false;
                }

                this.Level++;

                if (this.Level > maxZoomLevel)
                {
                    this.Level = maxZoomLevel;
                    return false;
                }

                width = Viewport.Width / 1.5f;
                height = Viewport.Height / 1.5f;
            }
            else
            {
                this.Level--;

                if (this.Level < 1)
                {
                    this.Level = 1;
                    this.Viewport = this.matrix;
                    return false;
                }

                width = Viewport.Width * 1.5f;
                height = Viewport.Height * 1.5f;
            }

            // Current focus on viewport
            float x = Viewport.X + Viewport.Width * focusX;
            float y = Viewport.Y + Viewport.Height * focusY;

            // New left and top positions
            x = x - width * focusX;
            y = y - height * focusY;

            // Check if outside of world
            if (x < matrix.X)
            {
                x = matrix.X;
            }
            else if (x + width > matrix.Right)
            {
                x = matrix.Right - width;
            }

            if (y < matrix.Y)
            {
                y = matrix.Y;
            }
            else if (y + height > matrix.Bottom)
            {
                y = matrix.Bottom - height;
            }

            // Set new viewport
            this.Viewport = new Rectangle((int)x, (int)y, (int)width, (int)height);

            this.ClipMargin = this.Viewport.Width * 0.05f;

            return true;
        }

        public void Reset()
        {
            this.Viewport = this.matrix;
            this.Level = 1;
        }

        public void WorldToScreen(ref int x, ref int y)
        {
            x = (int)((x - Viewport.X) / (1f * Viewport.Width) * screen.Width);
            y = (int)((y - Viewport.Y) / (1f * Viewport.Height) * screen.Height);
        }

        public void ScreenToWorld(ref int x, ref int y)
        {
            x = Viewport.X + Viewport.Width * x;
            y = Viewport.Y + Viewport.Height * y;
        }
    }
}
