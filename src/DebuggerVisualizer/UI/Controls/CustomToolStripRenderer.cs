// Windows 7 ToolStrip Renderer
//
// Andrea Martinelli
// http://at-my-window.blogspot.com/?page=windows7renderer
//
// Based on Office 2007 Renderer by Phil Wright
// http://www.componentfactory.com

using MathNet.MatrixDebuggerVisualizer.UI.Themes;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Windows.Forms;

namespace MathNet.MatrixDebuggerVisualizer.UI
{
    /// <summary>
    /// Draw ToolStrips using the Windows 7 themed appearance.
    /// </summary>
    public class CustomToolStripRenderer : ToolStripProfessionalRenderer
    {
        private static CustomToolStripRenderer _instance;

        public static CustomToolStripRenderer Instance
        {
            get
            {
                return _instance ?? (_instance = new CustomToolStripRenderer());
            }
        }

        #region Static Metrics

        private const int _gripOffset = 1;
        private const int _gripSquare = 2;
        private const int _gripSize = 3;
        private const int _gripMove = 4;
        private const int _gripLines = 3;
        private const int _checkInset = 0;
        private const int _marginInset = 2;
        private const int _separatorInset = 25;
        private const float _cutToolItemMenu = 0f;
        private const float _cutContextMenu = 0f;
        private const float _cutMenuItemBack = 0f;
        private const float _contextCheckTickThickness = 1.6f;
        private static Blend _statusStripBlend;

        #endregion

        #region Static Colors

        private static Color clrWindows7text = Color.FromArgb(30, 57, 91);
        private static Color clrDarkBottomBegin = BlueTheme.ToolStrip.BackgroundNormal;

        // Color scheme values
        private static Color _textDisabled = Color.FromArgb(167, 167, 167);
        private static Color _textStatusStripItem = clrWindows7text;
        private static Color _textContextMenuItem = clrWindows7text;

        private static Color _itemToolItemSelectedColors = Color.FromArgb(248, 251, 254);
        private static Color _itemToolItemPressedColors = Color.FromArgb(207, 220, 237);
        private static Color _itemToolItemCheckedColors = Color.FromArgb(207, 220, 237);
        private static Color _itemToolItemCheckPressColors = Color.FromArgb(207, 220, 237);

        private static Color _statusBarDark = Color.FromArgb(204, 217, 234);
        private static Color _statusBarLight = Color.FromArgb(241, 245, 251);

        private static Color _contextMenuBack = Color.FromArgb(240, 240, 240);

        private static Color _itemDisabledColors = Color.FromArgb(127, 236, 241, 247);

        #endregion

        static CustomToolStripRenderer()
        {
            // One time creation of the blend for the status strip gradient brush
            _statusStripBlend = new Blend();
            _statusStripBlend.Positions = new float[] { 0.0f, 0.25f, 0.25f, 0.57f, 0.86f, 1.0f };
            _statusStripBlend.Factors = new float[] { 0.1f, 0.6f, 1.0f, 0.4f, 0.0f, 0.95f };
        }

        /// <summary>
        /// Initialize a new instance of the Windows7Renderer class.
        /// </summary>
        public CustomToolStripRenderer()
            : base()
        {
            RoundedEdges = false;
        }

        protected override void InitializeItem(ToolStripItem item)
        {
            base.InitializeItem(item);

            if (item.DisplayStyle == ToolStripItemDisplayStyle.Image)
            {
                var m = item.Margin;
                m.Left += 4;
                m.Right += 4;
                item.Margin = m;
            }

            {
                var a = item as ToolStripSplitButton;
                if (a != null)
                {

                    a.DropDownButtonWidth = 17;

                    if (a.DisplayStyle == ToolStripItemDisplayStyle.Image)
                    {
                        a.Padding = new Padding(3, 0, 3, 0);
                    }
                    else if (a.DisplayStyle == ToolStripItemDisplayStyle.Text)
                    {
                        a.Padding = new Padding(14, 3, 15, 3);
                    }
                    else
                    {
                        a.Padding = new Padding(25, 3, 0, 3);
                        a.TextAlign = ContentAlignment.MiddleRight;
                    }
                }
            }

            {
                var a = item as ToolStripDropDownButton;
                if (a != null)
                {

                    if (a.DisplayStyle == ToolStripItemDisplayStyle.Image)
                    {
                        a.Padding = new Padding(7, 0, 7, 0);
                        // a.Margin = new Padding(4, 0, 4, 0);
                    }
                    else if (a.DisplayStyle == ToolStripItemDisplayStyle.Text)
                    {
                        a.Padding = new Padding(14, 3, 15, 3);
                    }
                    else
                    {
                        a.Padding = new Padding(25, 3, 0, 3);
                        a.TextAlign = ContentAlignment.MiddleRight;
                    }
                }
            }

            {
                var a = item as ToolStripButton;
                if (a != null)
                {

                    if (a.DisplayStyle == ToolStripItemDisplayStyle.Image)
                    {
                        a.Padding = new Padding(4, 0, 4, 0);
                        // a.Margin = new Padding(4, 0, 4, 0);
                    }
                    else if (a.DisplayStyle == ToolStripItemDisplayStyle.Text)
                    {
                        a.Padding = new Padding(15, 3, 15, 3);
                    }
                    else
                    {
                        a.Padding = new Padding(25, 3, 0, 3);
                        a.TextAlign = ContentAlignment.MiddleRight;
                    }
                }
            }

            if (item is ToolStripSeparator)
            {
                item.Height++;
            }

            if (item is ToolStripOverflowButton)
            {
                item.Width += 25;
            }
        }

        protected override void Initialize(ToolStrip toolStrip)
        {
            base.Initialize(toolStrip);

            if (toolStrip is MenuStrip)
            {
                toolStrip.CanOverflow = false;
            }
            else if (toolStrip is ContextMenuStrip)
            {
                // No nothing
            }
            else
            {

                toolStrip.AutoSize = false;
                toolStrip.Height = 32;

                toolStrip.Padding = new Padding(5, 2, 5, 2);
                toolStrip.CanOverflow = false;
                toolStrip.GripStyle = ToolStripGripStyle.Hidden;

            }

            //toolStrip.Font = new Font("Segoe UI", 9);
        }

        #region Render overrides

        /// <summary>
        /// Raises the RenderButtonBackground event. 
        /// </summary>
        /// <param name="e">An ToolStripItemRenderEventArgs containing the event data.</param>
        protected override void OnRenderButtonBackground(ToolStripItemRenderEventArgs e)
        {
            // Cast to correct type
            ToolStripButton button = (ToolStripButton)e.Item;

            if (button.Selected || button.Pressed || button.Checked)
                RenderToolButtonBackground(e.Graphics, button, e.ToolStrip);
        }

        /// <summary>
        /// Raises the RenderItemText event. 
        /// </summary>
        /// <param name="e">An ToolStripItemTextRenderEventArgs containing the event data.</param>
        protected override void OnRenderItemText(ToolStripItemTextRenderEventArgs e)
        {
            if ((e.ToolStrip is MenuStrip) ||
                (e.ToolStrip is ToolStrip) ||
                (e.ToolStrip is ContextMenuStrip) ||
                (e.ToolStrip is ToolStripDropDownMenu))
            {
                // We set the color depending on the enabled state
                if (!e.Item.Enabled)
                {
                    e.TextColor = _textDisabled;
                }
                else
                {
                    /*    if ((e.ToolStrip is MenuStrip) && !e.Item.Pressed && !e.Item.Selected)
                            e.TextColor = _textMenuStripItem;*/
                    if ((e.ToolStrip is StatusStrip) && !e.Item.Pressed && !e.Item.Selected)
                        e.TextColor = _textStatusStripItem;
                    else if ((e.ToolStrip is MenuStrip))
                        e.TextColor = _textStatusStripItem;
                    else if (e.ToolStrip is ToolStripDropDown)
                        e.TextColor = Color.Black;
                    else
                        e.TextColor = _textContextMenuItem;
                }

                e.TextRectangle = AdjustDrawRectangle(e.Item, e.TextRectangle);

                // All text is draw using the ClearTypeGridFit text rendering hint
                using (var clearTypeGridFit = new UseClearTypeGridFit(e.Graphics))
                {
                    base.OnRenderItemText(e);
                }
            }
            else
            {
                base.OnRenderItemText(e);
            }
        }

        /// <summary>
        /// Raises the RenderToolStripBackground event. 
        /// </summary>
        /// <param name="e">An ToolStripRenderEventArgs containing the event data.</param>
        protected override void OnRenderToolStripBackground(ToolStripRenderEventArgs e)
        {
            if ((e.ToolStrip is ContextMenuStrip) ||
                (e.ToolStrip is ToolStripDropDownMenu))
            {
                // Create border and clipping paths
                using (GraphicsPath borderPath = CreateBorderPath(e.AffectedBounds, _cutContextMenu),
                                      clipPath = CreateClipBorderPath(e.AffectedBounds, _cutContextMenu))
                {
                    // Clip all drawing to within the border path
                    using (var clipping = new UseClipping(e.Graphics, clipPath))
                    {
                        // Create the background brush
                        using (var backBrush = new SolidBrush(_contextMenuBack))
                            e.Graphics.FillPath(backBrush, borderPath);
                    }
                }
            }
            else if (e.ToolStrip is StatusStrip)
            {
                // We do not paint the top two pixel lines, so are drawn by the status strip border render method
                RectangleF backRect = new RectangleF(0, 1.5f, e.ToolStrip.Width, e.ToolStrip.Height - 2);

                // Cannot paint a zero sized area
                if ((backRect.Width > 0) && (backRect.Height > 0))
                {
                    using (var backBrush = new SolidBrush(_statusBarLight))
                    {
                        //backBrush.Blend = _statusStripBlend;
                        e.Graphics.FillRectangle(backBrush, backRect);
                    }
                    var topRect = new Rectangle(0, 0, (int)backRect.Width, 5);
                    using (var backBrush = new SolidBrush(_statusBarDark))
                    {
                        //backBrush.Blend = _statusStripBlend;
                        e.Graphics.FillRectangle(backBrush, topRect);
                    }
                }
            }
            else
            {
                //  base.OnRenderToolStripBackground(e);
                ToolStrip toolStrip = e.ToolStrip;

                if ((toolStrip.Height > 0) && (toolStrip.Width > 0))
                {
                    var height = toolStrip.Height;
                    var center = height / 2;
                    var width = toolStrip.Width;

                    var topRect = new Rectangle(0, 0, width, height);

                    //using (var topBrush = new SolidBrush(clrWindows7topBegin))
                    //    e.Graphics.FillRectangle(topBrush, topRect);

                    using (var bottomBrush = new SolidBrush(clrDarkBottomBegin))
                        e.Graphics.FillRectangle(bottomBrush, topRect);
                }
            }
        }

        #endregion

        #region Implementation

        private void RenderToolButtonBackground(Graphics g,
                                                ToolStripButton button,
                                                ToolStrip toolstrip)
        {
            // We only draw a background if the item is selected or being pressed
            if (button.Enabled)
            {
                if (button.Checked)
                {
                    if (button.Pressed)
                        DrawGradientToolItem(g, button, _itemToolItemPressedColors);
                    else if (button.Selected)
                        DrawGradientToolItem(g, button, _itemToolItemCheckPressColors);
                    else
                        DrawGradientToolItem(g, button, _itemToolItemCheckedColors);
                }
                else
                {
                    if (button.Pressed)
                        DrawGradientToolItem(g, button, _itemToolItemPressedColors);
                    else if (button.Selected)
                        DrawGradientToolItem(g, button, _itemToolItemSelectedColors);
                }
            }
            else
            {
                if (button.Selected)
                {
                    // Get the mouse position in tool strip coordinates
                    Point mousePos = toolstrip.PointToClient(Control.MousePosition);

                    // If the mouse is not in the item area, then draw disabled
                    if (!button.Bounds.Contains(mousePos))
                        DrawGradientToolItem(g, button, _itemDisabledColors);
                }
            }
        }

        private void DrawGradientToolItem(Graphics g,
                                          ToolStripItem item,
                                          Color colors)
        {
            // Perform drawing into the entire background of the item
            DrawGradientItem(g, new Rectangle(Point.Empty, item.Bounds.Size), colors);
        }

        private void DrawGradientItem(Graphics g,
                                      Rectangle backRect,
                                      Color colors)
        {
            // Cannot paint a zero sized area
            if ((backRect.Width > 0) && (backRect.Height > 0))
            {
                // Draw the background of the entire item
                DrawGradientBack(g, backRect, colors);

                // Draw the border of the entire item
                DrawGradientBorder(g, backRect, colors);
            }
        }

        private void DrawGradientBack(Graphics g,
                                      Rectangle backRect,
                                      Color colors)
        {
            // Reduce rect draw drawing inside the border
            backRect.Inflate(-1, -1);

            int y2 = backRect.Height / 2;
            Rectangle backRect1 = new Rectangle(backRect.X, backRect.Y, backRect.Width, y2);
            Rectangle backRect2 = new Rectangle(backRect.X, backRect.Y + y2, backRect.Width, backRect.Height - y2);
            Rectangle backRect1I = backRect1;
            Rectangle backRect2I = backRect2;
            backRect1I.Inflate(1, 1);
            backRect2I.Inflate(1, 1);

            using (SolidBrush insideBrush1 = new SolidBrush(colors))
            {
                // TODO: fix
                g.FillRectangle(insideBrush1, backRect1);
                g.FillRectangle(insideBrush1, backRect2);
            }

            y2 = backRect.Height / 2;
            backRect1 = new Rectangle(backRect.X, backRect.Y, backRect.Width, y2);
            backRect2 = new Rectangle(backRect.X, backRect.Y + y2, backRect.Width, backRect.Height - y2);
            backRect1I = backRect1;
            backRect2I = backRect2;
            backRect1I.Inflate(1, 1);
            backRect2I.Inflate(1, 1);

            using (SolidBrush fillBrush1 = new SolidBrush(colors))
            {
                // TODO: fix

                // Reduce rect one more time for the innermost drawing
                backRect.Inflate(-1, -1);

                y2 = (backRect.Height / 2) + 1;
                backRect1 = new Rectangle(backRect.X, backRect.Y, backRect.Width, y2);
                backRect2 = new Rectangle(backRect.X, backRect.Y + y2, backRect.Width, backRect.Height - y2);

                g.FillRectangle(fillBrush1, backRect1);
                g.FillRectangle(fillBrush1, backRect2);
            }
        }

        private void DrawGradientBorder(Graphics g,
                                        Rectangle backRect,
                                        Color colors)
        {
            // Drawing with anti aliasing to create smoother appearance
            using (UseAntiAlias uaa = new UseAntiAlias(g))
            {
                Rectangle backRectI = backRect;
                backRectI.Inflate(1, 1);

                // Finally draw the border around the menu item
                using (var borderBrush = new SolidBrush(colors))
                {
                    // Convert the brush to a pen for DrawPath call
                    using (var borderPen = new Pen(borderBrush))
                    {
                        // Create border path around the entire item
                        using (var borderPath = CreateBorderPath(backRect, _cutMenuItemBack))
                            g.DrawPath(borderPen, borderPath);
                    }
                }
            }
        }

        private GraphicsPath CreateBorderPath(Rectangle rect,
                                              Rectangle exclude,
                                              float cut)
        {
            // If nothing to exclude, then use quicker method
            if (exclude.IsEmpty)
                return CreateBorderPath(rect, cut);

            // Drawing lines requires we draw inside the area we want
            rect.Width--;
            rect.Height--;

            // Create an array of points to draw lines between
            List<PointF> pts = new List<PointF>();

            float l = rect.X;
            float t = rect.Y;
            float r = rect.Right;
            float b = rect.Bottom;
            float x0 = rect.X + cut;
            float x3 = rect.Right - cut;
            float y0 = rect.Y + cut;
            float y3 = rect.Bottom - cut;
            float cutBack = (cut == 0f ? 1 : cut);

            // Does the exclude intercept the top line
            if ((rect.Y >= exclude.Top) && (rect.Y <= exclude.Bottom))
            {
                float x1 = exclude.X - 1 - cut;
                float x2 = exclude.Right + cut;

                if (x0 <= x1)
                {
                    pts.Add(new PointF(x0, t));
                    pts.Add(new PointF(x1, t));
                    pts.Add(new PointF(x1 + cut, t - cutBack));
                }
                else
                {
                    x1 = exclude.X - 1;
                    pts.Add(new PointF(x1, t));
                    pts.Add(new PointF(x1, t - cutBack));
                }

                if (x3 > x2)
                {
                    pts.Add(new PointF(x2 - cut, t - cutBack));
                    pts.Add(new PointF(x2, t));
                    pts.Add(new PointF(x3, t));
                }
                else
                {
                    x2 = exclude.Right;
                    pts.Add(new PointF(x2, t - cutBack));
                    pts.Add(new PointF(x2, t));
                }
            }
            else
            {
                pts.Add(new PointF(x0, t));
                pts.Add(new PointF(x3, t));
            }

            pts.Add(new PointF(r, y0));
            pts.Add(new PointF(r, y3));
            pts.Add(new PointF(x3, b));
            pts.Add(new PointF(x0, b));
            pts.Add(new PointF(l, y3));
            pts.Add(new PointF(l, y0));

            // Create path using a simple set of lines that cut the corner
            GraphicsPath path = new GraphicsPath();

            // Add a line between each set of points
            for (int i = 1; i < pts.Count; i++)
                path.AddLine(pts[i - 1], pts[i]);

            // Add a line to join the last to the first
            path.AddLine(pts[pts.Count - 1], pts[0]);

            return path;
        }

        private GraphicsPath CreateBorderPath(Rectangle rect, float cut)
        {
            // Drawing lines requires we draw inside the area we want
            rect.Width--;
            rect.Height--;

            // Create path using a simple set of lines that cut the corner
            GraphicsPath path = new GraphicsPath();
            path.AddLine(rect.Left + cut, rect.Top, rect.Right - cut, rect.Top);
            path.AddLine(rect.Right - cut, rect.Top, rect.Right, rect.Top + cut);
            path.AddLine(rect.Right, rect.Top + cut, rect.Right, rect.Bottom - cut);
            path.AddLine(rect.Right, rect.Bottom - cut, rect.Right - cut, rect.Bottom);
            path.AddLine(rect.Right - cut, rect.Bottom, rect.Left + cut, rect.Bottom);
            path.AddLine(rect.Left + cut, rect.Bottom, rect.Left, rect.Bottom - cut);
            path.AddLine(rect.Left, rect.Bottom - cut, rect.Left, rect.Top + cut);
            path.AddLine(rect.Left, rect.Top + cut, rect.Left + cut, rect.Top);
            return path;
        }

        private GraphicsPath CreateClipBorderPath(Rectangle rect, float cut)
        {
            // Clipping happens inside the rect, so make 1 wider and taller
            rect.Width++;
            rect.Height++;

            // Now create a path based on this inner rectangle
            return CreateBorderPath(rect, cut);
        }

        private GraphicsPath CreateClipBorderPath(Rectangle rect,
                                                  Rectangle exclude,
                                                  float cut)
        {
            // Clipping happens inside the rect, so make 1 wider and taller
            rect.Width++;
            rect.Height++;
            //rect.Inflate(10, 10);
            // Now create a path based on this inner rectangle
            return CreateBorderPath(rect, exclude, cut);
        }

        #endregion

        private Rectangle AdjustDrawRectangle(ToolStripItem item, Rectangle rectangle)
        {
            if (LooksPressed(item))
            {
                var r = rectangle;
                r.X++;
                r.Y++;
                return r;
            }

            return rectangle;

        }

        private bool LooksPressed(ToolStripItem item)
        {
            var parent = item.GetCurrentParent();

            if (parent is MenuStrip) return false;
            if (parent is ContextMenuStrip) return false;
            if (parent is ToolStripDropDown) return false;

            if (item.Pressed) return true;

            {
                var a = item as ToolStripDropDownButton;
                if (a != null && a.IsOnDropDown) return true;
            }

            {
                var a = item as ToolStripSplitButton;
                if (a != null && a.IsOnDropDown) return true;
            }

            return false;
        }
    }
}
