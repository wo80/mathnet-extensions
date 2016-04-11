using System.Drawing;
using System.Windows.Forms;

namespace ExtensionsBenchmark
{
    public class CustomListBox : ListBox
    {
        TextFormatFlags TextFormat = TextFormatFlags.VerticalCenter | TextFormatFlags.SingleLine | TextFormatFlags.WordEllipsis;
        
        public CustomListBox()
        {
            base.DrawMode = DrawMode.OwnerDrawVariable;
            base.BorderStyle = BorderStyle.None;

            this.selectionColor = Color.DodgerBlue;
        }

        protected Color selectionColor;

        public Color SelectionColor
        {
            get { return selectionColor; }
            set { selectionColor = value; }
        }

        protected override void OnDrawItem(DrawItemEventArgs e)
        {
            if (this.Items.Count == 0)
            {
                return;
            }

            var g = e.Graphics;

            var background = e.State.HasFlag(DrawItemState.Selected) ? selectionColor : this.BackColor;

            using (var brush = new SolidBrush(background))
            {
                g.FillRectangle(brush, e.Bounds);
            }

            e.DrawFocusRectangle();

            string text = this.Items[e.Index].ToString();

            TextRenderer.DrawText((IDeviceContext)g, text, base.Font, e.Bounds, e.ForeColor, TextFormat);
        }
    }
}
