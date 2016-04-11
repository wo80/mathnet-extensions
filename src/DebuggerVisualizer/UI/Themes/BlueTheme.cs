using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.MatrixDebuggerVisualizer.UI.Themes
{
    class BlueTheme
    {
        public static class TabControl
        {
            public static readonly Color ForegroundNormal = Color.White;
            public static readonly Color ForegroundActive = Color.Black;
            public static readonly Color BackgroundNormal = Color.FromArgb(41, 57, 85);
            public static readonly Color BackgroundActive = Color.FromArgb(255, 240, 208);
        }

        public static class ToolStrip
        {
            public static readonly Color BackgroundNormal = Color.FromArgb(214, 219, 233);
            public static readonly Color BackgroundActive = Color.FromArgb(255, 240, 208);
        }

        public static class View
        {
            public static readonly Color ForegroundNormal = Color.Black;
            public static readonly Color BackgroundNormal = Color.WhiteSmoke;
        }

        public static class Panel
        {
            public static readonly Color ForegroundNormal = Color.Black;
            public static readonly Color BackgroundNormal = Color.FromArgb(242, 242, 242);
        }
    }
}
