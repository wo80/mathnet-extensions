
namespace DebuggerVisualizerTest
{
    using System.Text.RegularExpressions;

    static class Util
    {
        public static bool GetDenseRandomSize(string text, int value, out int rows, out int cols)
        {
            rows = value;
            cols = value;

            var match = Regex.Match(text, "(?<rows>\\d+)\\s?x\\s?(?<cols>\\d+)");

            if (match.Success)
            {
                int.TryParse(match.Groups["rows"].Value.Trim(), out rows);
                int.TryParse(match.Groups["cols"].Value.Trim(), out cols);

                return true;
            }

            return false;
        }

        public static bool GetSparseRandomSize(string text, out int rows, out int cols, out bool symmetric)
        {
            bool success = GetDenseRandomSize(text, 100, out rows, out cols);

            symmetric = text.ToLower().Contains("symmetric");

            return success;
        }

        public static bool GetSpecialSize(string text, out int nx, out int ny, out bool laplacian)
        {
            nx = 20;
            ny = 20;

            laplacian = true;

            var match = Regex.Match(text, "nx = (?<nx>\\d+), ny = (?<ny>\\d+)");

            if (match.Success)
            {
                int.TryParse(match.Groups["nx"].Value.Trim(), out nx);
                int.TryParse(match.Groups["ny"].Value.Trim(), out ny);

                laplacian = text.StartsWith("Laplacian");

                return true;
            }

            return false;
        }
    }
}
