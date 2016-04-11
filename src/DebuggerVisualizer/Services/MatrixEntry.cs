
namespace MathNet.MatrixDebuggerVisualizer.Services
{
    public struct MatrixEntry
    {
        public int Row;
        public int Column;
        public double Value;
        public double Z;
        public int Complex;
        public int Info;

        public MatrixEntry(int i, int j, double value)
        {
            this.Row = i;
            this.Column = j;
            this.Value = value;
            this.Z = 0.0;
            this.Complex = 0;
            this.Info = 0;
        }

        public MatrixEntry(int i, int j, double value, double z)
        {
            this.Row = i;
            this.Column = j;
            this.Value = value;
            this.Z = z;
            this.Complex = 1;
            this.Info = 0;
        }
    }
}
