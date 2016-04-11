
namespace MathNet.MatrixDebuggerVisualizer.Services
{
    /// <summary>
    /// Provides information on a sparse matrix.
    /// </summary>
    public class StorageInfo
    {
        /// <summary>Row count</summary>
        public int RowCount { get; set; }

        /// <summary>Column count</summary>
        public int ColumnCount { get; set; }

        /// <summary>Nonzeros count</summary>
        public int ValueCount { get; set; }

        /// <summary>Number of bytes allocated by the storage.</summary>
        public long TotalBytes { get; set; }


        /// <summary>number of nonzero elements in strict lower part</summary>
        public int NonZerosLower { get; set; }

        /// <summary>number of nonzero elements in strict upper part</summary>
        public int NonZerosUpper { get; set; }

        /// <summary>number of nonzero elements in main diagonal</summary>
        public int NonZerosDiagonal { get; set; }


        /// <summary>Average number of nonzero elements/row.</summary>
        public double NonZerosPerRow { get; set; }

        /// <summary>Standard deviation for the average number of nonzero elements/row.</summary>
        public double NonZerosPerRowDev { get; set; }

        /// <summary>Average number of nonzero elements/column.</summary>
        public double NonZerosPerColumn { get; set; }

        /// <summary>Standard deviation for the average number of nonzero elements/column.</summary>
        public double NonZerosPerColumnDev { get; set; }


        /// <summary>max length of rows</summary>
        public int MaxRowLength { get; set; }

        /// <summary>min length of rows</summary>
        public int MinRowLength { get; set; }

        /// <summary>number of zero rows</summary>
        public int ZeroRowsCount { get; set; }

        /// <summary>max length of columns</summary>
        public int MaxColumnLength { get; set; }

        /// <summary>min length of columns</summary>
        public int MinColumnLength { get; set; }

        /// <summary>number of zero columns</summary>
        public int ZeroColumnsCount { get; set; }


        /// <summary>Lower bandwidth</summary>
        public int LowerBandwidth { get; set; }

        /// <summary>Upper bandwidth</summary>
        public int UpperBandwidth { get; set; }

        /// <summary>Maximum bandwidth</summary>
        public int MaximumBandwidth { get; set; }

        /// <summary>Average bandwidth</summary>
        public double AverageBandwidth { get; set; }


        /// <summary>average distance of a(i,j) from diag.</summary>
        public double DiagonalDistanceAverage { get; set; }

        /// <summary>standard deviation for above average.</summary>
        public double DiagonalDistanceDeviation { get; set; }

        /// <summary>number of elements in each diagonal.</summary>
        /// <remarks>
        /// int array containing the numbers of elements in each of
        /// the nrow+ncol-1 diagonals of A. dist(k) contains the number of elements in
        /// diagonal '-nrow+k'. k ranges from 1 to (nrow+ncol-1).
        /// </remarks>
        public int[] DiagonalsCount { get; set; }


        /// <summary>number of weakly diagonally dominant rows.</summary>
        public int DiagonallyDominantRows { get; set; }

        /// <summary>number of weakly diagonally dominant columns.</summary>
        public int DiagonallyDominantColumns { get; set; }


        /// <summary>Max. absolute value.</summary>
        public double MaxAbsoluteValue { get; set; }

        /// <summary>One norm.</summary>
        public double OneNorm { get; set; }

        /// <summary>Infinity norm.</summary>
        public double InfNorm { get; set; }

        /// <summary>Frobenius norm.</summary>
        public double FrobeniusNorm { get; set; }


        /// <summary></summary>
        public int Band80 { get; set; }

        /// <summary>Max. absolute value.</summary>
        public int Band90 { get; set; }


        /// <summary>Number of matching elements in symmetry.</summary>
        public int ElementsInSymmetry { get; set; }

        /// <summary>Frobenius norm of symmetric part</summary>
        public double FrobeniusNormSymPart { get; set; }

        /// <summary>Frobenius norm of non-symmetric part</summary>
        public double FrobeniusNormNonSymPart { get; set; }
    }
}
