
namespace MathNet.MatrixDebuggerVisualizer.Services
{
    using System;
    using System.Collections.Generic;

    public interface IStorageAdapter
    {
        int RowCount { get; }
        int ColumnCount { get; }

        object GetMatrix();
        object GetStorage();

        string GetDataTypeName();

        void Save(string file);

        StorageInfo GetStorageInfo(int level);

        /// <summary>
        /// Get matrix value at index (i, j).
        /// </summary>
        /// <param name="i">Row index</param>
        /// <param name="j">Column index</param>
        /// <returns>Tuple: [real part, imaginary part]</returns>
        Tuple<double, double> GetValue(int i, int j);

        /// <summary>
        /// Enumerates all non-zero values of the matrix.
        /// </summary>
        /// <returns></returns>
        IEnumerable<MatrixEntry> Enumerate();

        /// <summary>
        /// Enumerates all non-zero values of the submatrix.
        /// </summary>
        /// <param name="rowStart"></param>
        /// <param name="rowEnd"></param>
        /// <param name="colStart"></param>
        /// <param name="colEnd"></param>
        /// <returns></returns>
        IEnumerable<MatrixEntry> EnumerateSubmatrix(int rowStart, int rowEnd, int colStart, int colEnd);

        /// <summary>
        /// Enumerates all non-zero values of the row with given index.
        /// </summary>
        /// <param name="i">Row index.</param>
        /// <returns></returns>
        IEnumerable<MatrixEntry> EnumerateRow(int i);
    }
}
