
namespace MathNet.MatrixDebuggerVisualizer.Complex
{
    using MathNet.MatrixDebuggerVisualizer.Services;
    using MathNet.Numerics.Data.Matlab;
    using MathNet.Numerics.Data.Text;
    using MathNet.Numerics.LinearAlgebra.Complex;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;
    using System.Numerics;
    using System.Collections.Generic;

    public class DenseStorageAdapter : IStorageAdapter
    {
        DenseMatrix matrix;
        DenseColumnMajorMatrixStorage<Complex> storage;

        DenseStorageInfoService<Complex> service;

        public DenseStorageAdapter(DenseMatrix matrix)
        {
            this.matrix = matrix;
            this.storage = matrix.Storage as DenseColumnMajorMatrixStorage<Complex>;

            this.service = new DenseStorageInfoService<Complex>(matrix, new MathService());
        }

        public int RowCount { get { return this.storage.RowCount; } }

        public int ColumnCount { get { return this.storage.ColumnCount; } }

        public object GetMatrix()
        {
            return matrix;
        }

        public object GetStorage()
        {
            return storage;
        }

        public string GetDataTypeName()
        {
            return "Complex";
        }

        public StorageInfo GetStorageInfo(int level)
        {
            return service.StorageInfo;
        }

        public void Save(string file)
        {
            if (file.EndsWith(".mat"))
            {
                MatlabWriter.Write<Complex>(file, matrix, "A");
            }
            else if (file.EndsWith(".mtx"))
            {
                MatrixMarketWriter.WriteMatrix<Complex>(file, matrix);
            }
            else if (file.EndsWith(".csv"))
            {
                DelimitedWriter.Write<Complex>(file, matrix, ";");
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        public Tuple<double, double> GetValue(int i, int j)
        {
            if (i >= 0 && j >= 0 && i < storage.RowCount && j < storage.ColumnCount)
            {
                var value = matrix.At(i, j);

                return new Tuple<double, double>(value.Real, value.Imaginary);
            }

            return new Tuple<double, double>(0.0, 0.0);
        }

        public IEnumerable<MatrixEntry> Enumerate()
        {
            int rows = storage.RowCount;
            int cols = storage.ColumnCount;

            Complex a;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    a = matrix[i, j];

                    yield return new MatrixEntry(i, j, a.Real, a.Imaginary);
                }
            }
        }

        public IEnumerable<MatrixEntry> EnumerateSubmatrix(int rowStart, int rowEnd, int colStart, int colEnd)
        {
            int rows = Math.Min(rowEnd, storage.RowCount);
            int cols = Math.Min(colEnd, storage.ColumnCount);

            Complex a;

            for (int i = rowStart; i < rows; i++)
            {
                for (int j = colStart; j < cols; j++)
                {
                    a = matrix[i, j];

                    yield return new MatrixEntry(i, j, a.Real, a.Imaginary);
                }
            }
        }

        public IEnumerable<MatrixEntry> EnumerateRow(int i)
        {
            int rows = storage.RowCount;
            int cols = storage.ColumnCount;

            Complex a;

            for (int j = 0; j < cols; j++)
            {
                a = matrix[i, j];

                yield return new MatrixEntry(i, j, a.Real, a.Imaginary);
            }
        }
    }
}
