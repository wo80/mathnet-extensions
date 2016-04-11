
namespace MathNet.MatrixDebuggerVisualizer.Double
{
    using MathNet.MatrixDebuggerVisualizer.Services;
    using MathNet.Numerics.Data.Matlab;
    using MathNet.Numerics.Data.Text;
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;
    using System.Collections.Generic;

    public class SparseStorageAdapter : IStorageAdapter
    {
        SparseMatrix matrix;
        SparseCompressedRowMatrixStorage<double> storage;

        SparseStorageInfoService<double> service;

        public SparseStorageAdapter(SparseMatrix matrix)
        {
            this.matrix = matrix;
            this.storage = matrix.Storage as SparseCompressedRowMatrixStorage<double>;

            this.service = new SparseStorageInfoService<double>(matrix, new MathService());
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
            return "Double";
        }

        public StorageInfo GetStorageInfo(int level)
        {
            service.Update(level);

            return service.StorageInfo;
        }

        public void Save(string file)
        {
            if (file.EndsWith(".mat"))
            {
                MatlabWriter.Write<double>(file, matrix, "A");
            }
            else if (file.EndsWith(".mtx"))
            {
                MatrixMarketWriter.WriteMatrix<double>(file, matrix);
            }
            else if (file.EndsWith(".csv"))
            {
                DelimitedWriter.Write<double>(file, matrix, ";");
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
                return new Tuple<double, double>(matrix.At(i, j), 0.0);
            }

            return new Tuple<double, double>(0.0, 0.0);
        }

        public IEnumerable<MatrixEntry> Enumerate()
        {
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int end, n = storage.RowCount;

            for (int i = 0; i < n; i++)
            {
                end = ap[i + 1];

                for (int j = ap[i]; j < end; j++)
                {
                    yield return new MatrixEntry(i, ai[j], ax[j]);
                }
            }
        }

        public IEnumerable<MatrixEntry> EnumerateSubmatrix(int rowStart, int rowEnd, int colStart, int colEnd)
        {
            var ap = storage.RowPointers;
            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            int col, end, n = Math.Min(rowEnd, storage.RowCount);

            for (int i = rowStart; i < n; i++)
            {
                end = ap[i + 1];

                for (int j = ap[i]; j < end; j++)
                {
                    col = ai[j];

                    if (col >= colStart)
                    {
                        yield return new MatrixEntry(i, ai[j], ax[j]);
                    }

                    if (col > colEnd)
                    {
                        break;
                    }
                }
            }
        }

        public IEnumerable<MatrixEntry> EnumerateRow(int i)
        {
            int start = storage.RowPointers[i];
            int end = storage.RowPointers[i + 1];

            var ai = storage.ColumnIndices;
            var ax = storage.Values;

            for (int j = start; j < end; j++)
            {
                yield return new MatrixEntry(i, ai[j], ax[j]);
            }
        }
    }
}
