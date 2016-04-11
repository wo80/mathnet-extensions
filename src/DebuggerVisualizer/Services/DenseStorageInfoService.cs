using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Storage;
using System;

namespace MathNet.MatrixDebuggerVisualizer.Services
{
    public class DenseStorageInfoService<T> : IStorageInfoService
        where T : struct, IEquatable<T>, IFormattable
    {
        Matrix<T> A;
        DenseColumnMajorMatrixStorage<T> S;

        IMathService<T> math;

        public StorageInfo StorageInfo { get; private set; }

        public DenseStorageInfoService(Matrix<T> matrix, IMathService<T> mathService)
        {
            A = matrix;
            S = matrix.Storage as DenseColumnMajorMatrixStorage<T>;

            this.math = mathService;

            this.StorageInfo = new StorageInfo();

            this.StorageInfo.RowCount = S.RowCount;
            this.StorageInfo.ColumnCount = S.ColumnCount;
            this.StorageInfo.ValueCount = S.RowCount * S.ColumnCount;

            int size = Helper.SizeOf<T>();

            this.StorageInfo.TotalBytes = size * S.RowCount * S.ColumnCount;
        }

        public void Update(int level)
        {
        }
    }
}
