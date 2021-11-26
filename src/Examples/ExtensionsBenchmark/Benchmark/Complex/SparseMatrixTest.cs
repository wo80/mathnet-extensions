
namespace ExtensionsBenchmark.Benchmark.Complex
{
    using MathNet.Numerics;
    using MathNet.Numerics.LinearAlgebra.Complex;
    using System;
    using System.Collections.Generic;
    using System.Numerics;

    using DVector = MathNet.Numerics.LinearAlgebra.Vector<double>;

    class SparseMatrixTest : IBenchmarkCollection
    {
        private const double EPS = 1e-12;
        private const int RAND_SEED = 230183;

        private static Random rand = new Random(RAND_SEED);

        List<IBenchmark> items;

        public SparseMatrixTest()
        {
            items = new List<IBenchmark>();

            items.Add(new Test_MULT_MV());
            items.Add(new Test_TMULT1());
            items.Add(new Test_TMULT2());
            items.Add(new Test_ADD());
            items.Add(new Test_KRON());
            //items.Add(new Test_MULT_MM());
            //items.Add(new Test_TRANS());
            items.Add(new Test_LOWER());
            items.Add(new Test_UPPER());
            items.Add(new Test_ROW_GET());
            items.Add(new Test_ROWS_GET());
            items.Add(new Test_COLS_GET());
            items.Add(new Test_ROWS_SET());
            //items.Add(new Test_SUBMAT());
            items.Add(new Test_PERM_R());
            items.Add(new Test_PERM_C());
            //items.Add(new Test_CLEAR_R());
            //items.Add(new Test_CLEAR_C());
            items.Add(new Test_NORM_R());
            items.Add(new Test_NORM_C());
            items.Add(new Test_SCAL_R());
            items.Add(new Test_SCAL_C());
            items.Add(new Test_NLZE_R());
            items.Add(new Test_NLZE_C());
            items.Add(new Test_DIAG_S());
            items.Add(new Test_DIAG_M());
        }

        public List<IBenchmark> Items
        {
            get { return items; }
        }

        public string Name
        {
            get { return "SparseMatrix (Double)"; }
        }

        public void Initialize()
        {
            foreach (BenchmarkBase item in items)
            {
                item.Initialize();
            }
        }

        #region Implementation

        class Test_MULT_MV : BenchmarkBase
        {
            public override string Name
            {
                get { return "Multiply (Matrix-Vector)"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var x = Create.Vector(cols, 1.0);
                var y = Create.Array(cols, 1.0);
                var z = (DenseVector)x.Clone();

                Util.Tic();

                for (int i = 0; i < repeat; i++)
                {
                    A.Multiply(2.0, Util.GetData(x), 0.1, y);
                }

                result.Time2 = Util.Toc();


                Util.Tic();

                var temp = x.Clone();

                for (int i = 0; i < repeat; i++)
                {
                    //z = (DenseVector)(2.0 * A * x + 0.1 * z);
                    temp.Multiply(2.0, x);
                    A.Multiply(x, x);
                    z.Multiply(0.1, z);
                    z.Add(x, z);
                }

                result.Time1 = Util.Toc();

                result.Success = Helper.Equals(z, y, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(200, 200, 10, 100);
            }
        }

        class Test_TMULT1 : BenchmarkBase
        {
            public override string Name
            {
                get { return "Transpose Multiply 1 (Matrix-Vector)"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var x = Create.Vector(cols, 1.0);
                var y = Create.Array(cols, 1.0);
                var z = (DenseVector)x.Clone();

                Util.Tic();

                for (int i = 0; i < repeat; i++)
                {
                    A.TransposeMultiply(2.0, Util.GetData(x), 0.1, y);
                }

                result.Time2 = Util.Toc();


                Util.Tic();

                var temp = x.Clone();

                for (int i = 0; i < repeat; i++)
                {
                    // z = 2.0 * A^t * x + 0.1 * z;
                    temp.Multiply(2.0, x);

                    // Completely unusable for large sparse matrices.
                    A.TransposeThisAndMultiply(x, x);

                    z.Multiply(0.1, z);
                    z.Add(x, z);
                }

                result.Time1 = Util.Toc();

                result.Success = Helper.Equals(z, y, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(200, 200, 10);
            }
        }

        // Difference to test above: repeated 100 times.
        class Test_TMULT2 : BenchmarkBase
        {
            public override string Name
            {
                get { return "Transpose Multiply 2 (Matrix-Vector)"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var x = Create.Vector(cols, 1.0);
                var y = Create.Array(cols, 1.0);
                var z = (DenseVector)x.Clone();

                Util.Tic();

                for (int i = 0; i < repeat; i++)
                {
                    A.TransposeMultiply(2.0, Util.GetData(x), 0.1, y);
                }

                result.Time2 = Util.Toc();


                Util.Tic();

                var temp = x.Clone();

                // We really don't like to compute/allocate the transpose, but ...
                A = (SparseMatrix)A.Transpose();

                for (int i = 0; i < repeat; i++)
                {
                    // z = 2.0 * A^t * x + 0.1 * z;
                    temp.Multiply(2.0, x);
                    A.Multiply(x, x);
                    z.Multiply(0.1, z);
                    z.Add(x, z);
                }

                result.Time1 = Util.Toc();

                result.Success = Helper.Equals(z, y, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(200, 200, 10, 100);
            }
        }

        class Test_ADD : BenchmarkBase
        {
            public override string Name
            {
                get { return "Add"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                var B = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                Util.Tic();
                var C = A.FastAdd(B);
                result.Time2 = Util.Toc();

                Util.Tic();
                //var D = A.Add(B);
                var D = A + B;
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(C, D, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(1000, 250, 10);
            }
        }

        class Test_KRON : BenchmarkBase
        {
            public override string Name
            {
                get { return "Kronecker Product"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                var S = CreateSparse.RandomSymmetric(4, 0.5);

                var result = new BenchmarkResult();

                Util.Tic();
                var B = S.FastKroneckerProduct(A);
                result.Time2 = Util.Toc();

                Util.Tic();
                var C = S.KroneckerProduct(A);
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(50, 50, 10);
            }
        }

        /*
        class Test_MULT_MM : BenchmarkBase
        {
            public override string Name
            {
                get { return "Multiply"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                var B = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                Util.Tic();
                var C = A.FastMultiply(B);
                result.Time2 = Util.Toc();

                Util.Tic();
                var D = A.Multiply(B);
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(C, D, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(100, 100, 10);
            }
        }

        class Test_TRANS : BenchmarkBase
        {
            public override string Name
            {
                get { return "Transpose"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                Util.Tic();
                var B = A.FastTranspose();
                result.Time2 = Util.Toc();

                Util.Tic();
                var C = A.ConjugateTranspose();
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }
        //*/

        class Test_LOWER : BenchmarkBase
        {
            public override string Name
            {
                get { return "Lower Triangle"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                var B = (SparseMatrix)A.Clone();

                var result = new BenchmarkResult();

                Util.Tic();
                B.FastLowerTriangle();
                result.Time2 = Util.Toc();

                Util.Tic();
                var C = (SparseMatrix)A.LowerTriangle();
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                // Strict lower.
                B = (SparseMatrix)A.Clone();

                Util.Tic();
                B.FastLowerTriangle(true);
                result.Time2 += Util.Toc();

                Util.Tic();
                C = (SparseMatrix)A.StrictlyLowerTriangle();
                result.Time1 += Util.Toc();

                result.Success &= MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }

        class Test_UPPER : BenchmarkBase
        {
            public override string Name
            {
                get { return "Upper Triangle"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                var B = (SparseMatrix)A.Clone();

                var result = new BenchmarkResult();

                Util.Tic();
                B.FastUpperTriangle();
                result.Time2 = Util.Toc();

                Util.Tic();
                var C = (SparseMatrix)A.UpperTriangle();
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                // Strict upper.
                B = (SparseMatrix)A.Clone();

                Util.Tic();
                B.FastUpperTriangle(true);
                result.Time2 += Util.Toc();

                Util.Tic();
                C = (SparseMatrix)A.StrictlyUpperTriangle();
                result.Time1 += Util.Toc();

                result.Success &= MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(100, 200, 10);
            }
        }

        class Test_ROW_GET : BenchmarkBase
        {
            public override string Name
            {
                get { return "Get Row"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var rowIndices = new int[repeat];

                for (int i = 0; i < repeat; i++)
                {
                    rowIndices[i] = rand.Next(0, rows - 1);
                }

                var x = new DenseVector(cols);
                var y = new DenseVector(cols);

                Util.Tic();
                for (int i = 0; i < repeat; i++)
                {
                    A.FastGetRow(rowIndices[i], x);
                }
                result.Time2 = Util.Toc();


                Util.Tic();
                for (int i = 0; i < repeat; i++)
                {
                    A.Row(rowIndices[i], y);
                }
                result.Time1 = Util.Toc();

                result.Success = Helper.Equals(x, (Complex[])y, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(600, 200, 10, 400);
            }
        }

        class Test_ROWS_GET : BenchmarkBase
        {
            public override string Name
            {
                get { return "Extract Rows"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var rowIndices = new int[repeat];

                for (int i = 0; i < repeat; i++)
                {
                    rowIndices[i] = rand.Next(0, rows - 1);
                }

                Util.Tic();
                var B = A.FastGetRows(rowIndices);
                result.Time2 = Util.Toc();


                Util.Tic();
                var C = new SparseMatrix(repeat, cols);
                for (int i = 0; i < repeat; i++)
                {
                    C.SetRow(i, A.Row(rowIndices[i]));
                }
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(600, 200, 10, 200);
            }
        }

        class Test_COLS_GET : BenchmarkBase
        {
            public override string Name
            {
                get { return "Extract Columns"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var columnIndices = new int[repeat];

                for (int i = 0; i < repeat; i++)
                {
                    columnIndices[i] = rand.Next(0, cols - 1);
                }

                Util.Tic();
                var B = A.FastGetColumns(columnIndices);
                result.Time2 = Util.Toc();


                Util.Tic();
                var C = new SparseMatrix(rows, repeat);
                for (int i = 0; i < repeat; i++)
                {
                    C.SetColumn(i, A.Column(columnIndices[i]));
                }
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(400, 200, 10, 100);
            }
        }

        class Test_ROWS_SET : BenchmarkBase
        {
            public override string Name
            {
                get { return "Replace Rows"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                var B = (SparseMatrix)A.Clone();

                var S = Create.SparseMatrix(repeat, cols, 2 * density, false);

                var result = new BenchmarkResult();

                var rowIndices = new int[repeat];

                for (int i = 0; i < repeat; i++)
                {
                    rowIndices[i] = rand.Next(0, rows - 1);
                }

                Util.Tic();
                A.FastSetRows(rowIndices, S);
                result.Time2 = Util.Toc();


                Util.Tic();
                for (int i = 0; i < repeat; i++)
                {
                    B.SetRow(rowIndices[i], S.Row(i));
                }
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(A, B, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 100, 10, 50);
            }
        }

        /*

        class Test_SUBMAT : BenchmarkBase
        {
            public override string Name
            {
                get { return "Sub-Matrix"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                int r2 = rows / 2;
                int c2 = cols / 2;

                Util.Tic();
                var B11 = A.FastSubMatrix(0, r2, 0, c2);
                var B12 = A.FastSubMatrix(r2, r2, 0, c2);
                var B21 = A.FastSubMatrix(r2, r2, 0, c2);
                var B22 = A.FastSubMatrix(r2, r2, c2, c2);
                result.Time2 = Util.Toc();

                Util.Tic();
                var C11 = (SparseMatrix)A.SubMatrix(0, r2, 0, c2);
                var C12 = (SparseMatrix)A.SubMatrix(r2, r2, 0, c2);
                var C21 = (SparseMatrix)A.SubMatrix(r2, r2, 0, c2);
                var C22 = (SparseMatrix)A.SubMatrix(r2, r2, c2, c2);
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B11, C11, EPS);
                result.Success &= MatrixComparer.Equals(B12, C12, EPS);
                result.Success &= MatrixComparer.Equals(B21, C21, EPS);
                result.Success &= MatrixComparer.Equals(B22, C22, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }

        //*/

        class Test_PERM_R : BenchmarkBase
        {
            public override string Name
            {
                get { return "Permute Rows"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var p = Util.Permutation(A.RowCount, rand);

                var B = (SparseMatrix)A.Clone();
                var C = (SparseMatrix)A.Clone();

                var result = new BenchmarkResult();

                Util.Tic();
                C.FastPermuteRows(p);
                result.Time2 = Util.Toc();

                Util.Tic();
                B.PermuteRows(new Permutation(p));
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(50, 50, 10);
            }
        }

        class Test_PERM_C : BenchmarkBase
        {
            public override string Name
            {
                get { return "Permute Columns"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var p = Util.Permutation(A.RowCount, rand);

                var B = (SparseMatrix)A.Clone();
                var C = (SparseMatrix)A.Clone();

                var result = new BenchmarkResult();

                Util.Tic();
                C.FastPermuteColumns(p);
                result.Time2 = Util.Toc();

                Util.Tic();
                B.PermuteColumns(new Permutation(p));
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(50, 50, 10);
            }
        }

        /*
        class Test_CLEAR_R : BenchmarkBase
        {
            public override string Name
            {
                get { return "Clear Rows"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);
                
                var indices = new int[A.RowCount / 5];

                for (int i = 0; i < indices.Length; i++)
                {
                    // Clear every 5th row.
                    indices[i] = 5 * i;
                }

                var result = new BenchmarkResult();

                var B = (SparseMatrix)A.Clone();
                var C = (SparseMatrix)A.Clone();

                Util.Tic();
                B.FastClearRows(indices);
                result.Time2 = Util.Toc();

                Util.Tic();
                C.ClearRows(indices);
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }

        class Test_CLEAR_C : BenchmarkBase
        {
            public override string Name
            {
                get { return "Clear Columns"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var indices = new int[A.ColumnCount / 5];

                for (int i = 0; i < indices.Length; i++)
                {
                    // Clear every 5th column.
                    indices[i] = 5 * i;
                }

                var result = new BenchmarkResult();

                var B = (SparseMatrix)A.Clone();
                var C = (SparseMatrix)A.Clone();

                Util.Tic();
                B.FastClearColumns(indices);
                result.Time2 = Util.Toc();

                Util.Tic();
                C.ClearColumns(indices);
                result.Time1 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }
        //*/

        class Test_NORM_R : BenchmarkBase
        {
            public override string Name
            {
                get { return "Row Norms"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                DVector b;
                var a = new double[A.RowCount];

                var result = new BenchmarkResult();

                Util.Tic();
                b = A.RowNorms(1);
                result.Time1 += Util.Toc();

                Util.Tic();
                A.FastRowNorms(1, a);
                result.Time2 += Util.Toc();

                result.Success = Helper.Equals(b, a, 0.0);

                Util.Tic();
                b = A.RowNorms(2);
                result.Time1 += Util.Toc();

                Util.Tic();
                A.FastRowNorms(2, a);
                result.Time2 += Util.Toc();

                result.Success &= Helper.Equals(b, a, 0.0);

                Util.Tic();
                b = A.RowNorms(double.PositiveInfinity);
                result.Time1 += Util.Toc();

                Util.Tic();
                A.FastRowNorms(0, a);
                result.Time2 += Util.Toc();

                result.Success &= Helper.Equals(b, a, 0.0);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }

        class Test_NORM_C : BenchmarkBase
        {
            public override string Name
            {
                get { return "Column Norms"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                DVector b;
                var a = new double[A.ColumnCount];

                var result = new BenchmarkResult();

                Util.Tic();
                b = A.ColumnNorms(1);
                result.Time1 += Util.Toc();

                Util.Tic();
                A.FastColumnNorms(1, a);
                result.Time2 += Util.Toc();

                result.Success = Helper.Equals(b, a, 0.0);

                Util.Tic();
                b = A.ColumnNorms(2);
                result.Time1 += Util.Toc();

                Util.Tic();
                A.FastColumnNorms(2, a);
                result.Time2 += Util.Toc();

                result.Success &= Helper.Equals(b, a, 0.0);

                Util.Tic();
                b = A.ColumnNorms(double.PositiveInfinity);
                result.Time1 += Util.Toc();

                Util.Tic();
                A.FastColumnNorms(0, a);
                result.Time2 += Util.Toc();

                result.Success &= Helper.Equals(b, a, 0.0);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(500, 250, 10);
            }
        }

        class Test_SCAL_R : BenchmarkBase
        {
            public override string Name
            {
                get { return "Scale Rows"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var temp = A.RowNorms(1);
                var a = DenseVector.Create(rows, i => temp[i]);
                var b = Util.GetData(a);
                var D = SparseMatrix.OfDiagonalVector(a);

                var result = new BenchmarkResult();

                Util.Tic();
                var B = D.Multiply(A);
                result.Time1 = Util.Toc();

                Util.Tic();
                A.FastScaleRows(b, A);
                result.Time2 = Util.Toc();

                result.Success = MatrixComparer.Equals(A, B, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(1000, 250, 10);
            }
        }

        class Test_SCAL_C : BenchmarkBase
        {
            public override string Name
            {
                get { return "Scale Columns"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var temp = A.RowNorms(1);
                var a = DenseVector.Create(cols, i => temp[i]);
                var b = Util.GetData(a);
                var D = SparseMatrix.OfDiagonalVector(a);

                var result = new BenchmarkResult();

                Util.Tic();
                var B = A.Multiply(D);
                result.Time1 = Util.Toc();

                Util.Tic();
                A.FastScaleColumns(b, A);
                result.Time2 = Util.Toc();

                result.Success = MatrixComparer.Equals(A, B, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(1000, 250, 10);
            }
        }

        class Test_NLZE_R : BenchmarkBase
        {
            public override string Name
            {
                get { return "Normalize Rows"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                Util.Tic();
                var B = A.NormalizeRows(1);
                result.Time1 = Util.Toc();

                Util.Tic();
                A.FastNormalizeRows(1);
                result.Time2 = Util.Toc();

                result.Success = MatrixComparer.Equals(A, B, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(1000, 250, 10);
            }
        }

        class Test_NLZE_C : BenchmarkBase
        {
            public override string Name
            {
                get { return "Normalize Columns"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                Util.Tic();
                var B = A.NormalizeColumns(1);
                result.Time1 = Util.Toc();

                Util.Tic();
                A.FastNormalizeColumns(1);
                result.Time2 = Util.Toc();

                result.Success = MatrixComparer.Equals(A, B, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(1000, 250, 10);
            }
        }

        class Test_DIAG_S : BenchmarkBase
        {
            public override string Name
            {
                get { return "Add Diagonal (Scalar)"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var C = (SparseMatrix)A.Clone();

                Util.Tic();
                var a = DenseVector.Create(rows, 1.0);
                var D = DiagonalMatrix.OfDiagonal(rows, cols, a);
                var B = A.Add(D);
                result.Time1 = Util.Toc();

                Util.Tic();
                C.FastAddDiagonalScalar(1.0);
                result.Time2 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(200, 200, 10);
            }
        }

        class Test_DIAG_M : BenchmarkBase
        {
            public override string Name
            {
                get { return "Add Diagonal (Matrix)"; }
            }

            public override BenchmarkResult Run(BenchmarkSetup config)
            {
                int rows = config.RowCount;
                int cols = config.ColumnCount;
                int repeat = config.Repeat;
                double density = config.Density;
                bool symmetric = config.Symmetric;

                var A = Create.SparseMatrix(rows, cols, density, symmetric);

                var result = new BenchmarkResult();

                var a = DenseVector.Create(rows, (i) => 1.0 + i);
                var D = DiagonalMatrix.OfDiagonal(rows, cols, a);
                var C = (SparseMatrix)A.Clone();

                Util.Tic();
                var B = A.Add(D);
                result.Time1 = Util.Toc();

                Util.Tic();
                C.FastAddDiagonalMatrix(Util.GetData(a), C);
                result.Time2 = Util.Toc();

                result.Success = MatrixComparer.Equals(B, C, EPS);

                return result;
            }

            public override List<BenchmarkSetup> GetConfigurations()
            {
                return SetupTests(200, 200, 10);
            }
        }

        #endregion
    }
}
