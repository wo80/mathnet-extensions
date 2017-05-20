
namespace MathNet.MatrixDebuggerVisualizer.Services
{
    using MathNet.Numerics;
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Storage;
    using System;

    public class SparseStorageInfoService<T> : IStorageInfoService
        where T : struct, IEquatable<T>, IFormattable
    {
        Matrix<T> A;
        SparseCompressedRowMatrixStorage<T> S;

        IMathService<T> math;

        int currentLevel;

        public StorageInfo StorageInfo { get; private set; }

        public SparseStorageInfoService(Matrix<T> matrix, IMathService<T> mathService)
        {
            A = matrix;
            S = matrix.Storage as SparseCompressedRowMatrixStorage<T>;

            this.math = mathService;

            this.StorageInfo = new StorageInfo();

            this.StorageInfo.RowCount = S.RowCount;
            this.StorageInfo.ColumnCount = S.ColumnCount;
            this.StorageInfo.ValueCount = S.ValueCount;

            int size = Helper.SizeOf<T>();
            int i = Constants.SizeOfInt;

            this.StorageInfo.TotalBytes = i * S.RowPointers.Length +
                i * S.ColumnIndices.Length + size * S.Values.Length;
        }

        public void Update(int level)
        {
            if (currentLevel < level && currentLevel < 1)
            {
                NonzerosInfo();
            }

            if (currentLevel < level && currentLevel < 2)
            {
                DiagonalInfo();
                Bandwidth();
            }

            if (currentLevel < level && currentLevel < 3)
            {
                DiagonalDominant();
            }

            if (currentLevel < level && currentLevel < 4)
            {
                SymmetryInfo();
            }

            currentLevel = level;
        }

        /// <summary>
        /// Computes the lower, upper, maximum, and average bandwidths.
        /// </summary>
        /// <remarks>
        /// Lower bandwidth: ml = Max ( i-j | all  a(i,j) != 0 )
        /// Upper bandwidth: mu = Max ( j-i | all  a(i,j) != 0 )
        /// Maximum bandwidth: Max ( Max [ j | a(i,j) != 0 ] - Min [ j | a(i,j) != 0 ] )
        /// </remarks>
        private void Bandwidth()
        {
            int n = S.RowCount;

            var ap = S.RowPointers;
            var ai = S.ColumnIndices;

            int minr, maxr, end;

            int ml = -n; // Lower bandwidth
            int mu = -n; // Upper bandwidth

            double bndav = 0.0; // Maximum bandwidth
            int iband = 0; // Average Bandwidth

            for (int i = 0; i < n; ++i)
            {
                end = ap[i + 1];

                // Make sure that row is not empty.
                if (end > ap[i])
                {
                    minr = ai[ap[i]];
                    maxr = ai[end - 1];

                    ml = Math.Max(ml, i - minr);
                    mu = Math.Max(mu, maxr - i);

                    iband = Math.Max(iband, maxr - minr + 1);
                    bndav += (double)(maxr - minr + 1);
                }
            }

            bndav /= n;

            var info = this.StorageInfo;

            info.LowerBandwidth = Math.Max(0, ml);
            info.UpperBandwidth = Math.Max(0, mu);
            info.MaximumBandwidth = iband;
            info.AverageBandwidth = bndav;
        }

        /// <summary>
        /// Computes maximum numbers of nonzero elements per row, minimum numbers of
        /// nonzero elements per row, and numbers of zero rows.
        /// Computes the number of nonzero elements in strict lower part, strict upper
        /// part, and main diagonal.
        /// Computes average number of nonzero elements/row and standard deviation for this average.
        /// </summary>
        /// <remarks>
        /// Standard deviation will not be correct for symmetric storage.
        /// </remarks>
        private void NonzerosInfo()
        {
            int rows = S.RowCount;
            int cols = S.ColumnCount;

            var ap = S.RowPointers;
            var ai = S.ColumnIndices;

            int i, j, k, start, end, length;

            int nz = ap[rows];

            int nlow = 0; // number of nonzero elements in strict lower part
            int ndiag = 0; // number of nonzero elements in main diagonal

            int max = 0; // max length of rows
            int min = rows; // min length of rows
            int zero = 0; // number of zero rows

            double av = nz / (double)rows; // average
            double st = 0.0; // standard deviation
            double t;

            var ccount = new int[cols]; // column counts

            for (i = 0; i < rows; i++)
            {
                start = ap[i];
                end = ap[i + 1];

                for (k = start; k < end; k++)
                {
                    j = ai[k];

                    ccount[j]++;

                    if (j < i)
                    {
                        nlow++;
                    }

                    if (j == i)
                    {
                        ndiag++;
                    }
                }

                length = end - start;

                if (length <= 0)
                {
                    zero++;
                }

                max = Math.Max(max, length);
                min = Math.Min(min, length);

                t = length - av;
                st += t * t;
            }

            var info = this.StorageInfo;

            info.NonZerosLower = nlow;
            info.NonZerosUpper = nz - nlow - ndiag;
            info.NonZerosDiagonal = ndiag;

            // Update row statistics

            info.NonZerosPerRow = av;
            info.NonZerosPerRowDev = Math.Sqrt(st / rows);

            info.MaxRowLength = max;
            info.MinRowLength = min;
            info.ZeroRowsCount = zero;

            // Now process the column counts

            max = 0; // max length of columns
            min = cols; // min length of columns
            zero = 0; // number of zero columns

            av = nz / (double)cols; // average
            st = 0.0; // standard deviation

            for (i = 0; i < cols; i++)
            {
                length = ccount[i];

                if (length <= 0)
                {
                    zero++;
                }

                max = Math.Max(max, length);
                min = Math.Min(min, length);

                t = length - av;
                st += t * t;
            }

            // Update column statistics

            info.NonZerosPerColumn = av;
            info.NonZerosPerColumnDev = Math.Sqrt(st / rows);

            info.MaxColumnLength = max;
            info.MinColumnLength = min;
            info.ZeroColumnsCount = zero;
        }

        /// <summary>
        /// Computes the numbers of elements in each diagonal, the average distance
        /// of a(i,j) from diag and standard deviation for this average.
        /// </summary>
        private void DiagonalInfo()
        {
            int n = S.RowCount;

            var ap = S.RowPointers;
            var ai = S.ColumnIndices;

            int i, j, k, end, nnz = ap[n];

            double t;
            double dist = 0.0; // average distance of a(i,j) from diag.
            double std = 0.0; // standard deviation for above average.

            var count = new int[n + S.ColumnCount];

            if (nnz > 0)
            {
                for (i = 0; i < n; ++i)
                {
                    end = ap[i + 1];

                    for (k = ap[i]; k < end; k++)
                    {
                        j = ai[k];
                        dist += (double)Math.Abs(j - i);
                        count[n + j - i]++;
                    }
                }

                dist /= nnz;

                for (i = 0; i < n; i++)
                {
                    end = ap[i + 1];

                    for (k = ap[i]; k < end; k++)
                    {
                        t = dist - (double)Math.Abs(ai[k] - i);
                        std += t * t;
                    }
                }

                std = Math.Sqrt(std / (double)(nnz));
            }

            var info = this.StorageInfo;

            info.DiagonalDistanceAverage = dist;
            info.DiagonalDistanceDeviation = std;

            info.DiagonalsCount = count;

            info.Band90 = BandPart(count, 90);
            info.Band80 = BandPart(count, 80);

            /*
            int ndiag = 0;
            int[] ioff;
            float[] dcount;
            ImportantDiagonals(n, nnz, count, 1, out ndiag, out ioff, out dcount);

            Console.WriteLine("The {0} most important diagonals (offsets)", ndiag);
            Console.WriteLine("and the accumulated percentages they represent:");
            for (j = 0; j < ndiag; ++j)
            {
                Console.WriteLine(" {0,8}  {0,5:0.0} ", ioff[j], dcount[j]);
            }
            //*/
        }

        /// <summary>
        /// Computes the bandwidth of the banded matrix, which contains 'nper' percent
        /// of the original matrix.
        /// </summary>
        /// <param name="dist">int array containing the numbers of elements in the matrix
        /// with different distance of row indices and column indices.</param>
        /// <param name="nper">percentage of matrix  within the bandwidth</param>
        /// <returns>the width of the band</returns>
        private int BandPart(int[] dist, int nper)
        {
            int n = S.ColumnCount;

            var ap = S.RowPointers;
            var ai = S.ColumnIndices;

            int nnz = ap[n] - ap[0];

            if (nnz == 0)
            {
                return 0;
            }

            int iacc = dist[n];

            int band = 0;
            int j = 1;

            while (true)
            {
                iacc = iacc + dist[n + j] + dist[n - j];
                j++;

                if (iacc * 100 > nnz * nper)
                {
                    break;
                }

                band++;
            }

            return band;
        }

        /// <summary>
        /// this routine computes the most important diagonals.
        /// </summary>
        /// <param name="n">column dimension of matrix</param>
        /// <param name="nnz">number of nonzero elements of matrix</param>
        /// <param name="dist">int array containing the numbers of elements in the
        /// matrix with different distance of row indices and column indices.</param>
        /// <param name="ipar1">percentage of nonzero elements of A that a diagonal
        /// should have in order to be an important diagonal</param>
        /// <param name="ndiag">number of the most important diagonals</param>
        /// <param name="ioff">the offsets with respect to the main diagonal</param>
        /// <param name="dcount">the accumulated percentages</param>
        private void ImportantDiagonals(int n, int nnz, int[] dist, int ipar1,
            out int ndiag, out int[] ioff, out float[] dcount)
        {
            ioff = new int[20];
            dcount = new float[20];

            // Local variables
            int i, j, k, n2, ii, jmax, itot;

            // Function Body
            n2 = n + n;
            ndiag = 10;
            ndiag = Math.Min(n2, ndiag);
            itot = 0;
            ii = 0;

            // sort diagonals by decreasing order of weights.
            while (ii < ndiag)
            {
                jmax = 0;
                i = 0;
                for (k = 0; k < n2; ++k)
                {
                    j = dist[k];
                    if (j >= jmax)
                    {
                        i = k;
                        jmax = j;
                    }
                }

                // permute
                // save offsets and accumulated count if diagonal is acceptable
                // (if it has at least ipar1*nnz/100 nonzero elements)
                // quit if no more acceptable diagonals --

                if (jmax * 100 < ipar1 * nnz)
                {
                    break;
                }

                ioff[ii] = i - n;
                dist[i] = -jmax;
                itot += jmax;
                dcount[ii] = (itot * 100) / (float)nnz;

                ii++;
            }

            ndiag = ii;
        }

        /// <summary>
        /// Computes the percentage of weakly diagonally dominant rows/columns
        /// </summary>
        private void DiagonalDominant()
        {
            int n = S.RowCount;

            var ax = S.Values;
            var ap = S.RowPointers;
            var ai = S.ColumnIndices;

            var rwork = new double[n];
            var cwork = new double[n];
            var dwork = new double[n];

            int i, j, k, end;

            double a, b, max = 0.0;
            double fnorm = 0.0; // Frobenius norm

            for (i = 0; i < n; i++)
            {
                end = ap[i + 1];

                for (k = ap[i]; k < end; k++)
                {
                    j = ai[k];
                    a = math.Abs(ax[k]);

                    fnorm += math.Square(ax[k]);

                    if (j == i)
                    {
                        dwork[i] = a;
                    }
                    else
                    {
                        rwork[i] += a;
                        cwork[j] += a;
                    }

                    max = Math.Max(max, a);
                }
            }

            int ddomc = 0; // number of weakly diagonally dominant columns
            int ddomr = 0; // number of weakly diagonally dominant rows

            a = b = 0.0;

            for (i = 0; i < n; i++)
            {
                if (cwork[i] <= dwork[i])
                {
                    ddomc++;
                }

                if (rwork[i] <= dwork[i])
                {
                    ddomr++;
                }

                a = Math.Max(rwork[i] + dwork[i], a); // Infinity norm
                b = Math.Max(cwork[i] + dwork[i], b); // 1-norm
            }

            var info = this.StorageInfo;

            info.DiagonallyDominantRows = ddomr;
            info.DiagonallyDominantColumns = ddomc;

            info.MaxAbsoluteValue = max;
            info.FrobeniusNorm = Math.Sqrt(fnorm);
            info.InfNorm = a;
            info.OneNorm = b;
        }

        /// <summary>
        /// Computes the Frobenius norm of the symmetric and non-symmetric parts of A,
        /// number of matching elements in symmetry and relative symmetry match.
        /// </summary>
        private void SymmetryInfo()
        {
            int n = S.ColumnCount;

            var ax = S.Values;
            var ap = S.RowPointers;
            var ai = S.ColumnIndices;

            var B = A.Transpose().Storage as SparseCompressedRowMatrixStorage<T>;

            var bx = B.Values;
            var bp = B.RowPointers;
            var bi = B.ColumnIndices;

            int i, j1, k1, k2, j2;
            int nnz, k1max, k2max;
            double st, fnorm;

            double fas = 0.0;
            double fan = 0.0;
            int imatch = 0;

            nnz = ap[n] - ap[0];
            st = 0.0;

            for (i = 0; i < n; ++i)
            {
                k1 = ap[i];
                k2 = bp[i];
                k1max = ap[i + 1];
                k2max = bp[i + 1];

                while (k1 < k1max && k2 < k2max)
                {

                    j1 = ai[k1];
                    j2 = bi[k2];

                    if (j1 == j2)
                    {
                        fas += math.Sum2(ax[k1], bx[k2]);
                        fan += math.Difference2(ax[k1], bx[k2]);

                        st += math.Square(ax[k1]);

                        imatch++;
                    }

                    k1++;
                    k2++;

                    if (j1 < j2)
                    {
                        k2--;
                    }

                    if (j1 > j2)
                    {
                        k1--;
                    }
                }
            }

            fas *= 0.25;
            fan *= 0.25;

            if (imatch == nnz)
            {
                st = 0.0;
            }
            else
            {
                fnorm = FrobeniusNorm(B);

                st = (fnorm * fnorm - st) * 0.5;
                if (st < 0.0)
                {
                    st = 0.0;
                }
            }

            var info = this.StorageInfo;

            info.ElementsInSymmetry = imatch;

            info.FrobeniusNormSymPart = Math.Sqrt(fas + st);
            info.FrobeniusNormNonSymPart = fan < 0.0 ? 0.0 : Math.Sqrt(fan + st);
        }

        /// <summary>
        /// Computes the Frobenius norm of A.
        /// </summary>
        /// <returns>Frobenius norm of A.</returns>
        private double FrobeniusNorm(SparseCompressedRowMatrixStorage<T> A)
        {
            int n = A.RowCount;

            var ax = A.Values;
            var ap = A.RowPointers;
            var ai = A.ColumnIndices;

            int i, k;
            double fnorm = 0.0;

            for (i = 0; i < n; ++i)
            {
                for (k = ap[i]; k < ap[i + 1]; ++k)
                {
                    fnorm += math.Square(ax[k]);
                }
            }

            return Math.Sqrt(fnorm);
        }
    }
}
