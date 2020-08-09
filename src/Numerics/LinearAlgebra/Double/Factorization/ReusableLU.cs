
// <copyright file="DenseLU.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
//
// Copyright (c) 2009-2013 Math.NET
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using System;
using MathNet.Numerics.Providers.LinearAlgebra;

namespace MathNet.Numerics.LinearAlgebra.Double.Factorization
{
    /// <summary>
    /// <para>A class which encapsulates the functionality of an LU factorization.</para>
    /// <para>For a matrix A, the LU factorization is a pair of lower triangular matrix L and
    /// upper triangular matrix U so that A = L*U.</para>
    /// </summary>
    /// <remarks>
    /// In contrast to the Math.NET LU implementation, <see cref="ReusableLU"/> will compute the
    /// factorization when requested and re-use allocated storage from previous factorizations.
    /// </remarks>
    internal sealed class ReusableLU
    {
        private readonly DenseMatrix factors;
        private readonly int[] pivots;

        /// <summary>
        /// Initializes a new instance of the <see cref="ReusableLU"/> class.
        /// </summary>
        /// <param name="n"></param>
        public ReusableLU(int n)
        {
            factors = new DenseMatrix(n);
            pivots = new int[n];
        }

        /// <summary>
        /// Gets the determinant of the matrix for which the LU factorization was computed.
        /// </summary>
        public double Determinant
        {
            get
            {
                int n = factors.RowCount;
                var det = 1.0;
                for (var i = 0; i < n; i++)
                {
                    if (pivots[i] != i)
                    {
                        det *= -factors.Values[i * (n + 1)];
                    }
                    else
                    {
                        det *= factors.Values[i * (n + 1)];
                    }
                }

                return det;
            }
        }

        /// <summary>
        /// Compute the LU factorization.
        /// </summary>
        /// <param name="matrix">The matrix to factor.</param>
        /// <exception cref="ArgumentNullException">If <paramref name="matrix"/> is <c>null</c>.</exception>
        /// <exception cref="ArgumentException">If <paramref name="matrix"/> is not a square matrix.</exception>
        public void Compute(DenseMatrix matrix)
        {
            if (matrix == null)
            {
                throw new ArgumentNullException("matrix");
            }

            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare);
            }

            if (factors.RowCount != matrix.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixDimensions);
            }

            matrix.CopyTo(factors);

            LinearAlgebraControl.Provider.LUFactor(factors.Values, factors.RowCount, pivots);
        }

        /// <summary>
        /// Solves a system of linear equations, <c>Ax = b</c>, with A LU factorized.
        /// </summary>
        /// <param name="input">The right hand side vector, <c>b</c>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <c>x</c>.</param>
        public void Solve(double[] input, double[] result)
        {
            // Check for proper arguments.
            if (input == null)
            {
                throw new ArgumentNullException("input");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            // Check for proper dimensions.
            if (input.Length != result.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            if (input.Length != factors.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixDimensions);
            }

            // Copy the contents of input to result.
            Buffer.BlockCopy(input, 0, result, 0, input.Length * Constants.SizeOfDouble);

            // LU solve by overwriting result.
            LinearAlgebraControl.Provider.LUSolveFactored(1, factors.Values, factors.RowCount, pivots, result);
        }
    }
}
