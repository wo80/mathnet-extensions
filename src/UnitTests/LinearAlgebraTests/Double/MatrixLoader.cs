using MathNet.Numerics.LinearAlgebra.Double;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace MathNet.Numerics.Extensions.UnitTests.LinearAlgebraTests.Double
{
    /// <summary>
    /// Base class for matrix tests.
    /// </summary>
    public static class MatrixLoader
    {
        /// <summary>
        /// Gets or sets test matrices instances to use.
        /// </summary>
        private static Dictionary<int, SparseMatrix> TestMatricesA { get; set; }

        /// <summary>
        /// Gets or sets test matrices instances to use.
        /// </summary>
        private static Dictionary<int, SparseMatrix> TestMatricesB { get; set; }

        public static SparseMatrix A(int size)
        {
            if (TestMatricesA == null)
            {
                SetupMatrices();
            }

            return (SparseMatrix)TestMatricesA[size].Clone();
        }

        public static SparseMatrix B(int size, bool transpose)
        {
            if (TestMatricesB == null)
            {
                SetupMatrices();
            }

            if (transpose)
            {
                if (size < 10 || size > 99)
                {
                    throw new Exception("Matrix size must be in range 10-99.");
                }

                size = 10 * (size % 10) + size / 10;
            }

            return (SparseMatrix)TestMatricesB[size].Clone();
        }

        /// <summary>
        /// Setup test matrices.
        /// </summary>
        private static void SetupMatrices()
        {
            var A55 = CreateSparse.Random(5, 5, 0.5);
            var A57 = CreateSparse.Random(5, 7, 0.5);
            var A75 = CreateSparse.Random(7, 5, 0.5);

            var B55 = CreateSparse.Random(5, 5, 0.5);
            var B57 = CreateSparse.Random(5, 7, 0.5);
            var B75 = CreateSparse.Random(7, 5, 0.5);

            // Set all diagonal entries to 1.0.
            for (int i = 0; i < 5; i++)
            {
                A55[i, i] = 1.0;
                A57[i, i] = 1.0;
                A75[i, i] = 1.0;
                B55[i, i] = 1.0;
                B57[i, i] = 1.0;
                B75[i, i] = 1.0;
            }

            // Ensure there are no zero rows/columns in A.
            A57[4, 5] = 0.1;
            A57[4, 6] = 0.2;
            A75[5, 4] = 0.1;
            A75[6, 4] = 0.2;

            TestMatricesA = new Dictionary<int, SparseMatrix>();
            TestMatricesB = new Dictionary<int, SparseMatrix>();

            TestMatricesA.Add(55, A55);
            TestMatricesA.Add(57, A57);
            TestMatricesA.Add(75, A75);
            TestMatricesB.Add(55, B55);
            TestMatricesB.Add(57, B57);
            TestMatricesB.Add(75, B75);
        }
    }
}
