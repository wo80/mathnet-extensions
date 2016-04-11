#if !MIT_ONLY

namespace MathNet.Numerics.LinearAlgebra.Complex.Solvers
{
    using MathNet.Numerics.LinearAlgebra;
    using MathNet.Numerics.LinearAlgebra.Solvers;
    using System;
    using System.Numerics;

    /// <summary>
    /// Preconditioned conjugate gradient solver.
    /// </summary>
    public sealed class PCG : IIterativeSolver<Complex>
    {
        private const double TINY = 1.0e-292;

        /// <summary>
        /// Solves the matrix equation Ax = b, where A is the coefficient matrix, b is the
        /// solution vector and x is the unknown vector.
        /// </summary>
        /// <param name="matrix">The coefficient <see cref="Matrix"/>, <c>A</c>.</param>
        /// <param name="input">The solution <see cref="Vector"/>, <c>b</c>.</param>
        /// <param name="result">The result <see cref="Vector"/>, <c>x</c>.</param>
        /// <param name="iterator">The iterator to use to control when to stop iterating.</param>
        /// <param name="preconditioner">The preconditioner to use for approximations.</param>
        public void Solve(Matrix<Complex> matrix, Vector<Complex> input, Vector<Complex> result,
            Iterator<Complex> iterator, IPreconditioner<Complex> preconditioner)
        {
            var A = (SparseMatrix)matrix;
            var M = preconditioner;

            var b = (DenseVector)input;
            var x = (DenseVector)result;

            double atolf = 0.0;
            double rtol_1 = 0.0;
            int recompute_residual_p = 0;

            var p = new DenseVector(b.Count);
            var s = new DenseVector(b.Count);
            var r = new DenseVector(b.Count);

            Complex alpha, beta;
            Complex gamma, gamma_old;

            Complex bi_prod;

            Complex sdotp;
            bool recompute_true_residual = false;

            int i = 0;

            M.Initialize(A);

            // Start pcg solve

            // bi_prod = <C*b,b>
            //VectorHelper.Clear(p);
            M.Approximate(input, p);
            bi_prod = VectorHelper.DotProduct(p, b);

            if (bi_prod.Real > 0.0)
            {
                if (atolf > 0)  // mixed relative and absolute tolerance
                    bi_prod += atolf;
            }
            else    // bi_prod==0.0: the rhs vector b is zero
            {
                // Set x equal to zero and return
                VectorHelper.Copy(b, x);
                return;
            }

            // r = b - Ax
            VectorHelper.Copy(b, r);
            A.Multiply(-1.0, result, 1.0, r);

            // p = C*r
            //VectorHelper.Clear(p);
            M.Approximate(r, p);

            // gamma = <r,p>
            gamma = VectorHelper.DotProduct(r, p);

            while (iterator.DetermineStatus(i, x, b, r) == IterationStatus.Continue)
            {
                // the core CG calculations...
                i++;

                // At user request, periodically recompute the residual from the formula
                // r = b - A x (instead of using the recursive definition). Note that this
                // is potentially expensive and can lead to degraded convergence (since it
                // essentially a "restarted CG").
                recompute_true_residual = (recompute_residual_p > 0) && !((i % recompute_residual_p) == 0);

                // s = A*p
                A.Multiply(1.0, p, 0.0, s);

                // alpha = gamma / <s,p>
                sdotp = VectorHelper.DotProduct(s, p);
                if (sdotp == 0.0)
                {
                    throw new NumericalBreakdownException();
                }
                alpha = gamma / sdotp;

                gamma_old = gamma;

                // x = x + alpha*p
                VectorHelper.Add(alpha, p, x, x);

                // r = r - alpha*s
                if (recompute_true_residual)
                {
                    //Recomputing the residual...
                    VectorHelper.Copy(b, r);
                    A.Multiply(-1.0, result, 1.0, r);
                }
                else
                {
                    VectorHelper.Add(-alpha, s, r, r);
                }

                // s = C*r
                VectorHelper.Clear(s);
                M.Approximate(r, s);

                // gamma = <r,s>
                gamma = VectorHelper.DotProduct(r, s);

                // residual-based stopping criteria: ||r_new-r_old||_C < rtol ||b||_C
                if (rtol_1 > 0)
                {
                    // use that ||r_new-r_old||_C^2 = (r_new ,C r_new) + (r_old, C r_old)
                    if (Complex.Abs((gamma + gamma_old) / bi_prod) < rtol_1 * rtol_1)
                    {
                        break;
                    }
                }

                // ... gamma should be >=0.  IEEE subnormal numbers are < 2**(-1022)=2.2e-308
                // (and >= 2**(-1074)=4.9e-324).  So a gamma this small means we're getting
                // dangerously close to subnormal or zero numbers (usually if gamma is small,
                // so will be other variables).  Thus further calculations risk a crash.
                // Such small gamma generally means no hope of progress anyway.
                if (Complex.Abs(gamma) < TINY)
                {
                    throw new NumericalBreakdownException();
                }

                // beta = gamma / gamma_old
                beta = gamma / gamma_old;

                // p = s + beta p
                if (recompute_true_residual)
                {
                    VectorHelper.Copy(s, p);
                }
                else
                {
                    p.Scale(beta);
                    VectorHelper.Add(1.0, s, p, p);
                }
            }
        }
    }
}

#endif