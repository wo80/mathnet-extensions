
namespace MathNet.Numerics.OdeSolvers.RK
{
    using System;

    public class ErrorController : IErrorController
    {
        double minscale;
        double maxscale;
        double hmax;
        double safe;
        double alpha;
        double beta;

        double hnext;
        double errold;

        bool rejected;

        int nrejected;
        int naccepted;

        public double NextStepSize { get { return hnext; } }

        /// <summary>
        /// Gets the number of rejected steps (due to error test).
        /// </summary>
        /// <remarks>
        /// Step rejections in the first step are not counted.
        /// </remarks>
        public int Rejected { get { return nrejected; } }

        /// <summary>
        /// Gets the number of accepted steps.
        /// </summary>
        public int Accepted { get { return naccepted; } }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="order">The order of the Runge-Kutta integrator.</param>
        /// <param name="minscale">Minimum scaling for step size selection (minscale &lt= hnew/hold).</param>
        /// <param name="maxscale">Maximum scaling for step size selection (hnew/hold &lt= maxscale).</param>
        /// <param name="hmax">Maximal step size (default tend - t).</param>
        /// <param name="safe">Safety factor in step size prediction (default 0.9).</param>
        /// <param name="beta">Stabilized step size control (default 0.0).</param>
        /// <remarks>
        /// The parameter "beta" is for stabilized step size control (see "Solving Ordinary Differential Equations",
        /// Hairer & Wanner, section IV.2, PI step size control).
        /// 
        /// DOPRI5: larger values of beta(&lt;= 0.1) make the step size control more stable. DOPRI5 needs
        ///    a larger beta than Higham & Hall (default 0.04).
        ///    
        /// DOPRI853: Positive values of beta (&lt;= 0.04) make the step size control more stable (default 0.0).
        /// </remarks>
        public ErrorController(int order, double minscale, double maxscale, double hmax, double safe, double beta)
        {
            if (safe >= 1.0 || safe <= 1e-4)
            {
                throw new ArgumentException("Curious input for safety factor.", nameof(safe));
            }

            if (beta > 0.2)
            {
                throw new ArgumentException("Curious input for beta.", nameof(beta));
            }

            if (beta < 0.0)
            {
                beta = 0.0;
            }

            this.errold = 1.0e-4;

            this.minscale = minscale;
            this.maxscale = maxscale;
            this.hmax = hmax;
            this.safe = safe;
            this.alpha = 1.0 / order - beta * 0.75;
            this.beta = beta;
        }

        public void Reset()
        {
            nrejected = 0;
            naccepted = 0;

            rejected = false;
        }

        public bool Success(double err, double posneg, ref double h)
        {
            // Computation of HNEW
            if (err <= 1.0)
            {
                // LUND-stabilization
                double scale = safe * Math.Pow(errold, beta) / Math.Pow(err, alpha);

                // We require  minscale <= hnext/h <= maxscale
                scale = Math.Min(maxscale, Math.Max(minscale, scale));

                hnext = scale * h;

                if (Math.Abs(hnext) > hmax)
                {
                    hnext = posneg * hmax;
                }

                if (rejected)
                {
                    hnext = posneg * Math.Min(Math.Abs(hnext), Math.Abs(h));
                }
                else
                {
                    hnext = h * scale;
                }

                errold = Math.Max(err, 1e-4);

                rejected = false;

                naccepted++;
            }
            else
            {
                // Step is rejected
                hnext = h * Math.Max(minscale, safe / Math.Pow(err, alpha));

                rejected = true;

                if (naccepted > 0)
                {
                    nrejected++;
                }
            }

            h = hnext;

            return !rejected;
        }
    }
}
