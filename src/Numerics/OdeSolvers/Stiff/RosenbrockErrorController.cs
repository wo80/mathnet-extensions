using System;

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    public class RosenbrockErrorController
    {
        bool predictive = true;

        double minscale;
        double maxscale;
        double hmax;
        double safe;

        bool rejected;

        int nrejected;
        int naccepted;

        double hnew;
        double hacc = 0.0;
        double erracc = 0.0;

        public double NextStepSize { get { return hnew; } }

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
        /// <param name="minscale">Minimum scaling for step size selection (minscale &lt= hnew/h) (default = 0.2).</param>
        /// <param name="maxscale">Maximum scaling for step size selection (hnew/h &lt= maxscale) (default = 6.0).</param>
        /// <param name="hmax">Maximal step size (default tend - t).</param>
        /// <param name="safe">Safety factor in step size prediction (default 0.9).</param>
        public RosenbrockErrorController(double minscale, double maxscale, double hmax, double safe, bool predictive)
        {
            if (safe >= 1.0 || safe <= 1e-4)
            {
                throw new ArgumentException("Curious input for safety factor.", nameof(safe));
            }

            if (minscale > 1.0)
            {
                throw new ArgumentException(" Curious input for minscale.", nameof(minscale));
            }

            if (maxscale < 1.0)
            {
                throw new ArgumentException(" Curious input for maxscale.", nameof(maxscale));
            }

            this.minscale = minscale;
            this.maxscale = maxscale;
            this.hmax = hmax;
            this.safe = safe;
            this.predictive = predictive;
        }

        public void Reset()
        {
            nrejected = 0;
            naccepted = 0;

            rejected = false;
        }
        
        public bool Success(double err, double posneg, ref double h)
        {
            // Computation of HNEW (we require 0.2 <= hnew/h <= 6.0)
            double scale = Math.Max(minscale, Math.Min(maxscale, safe / Math.Pow(err, 0.25)));

            hnew = h * scale;

            double facgus;

            // Is the error small enough?
            if (err <= 1.0)
            {
                // Step is accepted
                naccepted++;

                if (predictive)
                {
                    // Predictive controller of Gustafsson
                    if (naccepted > 1)
                    {
                        facgus = h / hacc * safe / Math.Pow(err * err / erracc, 0.25);
                        facgus = Math.Max(minscale, Math.Min(maxscale, facgus));
                        scale = Math.Min(scale, facgus);
                        hnew = h * scale;
                    }

                    hacc = h;
                    erracc = Math.Max(0.01, err);
                }

                if (Math.Abs(hnew) > hmax)
                {
                    hnew = posneg * hmax;
                }

                if (rejected)
                {
                    hnew = posneg * Math.Min(Math.Abs(hnew), Math.Abs(h));
                }

                rejected = false;
            }
            else
            {
                // Step is rejected
                rejected = true;
                
                if (naccepted > 0)
                {
                    nrejected++;
                }
            }

            h = hnew;

            return !rejected;
        }
    }
}
