using System;

namespace MathNet.Numerics.OdeSolvers
{
    public class ErrorControllerRos
    {
        //double maxscale = 5.0;
        //double minscale = 0.16666666666666666;
        //double safe = 0.9;

        bool pred = true;

        double minscale;
        double maxscale;
        double hmaxn;
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
        /// <param name="minscale">Minimum scaling for step size selection (minscale &lt= hnew/hold).</param>
        /// <param name="maxscale">Maximum scaling for step size selection (hnew/hold &lt= maxscale).</param>
        /// <param name="hmaxn">Maximal step size (default tend - t).</param>
        /// <param name="safe">Safety factor in step size prediction (default 0.9).</param>
        public ErrorControllerRos(double minscale, double maxscale, double hmaxn, double safe)
        {
            if (safe >= 1.0 || safe <= 1e-4)
            {
                throw new ArgumentException("Curious input for safety factor.", nameof(safe));
            }
            
            this.minscale = minscale;
            this.maxscale = maxscale;
            this.hmaxn = hmaxn;
            this.safe = safe;
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
            double fac = Math.Max(minscale, Math.Min(maxscale, Math.Pow(err, 0.25) / safe));

            hnew = h / fac;

            double facgus;

            // Is the error small enough?
            if (err <= 1.0)
            {
                // Step is accepted
                naccepted++;

                if (pred)
                {
                    // Predictive controller of Gustafsson
                    if (naccepted > 1)
                    {
                        facgus = hacc / h * Math.Pow(err * err / erracc, 0.25) / safe;
                        facgus = Math.Max(minscale, Math.Min(maxscale, facgus));
                        fac = Math.Max(fac, facgus);
                        hnew = h / fac;
                    }

                    hacc = h;
                    erracc = Math.Max(.01, err);
                }

                if (Math.Abs(hnew) > hmaxn)
                {
                    hnew = posneg * hmaxn;
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
