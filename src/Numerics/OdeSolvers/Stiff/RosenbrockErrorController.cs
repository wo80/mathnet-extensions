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

        // TODO: minscale <> maxscale

        /*
         * 
         *    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
         *              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
         *                 WORK(3) &lt;= HNEW/HOLD &lt;= WORK(4)
         *              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=6.D0
         *              

            //  FAC1,FAC2     Parameters for step size selection
            if (work[2] == 0.0)
            {
                maxscale = 5.0;
            }
            else
            {
                maxscale = 1.0 / work[2];
            }
            if (work[3] == 0.0)
            {
                minscale = 0.16666666666666666;
            }
            else
            {
                minscale = 1.0 / work[3];
            }
            if (maxscale < 1.0 || minscale > 1.0)
            {

                Console.WriteLine(" Curious input WORK(3,4)=", work[2], work[3]);

                arret = true;
            }

            // SAFE     Safety factor in step size prediction
            if (work[4] == 0.0)
            {
                safe = 0.9;
            }
            else
            {
                safe = work[4];
                if (safe <= 0.001 || safe >= 1.0)
                {

                    Console.WriteLine(" Curious input for WORK(5)=", work[4]);

                    arret = true;
                }
            }
        //*/

        /// <summary>
        /// 
        /// </summary>
        /// <param name="minscale">Minimum scaling for step size selection (minscale &lt= hnew/hold).</param>
        /// <param name="maxscale">Maximum scaling for step size selection (hnew/hold &lt= maxscale).</param>
        /// <param name="hmax">Maximal step size (default tend - t).</param>
        /// <param name="safe">Safety factor in step size prediction (default 0.9).</param>
        public RosenbrockErrorController(double minscale, double maxscale, double hmax, double safe, bool predictive)
        {
            if (safe >= 1.0 || safe <= 1e-4)
            {
                throw new ArgumentException("Curious input for safety factor.", nameof(safe));
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
            double fac = Math.Max(minscale, Math.Min(maxscale, Math.Pow(err, 0.25) / safe));

            hnew = h / fac;

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
                        facgus = hacc / h * Math.Pow(err * err / erracc, 0.25) / safe;
                        facgus = Math.Max(minscale, Math.Min(maxscale, facgus));
                        fac = Math.Max(fac, facgus);
                        hnew = h / fac;
                    }

                    hacc = h;
                    erracc = Math.Max(.01, err);
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
