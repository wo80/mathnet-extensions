using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.Numerics.OdeSolvers
{
    class ErrorController
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

        public int Rejected { get { return nrejected; } }

        public int Accepted { get { return naccepted; } }

        public ErrorController(int order, double minscale, double maxscale, double hmax, double safe, double beta)
        {
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
