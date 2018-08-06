
namespace MathNet.Numerics.OdeSolvers
{
    using System;

    class StiffnessChecker
    {
        double dist;
        int nstiff;
        int n;

        int nonsti = 0;
        int iasti = 0;

        bool enabled;

        public bool Enabled
        {
            get { return enabled; }
            set { enabled = value; }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dist">Distance of the boundary of the stability region to the origin.</param>
        /// <param name="nstiff"></param>
        /// <remarks>
        /// DormandPrince5    dist = 3.25
        /// DormandPrince853  dist = 6.1
        /// </remarks>
        public StiffnessChecker(double dist, int nstiff = 1000)
        {
            this.dist = dist;
            this.nstiff = nstiff;
        }

        public bool Check(int naccpt, double h, double[] ki, double[] kj, double[] yi, double[] yj)
        {
            // Stiffness detection
            if (enabled && (naccpt % nstiff == 0 || iasti > 0))
            {
                double hlamb = 0.0;

                double temp, stnum = 0.0, stden = 0.0;

                for (int i = 0; i < n; ++i)
                {
                    temp = ki[i] - kj[i];
                    stnum += temp * temp;

                    temp = yi[i] - yj[i];
                    stden += temp * temp;
                }
                if (stden > 0.0)
                {
                    hlamb = Math.Abs(h) * Math.Sqrt(stnum / stden);
                }
                if (hlamb > dist)
                {
                    nonsti = 0;
                    ++iasti;
                    if (iasti == 15)
                    {
                        return false;
                    }
                }
                else
                {
                    ++nonsti;
                    if (nonsti == 6)
                    {
                        iasti = 0;
                    }
                }
            }
            return true;
        }
    }
}
