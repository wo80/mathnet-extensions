
namespace MathNet.Numerics.OdeSolvers.RK
{
    using System;

    /// <summary>
    /// Check stiffness while integrating ODE's with Runge-Kutta
    /// </summary>
    public class StiffnessChecker
    {
        double dist;
        int nstiff;
        
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
        /// <param name="nstiff">Test for stiffness is activated after step number j*nstiff (default = 1000).</param>
        /// <remarks>
        /// DOPRI5    dist = 3.25
        /// DOPRI853  dist = 6.1
        /// </remarks>
        public StiffnessChecker(double dist, int nstiff = 1000)
        {
            this.dist = dist;
            this.nstiff = nstiff;
        }

        public bool Check(int naccpt, double h, double[] ki, double[] kj, double[] yi, double[] yj)
        {
            int n = yi.Length;

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
