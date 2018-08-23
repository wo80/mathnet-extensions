// Based on Fortran code RADAU5
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;

    using S_fp = System.Action<int, double, double[], double[]>;
    using J_fp = System.Action<int, double, double[], double[], int>;
    using M_fp = System.Action<int, double[], int>;

    /// <summary>
    /// 
    /// Authors: E. Hairer and G. Wanner
    ///
    /// This code is part of the book:
    ///         E. Hairer and G. Wanner
    ///         Solving Ordinary Differential Equations II.
    ///         Stiff and differential-algebraic problems. (2nd edition)
    ///         Springer-Verlag (1996)
    /// </summary>
    public class Radau5
    {

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        struct conra5
        {
            public int nn, nn2, nn3, nn4;
            public double xsol, hsol, c2m1, c1m1;
        }

        conra5 conra5_1;

        /* Subroutine */
        int radau5_(int n, S_fp fcn, double x, double[] y,
            double xend, double h, double[] rtol, double[] atol,
            int itol, J_fp jac, int ijac, int mljac, int mujac,
            M_fp mas, int imas, int mlmas, int mumas, S_fp solout,
            int iout, double[] work, int lwork, int[] iwork,
            int liwork, double[] rpar, int[] ipar, int idid)
        {
            /* System generated locals */
            int i1;

            /* Local variables */
            int i, m1, m2, nm1, nit, lde1;
            double facl;
            int ndec, njac;
            double facr, safe;
            int ijob, nfcn;
            bool pred;
            double hmax;
            int nmax;
            double thet, expm;
            int nsol;
            double quot;
            int nind1, nind2, nind3;
            double quot1, quot2;
            int ldjac;
            bool jband;
            int ldmas;
            bool arret;
            double fnewt;
            int nstep;
            double tolst;
            int ldmas2, naccpt;

            int nrejct;
            bool implct;
            int istore;
            bool startn;
            double uround;

            /**
             *     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
             *     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
             *                     M*Y'=F(X,Y).
             *     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I)
             *     OR EXPLICIT (M=I).
             *     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA)
             *     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT.
             *     CF. SECTION IV.8
             *
             *     AUTHORS: E. HAIRER AND G. WANNER
             *              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
             *              CH-1211 GENEVE 24, SWITZERLAND
             *              E-MAIL:  Ernst.Hairer@math.unige.ch
             *                       Gerhard.Wanner@math.unige.ch
             *
             *     THIS CODE IS PART OF THE BOOK:
             *         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
             *         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
             *         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
             *         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
             *
             *     VERSION OF JULY 9, 1996
             *     (latest small correction: January 18, 2002)
             *
             *     INPUT PARAMETERS
             *     ----------------
             *     N           DIMENSION OF THE SYSTEM
             *
             *     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
             *                 VALUE OF F(X,Y):
             *                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),F(N)
             *                    F(1)=...   ETC.
             *                 RPAR, IPAR (SEE BELOW)
             *
             *     X           INITIAL X-VALUE
             *
             *     Y(N)        INITIAL VALUES FOR Y
             *
             *     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
             *
             *     H           INITIAL STEP SIZE GUESS;
             *                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
             *                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
             *                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
             *                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
             *
             *     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
             *                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
             *
             *     ITOL        SWITCH FOR RTOL AND ATOL:
             *                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
             *                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
             *                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
             *                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
             *                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
             *                     RTOL(I)*ABS(Y(I))+ATOL(I).
             *
             *     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
             *                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
             *                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
             *                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
             *                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
             *                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N)
             *                    DFY(1,1)= 0...
             *                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
             *                 FURNISHED BY THE CALLING PROGRAM.
             *                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
             *                    BE FULL AND THE PARTIAL DERIVATIVES ARE
             *                    STORED IN DFY AS
             *                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
             *                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
             *                    THE PARTIAL DERIVATIVES ARE STORED
             *                    DIAGONAL-WISE AS
             *                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
             *
             *     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
             *                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
             *                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
             *                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
             *
             *     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
             *                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
             *                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
             *                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
             *                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
             *                       THE MAIN DIAGONAL).
             *
             *     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
             *                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
             *                 NEED NOT BE DEFINED IF MLJAC=N.
             *
             *     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
             *     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
             *
             *     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
             *                 MATRIX M.
             *                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
             *                 MATRIX AND NEEDS NOT TO BE DEFINED;
             *                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
             *                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
             *                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
             *                    DOUBLE PRECISION AM(LMAS,N)
             *                    AM(1,1)= 0....
             *                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
             *                    AS FULL MATRIX LIKE
             *                         AM(I,J) = M(I,J)
             *                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
             *                    DIAGONAL-WISE AS
             *                         AM(I-J+MUMAS+1,J) = M(I,J).
             *
             *     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
             *                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
             *                       MATRIX, MAS IS NEVER CALLED.
             *                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
             *
             *     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
             *                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
             *                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
             *                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
             *                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
             *                       THE MAIN DIAGONAL).
             *                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
             *
             *     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
             *                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
             *                 NEED NOT BE DEFINED IF MLMAS=N.
             *                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
             *
             *     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
             *                 NUMERICAL SOLUTION DURING INTEGRATION.
             *                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
             *                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
             *                 IT MUST HAVE THE FORM
             *                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,
             *                                       RPAR,IPAR,IRTRN)
             *                    DOUBLE PRECISION X,Y(N),CONT(LRC)
             *                    ....
             *                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
             *                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
             *                    THE FIRST GRID-POINT).
             *                 "XOLD" IS THE PRECEEDING GRID-POINT.
             *                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
             *                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM.
             *
             *          -----  CONTINUOUS OUTPUT: -----
             *                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
             *                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
             *                 THE FUNCTION
             *                        >>>   CONTR5(I,S,CONT,LRC)   <<<
             *                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
             *                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
             *                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
             *                 DO NOT CHANGE THE ENTRIES OF CONT(LRC), IF THE
             *                 DENSE OUTPUT FUNCTION IS USED.
             *
             *     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
             *                    IOUT=0: SUBROUTINE IS NEVER CALLED
             *                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
             *
             *     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
             *                 WORK(1), WORK(2),..0, WORK(20) SERVE AS PARAMETERS
             *                 FOR THE CODE. FOR STANDARD USE OF THE CODE
             *                 WORK(1),..0,WORK(20) MUST BE SET TO ZERO BEFORE
             *                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
             *                 WORK(21),..0,WORK(LWORK) SERVE AS WORKING SPACE
             *                 FOR ALL VECTORS AND MATRICES.
             *                 "LWORK" MUST BE AT LEAST
             *                             N*(LJAC+LMAS+3*LE+12)+20
             *                 WHERE
             *                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
             *                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.0)
             *                 AND
             *                    LMAS=0              IF IMAS=0
             *                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
             *                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.0)
             *                 AND
             *                    LE=N               IF MLJAC=N (FULL JACOBIAN)
             *                    LE=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.0)
             *
             *                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
             *                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
             *                 STORAGE REQUIREMENT IS
             *                             LWORK = 4*N*N+12*N+20.
             *                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST
             *                          N*(LJAC+12)+(N-M1)*(LMAS+3*LE)+20
             *                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE THE
             *                 NUMBER N CAN BE REPLACED BY N-M1.
             *
             *     LWORK       DECLARED LENGTH OF ARRAY "WORK".
             *
             *     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
             *                 IWORK(1),IWORK(2),...0,IWORK(20) SERVE AS PARAMETERS
             *                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..0,
             *                 IWORK(20) TO ZERO BEFORE CALLING.
             *                 IWORK(21),...0,IWORK(LIWORK) SERVE AS WORKING AREA.
             *                 "LIWORK" MUST BE AT LEAST 3*N+20.
             *
             *     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
             *
             *     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
             *                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
             *                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.
             *
             * ----------------------------------------------------------------------
             *
             *     SOPHISTICATED SETTING OF PARAMETERS
             *     -----------------------------------
             *              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
             *              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
             *              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
             *              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
             *
             *    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
             *              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
             *              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
             *              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N)
             *              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1).
             *
             *    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
             *              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
             *
             *    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE
             *              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP.
             *              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7.
             *
             *    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION
             *              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD.
             *              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED.
             *              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS
             *              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN
             *              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.0).
             *              DEFAULT IS IWORK(4)=0.
             *
             *       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
             *       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
             *       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
             *       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER.
             *       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
             *       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
             *
             *    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR
             *              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
             *              DEFAULT IWORK(5)=N.
             *
             *    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0.
             *
             *    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0.
             *
             *    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY
             *              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
             *              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL
             *              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1.
             *              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS;
             *              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES
             *              OFTEN SLIGHTLY FASTER RUNS
             *
             *       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT
             *            Y(I)' = Y(I+M2)   FOR  I=1,...0,M1,
             *       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME
             *       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10).
             *       E.G.0, FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE
             *       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2.
             *       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS:
             *       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE
             *              JACOBIAN HAVE TO BE STORED
             *              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL
             *                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J)
             *                FOR I=1,N-M1 AND J=1,N.
             *              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM )
             *                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2)
             *                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM.
             *       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL
             *                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM)
             *                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2
             *                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH
             *                    OF THESE MM+1 SUBMATRICES
             *       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES
             *                NEED NOT BE DEFINED IF MLJAC=N-M1
             *       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND
             *              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
             *              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF
             *              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX.
             *              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL
             *                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1.
             *              ELSE, THE MASS MATRIX IS BANDED
             *                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1)
             *       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL
             *                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX
             *       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX
             *                NEED NOT BE DEFINED IF MLMAS=N-M1
             *
             *    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0.
             *
             *    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1.
             *
             * ----------
             *
             *    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
             *
             *    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
             *              DEFAULT 0.9D0.
             *
             *    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
             *              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
             *              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER
             *              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO
             *              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.
             *              DEFAULT 0.001D0.
             *
             *    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
             *              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER.
             *              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0)
             *
             *    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE
             *              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A
             *              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR
             *              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE
             *              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS
             *              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD.
             *              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 .
             *
             *    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
             *
             *    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION
             *              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
             *                 WORK(8) <= HNEW/HOLD <= WORK(9)
             *              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0
             *
             * -----------------------------------------------------------------------
             *
             *     OUTPUT PARAMETERS
             *     -----------------
             *     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
             *                 (AFTER SUCCESSFUL RETURN X=XEND).
             *
             *     Y(N)        NUMERICAL SOLUTION AT X
             *
             *     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
             *
             *     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
             *                   IDID= 1  COMPUTATION SUCCESSFUL,
             *                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
             *                   IDID=-1  INPUT IS NOT CONSISTENT,
             *                   IDID=-2  LARGER NMAX IS NEEDED,
             *                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
             *                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
             *
             *   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
             *                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
             *   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
             *                      OR NUMERICALLY)
             *   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
             *   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
             *   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
             *                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
             *   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES
             *   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH
             *                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS,
             *                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED
             */

            /* Parameter adjustments */
            //--y;
            //--rtol;
            //--atol;
            //--work;
            //--iwork;
            //--rpar;
            //--ipar;

            /* Function Body */
            nfcn = 0;
            njac = 0;
            nstep = 0;
            naccpt = 0;
            nrejct = 0;
            ndec = 0;
            nsol = 0;
            arret = false;
            /* -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0 */
            if (work[1] == 0.0)
            {
                uround = 1e-16;
            }
            else
            {
                uround = work[1];
                if (uround <= 1e-19 || uround >= 1.0)
                {

                    Console.WriteLine(" COEFFICIENTS HAVE 20 DIGITS, UROUND=");
                    Console.WriteLine(work[1])
                        ;

                    arret = true;
                }
            }
            /* -------- CHECK AND CHANGE THE TOLERANCES */
            expm = 0.66666666666666663;
            if (itol == 0)
            {
                if (atol[1] <= 0.0 || rtol[1] <= uround * 10.0)
                {

                    Console.WriteLine(" TOLERANCES ARE TOO SMALL");

                    arret = true;
                }
                else
                {
                    quot = atol[1] / rtol[1];
                    rtol[1] = Math.Pow(rtol[1], expm) * 0.1;
                    atol[1] = rtol[1] * quot;
                }
            }
            else
            {
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    if (atol[i] <= 0.0 || rtol[i] <= uround * 10.0)
                    {

                        Console.WriteLine(" TOLERANCES(");
                        Console.WriteLine(i);
                        Console.WriteLine(") ARE TOO SMALL");

                        arret = true;
                    }
                    else
                    {
                        quot = atol[i] / rtol[i];
                        rtol[i] = Math.Pow(rtol[i], expm) * 0.1;
                        atol[i] = rtol[i] * quot;
                    }
                }
            }
            /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
            if (iwork[2] == 0)
            {
                nmax = 100000;
            }
            else
            {
                nmax = iwork[2];
                if (nmax <= 0)
                {

                    Console.WriteLine(" WRONG INPUT IWORK(2)=");
                    Console.WriteLine(iwork[2]);

                    arret = true;
                }
            }
            /* -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS */
            if (iwork[3] == 0)
            {
                nit = 7;
            }
            else
            {
                nit = iwork[3];
                if (nit <= 0)
                {

                    Console.WriteLine(" CURIOUS INPUT IWORK(3)=");
                    Console.WriteLine(iwork[3]);

                    arret = true;
                }
            }
            /* -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS */
            if (iwork[4] == 0)
            {
                startn = false;
            }
            else
            {
                startn = true;
            }
            /* -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS */
            nind1 = iwork[5];
            nind2 = iwork[6];
            nind3 = iwork[7];
            if (nind1 == 0)
            {
                nind1 = n;
            }
            if (nind1 + nind2 + nind3 != n)
            {

                Console.WriteLine(" CURIOUS INPUT FOR IWORK(5,6,7)=");
                Console.WriteLine(nind1);
                Console.WriteLine(nind2);
                Console.WriteLine(nind3);

                arret = true;
            }
            /* -------- PRED   STEP SIZE CONTROL */
            if (iwork[8] <= 1)
            {
                pred = true;
            }
            else
            {
                pred = false;
            }
            /* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
            m1 = iwork[9];
            m2 = iwork[10];
            nm1 = n - m1;
            if (m1 == 0)
            {
                m2 = n;
            }
            if (m2 == 0)
            {
                m2 = m1;
            }
            if (m1 < 0 || m2 < 0 || m1 + m2 > n)
            {

                Console.WriteLine(" CURIOUS INPUT FOR IWORK(9,10)=");
                Console.WriteLine(m1);
                Console.WriteLine(m2);

                arret = true;
            }
            /* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
            if (work[2] == 0.0)
            {
                safe = 0.9;
            }
            else
            {
                safe = work[2];
                if (safe <= 0.001 || safe >= 1.0)
                {

                    Console.WriteLine(" CURIOUS INPUT FOR WORK(2)=");
                    Console.WriteLine(work[2])
                        ;

                    arret = true;
                }
            }
            /* ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
            if (work[3] == 0.0)
            {
                thet = 0.001;
            }
            else
            {
                thet = work[3];
                if (thet >= 1.0)
                {

                    Console.WriteLine(" CURIOUS INPUT FOR WORK(3)=");
                    Console.WriteLine(work[3])
                        ;

                    arret = true;
                }
            }
            /* --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
            tolst = rtol[1];
            if (work[4] == 0.0)
            {
                fnewt = Math.Max(uround * 10 / tolst, Math.Min(0.03, Math.Pow(tolst, 0.5)));
            }
            else
            {
                fnewt = work[4];
                if (fnewt <= uround / tolst)
                {

                    Console.WriteLine(" CURIOUS INPUT FOR WORK(4)=");
                    Console.WriteLine(work[4]);

                    arret = true;
                }
            }
            /* --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST. */
            if (work[5] == 0.0)
            {
                quot1 = 1.0;
            }
            else
            {
                quot1 = work[5];
            }
            if (work[6] == 0.0)
            {
                quot2 = 1.2;
            }
            else
            {
                quot2 = work[6];
            }
            if (quot1 > 1.0 || quot2 < 1.0)
            {

                Console.WriteLine(" CURIOUS INPUT FOR WORK(5,6)=");
                Console.WriteLine(quot1);
                Console.WriteLine(quot2);

                arret = true;
            }
            /* -------- MAXIMAL STEP SIZE */
            if (work[7] == 0.0)
            {
                hmax = xend - x;
            }
            else
            {
                hmax = work[7];
            }
            /* -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION */
            if (work[8] == 0.0)
            {
                facl = 5.0;
            }
            else
            {
                facl = 1.0 / work[8];
            }
            if (work[9] == 0.0)
            {
                facr = 0.125;
            }
            else
            {
                facr = 1.0 / work[9];
            }
            if (facl < 1.0 || facr > 1.0)
            {

                Console.WriteLine(" CURIOUS INPUT WORK(8,9)=");
                Console.WriteLine(work[8]);
                Console.WriteLine(work[9]);

                arret = true;
            }
            /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
            /*         COMPUTATION OF ARRAY ENTRIES */
            /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
            /* ---- IMPLICIT, BANDED OR NOT ? */
            implct = imas != 0;
            jband = mljac < nm1;
            /* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
            /* -- JACOBIAN  AND  MATRICES E1, E2 */
            if (jband)
            {
                ldjac = mljac + mujac + 1;
                lde1 = mljac + ldjac;
            }
            else
            {

                mljac = nm1;

                mujac = nm1;
                ldjac = nm1;
                lde1 = nm1;
            }
            /* -- MASS MATRIX */
            if (implct)
            {
                if (mlmas != nm1)
                {
                    ldmas = mlmas + mumas + 1;
                    if (jband)
                    {
                        ijob = 4;
                    }
                    else
                    {
                        ijob = 3;
                    }
                }
                else
                {

                    mumas = nm1;
                    ldmas = nm1;
                    ijob = 5;
                }
                /* ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC" */
                if (mlmas > mljac || mumas > mujac)
                {

                    Console.WriteLine("BANDWITH OF \"MAS\" NOT SMALLER THAN BANDWITH OF \"JAC\"");

                    arret = true;
                }
            }
            else
            {
                ldmas = 0;
                if (jband)
                {
                    ijob = 2;
                }
                else
                {
                    ijob = 1;
                    if (n > 2 && iwork[1] != 0)
                    {
                        ijob = 7;
                    }
                }
            }
            ldmas2 = Math.Max(1, ldmas);
            /* ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN */
            if ((implct || jband) && ijob == 7)
            {

                Console.WriteLine(" HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH FULL JACOBIAN");

                arret = true;
            }
            /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
            var z1 = new double[n];
            var z2 = new double[n];
            var z3 = new double[n];
            var y0 = new double[n];
            var scal = new double[n];
            var f1 = new double[n];
            var f2 = new double[n];
            var f3 = new double[n];
            var con = new double[n << 2];
            var _jac = new double[n * ldjac];
            var _mas = new double[nm1 * ldmas];
            var e1 = new double[nm1 * lde1];
            var e2r = new double[nm1 * lde1];
            var e2i = new double[nm1 * lde1];

            /* ------ TOTAL STORAGE REQUIREMENT ----------- */
            //istore = iee2i + nm1 * lde1 - 1;
            //if (istore > lwork)
            //{
            //    Console.WriteLine(" INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", istore);
            //    arret = true;
            //}

            /* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
            var ip1 = new int[nm1];
            var ip2 = new int[nm1];
            var iph = new int[nm1];

            /* --------- TOTAL REQUIREMENT --------------- */
            //istore = ieiph + nm1 - 1;
            //if (istore > liwork)
            //{
            //    Console.WriteLine(" INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", istore);
            //    arret = true;
            //}

            /* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
            if (arret)
            {
                idid = -1;
                return 0;
            }

            /* -------- CALL TO CORE INTEGRATOR ------------ */
            radcor_(n, fcn, x, y, xend, hmax, h, rtol, atol,
                itol, jac, ijac, mljac, mujac, mas, mlmas, mumas, (
                S_fp)solout, iout, idid, nmax, uround, safe, thet, fnewt,
                quot1, quot2, nit, ijob, startn, nind1, nind2, nind3,
                pred, facl, facr, m1, m2, nm1, implct, jband, ldjac,
                lde1, ldmas2, z1, z2, z3, y0,
                 scal, f1, f2, f3,
                _jac, e1, e2r, e2i, _mas,
                ip1, ip2, iph, con, ref nfcn,
                ref njac, ref nstep, ref naccpt, ref nrejct, ref ndec, ref nsol);

            iwork[14] = nfcn;
            iwork[15] = njac;
            iwork[16] = nstep;
            iwork[17] = naccpt;
            iwork[18] = nrejct;
            iwork[19] = ndec;
            iwork[20] = nsol;

            /* -------- RESTORE TOLERANCES */
            expm = 1.0 / expm;
            if (itol == 0)
            {
                quot = atol[1] / rtol[1];
                rtol[1] = Math.Pow(rtol[1] * 10.0, expm);
                atol[1] = rtol[1] * quot;
            }
            else
            {
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    quot = atol[i] / rtol[i];
                    rtol[i] = Math.Pow(rtol[i] * 10.0, expm);
                    atol[i] = rtol[i] * quot;
                }
            }
            /* ----------- RETURN ----------- */
            return 0;
        } /* radau5_ */


        /* Subroutine */
        int radcor_(int n, S_fp fcn, double x, double[]
            y, double xend, double hmax, double h, double[]
            rtol, double[] atol, int itol, J_fp jac, int ijac,
            int mljac, int mujac, M_fp mas, int mlmas, int
            mumas, S_fp solout, int iout, int idid, int nmax,
            double uround, double safe, double thet, double
            fnewt, double quot1, double quot2, int nit, int
            ijob, bool startn, int nind1, int nind2, int nind3,
            bool pred, double facl, double facr, int m1,
            int m2, int nm1, bool implct, bool banded, int
            ldjac, int lde1, int ldmas, double[] z1, double[] z2,
            double[] z3, double[] y0, double[] scal, double[] f1,
            double[] f2, double[] f3, double[] fjac, double[] e1,
            double[] e2r, double[] e2i, double[] fmas, int[] ip1,
            int[] ip2, int[] iphes, double[] cont, ref int nfcn,
            ref int njac, ref int nstep, ref int naccpt, ref int nrejct,
            ref int ndec, ref int nsol)
        {
            /* Format strings */
            //static char fmt_979[] = "(  EXIT OF RADAU5 AT X= ,e18.4)";

            /* System generated locals */
            int i1, i2, i3, i4;
            double d1, d2, d3;


            /* Local variables */
            int i, j, k, l;
            double a1, a2, c1, c2, a3;
            int j1, n2, n3;
            double u1, ak;
            int md;
            double t11, t12, t13, t21, t22, t23, t31;
            int mm;
            double qt, dd1, dd2, dd3, ak1, ak2, ak3, f1i, f2i, f3i, c1q, c2q, c3q,
                 z1i, z2i, z3i, sq6, fac, ti11, cno;
            int lrc;
            double ti12, ti13, ti21, ti22, ti23, ti31, ti32, ti33;
            int ier;
            double xph, thq, err = 0, fac1, cfac, hacc = 0, c1mc2, beta;
            int lbeg;
            double alph, hold;
            int lend;
            double delt, hnew;
            bool last;
            double hopt, xold;
            int newt;
            double dyno, dyth, quot, hhfac, betan, alphn, denom, theta, ysafe,
                hmaxn;
            int nsing;
            bool first;
            int irtrn = 0, nrsol, nsolu;
            double qnewt, xosol, acont3;
            bool index1, index2, index3, caljac;
            double faccon;
            bool calhes;
            double erracc = 0;
            int mujacj;
            bool reject;
            double facgus;
            int mujacp;
            double dynold = 0, posneg;
            double thqold = 0;

            /* Fortran I/O blocks */


            /*     CORE INTEGRATOR FOR RADAU5 */
            /*     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED */


            /*  INITIALISATIONS */

            /* --------- DUPLIFY N FOR COMMON BLOCK CONT ----- */
            /* Parameter adjustments */
            //--cont;
            //--f3;
            //--f2;
            //--f1;
            //--scal;
            //--y0;
            //--z3;
            //--z2;
            //--z1;
            //--y;
            //--rtol;
            //--atol;
            //--iphes;
            //--ip2;
            //--ip1;
            //fjac_offset = 1 + ldjac;
            //fjac -= fjac_offset;
            //e2i_dim1 = lde1;
            //e2i_offset = 1 + e2i_dim1;
            //e2i -= e2i_offset;
            //e2r_dim1 = lde1;
            //e2r_offset = 1 + e2r_dim1;
            //e2r -= e2r_offset;
            //e1_dim1 = lde1;
            //e1_offset = 1 + e1_dim1;
            //e1 -= e1_offset;
            //fmas_dim1 = ldmas;
            //fmas_offset = 1 + fmas_dim1;
            //fmas -= fmas_offset;
            //--rpar;
            //--ipar;

            /* Function Body */
            conra5_1.nn = n;
            conra5_1.nn2 = n << 1;
            conra5_1.nn3 = n * 3;
            lrc = n << 2;
            /* -------- CHECK THE INDEX OF THE PROBLEM ----- */
            index1 = nind1 != 0;
            index2 = nind2 != 0;
            index3 = nind3 != 0;
            /* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
            if (implct)
            {
                mas(nm1, fmas, ldmas);
            }
            /* ---------- CONSTANTS --------- */
            sq6 = Math.Sqrt(6.0);
            c1 = (4.0 - sq6) / 10.0;
            c2 = (sq6 + 4.0) / 10.0;
            conra5_1.c1m1 = c1 - 1.0;
            conra5_1.c2m1 = c2 - 1.0;
            c1mc2 = c1 - c2;
            dd1 = -(sq6 * 7.0 + 13.0) / 3.0;
            dd2 = (sq6 * 7.0 - 13.0) / 3.0;
            dd3 = -0.33333333333333331;
            u1 = (Math.Pow(81.0, 0.33333333333333331) + 6.0 - Math.Pow(9.0, 0.33333333333333331)) / 30.0;
            alph = (12.0 - Math.Pow(81.0, 0.33333333333333331) + Math.Pow(9.0, 0.33333333333333331)) / 60.0;
            beta = (Math.Pow(81.0, 0.33333333333333331) + Math.Pow(9.0, 0.33333333333333331)) * Math.Sqrt(3.0) / 60.0;
            /* Computing 2nd power */
            d1 = alph;
            /* Computing 2nd power */
            d2 = beta;
            cno = d1 * d1 + d2 * d2;
            u1 = 1.0 / u1;
            alph /= cno;
            beta /= cno;
            t11 = 0.091232394870892942792;
            t12 = -0.14125529502095420843;
            t13 = -0.030029194105147424492;
            t21 = 0.24171793270710701896;
            t22 = 0.20412935229379993199;
            t23 = 0.38294211275726193779;
            t31 = 0.96604818261509293619;
            ti11 = 4.325579890063155351;
            ti12 = 0.33919925181580986954;
            ti13 = 0.54177053993587487119;
            ti21 = -4.1787185915519047273;
            ti22 = -0.32768282076106238708;
            ti23 = 0.47662355450055045196;
            ti31 = -0.50287263494578687595;
            ti32 = 2.5719269498556054292;
            ti33 = -0.59603920482822492497;
            if (m1 > 0)
            {
                ijob += 10;
            }
            posneg = d_sign(1.0, xend - x);
            hmaxn = Math.Min(Math.Abs(hmax), Math.Abs(xend - x));
            if (Math.Abs(h) <= uround * 10.0)
            {
                h = 1e-6;
            }
            h = Math.Min(Math.Abs(h), hmaxn);
            h = d_sign(h, posneg);
            hold = h;
            reject = false;
            first = true;
            last = false;
            if ((x + h * 1.0001 - xend) * posneg >= 0.0)
            {
                h = xend - x;
                last = true;
            }
            hopt = h;
            faccon = 1.0;
            cfac = safe * ((nit << 1) + 1);
            nsing = 0;
            xold = x;
            if (iout != 0)
            {
                irtrn = 1;
                nrsol = 1;
                xosol = xold;
                conra5_1.xsol = x;
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    cont[i] = y[i];
                }
                nsolu = n;
                conra5_1.hsol = hold;
                //solout(&nrsol, &xosol, &conra5_1.xsol, y, cont, &lrc, &nsolu, &irtrn);
                if (irtrn < 0)
                {
                    goto L179;
                }
            }
            linal_1.mle = mljac;
            linal_1.mue = mujac;
            linal_1.mbjac = mljac + mujac + 1;
            linal_1.mbb = mlmas + mumas + 1;
            linal_1.mdiag = linal_1.mle + linal_1.mue + 1;
            linal_1.mdiff = linal_1.mle + linal_1.mue - mumas;
            linal_1.mbdiag = mumas + 1;
            n2 = n << 1;
            n3 = n * 3;
            if (itol == 0)
            {
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    scal[i] = atol[1] + rtol[1] * Math.Abs(y[i]);
                }
            }
            else
            {
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    scal[i] = atol[i] + rtol[i] * Math.Abs(y[i]);
                }
            }
            hhfac = h;
            fcn(n, x, y, y0);
            ++(nfcn);
            /* --- BASIC INTEGRATION STEP */
            L10:
            /* *** *** *** *** *** *** *** */
            /*  COMPUTATION OF THE JACOBIAN */
            /* *** *** *** *** *** *** *** */
            ++(njac);
            if (ijac == 0)
            {
                /* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
                if (banded)
                {
                    /* --- JACOBIAN IS BANDED */
                    mujacp = mujac + 1;
                    md = Math.Min(linal_1.mbjac, m2);
                    i1 = m1 / m2 + 1;
                    for (mm = 1; mm <= i1; ++mm)
                    {
                        i2 = md;
                        for (k = 1; k <= i2; ++k)
                        {
                            j = k + (mm - 1) * m2;
                            L12:
                            f1[j] = y[j];
                            f2[j] = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(y[j])));
                            y[j] += f2[j];
                            j += md;
                            if (j <= mm * m2)
                            {
                                goto L12;
                            }
                            fcn(n, x, y, cont);
                            j = k + (mm - 1) * m2;
                            j1 = k;
                            /* Computing MAX */
                            i3 = 1, i4 = j1 - mujac;
                            lbeg = Math.Max(i3, i4) + m1;
                            L14:
                            /* Computing MIN */
                            i3 = m2, i4 = j1 + mljac;
                            lend = Math.Min(i3, i4) + m1;
                            y[j] = f1[j];
                            mujacj = mujacp - j1 - m1;
                            i3 = lend;
                            for (l = lbeg; l <= i3; ++l)
                            {
                                fjac[l + mujacj + j * ldjac] = (cont[l] - y0[l]) / f2[j];
                            }
                            j += md;
                            j1 += md;
                            lbeg = lend + 1;
                            if (j <= mm * m2)
                            {
                                goto L14;
                            }
                        }
                    }
                }
                else
                {
                    /* --- JACOBIAN IS FULL */
                    i1 = n;
                    for (i = 1; i <= i1; ++i)
                    {
                        ysafe = y[i];
                        delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                        y[i] = ysafe + delt;
                        fcn(n, x, y, cont);
                        i2 = n;
                        for (j = m1 + 1; j <= i2; ++j)
                        {
                            fjac[j - m1 + i * ldjac] = (cont[j] - y0[j]) /
                                delt;
                        }
                        y[i] = ysafe;
                    }
                }
            }
            else
            {
                /* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
                jac(n, x, y, fjac, ldjac);
            }
            caljac = true;
            calhes = true;
            L20:
            /* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
            fac1 = u1 / h;
            alphn = alph / h;
            betan = beta / h;
            decomr_(n, fjac, ldjac, fmas, ldmas, mlmas,
                mumas, m1, m2, nm1, fac1, e1, lde1, ip1, ier, ijob, calhes, iphes);
            if (ier != 0)
            {
                goto L78;
            }
            decomc_(n, fjac, ldjac, fmas, ldmas, mlmas,
                mumas, m1, m2, nm1, alphn, betan, e2r, e2i, lde1, ip2, ier, ijob);
            if (ier != 0)
            {
                goto L78;
            }
            ++(ndec);
            L30:
            ++(nstep);
            if (nstep > nmax)
            {
                goto L178;
            }
            if (Math.Abs(h) * 0.1 <= Math.Abs(x) * uround)
            {
                goto L177;
            }
            if (index2)
            {
                i1 = nind1 + nind2;
                for (i = nind1 + 1; i <= i1; ++i)
                {
                    scal[i] /= hhfac;
                }
            }
            if (index3)
            {
                i1 = nind1 + nind2 + nind3;
                for (i = nind1 + nind2 + 1; i <= i1; ++i)
                {
                    scal[i] /= hhfac * hhfac;
                }
            }
            xph = x + h;
            /* *** *** *** *** *** *** *** */
            /*  STARTING VALUES FOR NEWTON ITERATION */
            /* *** *** *** *** *** *** *** */
            if (first || startn)
            {
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    z1[i] = 0.0;
                    z2[i] = 0.0;
                    z3[i] = 0.0;
                    f1[i] = 0.0;
                    f2[i] = 0.0;
                    f3[i] = 0.0;
                }
            }
            else
            {
                c3q = h / hold;
                c1q = c1 * c3q;
                c2q = c2 * c3q;
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    ak1 = cont[i + n];
                    ak2 = cont[i + n2];
                    ak3 = cont[i + n3];
                    z1i = c1q * (ak1 + (c1q - conra5_1.c2m1) * (ak2 + (c1q -
                        conra5_1.c1m1) * ak3));
                    z2i = c2q * (ak1 + (c2q - conra5_1.c2m1) * (ak2 + (c2q -
                        conra5_1.c1m1) * ak3));
                    z3i = c3q * (ak1 + (c3q - conra5_1.c2m1) * (ak2 + (c3q -
                        conra5_1.c1m1) * ak3));
                    z1[i] = z1i;
                    z2[i] = z2i;
                    z3[i] = z3i;
                    f1[i] = ti11 * z1i + ti12 * z2i + ti13 * z3i;
                    f2[i] = ti21 * z1i + ti22 * z2i + ti23 * z3i;
                    f3[i] = ti31 * z1i + ti32 * z2i + ti33 * z3i;
                }
            }
            /* *** *** *** *** *** *** *** */
            /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
            /* *** *** *** *** *** *** *** */
            newt = 0;
            faccon = Math.Pow(Math.Max(faccon, uround), 0.8);
            theta = Math.Abs(thet);
            L40:
            if (newt >= nit)
            {
                goto L78;
            }
            /* ---     COMPUTE THE RIGHT-HAND SIDE */
            i1 = n;
            for (i = 1; i <= i1; ++i)
            {
                cont[i] = y[i] + z1[i];
            }
            fcn(n, x + c1 * h, cont, z1);
            i1 = n;
            for (i = 1; i <= i1; ++i)
            {
                cont[i] = y[i] + z2[i];
            }
            fcn(n, x + c2 * h, cont, z2);
            i1 = n;
            for (i = 1; i <= i1; ++i)
            {
                cont[i] = y[i] + z3[i];
            }
            fcn(n, xph, cont, z3);
            nfcn += 3;
            /* ---     SOLVE THE LINEAR SYSTEMS */
            i1 = n;
            for (i = 1; i <= i1; ++i)
            {
                a1 = z1[i];
                a2 = z2[i];
                a3 = z3[i];
                z1[i] = ti11 * a1 + ti12 * a2 + ti13 * a3;
                z2[i] = ti21 * a1 + ti22 * a2 + ti23 * a3;
                z3[i] = ti31 * a1 + ti32 * a2 + ti33 * a3;
            }
            slvrad_(n, fjac, ldjac, mljac, mujac, fmas,
                ldmas, mlmas, mumas, m1, m2, nm1, fac1, alphn, betan, e1,
                 e2r, e2i, lde1, z1, z2, z3, f1, f2, f3, cont, ip1, ip2, iphes, ier, ijob);
            ++(nsol);
            ++newt;
            dyno = 0.0;
            i1 = n;
            for (i = 1; i <= i1; ++i)
            {
                denom = scal[i];
                /* Computing 2nd power */
                d1 = z1[i] / denom;
                /* Computing 2nd power */
                d2 = z2[i] / denom;
                /* Computing 2nd power */
                d3 = z3[i] / denom;
                dyno = dyno + d1 * d1 + d2 * d2 + d3 * d3;
            }
            dyno = Math.Sqrt(dyno / n3);
            /* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
            if (newt > 1 && newt < nit)
            {
                thq = dyno / dynold;
                if (newt == 2)
                {
                    theta = thq;
                }
                else
                {
                    theta = Math.Sqrt(thq * thqold);
                }
                thqold = thq;
                if (theta < .99)
                {
                    faccon = theta / (1.0 - theta);
                    dyth = faccon * dyno * Math.Pow(theta, nit - 1 - newt) / fnewt; // TODO: pow_di
                    if (dyth >= 1.0)
                    {
                        qnewt = Math.Max(1e-4, Math.Min(20.0, dyth));
                        hhfac = Math.Pow(qnewt, -1.0 / (nit + 4.0 - 1 - newt)) * 0.8;
                        h = hhfac * h;
                        reject = true;
                        last = false;
                        if (caljac)
                        {
                            goto L20;
                        }
                        goto L10;
                    }
                }
                else
                {
                    goto L78;
                }
            }
            dynold = Math.Max(dyno, uround);
            i1 = n;
            for (i = 1; i <= i1; ++i)
            {
                f1i = f1[i] + z1[i];
                f2i = f2[i] + z2[i];
                f3i = f3[i] + z3[i];
                f1[i] = f1i;
                f2[i] = f2i;
                f3[i] = f3i;
                z1[i] = t11 * f1i + t12 * f2i + t13 * f3i;
                z2[i] = t21 * f1i + t22 * f2i + t23 * f3i;
                z3[i] = t31 * f1i + f2i;
            }
            if (faccon * dyno > fnewt)
            {
                goto L40;
            }
            /* --- ERROR ESTIMATION */
            estrad_(n, fjac, ldjac, mljac, mujac, fmas,
                ldmas, mlmas, mumas, h, dd1, dd2, dd3, fcn, nfcn, y0,
                y, ijob, x, m1, m2, nm1, e1, lde1, z1, z2, z3,
                cont, f1, f2, ip1, iphes, scal, err, first, reject, fac1);
            /* --- COMPUTATION OF HNEW */
            /* --- WE REQUIRE .2<=HNEW/H<=8. */
            fac = Math.Min(safe, cfac / (newt + (nit << 1)));
            quot = Math.Max(facr, Math.Min(facl, Math.Pow(err, 0.25) / fac));
            hnew = h / quot;
            /* *** *** *** *** *** *** *** */
            /*  IS THE ERROR SMALL ENOUGH ? */
            /* *** *** *** *** *** *** *** */
            if (err < 1.0)
            {
                /* --- STEP IS ACCEPTED */
                first = false;
                ++(naccpt);
                if (pred)
                {
                    /*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
                    if (naccpt > 1)
                    {
                        facgus = hacc / h * Math.Pow(err * err / erracc, 0.25) / safe;
                        facgus = Math.Max(facr, Math.Min(facl, facgus));
                        quot = Math.Max(quot, facgus);
                        hnew = h / quot;
                    }
                    hacc = h;
                    erracc = Math.Max(.01, err);
                }
                xold = x;
                hold = h;
                x = xph;
                i1 = n;
                for (i = 1; i <= i1; ++i)
                {
                    y[i] += z3[i];
                    z2i = z2[i];
                    z1i = z1[i];
                    cont[i + n] = (z2i - z3[i]) / conra5_1.c2m1;
                    ak = (z1i - z2i) / c1mc2;
                    acont3 = z1i / c1;
                    acont3 = (ak - acont3) / c2;
                    cont[i + n2] = (ak - cont[i + n]) / conra5_1.c1m1;
                    cont[i + n3] = cont[i + n2] - acont3;
                }
                if (itol == 0)
                {
                    i1 = n;
                    for (i = 1; i <= i1; ++i)
                    {
                        scal[i] = atol[1] + rtol[1] * Math.Abs(y[i]);
                    }
                }
                else
                {
                    i1 = n;
                    for (i = 1; i <= i1; ++i)
                    {
                        scal[i] = atol[i] + rtol[i] * Math.Abs(y[i]);
                    }
                }
                if (iout != 0)
                {
                    nrsol = naccpt + 1;
                    conra5_1.xsol = x;
                    xosol = xold;
                    i1 = n;
                    for (i = 1; i <= i1; ++i)
                    {
                        cont[i] = y[i];
                    }
                    nsolu = n;
                    conra5_1.hsol = hold;
                    //solout(&nrsol, &xosol, &conra5_1.xsol, y, cont, &lrc, &nsolu, &irtrn);
                    if (irtrn < 0)
                    {
                        goto L179;
                    }
                }
                caljac = false;
                if (last)
                {
                    h = hopt;
                    idid = 1;
                    return 0;
                }
                fcn(n, x, y, y0);
                ++(nfcn);
                hnew = posneg * Math.Min(Math.Abs(hnew), hmaxn);
                hopt = hnew;
                hopt = Math.Min(h, hnew);
                if (reject)
                {
                    hnew = posneg * Math.Min(Math.Abs(hnew), Math.Abs(h));
                }
                reject = false;
                if ((x + hnew / quot1 - xend) * posneg >= 0.0)
                {
                    h = xend - x;
                    last = true;
                }
                else
                {
                    qt = hnew / h;
                    hhfac = h;
                    if (theta <= thet && qt >= quot1 && qt <= quot2)
                    {
                        goto L30;
                    }
                    h = hnew;
                }
                hhfac = h;
                if (theta <= thet)
                {
                    goto L20;
                }
                goto L10;
            }
            else
            {
                /* --- STEP IS REJECTED */
                reject = true;
                last = false;
                if (first)
                {
                    h *= 0.1;
                    hhfac = 0.1;
                }
                else
                {
                    hhfac = hnew / h;
                    h = hnew;
                }
                if (naccpt >= 1)
                {
                    ++(nrejct);
                }
                if (caljac)
                {
                    goto L20;
                }
                goto L10;
            }
            /* --- UNEXPECTED STEP-REJECTION */
            L78:
            if (ier != 0)
            {
                ++nsing;
                if (nsing >= 5)
                {
                    goto L176;
                }
            }
            h *= 0.5;
            hhfac = 0.5;
            reject = true;
            last = false;
            if (caljac)
            {
                goto L20;
            }
            goto L10;
            /* --- FAIL EXIT */
            L176:

            Console.WriteLine((x));


            Console.WriteLine(" MATRIX IS REPEATEDLY SINGULAR, IER=");
            Console.WriteLine(ier);

            idid = -4;
            return 0;
            L177:

            Console.WriteLine((x));


            Console.WriteLine(" STEP SIZE T0O SMALL, H=");
            Console.WriteLine((h));

            idid = -3;
            return 0;
            L178:

            Console.WriteLine((x));


            Console.WriteLine(" MORE THAN NMAX =");
            Console.WriteLine((nmax));
            Console.WriteLine("STEPS ARE NEEDED");

            idid = -2;
            return 0;
            /* --- EXIT CAUSED BY SOLOUT */
            L179:

            Console.WriteLine((x));

            idid = 2;
            return 0;
        } /* radcor_ */


        /*     END OF SUBROUTINE RADCOR */

        /* *********************************************************** */

        double contr5_(int i, double x, double[] cont, int lrc)
        {
            /* System generated locals */
            double ret_val;

            /* Local variables */
            double s;

            /* ---------------------------------------------------------- */
            /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
            /*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
            /*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
            /*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5). */
            /* ---------------------------------------------------------- */
            /* Parameter adjustments */
            //--cont;

            /* Function Body */
            s = (x - conra5_1.xsol) / conra5_1.hsol;
            ret_val = cont[i] + s * (cont[i + conra5_1.nn] + (s - conra5_1.c2m1)
                 * (cont[i + conra5_1.nn2] + (s - conra5_1.c1m1) * cont[i + conra5_1.nn3]));
            return ret_val;
        } /* contr5_ */
    }
}