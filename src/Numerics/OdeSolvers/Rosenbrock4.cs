﻿using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using System;

namespace MathNet.Numerics.OdeSolvers
{
    using S_fp = Action<int, double, double[], double[]>;
    using J_fp = Action<int, double, double[], double[], int>;
    using M_fp = Action<int, double[], int>;
    using P_fp = Action<int, double, double, double[], double[], int, int, int>;

    public class Rosenbrock4_
    {
        struct linal
        {
            public int mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
        };

        linal linal_1 = new linal();

        struct conros
        {
            public double xold, h;
            public int n;
        };

        conros conros_ = new conros();

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c21, c31,
            c32, c41, c42, c43, c51, c52, c53, c54, c61, c62, c63, c64, c65,
            d21, d22, d23, d24, d25, d31, d32, d33, d34, d35;
        double c2, c3, c4, d1, d2, d3, d4, gamma;

        /* Subroutine */
        public int rodas_(int n, S_fp fcn, int ifcn, double x, double[] y, double xend, double h, double[] rtol,
        double[] atol, int itol, J_fp jac, int ijac, int mljac, int mujac, S_fp dfx, int idfx, M_fp mas,
        int imas, int mlmas, int mumas, P_fp solout, int iout,
        double[] work, int lwork, int[] iwork, int liwork, int idid)
        {
            /* Local variables */
            int i, m1, m2, nm1, lde;
            double fac1, fac2;
            int ndec, njac;
            double safe;
            int ijob, nfcn;
            bool pred;
            int meth;
            double hmax;
            int nmax, nsol, ldjac;
            bool jband;
            int ldmas;
            bool arret;
            int nstep, ldmas2, naccpt, nrejct;
            bool implct;
            bool autnms;
            double uround;

            /** ----------------------------------------------------------
             *     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
             *     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
             *     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4
             *     (WITH STEP SIZE CONTROL).
             *     C.F. SECTIONS IV.7  AND VI.3
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
             *     VERSION OF OCTOBER 28, 1996
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
             *     IFCN        GIVES INFORMATION ON FCN:
             *                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
             *                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
             *
             *     X           INITIAL X-VALUE
             *
             *     Y(N)        INITIAL VALUES FOR Y
             *
             *     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
             *
             *     H           INITIAL STEP SIZE GUESS;
             *                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
             *                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
             *                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
             *                 ADAPTS ITS STEP SIZE (IF H=0.D0, THE CODE PUTS H=1.D-6).
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
             *                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS
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
             *     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
             *                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X
             *                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1;
             *                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0).
             *                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM
             *                    SUBROUTINE DFX(N,X,Y,FX,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),FX(N)
             *                    FX(1)= 0...
             *
             *     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX:
             *                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE
             *                       DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED.
             *                    IDFX=1: DF/DX IS SUPPLIED BY SUBROUTINE DFX.
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
             *                    IS SET &lt;0, RODAS RETURNS TO THE CALLING PROGRAM.
             *
             *          -----  CONTINUOUS OUTPUT: -----
             *                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
             *                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
             *                 THE FUNCTION
             *                        >>>   CONTRO(I,S,CONT,LRC)
             *                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
             *                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
             *                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
             *
             *     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT:
             *                    IOUT=0: SUBROUTINE IS NEVER CALLED
             *                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
             *
             *     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
             *                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
             *                 "LWORK" MUST BE AT LEAST
             *                             N*(LJAC+LMAS+LE1+14)+20
             *                 WHERE
             *                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
             *                    LJAC=MLJAC+MUJAC+1  IF MLJAC&lt;N (BANDED JAC.0)
             *                 AND
             *                    LMAS=0              IF IMAS=0
             *                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
             *                    LMAS=MLMAS+MUMAS+1  IF MLMAS&lt;N (BANDED MASS-M.0)
             *                 AND
             *                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
             *                    LE1=2*MLJAC+MUJAC+1 IF MLJAC&lt;N (BANDED JAC.0).
             *                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
             *                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
             *                 STORAGE REQUIREMENT IS
             *                             LWORK = 2*N*N+14*N+20.
             *                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST
             *                          N*(LJAC+14)+(N-M1)*(LMAS+LE1)+20
             *                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE1 THE
             *                 NUMBER N CAN BE REPLACED BY N-M1.
             *
             *     LWORK       DECLARED LENGTH OF ARRAY "WORK".
             *
             *     IWORK       int WORKING SPACE OF LENGTH "LIWORK".
             *                 "LIWORK" MUST BE AT LEAST N+20.
             *
             *     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
             *
             *     RPAR, IPAR  REAL AND int PARAMETERS (OR PARAMETER ARRAYS) WHICH
             *                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
             *                 PROGRAM AND THE FCN, DFX, JAC, MAS, SOLOUT SUBROUTINES.
             *
             * ----------------------------------------------------------------------
             *
             *     SOPHISTICATED SETTING OF PARAMETERS
             *     -----------------------------------
             *              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
             *              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(4)
             *              AS WELL AS IWORK(1),IWORK(2) DIFFERENT FROM ZERO.
             *              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
             *
             *    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
             *              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
             *
             *    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
             *              IF IWORK(2).EQ.1  METHOD (SEE BOOK, PAGE 452)
             *              IF IWORK(2).EQ.2  SAME METHOD WITH DIFFERENT PARAMETERS
             *              IF IWORK(2).EQ.3  METHOD WITH COEFF. OF GERD STEINEBACH
             *              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1.
             *
             *    IWORK(3)  SWITCH FOR STEP SIZE STRATEGY
             *              IF IWORK(3).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
             *              IF IWORK(3).EQ.2  CLASSICAL APPROACH
             *              THE DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=1.
             *
             *       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT
             *            Y(I)' = Y(I+M2)   FOR  I=1,...,M1,
             *       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME
             *       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10).
             *       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE
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
             *                0&lt;=MLJAC&lt;N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM)
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
             *                0&lt;=MLMAS&lt;N-M1: LOWER BANDWIDTH OF THE MASS MATRIX
             *       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX
             *                NEED NOT BE DEFINED IF MLMAS=N-M1
             *
             *    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0.
             *
             *    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1.
             *
             *    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
             *
             *    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
             *
             *    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
             *              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
             *                 WORK(3) &lt;= HNEW/HOLD &lt;= WORK(4)
             *              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=6.D0
             *
             *    WORK(5)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
             *              DEFAULT 0.9D0.
             *
             * -----------------------------------------------------------------------
             *
             *     OUTPUT PARAMETERS
             *     -----------------
             *     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
             *                 (AFTER SUCCESSFUL RETURN X=XEND)
             *
             *     Y(N)        SOLUTION AT X
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
             *   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
             *   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
             * --------------------------------------------------------- */

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
            naccpt = 0;
            nrejct = 0;
            nstep = 0;
            njac = 0;
            ndec = 0;
            nsol = 0;
            arret = false;
            /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
            if (iwork[0] == 0)
            {
                nmax = 100000;
            }
            else
            {
                nmax = iwork[0];
                if (nmax <= 0)
                {

                    Console.WriteLine(" WRONG INPUT IWORK(1)=", iwork[0]);

                    arret = true;
                }
            }
            /* -------- METH   COEFFICIENTS OF THE METHOD */
            if (iwork[1] == 0)
            {
                meth = 1;
            }
            else
            {
                meth = iwork[1];
                if (meth <= 0 || meth >= 4)
                {

                    Console.WriteLine(" CURIOUS INPUT IWORK(2)=", iwork[1]);

                    arret = true;
                }
            }
            /* -------- PRED   STEP SIZE CONTROL */
            if (iwork[2] <= 1)
            {
                pred = true;
            }
            else
            {
                pred = false;
            }
            /* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
            m1 = iwork[8];
            m2 = iwork[9];
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

                Console.WriteLine(" CURIOUS INPUT FOR IWORK(9,10)=", m1, m2);

                arret = true;
            }
            /* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
            if (work[0] == 0.0)
            {
                uround = 1e-16;
            }
            else
            {
                uround = work[0];
                if (uround < 1e-16 || uround >= 1.0)
                {

                    Console.WriteLine(" COEFFICIENTS HAVE 16 DIGITS, UROUND=", work[0]);

                    arret = true;
                }
            }
            /* -------- MAXIMAL STEP SIZE */
            if (work[1] == 0.0)
            {
                hmax = xend - x;
            }
            else
            {
                hmax = work[1];
            }
            /* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
            if (work[2] == 0.0)
            {
                fac1 = 5.0;
            }
            else
            {
                fac1 = 1.0 / work[2];
            }
            if (work[3] == 0.0)
            {
                fac2 = 0.16666666666666666;
            }
            else
            {
                fac2 = 1.0 / work[3];
            }
            if (fac1 < 1.0 || fac2 > 1.0)
            {

                Console.WriteLine(" CURIOUS INPUT WORK(3,4)=", work[2], work[3]);

                arret = true;
            }
            /* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
            if (work[4] == 0.0)
            {
                safe = 0.9;
            }
            else
            {
                safe = work[4];
                if (safe <= 0.001 || safe >= 1.0)
                {

                    Console.WriteLine(" CURIOUS INPUT FOR WORK(5)=", work[4]);

                    arret = true;
                }
            }
            /* --------- CHECK IF TOLERANCES ARE O.K. */
            if (itol == 0)
            {
                if (atol[0] <= 0.0 || rtol[0] <= uround * 10.0)
                {

                    Console.WriteLine(" TOLERANCES ARE TOO SMALL");

                    arret = true;
                }
            }
            else
            {
                for (i = 0; i < n; ++i)
                {
                    if (atol[i] <= 0.0 || rtol[i] <= uround * 10.0)
                    {

                        Console.WriteLine(" TOLERANCES(%d) ARE TOO SMALL", i);

                        arret = true;
                    }
                }
            }
            /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
            /*         COMPUTATION OF ARRAY ENTRIES */
            /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
            /* ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ? */
            autnms = ifcn == 0;
            implct = imas != 0;
            jband = mljac < nm1;
            /* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
            /* -- JACOBIAN AND MATRIX E */
            if (jband)
            {
                ldjac = mljac + mujac + 1;
                lde = mljac + ldjac;
            }
            else
            {

                mljac = nm1;

                mujac = nm1;
                ldjac = nm1;
                lde = nm1;
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
                    ldmas = nm1;
                    ijob = 5;
                }
                /* ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC" */
                if (mlmas > mljac || mumas > mujac)
                {

                    Console.WriteLine("BANDWITH OF \"MAS\" NOT LARGER THAN BANDWITH OF \"JAC\"");

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
                }
            }
            ldmas2 = Math.Max(1, ldmas);
            /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
            var ynew = new double[n];
            var dy1 = new double[n];
            var dy = new double[n];
            var ak1 = new double[n];
            var ak2 = new double[n];
            var ak3 = new double[n];
            var ak4 = new double[n];
            var ak5 = new double[n];
            var ak6 = new double[n];
            var fx = new double[n];
            var con = new double[n << 2];
            var _jac = new double[n * ldjac];
            var _mas = new double[nm1 * ldmas];
            var e = new double[nm1 * lde];
            /* ------ TOTAL STORAGE REQUIREMENT ----------- */
            //istore = e + nm1 * lde - 1;
            //if (istore > lwork)
            //{
            //    Console.WriteLine(" INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", istore);
            //    arret = true;
            //}

            /* ------- ENTRY POINTS FOR int WORKSPACE ----- */
            var ip = new int[nm1];

            //istore = ip + nm1 - 1;
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
            roscor_(n, fcn, x, y, xend, hmax, h, rtol, atol,
                itol, jac, ijac, mljac, mujac, dfx, idfx, mas,
                mlmas, mumas, solout, iout, idid, nmax, uround, meth, 
                ijob, fac1, fac2, safe, autnms, implct, jband, pred, 
                ldjac, lde, ldmas2, ynew, dy1, dy, 
                ak1, ak2, ak3, ak4, ak5,
                     ak6, fx, _jac, e, 
                        _mas, ip, con, m1, m2, nm1, nfcn, 
                        njac, nstep, naccpt, nrejct, ndec, nsol);
            iwork[13] = nfcn;
            iwork[14] = njac;
            iwork[15] = nstep;
            iwork[16] = naccpt;
            iwork[17] = nrejct;
            iwork[18] = ndec;
            iwork[19] = nsol;
            /* ----------- RETURN ----------- */
            return 0;
        } /* rodas_ */


        /*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

        int roscor_(int n, S_fp fcn, double x, double[] y,
            double xend, double hmax, double h, double[]
            rtol, double[] atol, int itol, J_fp jac, int ijac,
            int mljac, int mujac, S_fp dfx, int idfx, M_fp mas,
            int mlmas, int mumas, P_fp solout, int iout, int
            idid, int nmax, double uround, int meth, int ijob,
            double fac1, double fac2, double safe, bool
            autnms, bool implct, bool banded, bool pred, int
            ldjac, int lde, int ldmas, double[] ynew, double[]
            dy1, double[] dy, double[] ak1, double[] ak2, double[]
            ak3, double[] ak4, double[] ak5, double[] ak6, double[]
            fx, double[] fjac, double[] e, double[] fmas, int[] ip,
            double[] cont, int m1, int m2, int nm1, int
            nfcn, int njac, int nstep, int naccpt, int nrejct,
            int ndec, int nsol)
        {
            var lu = new ReusableLU(n);

            /* System generated locals */
            double d__1;
            
            /* Local variables */
            int i, j, k, l;
            int j1;
            int md, mm;
            double sk, hd1 = 0, hd2 = 0, hd3 = 0, hd4 = 0;
            int nn2, nn3;
            double fac, hc21, hc31, hc32, hc41, hc42, hc43, hc51, hc52, hc53,
                hc54, hc61, hc62;
            int ier = 0, lrc;
            double hc63, hc64, hc65, err, hacc = 0;
            int lbeg, lend;
            double delt, hnew;
            bool last;
            double hopt = 0;

            double ysafe, hmaxn;
            int nsing;
            double xdelt;
            int irtrn;
            double erracc = 0;
            int mujacj;

            double facgus;
            bool reject;
            int mujacp;
            double posneg;


            /* ---------------------------------------------------------- */
            /*     CORE INTEGRATOR FOR RODAS */
            /*     PARAMETERS SAME AS IN RODAS WITH WORKSPACE ADDED */
            /* ---------------------------------------------------------- */
            /*         DECLARATIONS */
            /* ---------------------------------------------------------- */
            /* *** *** *** *** *** *** *** */
            /*  INITIALISATIONS */
            /* *** *** *** *** *** *** *** */
            /* Parameter adjustments */
            //--cont;
            //--fx;
            //--ak6;
            //--ak5;
            //--ak4;
            //--ak3;
            //--ak2;
            //--ak1;
            //--dy;
            //--dy1;
            //--ynew;
            //--y;
            //--rtol;
            //--atol;
            //fjac_offset = 1 + ldjac;
            //fjac -= 1 + ldjac;
            //--ip;
            //fmas_offset = 1 + ldmas;
            //fmas -= 1 + ldmas;
            //e_offset = 1 + lde;
            //e -= 1 + lde;
            //--rpar;
            //--ipar;

            /* Function Body */
            conros_.n = n;
            nn2 = n << 1;
            nn3 = n * 3;
            lrc = n << 2;
            /* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
            if (implct)
            {
                mas(nm1, fmas, ldmas);
            }
            /* ------ SET THE PARAMETERS OF THE METHOD ----- */
            rocoe_(meth);
            /* --- INITIAL PREPARATIONS */
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
            reject = false;
            last = false;
            nsing = 0;
            irtrn = 1;
            if (autnms)
            {
                hd1 = 0.0;
                hd2 = 0.0;
                hd3 = 0.0;
                hd4 = 0.0;
            }
            /* -------- PREPARE BAND-WIDTHS -------- */
            linal_1.mbdiag = mumas + 1;
            if (banded)
            {
                linal_1.mle = mljac;
                linal_1.mue = mujac;
                linal_1.mbjac = mljac + mujac + 1;
                linal_1.mbb = mlmas + mumas + 1;
                linal_1.mdiag = linal_1.mle + linal_1.mue + 1;
                linal_1.mdiff = linal_1.mle + linal_1.mue - mumas;
            }
            if (iout != 0)
            {
                conros_.xold = x;
                irtrn = 1;
                conros_.h = h;
                //solout(naccpt + 1, conros_.xold, x, y, cont, lrc, n, irtrn);
                if (irtrn < 0)
                {
                    goto L179;
                }
            }
            /* --- BASIC INTEGRATION STEP */
            L1:
            if (nstep > nmax)
            {
                goto L178;
            }
            if (Math.Abs(h) * 0.1 <= Math.Abs(x) * uround)
            {
                goto L177;
            }
            if (last)
            {
                h = hopt;

                idid = 1;
                return 0;
            }
            hopt = h;
            if ((x + h * 1.0001 - xend) * posneg >= 0.0)
            {

                h = xend - x;
                last = true;
            }
            /* *** *** *** *** *** *** *** */
            /*  COMPUTATION OF THE JACOBIAN */
            /* *** *** *** *** *** *** *** */
            fcn(n, x, y, dy1);
            ++(nfcn);
            ++(njac);
            if (ijac == 0)
            {
                /* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
                if (banded)
                {
                    /* --- JACOBIAN IS BANDED */
                    mujacp = mujac + 1;
                    md = Math.Min(linal_1.mbjac, n);
                    for (mm = 0; mm < m1 / m2 + 1; ++mm)
                    {
                        for (k = 0; k < md; ++k)
                        {
                            j = k + (mm - 1) * m2;
                            L12:
                            ak2[j] = y[j];
                            ak3[j] = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(y[j])));
                            y[j] += ak3[j];
                            j += md;
                            if (j <= mm * m2)
                            {
                                goto L12;
                            }
                            fcn(n, x, y, ak1);
                            j = k + (mm - 1) * m2;
                            j1 = k;
                            lbeg = Math.Max(1, j1 - mujac) + m1;
                            L14:
                            lend = Math.Min(m2, j1 + mljac) + m1;
                            y[j] = ak2[j];
                            mujacj = mujacp - j1 - m1;
                            for (l = lbeg - 1; l < lend; ++l)
                            {
                                fjac[l + mujacj + j * ldjac] = (ak1[l] - dy1[l]) / ak3[j];
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
                    for (i = 0; i < n; ++i)
                    {
                        ysafe = y[i];
                        delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                        y[i] = ysafe + delt;
                        fcn(n, x, y, ak1);
                        for (j = m1; j < n; ++j)
                        {
                            fjac[j - m1 + i * ldjac] = (ak1[j] - dy1[j]) / delt;
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
            if (!(autnms))
            {
                if (idfx == 0)
                {
                    /* --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X */
                    delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(x)));
                    xdelt = x + delt;
                    fcn(n, xdelt, y, ak1);
                    for (j = 0; j < n; ++j)
                    {
                        fx[j] = (ak1[j] - dy1[j]) / delt;
                    }
                }
                else
                {
                    /* --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X */
                    dfx(n, x, y, fx);
                }
            }
            L2:
            /* *** *** *** *** *** *** *** */
            /*  COMPUTE THE STAGES */
            /* *** *** *** *** *** *** *** */
            fac = 1.0 / (h * gamma);
            dc_decsol.decomr_(n, fjac/*[1 + ldjac]*/, ldjac, fmas/*[1 + ldmas]*/, ldmas, mlmas,
                mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip, ref ier, ijob,
                implct, ip);
            if (ier != 0)
            {
                goto L80;
            }
            ++(ndec);
            /* --- PREPARE FOR THE COMPUTATION OF THE 6 STAGES */
            hc21 = c21 / h;
            hc31 = c31 / h;
            hc32 = c32 / h;
            hc41 = c41 / h;
            hc42 = c42 / h;
            hc43 = c43 / h;
            hc51 = c51 / h;
            hc52 = c52 / h;
            hc53 = c53 / h;
            hc54 = c54 / h;
            hc61 = c61 / h;
            hc62 = c62 / h;
            hc63 = c63 / h;
            hc64 = c64 / h;
            hc65 = c65 / h;
            if (!(autnms))
            {
                hd1 = h * d1;
                hd2 = h * d2;
                hd3 = h * d3;
                hd4 = h * d4;
            }
            /* --- THE STAGES */
            dc_decsol.slvrod_(n, fjac/*[1 + ldjac]*/, ldjac, mljac, mujac, fmas/*[1 + ldmas]*/,
                ldmas, mlmas, mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip,
                dy1, ak1, fx, ynew, hd1, ijob, false);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a21 * ak1[i];
            }

            fcn(n, x + c2 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc21 * ak1[i];
            }
            dc_decsol.slvrod_(n, fjac/*[1 + ldjac]*/, ldjac, mljac, mujac, fmas/*[1 + ldmas]*/,
                ldmas, mlmas, mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip,
                dy, ak2, fx, ynew, hd2, ijob, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a31 * ak1[i] + a32 * ak2[i];
            }

            fcn(n, x + c3 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc31 * ak1[i] + hc32 * ak2[i];
            }
            dc_decsol.slvrod_(n, fjac/*[1 + ldjac]*/, ldjac, mljac, mujac, fmas/*[1 + ldmas]*/,
                ldmas, mlmas, mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip,
                dy, ak3, fx, ynew, hd3, ijob, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a41 * ak1[i] + a42 * ak2[i] + a43 * ak3[i];
            }

            fcn(n, x + c4 * h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ynew[i] = hc41 * ak1[i] + hc42 * ak2[i] + hc43 * ak3[i];
            }
            dc_decsol.slvrod_(n, fjac/*[1 + ldjac]*/, ldjac, mljac, mujac, fmas/*[1 + ldmas]*/,
                ldmas, mlmas, mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip,
                dy, ak4, fx, ynew, hd4, ijob, true);
            for (i = 0; i < n; i++)
            {
                ynew[i] = y[i] + a51 * ak1[i] + a52 * ak2[i] + a53 * ak3[i]
                    + a54 * ak4[i];
            }

            fcn(n, x + h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                ak6[i] = hc52 * ak2[i] + hc54 * ak4[i] + hc51 * ak1[i] + hc53
                 * ak3[i];
            }
            dc_decsol.slvrod_(n, fjac/*[1 + ldjac]*/, ldjac, mljac, mujac, fmas/*[1 + ldmas]*/,
                ldmas, mlmas, mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip,
                dy, ak5, fx, ak6, 0.0, ijob, true);

            /* ------------ EMBEDDED SOLUTION --------------- */
            for (i = 0; i < n; i++)
            {
                ynew[i] += ak5[i];
            }

            fcn(n, x + h, ynew, dy);
            for (i = 0; i < n; i++)
            {
                cont[i] = hc61 * ak1[i] + hc62 * ak2[i] + hc65 * ak5[i] +
                    hc64 * ak4[i] + hc63 * ak3[i];
            }
            dc_decsol.slvrod_(n, fjac/*[1 + ldjac]*/, ldjac, mljac, mujac, fmas/*[1 + ldmas]*/,
                ldmas, mlmas, mumas, m1, m2, nm1, fac, e/*[1 + lde]*/, lde, ip,
                dy, ak6, fx, cont, 0.0, ijob, true);

            /* ------------ NEW SOLUTION --------------- */
            for (i = 0; i < n; i++)
            {
                ynew[i] += ak6[i];
            }

            nsol += 6;

            nfcn += 5;
            /* ------------ DENSE OUTPUT ---------- */
            if (iout != 0)
            {
                for (i = 0; i < n; ++i)
                {
                    cont[i] = y[i];
                    cont[i + nn2] = d21 * ak1[i] + d22 * ak2[i] + d23 * ak3[i] + d24 * ak4[i] + d25 * ak5[i];
                    cont[i + nn3] = d31 * ak1[i] + d32 * ak2[i] + d33 * ak3[i] + d34 * ak4[i] + d35 * ak5[i];
                }
            }
            /* *** *** *** *** *** *** *** */
            /*  ERROR ESTIMATION */
            /* *** *** *** *** *** *** *** */
            ++(nstep);
            /* ------------ COMPUTE ERROR ESTIMATION ---------------- */
            err = 0.0;
            for (i = 0; i < n; i++)
            {
                if (itol == 0)
                {
                    sk = atol[0] + rtol[0] * Math.Max(Math.Abs(y[i]), Math.Abs(ynew[i]));
                }
                else
                {
                    sk = atol[i] + rtol[i] * Math.Max(Math.Abs(y[i]), Math.Abs(ynew[i]));
                }
                /* Computing 2nd power */
                d__1 = ak6[i] / sk;
                err += d__1 * d__1;
            }
            err = Math.Sqrt(err / n);
            /* --- COMPUTATION OF HNEW */
            /* --- WE REQUIRE .2<=HNEW/H<=6. */
            fac = Math.Max(fac2, Math.Min(fac1, Math.Pow(err, 0.25) / safe));
            hnew = h / fac;
            /* *** *** *** *** *** *** *** */
            /*  IS THE ERROR SMALL ENOUGH ? */
            /* *** *** *** *** *** *** *** */
            if (err <= 1.0)
            {
                /* --- STEP IS ACCEPTED */
                ++(naccpt);
                if (pred)
                {
                    /*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
                    if (naccpt > 1)
                    {
                        /* Computing 2nd power */
                        facgus = hacc / h * Math.Pow(err * err / erracc, 0.25) / safe;
                        facgus = Math.Max(fac2, Math.Min(fac1, facgus));
                        fac = Math.Max(fac, facgus);
                        hnew = h / fac;
                    }
                    hacc = h;
                    erracc = Math.Max(.01, err);
                }
                for (i = 0; i < n; ++i)
                {
                    y[i] = ynew[i];
                }
                conros_.xold = x;

                x += h;
                if (iout != 0)
                {
                    for (i = 0; i < n; ++i)
                    {
                        cont[conros_.n + i] = y[i];
                    }
                    irtrn = 1;
                    conros_.h = h;
                    //solout(naccpt + 1, conros_.xold, x, y, cont, lrc, n, irtrn);
                    if (irtrn < 0)
                    {
                        goto L179;
                    }
                }
                if (Math.Abs(hnew) > hmaxn)
                {
                    hnew = posneg * hmaxn;
                }
                if (reject)
                {
                    hnew = posneg * Math.Min(Math.Abs(hnew), Math.Abs(h));
                }
                reject = false;

                h = hnew;
                goto L1;
            }
            else
            {
                /* --- STEP IS REJECTED */
                reject = true;
                last = false;

                h = hnew;
                if (naccpt >= 1)
                {
                    ++(nrejct);
                }
                goto L2;
            }
            /* --- SINGULAR MATRIX */
            L80:
            ++nsing;
            if (nsing >= 5)
            {
                goto L176;
            }

            h *= 0.5;
            reject = true;
            last = false;
            goto L2;

            /* --- FAIL EXIT */
            L176:
            //do_fio("  EXIT OF RODAS AT X= ,e18.4", (x));
            Console.WriteLine(" MATRIX IS REPEATEDLY SINGULAR, IER=" + ier);
            return -4;

            L177:
            //do_fio("  EXIT OF RODAS AT X= ,e18.4", (x));
            Console.WriteLine(" STEP SIZE T0O SMALL, H=" + h);
            return -3;

            L178:
            //do_fio("  EXIT OF RODAS AT X= ,e18.4", (x));
            Console.WriteLine(" MORE THAN NMAX =" + nmax + "STEPS ARE NEEDED");
            return -2;

            /* --- EXIT CAUSED BY SOLOUT */
            L179:
            //do_fio("  EXIT OF RODAS AT X= ,e18.4", (x));
            return 2;
        } /* roscor_ */


        double contro_(int i, double x, double[] cont, int lrc)
        {
            /* System generated locals */
            double ret_val;

            /* Local variables */
            double s;

            /* ---------------------------------------------------------- */
            /*     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION */
            /*     WITH THE OUTPUT-SUBROUTINE FOR RODAS. IT PROVIDES AN */
            /*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
            /* ---------------------------------------------------------- */
            /* Parameter adjustments */
            //--cont;

            /* Function Body */
            s = (x - conros_.xold) / conros_.h;
            ret_val = cont[i] * (1 - s) + s * (cont[i + conros_.n] + (1 - s) * (cont[i + (conros_.n << 1)] + s * cont[i + conros_.n * 3]));
            return ret_val;
        } /* contro_ */


        /* Subroutine */
        int rocoe_(int meth)
        {
            double bet2p, bet3p, bet4p;

            switch (meth)
            {
                case 1:
                    c2 = 0.386;
                    c3 = 0.21;
                    c4 = 0.63;
                    bet2p = 0.0317;
                    bet3p = 0.0635;
                    bet4p = 0.3438;
                    d1 = 0.25;
                    d2 = -0.1043;
                    d3 = 0.1035;
                    d4 = -0.03620000000000023;
                    a21 = 1.544;
                    a31 = 0.9466785280815826;
                    a32 = 0.2557011698983284;
                    a41 = 3.314825187068521;
                    a42 = 2.896124015972201;
                    a43 = 0.9986419139977817;
                    a51 = 1.221224509226641;
                    a52 = 6.019134481288629;
                    a53 = 12.53708332932087;
                    a54 = -0.687886036105895;
                    c21 = -5.6688;
                    c31 = -2.430093356833875;
                    c32 = -0.2063599157091915;
                    c41 = -0.1073529058151375;
                    c42 = -9.594562251023355;
                    c43 = -20.47028614809616;
                    c51 = 7.496443313967647;
                    c52 = -10.24680431464352;
                    c53 = -33.99990352819905;
                    c54 = 11.7089089320616;
                    c61 = 8.083246795921522;
                    c62 = -7.981132988064893;
                    c63 = -31.52159432874371;
                    c64 = 16.31930543123136;
                    c65 = -6.058818238834054;
                    gamma = 0.25;
                    d21 = 10.12623508344586;
                    d22 = -7.487995877610167;
                    d23 = -34.80091861555747;
                    d24 = -7.992771707568823;
                    d25 = 1.025137723295662;
                    d31 = -0.6762803392801253;
                    d32 = 6.087714651680015;
                    d33 = 16.43084320892478;
                    d34 = 24.76722511418386;
                    d35 = -6.594389125716872;
                    break;
                case 2:
                    c2 = 0.3507221;
                    c3 = 0.2557041;
                    c4 = 0.681779;
                    bet2p = 0.0317;
                    bet3p = 0.0047369;
                    bet4p = 0.3438;
                    d1 = 0.25;
                    d2 = -0.06902209999999998;
                    d3 = -9.671999999999459e-4;
                    d4 = -0.08797900000000025;
                    a21 = 1.4028884;
                    a31 = 0.6581212688557198;
                    a32 = -1.320936088384301;
                    a41 = 7.131197445744498;
                    a42 = 16.02964143958207;
                    a43 = -5.561572550509766;
                    a51 = 22.73885722420363;
                    a52 = 67.38147284535289;
                    a53 = -31.2187749303856;
                    a54 = 0.7285641833203814;
                    c21 = -5.1043536;
                    c31 = -2.899967805418783;
                    c32 = 4.040399359702244;
                    c41 = -32.64449927841361;
                    c42 = -99.35311008728094;
                    c43 = 49.99119122405989;
                    c51 = -76.46023087151691;
                    c52 = -278.5942120829058;
                    c53 = 153.9294840910643;
                    c54 = 10.97101866258358;
                    c61 = -76.29701586804983;
                    c62 = -294.2795630511232;
                    c63 = 162.0029695867566;
                    c64 = 23.6516690309527;
                    c65 = -7.652977706771382;
                    gamma = 0.25;
                    d21 = -38.71940424117216;
                    d22 = -135.8025833007622;
                    d23 = 64.51068857505875;
                    d24 = -4.192663174613162;
                    d25 = -2.53193205033506;
                    d31 = -14.99268484949843;
                    d32 = -76.30242396627033;
                    d33 = 58.65928432851416;
                    d34 = 16.61359034616402;
                    d35 = -0.6758691794084156;
                    break;
                case 3:
                    /* Coefficients for RODAS with order 4 for linear parabolic problems */
                    /* Gerd Steinebach (1993) */
                    gamma = 0.25;
                    c2 = gamma * 3.0;
                    c3 = 0.21;
                    c4 = 0.63;
                    bet2p = 0.0;
                    bet3p = c3 * c3 * (c3 / 6.0 - gamma / 2.0) / (gamma * gamma);
                    bet4p = 0.3438;
                    d1 = 0.25;
                    d2 = -0.5;
                    d3 = -0.023504;
                    d4 = -0.0362;
                    a21 = 3.0;
                    a31 = 1.831036793486759;
                    a32 = 0.4955183967433795;
                    a41 = 2.304376582692669;
                    a42 = -0.05249275245743001;
                    a43 = -1.176798761832782;
                    a51 = -7.170454962423024;
                    a52 = -4.741636671481785;
                    a53 = -16.31002631330971;
                    a54 = -1.062004044111401;
                    c21 = -12.0;
                    c31 = -8.791795173947035;
                    c32 = -2.207865586973518;
                    c41 = 10.81793056857153;
                    c42 = 6.780270611428266;
                    c43 = 19.5348594464241;
                    c51 = 34.19095006749676;
                    c52 = 15.49671153725963;
                    c53 = 54.7476087596413;
                    c54 = 14.16005392148534;
                    c61 = 34.62605830930532;
                    c62 = 15.30084976114473;
                    c63 = 56.99955578662667;
                    c64 = 18.40807009793095;
                    c65 = -5.714285714285717;

                    d21 = 25.09876703708589;
                    d22 = 11.62013104361867;
                    d23 = 28.49148307714626;
                    d24 = -5.664021568594133;
                    d25 = 0.0;
                    d31 = 1.638054557396973;
                    d32 = -0.7373619806678748;
                    d33 = 8.47791821923899;
                    d34 = 15.9925314877952;
                    d35 = -1.882352941176471;
                    break;
            }
            return 0;
        }

    }
}