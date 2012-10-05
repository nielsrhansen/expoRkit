* Subroutine wrappers for the Expokit subroutines to be called from R. 
* The primary purpose of the wrapper is to select the subroutine for 
* the matrix-vector product appropriate for the storage format of the 
* sparse matrix. The required information is passed from R via the 
* sflag argument. 
*
*     Copyright (C) 2011 Niels Richard Hansen.
*
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the
* Free Software Foundation; either version 2, or (at your option) any
* later version.
*
* The subroutine is distributed in the hope that they will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, a copy is available at
* http://www.r-project.org/Licenses/
*
* --- General arguments to all the subroutines ------------------------*
*
*     n      : (input) order of the principal matrix A.
*
*     ia, ja, a (input): 
*           sparse matrix data stored in an appropriate
*           format. Currently CCS, CRS and COO supported with
*           implemented external matrix-vector product functions.
*
*     nz     : (input) the actual number of non zero entries.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*                      
*     v(n)   : (input) given operand vector.
*                      
*     w(n) : (output) computed approximation of exp(t*A)*v or 
*            exp(t*A)*v + t*phi(t*A)*u.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     itrace : (input) running mode. 0=silent, 1=print happy breakdown,
*              2>=print step-by-step info.  
*                   
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*     sflag  :  (input) storage flag. Selects the matrix-vector routine.
*               1 - Compressed Column Storage (CCS)
*               2 - Compressed Row Storage (CRS)
*               3 - COOrdinates storage format

      subroutine R_DGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,
     .     anorm, itrace,iflag,sflag )
      implicit none
      integer n, nz, ia(*), ja(*), m, lwsp, iflag, iwsp(m+2), mxstep,
     .     itrace, sflag
      double precision a(*), t, tol, anorm, v(n), w(n),
     .     wsp( n*(m+2)+5*(m+2)*(m+2)+6+1 )
      external dgcoov, dgccsv, dgcrsv

      lwsp = n*(m+2)+5*(m+2)*(m+2)+6+1 
      if ( sflag.eq.1) then
         call DGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, dgccsv, itrace,iflag )
      endif
      if ( sflag.eq.2) then
         call DGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, dgcrsv, itrace,iflag )
      endif
      if ( sflag.eq.3) then
         call DGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, dgcoov, itrace,iflag )
      endif
      END

      subroutine R_DMEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,
     .     anorm, itrace,iflag,sflag )
      implicit none
      integer n, nz, ia(*), ja(*), m, lwsp, iflag, iwsp(m+2), mxstep,
     .     itrace, sflag
      double precision a(*), t, tol, anorm, v(n), w(n),
     .     wsp( n*(m+2)+5*(m+2)*(m+2)+6+1  )
      external dgcoov, dgccsv, dgcrsv
      
      lwsp = n*(m+2)+5*(m+2)*(m+2)+6+1 
      if ( sflag.eq.1) then
         call DMEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, dgccsv, itrace,iflag )
      endif
      if ( sflag.eq.2) then
         call DMEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, dgcrsv, itrace,iflag )
      endif
      if ( sflag.eq.3) then
         call DMEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, dgcoov, itrace,iflag )
      endif
      END
      
      subroutine R_ZGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,
     .     anorm, itrace,iflag,sflag )
      implicit none
      integer n, nz, ia(*), ja(*), m, lwsp, iflag, iwsp(m+2), mxstep,
     .     itrace, sflag
      double precision t, tol, anorm
      complex*16 a(*), v(n), w(n), 
     .     wsp( n*(m+2)+5*(m+2)*(m+2)+6+1  )
      external zgcoov, zgccsv, zgcrsv

      lwsp = n*(m+2)+5*(m+2)*(m+2)+6+1 
      if ( sflag.eq.1) then
         call ZGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, zgccsv, itrace,iflag )
      endif
      if ( sflag.eq.2) then
         call ZGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, zgcrsv, itrace,iflag )
      endif
      if ( sflag.eq.3) then
         call ZGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+2, zgcoov, itrace,iflag )
      endif
      END

      subroutine R_DGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,
     .     anorm, itrace,iflag,sflag )
      implicit none
      integer n, nz, ia(*), ja(*), m, lwsp, iflag, iwsp(m+3), mxstep, 
     .     itrace, sflag
      double precision a(*), t, tol, anorm, u(n), v(n), w(n),
     .     wsp(n*(m+3)+5*(m+3)*(m+3)+6+1  )
      external dgcoov, dgccsv, dgcrsv

      lwsp = n*(m+3)+5*(m+3)*(m+3)+6+1 
      if ( sflag.eq.1) then
         call DGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+3, dgccsv, itrace,iflag )
      endif
      if ( sflag.eq.2) then
         call DGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+3, dgcrsv, itrace,iflag )
      endif
      if ( sflag.eq.3) then
         call DGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+3, dgcoov, itrace,iflag )
      endif
      END

      subroutine R_ZGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,
     .     anorm, itrace,iflag,sflag )
      implicit none
      integer n, nz, ia(*), ja(*), m, lwsp, iflag, iwsp(m+3), mxstep,
     .     itrace, sflag
      double precision t, tol, anorm
      complex*16 a(*), u(n), v(n), w(n),
     .     wsp( n*(m+3)+5*(m+3)*(m+3)+6+1 )
      external zgcoov, zgccsv, zgcrsv

      lwsp = n*(m+3)+5*(m+3)*(m+3)+6+1 
      if ( sflag.eq.1) then
         call ZGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+3, zgccsv, itrace,iflag )
      endif
      if ( sflag.eq.2) then
         call ZGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+3, zgcrsv, itrace,iflag )
      endif
      if ( sflag.eq.3) then
         call ZGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,anorm,
     .        wsp,lwsp, iwsp,m+3, zgcoov, itrace,iflag )
      endif
      END
