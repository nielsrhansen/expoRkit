*----------------------------------------------------------------------*
*     Subroutines adapted from the Expokit Fortran package for         *
*     computing the matrix exponential. See                            *
*     http://www.maths.uq.edu.au/expokit/                              *
*                                                                      *
*     Copyright (C) 1996-2006 Roger B. Sidje <rbs@maths.uq.edu.au>     *
*                                                                      *
*     The subroutines were adapted with the purpose of calling the     *
*     routines from R. The main adaptation is in the way a sparse      *
*     matrix is passed between subroutines where Expokit relies        *
*     on a common block statement. This is modified so that sparse     *
*     matrices are passed as arguments.                                *
*                                                                      *
*     This program is free software; you can redistribute it and/or    *
*     modify it under the terms of the GNU General Public License as   *
*     published by the Free Software Foundation; either version 2, or  *
*     (at your option) any later version.                              *
*                                                                      *
*     The subroutine is distributed in the hope that they will be      *
*     useful, but WITHOUT ANY WARRANTY; without even the implied       *
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *	 
*     See the GNU General Public License for more details.             *
*                                                                      *
*     You should have received a copy of the GNU General Public        *
*     License along with this program; if not, a copy is available     *
*     at http://www.r-project.org/Licenses/                            *
*                                                                      *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     The first section of subroutines is from the Expokit mataid.f    *
*     file with matrix-vector support routines. These have been        *
*     adapted to support passing the sparse matrix via arguments       *
*     instead of via a common block. The following note from the       *
*     mataid.f files still applies.                                    *
*----------------------------------------------------------------------*

*-------------------------------NOTE-----------------------------------*
*     This is an accessory to Expokit and it is not intended to be     *
*     complete. It is supplied primarily to ensure an unconstrained    *
*     distribution and portability of the package. The matrix-vector   *
*     multiplication routines supplied here fit the non symmetric      *
*     storage and for a symmetric matrix, the entire (not half) matrix *
*     is required.  If the sparsity pattern is known a priori, it is   *
*     recommended to use the most advantageous format and to devise    *
*     the most advantageous matrix-vector multiplication routine.      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine dgcoov ( x, y, a, ia, ja, n, nz )
      implicit none
      integer i, j, n, nz, ia(nz), ja(nz)
      double precision x(n), y(n), a(nz)
*
*---  Computes y = A*x. 
*---  A is in the COOrdinates storage format.
*---  In the original Expokit code, A is stored in a 'common block'.
*---  For array argument passing of the matrix from R this was rewritten.
*
      do j = 1,n
         y(j) = 0.0d0
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine dgcrsv ( x, y, a, ia, ja, n, nz )
      implicit none
      integer i, j, n, nz, ia(n+1), ja(nz)
      double precision x(n), y(n), a(nz)
*
*---  Computes y = A*x. 
*---  A is in the Compressed Row Storage (CRS) format.
*---  In the original Expokit code, A is stored in a 'common block'.
*---  For array argument passing of the matrix from R this was rewritten.
*

      do i = 1,n
         y(i) = 0.0d0
         do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
         enddo
      enddo
      END
*----------------------------------------------------------------------|
      subroutine dgccsv( x, y, a, ia, ja, n, nz )
      implicit none
      integer i, j, n, nz, ia(nz), ja(n+1)
      double precision x(n), y(n), a(nz)
*
*---  Computes y = A*x. 
*---  A is in the Compressed Column Storage (CCS) format.
*---  In the original Expokit code, A is stored in a 'common block'.
*---  For array argument passing of the matrix from R this was rewritten.
*
      do i = 1,n
         y(i) = 0.0d0
      enddo
      do j = 1,n
         do i = ja(j),ja(j+1)-1
            y(ia(i)) = y(ia(i)) + a(i)*x(j)
         enddo
      enddo
      end
*----------------------------------------------------------------------|
      subroutine zgcoov ( x, y, a, ia, ja, n, nz )
      implicit none
      integer i, j, n, nz, ia(nz), ja(nz)
      complex*16 x(n), y(n), a(nz)
*
*---  Computes y = A*x. Complex version.
*---  A is in the COOrdinates storage format.
*---  In the original Expokit code, A is stored in a 'common block'.
*---  For array argument passing of the matrix from R this was rewritten.
*
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )
 
      do j = 1,n
         y(j) = ZERO
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END
*----------------------------------------------------------------------|
      subroutine zgccsv( x, y, a, ia, ja, n, nz )
      implicit none
      integer i, j, n, nz, ia(nz), ja(n+1)
      complex*16 x(n), y(n), a(nz)
*
*---  Computes y = A*x. Complex version.
*---  A is in the Compressed Column Storage (CCS) format.
*---  In the original Expokit code, A is stored in a 'common block'.
*---  For array argument passing of the matrix from R this was rewritten.
*
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )
      do i = 1,n
         y(i) = ZERO
      enddo
      do j = 1,n
         do i = ja(j),ja(j+1)-1
            y(ia(i)) = y(ia(i)) + a(i)*x(j)
         enddo
      enddo
      end

*----------------------------------------------------------------------|
      subroutine zgcrsv ( x, y, a, ia, ja, n, nz )
      implicit none
      integer i, j, n, nz, ia(n+1), ja(nz)
      complex*16 x(n), y(n), a(nz)
*
*---  Computes y = A*x. Complex version.
*---  A is in the Compressed Row Storage (CRS) format.
*---  In the original Expokit code, A is stored in a 'common block'.
*---  For array argument passing of the matrix from R this was rewritten.
*
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )

      do i = 1,n
         y(i) = ZERO
         do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
         enddo
      enddo
      END
*----------------------------------------------------------------------|

*----------------------------------------------------------------------*
*     The second section contains the different routines for computing *
*     the matrix exponential. They depend on the external 'matvec'     *
*     function for the matrix-vector product (e.g. one of the          *
*     functions above).                                                *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------|
      subroutine DMEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep, 
     .     anorm, wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .     itrace,iflag )

      implicit none
      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp), mxstep,
     .     nz, ia(*), ja(*)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp), a(*),
     .     delta, gamma
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DMEXPV computes w = exp(t*A)*v - Customised for MARKOV CHAINS.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts 
*     via the external routine `matvec' performing the matrix-vector 
*     product.
*
*     This is a customised version for Markov Chains. This means that a
*     check is done within this code to ensure that the resulting vector 
*     w is a probability vector, i.e., w must have all its components 
*     in [0,1], with sum equal to 1. This check is done at some expense
*     and the user may try DGEXPV which is cheaper since it ignores 
*     probability constraints.
*
*     IMPORTANT: The check assumes that the transition rate matrix Q
*                satisfies Qe = 0, where e=(1,...,1)'. Don't use DMEXPV
*                if this condition does not hold. Use DGEXPV instead.
*                DMEXPV/DGEXPV require the matrix-vector product 
*                y = A*x = Q'*x, i.e, the TRANSPOSE of Q times a vector.
*                Failure to remember this leads to wrong results.
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
*
*     nz     : (input) the actual number of non zero entries.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested acurracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H     wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.  
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*     Markov chains are usually well-conditioned problems.
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1,
     .                 roundoff, s_round, x_round

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2, DASUM

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      sgn      = SIGN( 1.0d0,t )
      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      s_round  = 0.0d0
      x_error  = 0.0d0
      x_round  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
               call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ..
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  Markov model constraints ...
*
      j = 0
      do i = 1,n
         if ( w(i).lt.0.0d0 ) then
            w(i) = 0.0d0
            j = j + 1
         endif
      enddo
      p1 = DASUM( n, w,1 )
      if ( j.gt.0 ) call DSCAL( n, 1.0d0/p1, w,1 )
      roundoff = DABS( 1.0d0-p1 ) / DBLE( n )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc, roundoff )
      err_loc = MAX( err_loc, rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('wnorm', -1, beta, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('roundoff', -1, roundoff, 1)
         call dblepr('next_step', -1, t_new, 1)        
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      s_round = s_round + roundoff
      x_error = MAX( x_error, err_loc )
      x_round = MAX( x_round, roundoff )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = x_round
      wsp(4)  = s_round
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

      implicit none
      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a general matrix in
*     full, using the irreducible rational Pade approximation to the 
*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degree of the diagonal Pade to be used.
*                 A value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*                  
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr), 
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) goto 600
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = 0.0d0
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      enddo
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) then
*--   'Error - null H in input of DGPADM.'
         iflag = 3
         goto 600
      endif
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
*
*---  H2 = scale2*H*H ...
*
      call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) goto 600
      call DSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

      implicit none
      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a symmetric matrix
*     in full, using the irreducible rational Pade approximation to the 
*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix (both lower and upper parts).
*
*     t         : (input) time-scale (can be < 0).
*                  
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr), 
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) goto 600 
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = 0.0d0
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      enddo
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) then
*--   'Error - null H in input of DSPADM.'
         iflag = 3
         goto 600
      endif
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
*
*---  H2 = scale2*H*H ...
*
      call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DSYSV( 'U',m,m,wsp(iq),m,ipiv,wsp(ip),m,wsp(ih2),mm,iflag )
      if ( iflag.ne.0 ) goto 600
      call DSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine ZGPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

      implicit none
      double precision t
      integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      complex*16       H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a general complex 
*     matrix in full, using the irreducible rational Pade approximation
*     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*                  
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr), 
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                       0 - no problem
*                      <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
      double precision hnorm
      complex*16 cp, cq, scale, scale2, ZERO, ONE

      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) goto 600
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = ZERO
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,DBLE(wsp(i)) )
      enddo
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) then
*--   'Error - null H in input of ZGPADM.'
         iflag = 3
         goto 600
      endif
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = ONE
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
*
*---  H2 = scale2*H*H ...
*
      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
*
*---  initialise p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = ZERO
            wsp(iq + (j-1)*m + i-1) = ZERO
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,
     .             wsp(ih2),m, ZERO,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd.ne.0 ) then
         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, ZERO,wsp(ifree),m )
         iq = ifree
      else
         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, ZERO,wsp(ifree),m )
         ip = ifree
      endif
      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
      call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) goto 600
      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.ne.0 ) then
         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,
     .                ZERO,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
 600  END
c$$$ ZHPADM is currently not accessible from R. It relies on ZHESV,
c$$$ which is not available from the R distribution of LAPACK. 
c$$$ *----------------------------------------------------------------------|
c$$$ subroutine ZHPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
c$$$
c$$$      implicit none
c$$$      double precision t
c$$$      integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
c$$$      complex*16       H(ldh,m), wsp(lwsp)
c$$$
c$$$*-----Purpose----------------------------------------------------------|
c$$$*
c$$$*     Computes exp(t*H), the matrix exponential of an Hermitian matrix
c$$$*     in full, using the irreducible rational Pade approximation to the
c$$$*     exponential function exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
c$$$*     combined with scaling-and-squaring.
c$$$*
c$$$*-----Arguments--------------------------------------------------------|
c$$$*
c$$$*     ideg      : (input) the degre of the diagonal Pade to be used.
c$$$*                 a value of 6 is generally satisfactory.
c$$$*
c$$$*     m         : (input) order of H.
c$$$*
c$$$*     H(ldh,m)  : (input) argument matrix (both lower and upper parts).
c$$$*
c$$$*     t         : (input) time-scale (can be < 0).
c$$$*                  
c$$$*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
c$$$*
c$$$*     ipiv(m)   : (workspace)
c$$$*
c$$$*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
c$$$*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
c$$$*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c$$$*                 NOTE: if the routine was called with wsp(iptr), 
c$$$*                       then exp(tH) will start at wsp(iptr+iexph-1).
c$$$*
c$$$*     ns        : (output) number of scaling-squaring used.
c$$$*
c$$$*     iflag     : (output) exit flag.
c$$$*                       0 - no problem
c$$$*                      <0 - problem
c$$$*
c$$$*----------------------------------------------------------------------|
c$$$*     Roger B. Sidje (rbs@maths.uq.edu.au)
c$$$*     EXPOKIT: Software Package for Computing Matrix Exponentials.
c$$$*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
c$$$*----------------------------------------------------------------------|
c$$$*
c$$$      integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
c$$$      double precision hnorm
c$$$      complex*16 cp, cq, scale, scale2, ZERO, ONE
c$$$
c$$$      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
c$$$      intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX
c$$$
c$$$*---  check restrictions on input parameters ...
c$$$      mm = m*m
c$$$      iflag = 0
c$$$      if ( ldh.lt.m ) iflag = -1
c$$$      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
c$$$      if ( iflag.ne.0 ) goto 600
c$$$*
c$$$*---  initialise pointers ...
c$$$*
c$$$      icoef = 1
c$$$      ih2 = icoef + (ideg+1)
c$$$      ip  = ih2 + mm
c$$$      iq  = ip + mm
c$$$      ifree = iq + mm
c$$$*
c$$$*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
c$$$*     and set scale = t/2^ns ...
c$$$*
c$$$      do i = 1,m
c$$$         wsp(i) = ZERO
c$$$      enddo
c$$$      do j = 1,m
c$$$         do i = 1,m
c$$$            wsp(i) = wsp(i) + ABS( H(i,j) )
c$$$         enddo
c$$$      enddo
c$$$      hnorm = 0.0d0
c$$$      do i = 1,m
c$$$         hnorm = MAX( hnorm,DBLE(wsp(i)) )
c$$$      enddo
c$$$      hnorm = ABS( t*hnorm )
c$$$      if ( hnorm.eq.0.0d0 ) then
c$$$*---  'Error - null H in input of ZHPADM.'
c$$$         iflag = 3
c$$$         goto 600
c$$$      endif
c$$$      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
c$$$      scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
c$$$      scale2 = scale*scale
c$$$*
c$$$*---  compute Pade coefficients ...
c$$$*
c$$$      i = ideg+1
c$$$      j = 2*ideg+1
c$$$      wsp(icoef) = ONE
c$$$      do k = 1,ideg
c$$$         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
c$$$      enddo
c$$$*
c$$$*---  H2 = scale2*H*H ...
c$$$*
c$$$      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
c$$$*
c$$$*---  initialise p (numerator) and q (denominator) ...
c$$$*
c$$$      cp = wsp(icoef+ideg-1)
c$$$      cq = wsp(icoef+ideg)
c$$$      do j = 1,m
c$$$         do i = 1,m
c$$$            wsp(ip + (j-1)*m + i-1) = ZERO
c$$$            wsp(iq + (j-1)*m + i-1) = ZERO
c$$$         enddo
c$$$         wsp(ip + (j-1)*(m+1)) = cp
c$$$         wsp(iq + (j-1)*(m+1)) = cq
c$$$      enddo
c$$$*
c$$$*---  Apply Horner rule ...
c$$$*
c$$$      iodd = 1
c$$$      k = ideg - 1
c$$$ 100  continue
c$$$      iused = iodd*iq + (1-iodd)*ip
c$$$      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,
c$$$     .             wsp(ih2),m, ZERO,wsp(ifree),m )
c$$$      do j = 1,m
c$$$         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
c$$$      enddo
c$$$      ip = (1-iodd)*ifree + iodd*ip
c$$$      iq = iodd*ifree + (1-iodd)*iq
c$$$      ifree = iused
c$$$      iodd = 1-iodd
c$$$      k = k-1
c$$$      if ( k.gt.0 )  goto 100
c$$$*
c$$$*---  Obtain (+/-)(I + 2*(p\q)) ...
c$$$*
c$$$      if ( iodd.ne.0 ) then
c$$$         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
c$$$     .                H,ldh, ZERO,wsp(ifree),m )
c$$$         iq = ifree
c$$$      else
c$$$         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
c$$$     .                H,ldh, ZERO,wsp(ifree),m )
c$$$         ip = ifree
c$$$      endif
c$$$      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
c$$$      call ZHESV( 'U',m,m,wsp(iq),m,ipiv,wsp(ip),m,wsp(ih2),mm,iflag )
c$$$      if ( iflag.ne.0 ) goto 600
c$$$      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
c$$$      do j = 1,m
c$$$         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
c$$$      enddo
c$$$      iput = ip
c$$$      if ( ns.eq.0 .and. iodd.ne.0 ) then
c$$$         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
c$$$         goto 200
c$$$      endif
c$$$*
c$$$*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
c$$$*
c$$$      iodd = 1
c$$$      do k = 1,ns
c$$$         iget = iodd*ip + (1-iodd)*iq
c$$$         iput = (1-iodd)*ip + iodd*iq
c$$$         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,
c$$$     .                ZERO,wsp(iput),m )
c$$$         iodd = 1-iodd
c$$$      enddo
c$$$ 200  continue
c$$$      iexph = iput
c$$$ 600  END
*----------------------------------------------------------------------|
      subroutine DGCHBV( m, t, H,ldh, y, wsp, iwsp, iflag )

      implicit none
      integer          m, ldh, iflag, iwsp(m)
      double precision t, H(ldh,m), y(m)
      complex*16       wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DGCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is a General matrix.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) argument matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg)

      intrinsic DBLE
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         enddo
         call ZGESV( M, 1, WSP(iH),M, IWSP, WSP(iY),M, IFLAG )
         if ( IFLAG.ne.0 ) goto 600
*---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         enddo
      enddo
 600  END
*----------------------------------------------------------------------|
c$$$ DSCHBV is currently not accessible from R. It relies on ZSYSV,
c$$$ which is not available from the R distribution of LAPACK. 
c$$$*----------------------------------------------------------------------|
c$$$      subroutine DSCHBV( m, t, H,ldh, y, wsp, iwsp, iflag )
c$$$
c$$$      implicit none
c$$$      integer          m, ldh, iflag, iwsp(m)
c$$$      double precision t, H(ldh,m), y(m)
c$$$      complex*16       wsp(m*(m+2))
c$$$
c$$$*-----Purpose----------------------------------------------------------|
c$$$*
c$$$*---  DSCHBV computes y = exp(t*H)*y using the partial fraction
c$$$*     expansion of the uniform rational Chebyshev approximation
c$$$*     to exp(-x) of type (14,14). H is assumed to be symmetric.
c$$$*     About 14-digit accuracy is expected if the matrix H is negative
c$$$*     definite. The algorithm may behave poorly otherwise. 
c$$$*
c$$$*-----Arguments--------------------------------------------------------|
c$$$*
c$$$*     m       : (input) order of matrix H
c$$$*
c$$$*     t       : (input) time-scaling factor (can be < 0).
c$$$*
c$$$*     H(ldh,m): (input) symmetric matrix.
c$$$*
c$$$*     y(m)    : (input/output) on input the operand vector,
c$$$*               on output the resulting vector exp(t*H)*y.
c$$$*
c$$$*     iwsp(m) : (workspace)
c$$$*
c$$$*     wsp     : (workspace). Observe that a double precision vector of
c$$$*               length 2*m*(m+2) can be used as well when calling this
c$$$*               routine (thus avoiding an idle complex array elsewhere)
c$$$*
c$$$*----------------------------------------------------------------------|
c$$$*     Roger B. Sidje (rbs@maths.uq.edu.au)
c$$$*     EXPOKIT: Software Package for Computing Matrix Exponentials.
c$$$*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
c$$$*----------------------------------------------------------------------|
c$$$*
c$$$      integer ndeg, i, j, ip, ih, iy, iz
c$$$      parameter ( ndeg=7 )
c$$$      double precision alpha0
c$$$      complex*16 alpha(ndeg), theta(ndeg), w
c$$$
c$$$      intrinsic ABS,CMPLX,DBLE,MIN
c$$$      
c$$$*---  Pointers ...
c$$$
c$$$      ih = 1
c$$$      iy = ih + m*m
c$$$      iz = iy + m
c$$$
c$$$*---  Coefficients and poles of the partial fraction expansion ...
c$$$
c$$$      alpha0  =  0.183216998528140087D-11
c$$$      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
c$$$      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
c$$$      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
c$$$      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
c$$$      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
c$$$      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
c$$$      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)
c$$$
c$$$      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
c$$$      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
c$$$      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
c$$$      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
c$$$      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
c$$$      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
c$$$      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
c$$$*     
c$$$*---  Accumulation of the contribution of each pole ...
c$$$*
c$$$      do j = 1,m
c$$$         wsp(iz+j-1) = y(j)
c$$$         y(j) = y(j)*alpha0
c$$$      enddo
c$$$      do ip = 1,ndeg
c$$$*---     Solve each fraction using Gaussian elimination with pivoting...
c$$$         do j = 1,m
c$$$            do i = 1,m
c$$$               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
c$$$            enddo
c$$$            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
c$$$            wsp(iy+j-1) = wsp(iz+j-1)
c$$$         enddo
c$$$         call ZSYSV('U', M, 1, WSP(iH),M, IWSP, WSP(iY),M, W,1, IFLAG )
c$$$         if ( IFLAG.ne.0 ) goto 600
c$$$*---     Accumulate the partial result in y ...     
c$$$         do i = 1,m
c$$$            y(i) = y(i) + DBLE( alpha(ip)*wsp(iy+i-1) )
c$$$         enddo
c$$$      enddo
c$$$ 600  END
c$$$*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine ZGCHBV( m, t, H,ldh, y, wsp, iwsp, iflag )

      implicit none
      integer          m, ldh, iflag, iwsp(m)
      double precision t
      complex*16       H(ldh,m), y(m), wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  ZGCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is a General matrix.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the matrix H.
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) argument matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine.
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer     ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double      precision alpha0
      complex*16  alpha(2*ndeg), theta(2*ndeg)
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
      do ip = 1,ndeg
         theta(ndeg+ip) = CONJG( theta(ip) )
         alpha(ndeg+ip) = CONJG( alpha(ip) )
      enddo
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,2*ndeg
         alpha(ip) = 0.5d0*alpha(ip)
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         enddo
         call ZGESV( M, 1, WSP(iH),M, IWSP, WSP(iY),M, IFLAG )
         if ( IFLAG.ne.0 ) goto 600
*---     Accumulate the partial result in y ...     
         do i = 1,m
            y(i) = y(i) + alpha(ip)*wsp(iy+i-1)
         enddo
      enddo
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DNCHBV( m, t, H,ldh, y, wsp )

      implicit none
      integer          m, ldh
      double precision t, H(ldh,m), y(m), wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DNCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO
      integer ndeg, i, j, k, ip, ih, iy, iz
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), tmpc

      intrinsic ABS,DBLE,MIN
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            enddo
         enddo
         do i = 1,m-1
*---        Get pivot and exchange rows ...
            if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
               call ZSWAP( m-i+1, wsp(ih+(i-1)*m+i-1),m, 
     .                     wsp(ih+(i-1)*m+i),m )
               call ZSWAP( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
            endif
*---        Forward eliminiation ... 
            tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
            call ZAXPY( m-i, -tmpc, wsp(ih+i*m+i-1),m, wsp(ih+i*m+i),m )
            wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
         enddo
*---     Backward substitution ...    
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            enddo
            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
         enddo
*---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine ZNCHBV( m, t, H,ldh, y, wsp )

      implicit none
      integer          m, ldh
      double precision t
      complex*16       H(ldh,m), y(m), wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  ZNCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine.
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO
      integer ndeg, i, j, k, ip, ih, iy, iz
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), tmpc

      intrinsic ABS,DBLE,CONJG,MIN
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
      do ip = 1,ndeg
         theta(ndeg+ip) = CONJG( theta(ip) )
         alpha(ndeg+ip) = CONJG( alpha(ip) )
      enddo
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,2*ndeg
         alpha(ip) = 0.5d0*alpha(ip)
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            enddo
         enddo
         do i = 1,m-1
*---        Get pivot and exchange rows ...
            if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
               call ZSWAP( m-i+1, wsp(ih+(i-1)*m+i-1),m, 
     .                     wsp(ih+(i-1)*m+i),m )
               call ZSWAP( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
            endif
*---        Forward eliminiation ... 
            tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
            call ZAXPY( m-i, -tmpc, wsp(ih+i*m+i-1),m, wsp(ih+i*m+i),m )
            wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
         enddo
*---     Backward substitution ...    
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            enddo
            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
         enddo
*---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + alpha(ip)*wsp(iy+j-1)
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .                   itrace,iflag )

      implicit none
      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp), mxstep,
     .     nz, ia(*), ja(*)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp), a(*), 
     .     delta, gamma
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DGEXPV computes w = exp(t*A)*v - for a General matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts 
*     via the external routine `matvec' performing the matrix-vector 
*     product.
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the ccs format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
*
*     nz     : (input) the actual number of non zero entries.
*                               
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*                      
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in wsp and iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )
*
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif


* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .                   itrace,iflag )

      implicit none
      integer n, nz, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp), mxstep,
     .     ia(*), ja(*)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp), a(*),
     .     delta, gamma
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DSEXPV computes w = exp(t*A)*v - for a Symmetric matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
*
*     nz     : (input) the actual number of non zero entries.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, nbr of Tridiagonal matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, nbr of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )

*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hjj, hump, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         if ( j.gt.1 )
     .     call DAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = DDOT( n, wsp(j1v-n),1, wsp(j1v),1 )
         call DAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         wsp(ih+j*mh+j-1) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m-1) = 0.0d0
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
 600  END
*----------------------------------------------------------------------|

*----------------------------------------------------------------------|
      subroutine ZGEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .                   itrace,iflag )

      implicit none
      integer          n, nz, m, lwsp, liwsp, itrace, iflag, mxstep,
     .     iwsp(liwsp), ia(*), ja(*)
      double precision t, tol, anorm, delta, gamma
      complex*16       v(n), w(n), wsp(lwsp), a(*)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZGEXPV computes w = exp(t*A)*v
*     for a Zomplex (i.e., complex double precision) matrix A 
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
*
*     nz     : (input) the actual number of non zero entries.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        complex*16 x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )

*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|

      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hij

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DZNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = ONE
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         enddo
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hij = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      wsp(9)  = CMPLX( hump/vnorm )
      wsp(10) = CMPLX( beta/vnorm )
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine ZHEXPV( a, ia, ja, n, nz, m, t, v, w, tol,mxstep,anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .                   itrace,iflag )

      implicit none
      integer          n, nz, m, lwsp, liwsp, itrace, iflag, mxstep,
     .     iwsp(liwsp), ia(*), ja(*)
      double precision t, tol, anorm, delta, gamma
      complex*16       v(n), w(n), wsp(lwsp), a(*)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZHEXPV computes w = exp(t*A)*v for a Hermitian matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
*
*     nz     : (input) the actual number of non zero entries.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        complex*16 x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )

*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hjj

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      beta = DZNRM2( n, w,1 )
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         if ( j.gt.1 )
     .     call ZAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = ZDOTC( n, wsp(j1v-n),1, wsp(j1v),1 )
         call ZAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DZNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         wsp(ih+j*mh+j-1) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m-1) = ZERO
      wsp(ih+m*mh+m+1) = ONE
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         enddo
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif


* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hjj = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hjj,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      wsp(9)  = CMPLX( hump/vnorm )
      wsp(10) = CMPLX( beta/vnorm )
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,
     .     anorm, wsp,lwsp, iwsp,liwsp, matvec, delta, gamma, 
     .     itrace,iflag ) 

      implicit none
      integer n, nz, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp), mxstep,
     .     ia(*), ja(*)
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp), a(*),
     .     delta, gamma 
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a General matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
*
*     nz     : (input) the actual number of non zero entries.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at which the solution is needed (can be < 0).
*   
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+2)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )
  
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hij, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm
 
      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call DCOPY( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv), a, ia, ja, n, nz )
      call DAXPY( n, 1.0d0, u,1, wsp(iv),1 )
      beta = DNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call DSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = 1.0d0
      wsp(ih+(m-1)*mh+m) = 0.0d0
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = 1.0d0
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
*---  irreducible rational Pade approximation ...
      mx = mbrkdwn + MAX( 1,k1 )
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .             wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)
 

*---  error estimate ...
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif

*---  reject the step-size if the error is not acceptable ...   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and. 
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif 
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1 

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step 
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif
 
      step_min = MIN( step_min, t_step ) 
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )
 
      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1
 
 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep, 
     .     anorm, wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .     itrace,iflag ) 

      implicit none
      integer n, nz, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp), mxstep,
     .     ia(*), ja(*)
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp), a(*),
     .     delta, gamma
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DSPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a Symmetric matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
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
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+2)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 ) 
  
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hjj, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm
 
      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call DCOPY( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv), a, ia, ja, n, nz )
      call DAXPY( n, 1.0d0, u,1, wsp(iv),1 )
      beta = DNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call DSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         if ( j.gt.1 )
     .     call DAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = DDOT( n, wsp(j1v-n),1, wsp(j1v),1 )
         call DAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         wsp(ih+j*mh+j-1) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = 1.0d0
      wsp(ih+m*mh+m-1)   = 0.0d0
      wsp(ih+(m-1)*mh+m) = 0.0d0
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = 1.0d0
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .              wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)
 

* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and. 
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif 
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1 

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step 
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step ) 
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )
 
      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1
 
 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine ZGPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep,
     .     anorm, wsp,lwsp, iwsp,liwsp, matvec, delta,gamma, 
     .     itrace,iflag )

      implicit none
      integer          n, nz, m, lwsp, liwsp, itrace, iflag, mxstep,
     .     iwsp(liwsp), ia(*), ja(*)
      double precision t, tol, anorm, delta, gamma
      complex*16       u(n), v(n), w(n), wsp(lwsp), a(*)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a Zomplex (i.e., complex double 
*     precision matrix).
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
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
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+2)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 ) 

*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, SQR1
      complex*16 hij

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call ZCOPY( n, v,1, w,1 )
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv), a, ia, ja, n, nz )
      call ZAXPY( n, ONE, u,1, wsp(iv),1 )
      beta = DZNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call ZDSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DZNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            if ( itrace.ge.1 ) then
                call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
               call dblepr('h', -1, hj1j, 1)
            endif           
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = ONE
      wsp(ih+(m-1)*mh+m) = ZERO
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = ONE
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .             wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)


* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      hij = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iphih),1,ONE,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
 600  END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine ZHPHIV( a, ia, ja, n, nz, m, t, u, v, w, tol,mxstep, 
     .     anorm, wsp,lwsp, iwsp,liwsp, matvec, delta,gamma,
     .     itrace,iflag )

      implicit none
      integer          n, nz, m, lwsp, liwsp, itrace, iflag, mxstep,
     .     iwsp(liwsp), ia(*), ja(*)
      double precision t, tol, anorm, delta, gamma
      complex*16       u(n), v(n), w(n), wsp(lwsp), a(*)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZHPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is an Hermitian matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     ia,ja,a (input):
*           sparse matrix data stored in the appropriate storage format 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly.
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
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+2)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y, a, ia, ja, n, nz )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*     itrace : (input) running mode. 0=silent, 1>=print happy breakdown,
*              2>=print step-by-step info.
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxreject, ideg
      parameter( mxreject = 0,
     .           ideg     = 6 )

*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, SQR1
      complex*16 hjj

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) goto 600
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call ZCOPY( n, v,1, w,1 )
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv), a, ia, ja, n, nz )
      call ZAXPY( n, ONE, u,1, wsp(iv),1 )
      beta = DZNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call ZDSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
         if ( j.gt.1 )
     .     call ZAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = ZDOTC( n, wsp(j1v-n),1, wsp(j1v),1 )
         call ZAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DZNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
         if ( hj1j.le.break_tol ) then
            call intpr('happy breakdown:\n mbrkdwn', -1, j, 1)
            call dblepr('h', -1, hj1j, 1)
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         wsp(ih+j*mh+j-1) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, n, nz )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = ONE
      wsp(ih+m*mh+m-1)   = ZERO
      wsp(ih+(m-1)*mh+m) = ZERO
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = ONE
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .              wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)


* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ge.2 ) then
            call dblepr('t_step', -1, t_old, 1)
            call dblepr('err_loc',-1, err_loc, 1)
            call dblepr('err_required', -1, delta*t_old*tol, 1)
            call dblepr('stepsize rejected, stepping down to:', -1, 
     .           t_step, 1)
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            call intpr('The requested tolerance is too high.', 
     .           -1, 1, 0)
            call intpr('Rerun with a smaller value.', 
     .           -1, 1, 0)
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      hjj = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hjj,wsp(iv),n,wsp(iphih),1,ONE,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ge.2 ) then
         call dblepr('integration', -1, nstep, 1)
         call dblepr('scale-square', -1, ns, 1)
         call dblepr('step_size', -1, t_step, 1)
         call dblepr('err_loc', -1, err_loc, 1)
         call dblepr('next_step', -1, t_new, 1)
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
 600  END
*----------------------------------------------------------------------|
