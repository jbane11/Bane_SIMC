ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c F1F209.f
c Package of FORTRAN subroutine describing fits to inclusive inelastic
c electron scattering from proton, neutron, deuteron, or heaviere nuclei
c Proton fit is described in:
c   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision 
c    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
c    (arXiv:0712.3731). Submitted to Phys. Rev. C.
c Deuteron (and netron) fit is described in:
c    P.E. Bosted and M.E. Christy, ``Empirical Fit to Inelastic 
c    Electron-Deuteron and Electron-Neutron
c    Resonance Region Transverse Cross Sections, 
c    (arXiv:0711.0159), publichished in Phys. Rev. C 77, 065206 (2008). (
c New fits for A>2 by Vahe M. and Peter B. (to be publsihed).
c
c routine in this package:
c Main program test. Compiling and linking all code
c    in this package will make a "test" program.
c    Upon execution, the output should be:
c  1.000  1.500  0.617  0.204  0.231  0.312  0.410  1.644  0.072  0.102
c  1.000  2.000  0.472  0.129  0.159  0.292  0.250  2.329  0.199  0.191
c  1.000  2.500  0.382  0.169  0.227  0.417  0.360  2.973  0.241  0.292
c  1.000  3.000  0.321  0.224  0.315  0.511  0.500  3.562  0.327  0.359
c  1.000  3.500  0.276  0.218  0.309  0.517  0.521  4.060  0.348  0.387
c  2.000  1.500  0.763  0.079  0.091  0.110  0.104  0.400  0.154  0.188
c  2.000  2.000  0.641  0.058  0.091  0.157  0.104  0.764  0.184  0.216
c  2.000  2.500  0.553  0.090  0.139  0.241  0.162  1.174  0.215  0.246
c  2.000  3.000  0.485  0.137  0.214  0.324  0.252  1.558  0.233  0.270
c  2.000  3.500  0.433  0.136  0.227  0.354  0.272  1.919  0.267  0.283
c   1.000   0.747   0.432
c   2.000   0.157   0.132
c   3.000   0.048   0.049
c
c
c subroutine F1F2IN09. Returns inelastic F1, F2, and R for
c      nucleus with charge Z and atomic number A
c      for given value of Q2 and W**2. F1 and F2
c      are per nucleus (not per nucleon).
c subroutine pind. Returns F1, R, Sigma_L, and Sigma_T
c      for a proton with Fermi motion of a detueron
c subroutine resd. Returns F1 for a deuteron.
c subroutine resder. Returns error on F1 from fit.
c      Requires auxillary file "F1F207D2emat.dat"
c subroutine resmodd. Returns F1 for average of a free
c      proton and neutron
c subroutine christy507. Returns F1, R, sigma_L, and
c      sigma_T for the free proton.
c subroutine resmod507. Returns sigma_L or sigma_T
c      (called by christy507)
c subroutine mec. Called by resd to get extra terms
c      for dip region between quasi-elastic and Delta
c funtion fitemc. Fit to "EMC effect" used to 
c      get F1 and F2 for A>2
c 
c subroutine F1F2QE09. Returns quasi-elastic F1, F2 for
c      nucleus with charge Z and atomic number A
c      for given value of Q2 and W**2
c 
      implicit none
      integer iq,iw
      real*8 q2,w2,F1n,F2n,r,F1p,F2p,F1d,F2d,F1c,F2c
      real*8 f1dqe, f2dqe, nu, x, am/0.9383/, F1be, F2be,rbe

      do iq=1,2
       q2 = 1.0 * iq
       do iw=1,5
        w2 = 1. + 0.5*iw
        nu = (w2 - am**2 + q2) / 2. / am
        x = q2 / 2. / am / nu
        call F1F2IN09(0.D0, 1.D0, q2, w2, F1n, F2n,r)
        call F1F2IN09(1.D0, 1.D0, q2, w2, F1p, F2p,r)
        call F1F2IN09(1.D0, 2.D0, q2, w2, F1d, F2d,r)
        call F1F2IN09(4.D0, 9.D0, q2, w2, F1be, F2be,rbe)
        write(6,'(10f7.3)') q2,w2,x,f2n,f2p,f2d,
     >    f1p,f1be,r,rbe
       enddo
      enddo
      do iq=1,3
        q2 = 1.0 * iq
        w2 = am**2
        call F1F2QE09(1.D0, 2.D0, q2, w2, F1dqe, F2dqe)
        write(6,'(3f8.3)') q2, F1dqe, F2dqe
      enddo

      stop
      end

