#ifndef kernel_streching_h
#define kernel_streching_h
#include "exafmm.h"

namespace exafmm
{
  const complex_t I(0., 1.);     //!< Imaginary unit
  std::vector<real_t> prefactor; // sqrt( (n - |m|)! / (n + |m|)! )
  std::vector<real_t> Anm;       // (-1)^n / sqrt( (n + m)! / (n - m)! )
  std::vector<complex_t> Cnm;    // M2L translation matrix Cjknm
  double EPS = 1e-16;
  //!< L2 norm of vector X
  inline real_t norm(real_t *X)
  {
    return X[0] * X[0] + X[1] * X[1] + X[2] * X[2]; // L2 norm
  }
  //! Odd or even
  inline int oddOrEven(int n)
  {
    return (((n) & 1) == 1) ? -1 : 1;
  }

  //! Get r,theta,phi from x,y,z
  void cart2sph(real_t *dX, real_t &r, real_t &theta, real_t &phi)
  {
    r = sqrt(norm(dX));                   // r = sqrt(x^2 + y^2 + z^2)
    theta = r == 0 ? 0 : acos(dX[2] / r); // theta = acos(z / r)
    phi = atan2(dX[1], dX[0]);            // phi = atan(y / x)
  }

  //! Spherical to cartesian coordinates
  template <class T>
  void sph2cart(real_t r, real_t theta, real_t phi, T *spherical, T *cartesian)
  {
    cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0] // x component (not x itself)
                   + std::cos(theta) * std::cos(phi) / r * spherical[1] - std::sin(phi) / r / std::sin(theta) * spherical[2];
    cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0] // y component (not y itself)
                   + std::cos(theta) * std::sin(phi) / r * spherical[1] + std::cos(phi) / r / std::sin(theta) * spherical[2];
    cartesian[2] = std::cos(theta) * spherical[0] // z component (not z itself)
                   - std::sin(theta) / r * spherical[1];
  }

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t *Ynm, complex_t *YnmTheta)
  {
    real_t x = std::cos(alpha); // x = cos(alpha)
    real_t y = std::sin(alpha); // y = sin(alpha)
    real_t fact = 1;            // Initialize 2 * m + 1
    real_t pn = 1;              // Initialize Legendre polynomial Pn
    real_t rhom = 1;            // Initialize rho^m
    for (int m = 0; m < P; m++)
    {                                                                           // Loop over m in Ynm
      complex_t eim = std::exp(I * real_t(m * beta));                           //  exp(i * m * beta)
      real_t p = pn;                                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;                               //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim; // theta derivative of r^n * Ynm
      rhom *= rho;                                                              //  rho^m
      real_t rhon = rhom;                                                       //  rho^n
      for (int n = m + 1; n < P; n++)
      {                                                                                         //  Loop over n in Ynm
        int npm = n * n + n + m;                                                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;                                             //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                                                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                                                         //   Pnm-2
        p1 = p;                                                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);                                //   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim; // theta derivative
        rhon *= rho;                                                                            //   Update rho^n
      }                                                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                                                      //  Pn
      fact += 2;                                                                                //  2 * m + 1
    }                                                                                           // End loop over m in Ynm
  }

  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t *Ynm2)
  {
    real_t x = std::cos(alpha); // x = cos(alpha)
    real_t y = std::sin(alpha); // y = sin(alpha)
    real_t fact = 1;            // Initialize 2 * m + 1
    real_t pn = 1;              // Initialize Legendre polynomial Pn
    real_t rhom = 1.0 / rho;    // Initialize rho^(-m-1)
    for (int m = 0; m < 2 * P; m++)
    {                                                 // Loop over m in Ynm
      complex_t eim = std::exp(I * real_t(m * beta)); //  exp(i * m * beta)
      real_t p = pn;                                  //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                        //  Index of Ynm for m > 0
      int nmn = m * m;                                //  Index of Ynm for m < 0
      Ynm2[npn] = rhom * p * prefactor[npn] * eim;    //  rho^(-m-1) * Ynm for m > 0
      Ynm2[nmn] = std::conj(Ynm2[npn]);               //  Use conjugate relation for m < 0
      real_t p1 = p;                                  //  Pnm-1
      p = x * (2 * m + 1) * p1;                       //  Pnm using recurrence relation
      rhom /= rho;                                    //  rho^(-m-1)
      real_t rhon = rhom;                             //  rho^(-n-1)
      for (int n = m + 1; n < 2 * P; n++)
      {                                                          //  Loop over n in Ynm
        int npm = n * n + n + m;                                 //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                 //   Index of Ynm for m < 0
        Ynm2[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
        Ynm2[nmm] = std::conj(Ynm2[npm]);                        //   Use conjugate relation for m < 0
        real_t p2 = p1;                                          //   Pnm-2
        p1 = p;                                                  //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
        rhon /= rho;                                             //   rho^(-n-1)
      }                                                          //  End loop over n in Ynm
      pn = -pn * fact * y;                                       //  Pn
      fact += 2;                                                 //  2 * m + 1
    }                                                            // End loop over m in Ynm
  }

  void initKernel()
  {
    NTERM = 3 * (P * (P + 1) / 2);
    prefactor.resize(4 * P * P);
    Anm.resize(4 * P * P);
    Cnm.resize(P * P * P * P);
    for (int n = 0; n < 2 * P; n++)
    { // Loop over n in Anm
      for (int m = -n; m <= n; m++)
      {                         //  Loop over m in Anm
        int nm = n * n + n + m; //   Index of Anm
        int nabsm = abs(m);     //   |m|
        real_t fnmm = EPS;      //   Initialize (n - m)!
        for (int i = 1; i <= n - m; i++)
          fnmm *= i;       //   (n - m)!
        real_t fnpm = EPS; //   Initialize (n + m)!
        for (int i = 1; i <= n + m; i++)
          fnpm *= i;       //   (n + m)!
        real_t fnma = 1.0; //   Initialize (n - |m|)!
        for (int i = 1; i <= n - nabsm; i++)
          fnma *= i;       //   (n - |m|)!
        real_t fnpa = 1.0; //   Initialize (n + |m|)!
        for (int i = 1; i <= n + nabsm; i++)
          fnpa *= i;                                     //   (n + |m|)!
        prefactor[nm] = std::sqrt(fnma / fnpa);          //   sqrt( (n - |m|)! / (n + |m|)! )
        Anm[nm] = oddOrEven(n) / std::sqrt(fnmm * fnpm); //   (-1)^n / sqrt( (n + m)! / (n - m)! )
      }                                                  //  End loop over m in Anm
    }                                                    // End loop over n in Anm
    for (int j = 0, jk = 0, jknm = 0; j < P; j++)
    { // Loop over j in Cjknm
      for (int k = -j; k <= j; k++, jk++)
      { //  Loop over k in Cjknm
        for (int n = 0, nm = 0; n < P; n++)
        { //   Loop over n in Cjknm
          for (int m = -n; m <= n; m++, nm++, jknm++)
          {                                                               //    Loop over m in Cjknm
            const int jnkm = (j + n) * (j + n) + j + n + m - k;           //     Index C_{j+n}^{m-k}
            Cnm[jknm] = std::pow(I, real_t(abs(k - m) - abs(k) - abs(m))) //     Cjknm
                        * real_t(oddOrEven(j) * Anm[nm] * Anm[jk] / Anm[jnkm]) * EPS;
          } //    End loop over m in Cjknm
        }   //   End loop over n in Cjknm
      }     //  End loop over in k in Cjknm
    }       // End loop over in j in Cjknm
  }

  void P2P(Cell *Ci, Cell *Cj)
  {
    Body *Bi = Ci->BODY;
    Body *Bj = Cj->BODY;
    int ni = Ci->NBODY;
    int nj = Cj->NBODY;
    for (int i = 0; i < ni; i++)
    {
      real_t pot = 0;
      real_t ax = 0;
      real_t ay = 0;
      real_t az = 0;
      for (int j = 0; j < nj; j++)
      {
        for (int d = 0; d < 3; d++)
          dX[d] = Bi[i].X[d] - Bj[j].X[d];
        real_t R2 = norm(dX);
        if (R2 != 0)
        {
          real_t S2 = 2;                                                         //   2 * simga^2
          real_t RS = R2 / S2;                                                   //   R^2 / (2 * sigma^2)
          real_t cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf(std::sqrt(RS)) // cutoff function for first term
                                                              - std::sqrt(4 / M_PI * RS) * exp(-RS));

          ax += (Bi[i].alpha[1] * Bj[j].alpha[2] - Bi[i].alpha[2] * Bj[j].alpha[1]) * cutoff; // x component of first term

          ay += (Bi[i].alpha[2] * Bj[j].alpha[0] - Bi[i].alpha[0] * Bj[j].alpha[2]) * cutoff; // y component of first term

          az += (Bi[i].alpha[0] * Bj[j].alpha[1] - Bi[i].alpha[1] * Bj[j].alpha[0]) * cutoff; // z component of first term
          cutoff = 0.25 / M_PI / R2 / R2 / std::sqrt(R2) * (3 * erf(std::sqrt(RS))            // cutoff function for second term
                                                            - (2 * RS + 3) * std::sqrt(4 / M_PI * RS) * exp(-RS)) *
                   (Bi[i].alpha[0] * dX[0] + Bi[i].alpha[1] * dX[1] + Bi[i].alpha[2] * dX[2]);

          ax += (Bj[j].alpha[1] * dX[2] - Bj[j].alpha[2] * dX[1]) * cutoff; // x component of second term

          ay += (Bj[j].alpha[2] * dX[0] - Bj[j].alpha[0] * dX[2]) * cutoff; // y component of second term

          az += (Bj[j].alpha[0] * dX[1] - Bj[j].alpha[1] * dX[0]) * cutoff; // z component of second term
        }
      }
      Bi[i].p += pot;
      Bi[i].velocity[0] -= ax;
      Bi[i].velocity[1] -= ay;
      Bi[i].velocity[2] -= az;
    }
  }

  void P2M(Cell *C)
  {
    complex_t Ynm[P * P], YnmTheta[P * P];
    for (Body *B = C->BODY; B != C->BODY + C->NBODY; B++)
    {
      {
        for (int d = 0; d < 3; d++)
          dX[d] = B->X[d] - C->X[d];
        complex_t spherical[3];
        complex_t cartesian[3];
        real_t rho, alpha, beta;
        cart2sph(dX, rho, alpha, beta);
        evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
        for (int n = 0; n < P; n++)
        {
          for (int m = 0; m <= n; m++)
          {
            int nm = n * n + n + m;
            int nms = n * (n + 1) / 2 + m;
            spherical[0] = Ynm[nm] * double(n / rho);
            spherical[1] = YnmTheta[nm];
            spherical[2] = -Ynm[nm] * I * double(m);
            sph2cart(rho, alpha, beta, spherical, cartesian);
            C->M[3 * nms + 0] += cartesian[2] * double(B->alpha[1]) - cartesian[1] * double(B->alpha[2]);
            C->M[3 * nms + 1] += cartesian[0] * double(B->alpha[2]) - cartesian[2] * double(B->alpha[0]);
            C->M[3 * nms + 2] += cartesian[1] * double(B->alpha[0]) - cartesian[0] * double(B->alpha[1]);
          }
        }
      }
    }
  }

  void M2M(Cell *Ci)
  {
    complex_t Ynm[P * P], YnmTheta[P * P];
    for (Cell *Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
    {
      {
        for (int d = 0; d < 3; d++)
          dX[d] = Ci->X[d] - Cj->X[d];
        real_t rho, alpha, beta;
        cart2sph(dX, rho, alpha, beta);
        evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
        for (int j = 0; j < P; j++)
        {
          for (int k = 0; k <= j; k++)
          {
            int jk = j * j + j + k;
            int jks = j * (j + 1) / 2 + k;
            complex_t M[3]{0, 0, 0};
            for (int n = 0; n <= j; n++)
            {
              for (int m = -n; m <= std::min(k - 1, n); m++)
              {
                if (j - n >= k - m)
                {
                  int jnkm = (j - n) * (j - n) + j - n + k - m;
                  int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                  int nm = n * n + n + m;
                  for (int d = 0; d < 3; d++)
                  {
                    M[d] += Cj->M[3 * jnkms + d] * std::pow(I, real_t(m - abs(m))) * Ynm[nm] * real_t(oddOrEven(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
                  }
                }
              }
              for (int m = k; m <= n; m++)
              {
                if (j - n >= m - k)
                {
                  int jnkm = (j - n) * (j - n) + j - n + k - m;
                  int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                  int nm = n * n + n + m;
                  for (int d = 0; d < 3; d++)
                  {
                    M[d] += std::conj(Cj->M[3 * jnkms + d]) * Ynm[nm] * real_t(oddOrEven(k + n + m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
                  }
                }
              }
            }
            for (int d = 0; d < 3; d++)
            {
              Ci->M[3 * jks + d] += M[d] * EPS;
            }
          }
        }
      }
    }
  }

  void M2L(Cell *Ci, Cell *Cj)
  {
    complex_t Ynm2[4 * P * P];
    for (int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];
    real_t rho, alpha, beta;
    cart2sph(dX, rho, alpha, beta);
    evalLocal(rho, alpha, beta, Ynm2);
    for (int j = 0; j < P; j++)
    {
      for (int k = 0; k <= j; k++)
      {
        int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        complex_t L[3]{0, 0, 0};
        for (int n = 0; n < P; n++)
        {
          for (int m = -n; m < 0; m++)
          {
            int nm = n * n + n + m;
            int nms = n * (n + 1) / 2 - m;
            int jknm = jk * P * P + nm;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            for (int d = 0; d < 3; d++)
            {
              L[d] += std::conj(Cj->M[3 * nms + d]) * Cnm[jknm] * Ynm2[jnkm];
            }
          }
          for (int m = 0; m <= n; m++)
          {
            int nm = n * n + n + m;
            int nms = n * (n + 1) / 2 + m;
            int jknm = jk * P * P + nm;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            for (int d = 0; d < 3; d++)
            {
              L[d] += Cj->M[3 * nms + d] * Cnm[jknm] * Ynm2[jnkm];
            }
          }
        }
        for (int d = 0; d < 3; d++)
        {
          Ci->L[3 * jks + d] += L[d];
        }
      }
    }
  }

  void L2L(Cell *Cj)
  {
    complex_t Ynm[P * P], YnmTheta[P * P];
    for (Cell *Ci = Cj->CHILD; Ci != Cj->CHILD + Cj->NCHILD; Ci++)
    {
      for (int d = 0; d < 3; d++)
        dX[d] = Ci->X[d] - Cj->X[d];
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
      for (int j = 0; j < P; j++)
      {
        for (int k = 0; k <= j; k++)
        {
          int jk = j * j + j + k;
          int jks = j * (j + 1) / 2 + k;
          complex_t L[3]{0, 0, 0};
          for (int n = j; n < P; n++)
          {
            for (int m = j + k - n; m < 0; m++)
            {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nm = n * n + n - m;
              int nms = n * (n + 1) / 2 - m;
              for (int d = 0; d < 3; d++)
              {
                L[d] += std::conj(Cj->L[3 * nms + d]) * Ynm[jnkm] * real_t(oddOrEven(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
              }
            }
            for (int m = 0; m <= n; m++)
            {
              if (n - j >= abs(m - k))
              {
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                int nm = n * n + n + m;
                int nms = n * (n + 1) / 2 + m;
                for (int d = 0; d < 3; d++)
                {
                  L[d] += Cj->L[3 * nms + d] * std::pow(I, real_t(m - k - abs(m - k))) * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
                }
              }
            }
          }
          for (int d = 0; d < 3; d++)
          {
            Ci->L[3 * jks + d] += L[d] * EPS;
          }
        }
      }
    }
  }

  void L2P(Cell *Ci)
  {
    {
      complex_t Ynm[P * P], YnmTheta[P * P];
      for (Body *B = Ci->BODY; B != Ci->BODY + Ci->NBODY; B++)
      {
        for (int d = 0; d < 3; d++)
          dX[d] = B->X[d] - Ci->X[d];
        real_t spherical[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        real_t cartesian[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        real_t r, theta, phi;
        cart2sph(dX, r, theta, phi);
        evalMultipole(r, theta, phi, Ynm, YnmTheta);
        for (int n = 0; n < P; n++)
        {
          int nm = n * n + n;
          int nms = n * (n + 1) / 2;
          for (int d = 0; d != 3; ++d)
          {
            spherical[d][0] += (Ci->L[3 * nms + d] * Ynm[nm]).real() / r * n;
            spherical[d][1] += (Ci->L[3 * nms + d] * YnmTheta[nm]).real();
          }
          for (int m = 1; m <= n; m++)
          {
            nm = n * n + n + m;
            nms = n * (n + 1) / 2 + m;
            for (int d = 0; d != 3; ++d)
            {
              spherical[d][0] += 2 * std::real(Ci->L[3 * nms + d] * Ynm[nm]) / r * n;
              spherical[d][1] += 2 * std::real(Ci->L[3 * nms + d] * YnmTheta[nm]);
              spherical[d][2] += 2 * std::real(Ci->L[3 * nms + d] * Ynm[nm] * I) * m;
            }
          }
        }
        for (int d = 0; d != 3; ++d)
        {
          sph2cart(r, theta, phi, spherical[d], cartesian[d]);
        }
        B->velocity[0] += 0.25 / M_PI * (B->alpha[0] * cartesian[0][0] + B->alpha[1] * cartesian[0][1] + B->alpha[2] * cartesian[0][2]);
        B->velocity[1] += 0.25 / M_PI * (B->alpha[0] * cartesian[1][0] + B->alpha[1] * cartesian[1][1] + B->alpha[2] * cartesian[1][2]);
        B->velocity[2] += 0.25 / M_PI * (B->alpha[0] * cartesian[2][0] + B->alpha[1] * cartesian[2][1] + B->alpha[2] * cartesian[2][2]);
      }
    }
  }
}
#endif