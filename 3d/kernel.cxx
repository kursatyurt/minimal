#include "kernel.h"
using namespace exafmm;

int main(int argc, char **argv)
{
  P = atoi(argv[1]);
  initKernel();

  auto testVelocity = []()
  {
    // P2M
    Bodies jbodies(1);
    for (int d = 0; d < 3; d++)
      jbodies[0].X[d] = 2;
    for (int d = 0; d < 3; d++)
      jbodies[0].alpha[d] = 1;
    jbodies[0].radius = 0.1;
    Cells cells(4);
    Cell *Cj = &cells[0];
    Cj->X[0] = 3;
    Cj->X[1] = 1;
    Cj->X[2] = 1;
    Cj->R = 1;
    Cj->BODY = &jbodies[0];
    Cj->NBODY = jbodies.size();
    Cj->M.resize(NTERM, 0.0);
    P2M(Cj);

    // M2M
    Cell *CJ = &cells[1];
    CJ->CHILD = Cj;
    CJ->NCHILD = 1;
    CJ->X[0] = 4;
    CJ->X[1] = 0;
    CJ->X[2] = 0;
    CJ->R = 2;
    CJ->M.resize(NTERM, 0.0);
    M2M(CJ);

    // M2L
    Cell *CI = &cells[2];
    CI->X[0] = -4;
    CI->X[1] = 0;
    CI->X[2] = 0;
    CI->R = 2;
    CI->L.resize(NTERM, 0.0);
    M2L(CI, CJ);

    // L2L
    Cell *Ci = &cells[3];
    CI->CHILD = Ci;
    CI->NCHILD = 1;
    Ci->X[0] = -3;
    Ci->X[1] = 1;
    Ci->X[2] = 1;
    Ci->R = 1;
    Ci->L.resize(NTERM, 0.0);
    L2L(CI);

    // L2P
    Bodies bodies(1);
    bodies[0].X[0] = -2;
    bodies[0].X[1] = 2;
    bodies[0].X[2] = 2;
    bodies[0].alpha[0] = 1;
    bodies[0].alpha[1] = 1;
    bodies[0].alpha[2] = 1;
    bodies[0].radius = 0.1;
    for (int d = 0; d < 3; d++)
      bodies[0].velocity[d] = 0;
    Ci->BODY = &bodies[0];
    Ci->NBODY = bodies.size();
    L2P(Ci);

    // P2P
    Bodies bodies2(1);
    for (size_t b = 0; b < bodies2.size(); b++)
    {
      bodies2[b] = bodies[b];
      bodies2[b].radius = 0.1;
      for (int d = 0; d < 3; d++)
        bodies2[b].velocity[d] = 0;
    }
    Cj->NBODY = jbodies.size();
    Ci->NBODY = bodies2.size();
    Ci->BODY = &bodies2[0];
    P2P(Ci, Cj);

    // Verify results
    real_t FDif = 0, FNrm = 0;
    for (size_t b = 0; b < bodies.size(); b++)
    {
      FDif += (bodies[b].velocity[0] - bodies2[b].velocity[0]) * (bodies[b].velocity[0] - bodies2[b].velocity[0]) +
              (bodies[b].velocity[1] - bodies2[b].velocity[1]) * (bodies[b].velocity[1] - bodies2[b].velocity[1]) +
              (bodies[b].velocity[2] - bodies2[b].velocity[2]) * (bodies[b].velocity[2] - bodies2[b].velocity[2]);
      FNrm += bodies[b].velocity[0] * bodies[b].velocity[0] + bodies[b].velocity[1] * bodies[b].velocity[1] +
              bodies[b].velocity[2] * bodies[b].velocity[2];
    }
    printf("%-20s : %8.5e s\n", "Rel. L2 Error (Velocity)", sqrt(FDif / FNrm));
  };

  auto testStreching = []()
  {
    // P2M
    Bodies jbodies(1);
    for (int d = 0; d < 3; d++)
      jbodies[0].X[d] = 2;
    // for (int d = 0; d < 3; d++)
    jbodies[0].alpha[0] = 2.2;
    jbodies[0].alpha[1] = 1.2;
    jbodies[0].alpha[2] = 3.2;
    jbodies[0].radius = 0.125;
    Cells cells(4);
    Cell *Cj = &cells[0];
    Cj->X[0] = 3;
    Cj->X[1] = 1;
    Cj->X[2] = 1;
    Cj->R = 1;
    Cj->BODY = &jbodies[0];
    Cj->NBODY = jbodies.size();
    Cj->M.resize(NTERM, 0.0);
    P2MS(Cj);

    // M2M
    Cell *CJ = &cells[1];
    CJ->CHILD = Cj;
    CJ->NCHILD = 1;
    CJ->X[0] = 4;
    CJ->X[1] = 0;
    CJ->X[2] = 0;
    CJ->R = 2;
    CJ->M.resize(NTERM, 0.0);
    M2M(CJ);

    // M2L
    Cell *CI = &cells[2];
    CI->X[0] = -4;
    CI->X[1] = 0;
    CI->X[2] = 0;
    CI->R = 2;
    CI->L.resize(NTERM, 0.0);
    M2L(CI, CJ);

    // L2L
    Cell *Ci = &cells[3];
    CI->CHILD = Ci;
    CI->NCHILD = 1;
    Ci->X[0] = -3;
    Ci->X[1] = 1;
    Ci->X[2] = 1;
    Ci->R = 1;
    Ci->L.resize(NTERM, 0.0);
    L2L(CI);

    // L2P
    Bodies bodies(1);
    bodies[0].X[0] = -2.1;
    bodies[0].X[1] = 1.9;
    bodies[0].X[2] = 2.2;
    bodies[0].alpha[0] = 1.259;
    bodies[0].alpha[1] = 0.1;
    bodies[0].alpha[2] = 1.9;
    bodies[0].radius = 0.225;
    for (int d = 0; d < 3; d++)
      bodies[0].dadt[d] = 0;
    Ci->BODY = &bodies[0];
    Ci->NBODY = bodies.size();
    L2PS(Ci);

    // P2P
    Bodies bodies2(1);
    for (size_t b = 0; b < bodies2.size(); b++)
    {
      bodies2[b] = bodies[b];
      for (int d = 0; d < 3; d++)
        bodies2[b].dadt[d] = 0;
    }
    Cj->NBODY = jbodies.size();
    Ci->NBODY = bodies2.size();
    Ci->BODY = &bodies2[0];
    P2PS(Ci, Cj);

    // Verify results
    real_t FDif = 0, FNrm = 0;
    for (size_t b = 0; b < bodies.size(); b++)
    {
      // Write each particle data to screen
      FDif += (bodies[b].dadt[0] - bodies2[b].dadt[0]) * (bodies[b].dadt[0] - bodies2[b].dadt[0]) +
              (bodies[b].dadt[1] - bodies2[b].dadt[1]) * (bodies[b].dadt[1] - bodies2[b].dadt[1]) +
              (bodies[b].dadt[2] - bodies2[b].dadt[2]) * (bodies[b].dadt[2] - bodies2[b].dadt[2]);
      FNrm += bodies[b].dadt[0] * bodies[b].dadt[0] + bodies[b].dadt[1] * bodies[b].dadt[1] +
              bodies[b].dadt[2] * bodies[b].dadt[2];
    }
    printf("%-20s : %8.5e \n", "Rel. L2 Error (Streching)", sqrt(FDif / FNrm));
  };

  testVelocity();
  testStreching();

  return 0;
}
