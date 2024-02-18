#include "build_tree.h"
#include "kernel.h"
#include "timer.h"
#if EXAFMM_EAGER
#include "traverse_eager.h"
#elif EXAFMM_LAZY
#include "traverse_lazy.h"
#endif
using namespace exafmm;
#include <fstream>

int main(int argc, char **argv)
{
  const int numBodies = 1'800; // Number of bodies
  const int numTargets = 1'000; // Number of targets for checking answer
  P = 10;                       // Order of expansions
  ncrit = 64;                   // Number of bodies per leaf cell
  theta = 0.4;                  // Multipole acceptance criterion

  printf("--- %-16s ------------\n", "FMM Profiling"); // Start profiling
  //! Initialize bodies

  Bodies bodies(numBodies); // Initialize bodies

  auto initRandom = [&bodies]()
  {
    srand48(0); // Set seed for random number generator
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      for (int d = 0; d < 3; d++)
      {                                               //  Loop over dimension
        bodies[b].X[d] = drand48() * 2 * M_PI - M_PI; //   Initialize positions
        bodies[b].alpha[d] = drand48();               //  Initialize charge
      }                                               //  End loop over dimension
      bodies[b].radius = drand48();                   //  Initialize radius
      for (int d = 0; d < 3; d++)
      {
        bodies[b].velocity[d] = 0;
        bodies[b].dadt[d] = 0;
      }
    } // End loop over bodies
  };  // End initRandom

  auto initFromFile = [&bodies]()
  {
    std::ifstream file;
    file.open("init.txt");
    if (file.is_open())
    {
      std::string line;
      int i = 0;
      while (std::getline(file, line))
      {
        std::istringstream iss(line);
        double x, y, z, sigma, alpha_x, alpha_y, alpha_z;
        if (!(iss >> alpha_x >> alpha_y >> alpha_z >> sigma >> x >> y >> z))
        {
          break;
        } // error
        bodies[i].X[0] = x;
        bodies[i].X[1] = y;
        bodies[i].X[2] = z;
        bodies[i].radius = sigma;
        bodies[i].alpha[0] = alpha_x;
        bodies[i].alpha[1] = alpha_y;
        bodies[i].alpha[2] = alpha_z;
        bodies[i].velocity[0] = 0;
        bodies[i].velocity[1] = 0;
        bodies[i].velocity[2] = 0;
        bodies[i].dadt[0] = 0;
        bodies[i].dadt[1] = 0;
        bodies[i].dadt[2] = 0;
        i++;
      }
      file.close();
    }
    else
    {
      return;
    }
  };

  // initRandom();
  initFromFile();

  auto writeTovtk = [&bodies](int step)
  {
    std::ofstream file;
    file.open("output" + std::to_string(step) + ".vtk");
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << bodies.size() << " float\n";
    for (size_t b = 0; b < bodies.size(); b++)
    {
      file << bodies[b].X[0] << " " << bodies[b].X[1] << " " << bodies[b].X[2] << "\n";
    }
    // Add points as vertices
    file << "VERTICES " << bodies.size() << " " << 2 * bodies.size() << "\n";
    for (size_t b = 0; b < bodies.size(); b++)
    {
      file << "1 " << b << "\n";
    }
    file.close();
  };

  auto advance_euler = [&bodies](real_t dt)
  {
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      for (int d = 0; d < 3; d++)
      {
        bodies[b].X[d] += bodies[b].velocity[d] * dt; //  Update position
        bodies[b].velocity[d] = 0;                    //  Clear velocity
        bodies[b].alpha[d] += bodies[b].dadt[d] * dt;
        bodies[b].dadt[d] = 0;
      }
    } // End loop over bodies
  };

  //! Build tree

  //! FMM evaluation
  auto BiotSavart = [&bodies]()
  {
    // start("Build tree");             // Start timer
    Cells cells = buildTree(bodies); // Build tree
    // stop("Build tree");              // Stop timer

    // start("P2M & M2M");           // Start timer
    initKernel();                 // Initialize kernel
    upwardPass(cells);            // Upward pass for P2M, M2M
    // stop("P2M & M2M");            // Stop timer
    // start("M2L & P2P");           // Start timer
    horizontalPass(cells, cells); // Horizontal pass for M2L, P2P
    // stop("M2L & P2P");            // Stop timer
    // start("L2L & L2P");           // Start timer
    downwardPass(cells);          // Downward pass for L2L, L2P
    // stop("L2L & L2P");            // Stop timer
  };

  auto Streching = [&bodies]()
  {
    // start("Build tree");             // Start timer
    Cells cells = buildTree(bodies); // Build tree
    // stop("Build tree");              // Stop timer

    // start("P2M & M2M");            // Start timer
    initKernel();                  // Initialize kernel
    upwardPassS(cells);            // Upward pass for P2M, M2M
    // stop("P2M & M2M");             // Stop timer
    // start("M2L & P2P");            // Start timer
    horizontalPassS(cells, cells); // Horizontal pass for M2L, P2P
    // stop("M2L & P2P");             // Stop timer
    // start("L2L & L2P");            // Start timer
    downwardPassS(cells);          // Downward pass for L2L, L2P
    // stop("L2L & L2P");             // Stop timer
  };

  // writeTovtk(0);

  // for (int i = 1; i < 5'000'000; i++)
  // {
  //   start("Step");
  //   BiotSavart();
  //   Streching();
  //   advance_euler(1e-4);
  //   if (i % 1000 == 0)
  //       writeTovtk(i);
  //   stop("Step");

  // }


    BiotSavart();
    Streching();
  //! Direct N-Body
  start("Direct N-Body");                  // Start timer
  Bodies jbodies = bodies;                 // Save bodies in jbodies
  int stride = bodies.size() / numTargets; // Stride of sampling
  for (int b = 0; b < numTargets; b++)
  {                                 // Loop over target samples
    bodies[b] = bodies[b * stride]; //  Sample targets
  }                                 // End loop over target samples
  bodies.resize(numTargets);        // Resize bodies
  Bodies bodies2 = bodies;          // Backup bodies
  for (size_t b = 0; b < bodies.size(); b++)
  { // Loop over bodies
    for (int d = 0; d < 3; d++)
    {
      bodies[b].velocity[d] = 0;
      bodies[b].dadt[d] = 0;
    }
  }                        // End loop over bodies
  direct(bodies, jbodies); // Direct N-Body
  stop("Direct N-Body");   // Stop timer

  //! Verify result
  // real_t pDif = 0, pNrm = 0;
  real_t VelDiff = 0, VelNorm = 0;
  real_t StrDiff = 0, StrNorm = 0;
  for (size_t b = 0; b < bodies.size(); b++)
  {
    StrDiff += (bodies[b].dadt[0] - bodies2[b].dadt[0]) * (bodies[b].dadt[0] - bodies2[b].dadt[0]) +
               (bodies[b].dadt[1] - bodies2[b].dadt[1]) * (bodies[b].dadt[1] - bodies2[b].dadt[1]) +
               (bodies[b].dadt[2] - bodies2[b].dadt[2]) * (bodies[b].dadt[2] - bodies2[b].dadt[2]);
    StrNorm += bodies[b].dadt[0] * bodies[b].dadt[0] + bodies[b].dadt[1] * bodies[b].dadt[1] +
               bodies[b].dadt[2] * bodies[b].dadt[2];
    VelDiff += (bodies[b].velocity[0] - bodies2[b].velocity[0]) * (bodies[b].velocity[0] - bodies2[b].velocity[0]) +
               (bodies[b].velocity[1] - bodies2[b].velocity[1]) * (bodies[b].velocity[1] - bodies2[b].velocity[1]) +
               (bodies[b].velocity[2] - bodies2[b].velocity[2]) * (bodies[b].velocity[2] - bodies2[b].velocity[2]);
    VelNorm += bodies[b].velocity[0] * bodies[b].velocity[0] + bodies[b].velocity[1] * bodies[b].velocity[1] +
               bodies[b].velocity[2] * bodies[b].velocity[2];
  }                                                                            // End loop over bodies & bodies2
  printf("--- %-16s ------------\n", "FMM vs. direct");                        // Print message
                                                                               // printf("%-20s : %8.5e s\n", "Rel. L2 Error (p)", sqrt(pDif / pNrm)); // Print potential error
  printf("%-20s : %8.5e s\n", "Rel. L2 Error (Str)", sqrt(StrDiff / StrNorm)); // Print force error
  printf("%-20s : %8.5e s\n", "Rel. L2 Error (Vel)", sqrt(VelDiff / VelNorm)); // Print force error
  return 0;
}
