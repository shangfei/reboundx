/** * @file outgassing_force.c
 * @brief   A general central force.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Central Force$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_central_force`
 * Python Example          `CentralForce.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CentralForce.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds a general central acceleration of the form a=Acentral*r^gammacentral, outward along the direction from a central particle to the body.
 * Effect is turned on by adding Acentral and gammacentral parameters to a particle, which will act as the central body for the effect,
 * and will act on all other particles.
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * Acentral (double)             Yes         Normalization for central acceleration.
 * gammacentral (double)         Yes         Power index for central acceleration.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"
#include "tools.h"
#define DEBUG
// #define VARIATION

#define NEWFORM

#define NOD 10
double odarr[NOD];

// static void rebx_calculate_outgassing_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double DT, const int source_index){
//     for (int i=1; i<N; i++){
//         const struct reb_particle p = particles[i];
//         const struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, p, particles[0]);
//         const double dx = p.x - source.x;
//         const double dy = p.y - source.y;
//         const double dz = p.z - source.z;
//         const double r2 = dx*dx + dy*dy + dz*dz;
//         const double prefac = A*pow(r2, (gamma-1.)/2.);

//         particles[i].ax += prefac*dx;
//         particles[i].ay += prefac*dy;
//         particles[i].az += prefac*dz;
//         particles[source_index].ax -= p.m/source.m*prefac*dx;
//         particles[source_index].ay -= p.m/source.m*prefac*dy;
//         particles[source_index].az -= p.m/source.m*prefac*dz;
//     }
// }

    // const double yr2s = 31557600.;
    // const double au2m = 149597870700.;
    // const double a1 = 1.e-4/au2m*yr2s*yr2s/4./M_PI/M_PI;// 1.686e-5;

void rebx_outgassing_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    const double* const DeltaT  = rebx_get_param(sim->extras, force->ap, "DeltaT");
    const double* const alpha   = rebx_get_param(sim->extras, force->ap, "alpha");
    const double* const r0      = rebx_get_param(sim->extras, force->ap, "r0");
    const double* const a1      = rebx_get_param(sim->extras, force->ap, "a1");
    const double* const a2      = rebx_get_param(sim->extras, force->ap, "a2");
    // const double r0 = 2.808;
    // const double alpha = 0.111262;
    const double m = 2.15;
    const double n = 5.093;
    const double k = 4.6142;
    // const double a1 = 1.e-4*168.625;// 1.686e-5;
    // const double a2 = 0*1.e-6;
    double od, nd, oma, ota, esq, ma;
#ifdef DEBUG
    const char filename[30] = "acc.txt";
    const FILE* of = fopen(filename, "a");
#endif
#ifdef NEWFORM
    const double omega_s = -453.037;
#endif
    if (DeltaT != NULL){
        const struct reb_particle primary = particles[0];
        for (int i=1; i<2; i++){
            if(particles[i].hash == "jupiter") return;
            const struct reb_particle p = particles[i];
            const struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, p, primary);
/*             if(o.M < 0.) {
                oma = o.M + 2*M_PI - o.n*(*DeltaT);
            }
            else{
                oma = o.M - o.n*(*DeltaT);
            }
            if(oma < 0.) oma += 2*M_PI;
            // if(oma < -M_PI) oma += M_PI;
            esq = o.e*o.e;
            ota = oma + (2.*o.e - 0.25*o.e*esq)*sin(oma) + 1.25*esq*sin(2.*oma)+13./12.*o.e*esq*sin(3.*oma);
            // ota = o.f;
            od = o.a*(1-esq)/(1+o.e*cos(ota));
            // od = o.d */;

            // const double prefac = 0;
            // const double prefac = 2.11;
            if(*DeltaT == 0.){
                od = o.d;
            } else if(sim->steps_done < NOD){
                odarr[sim->steps_done] = o.d;
                od = o.d;
            } else {
                od = odarr[sim->steps_done % NOD];
                odarr[sim->steps_done % NOD] = o.d;
            }

#ifdef NEWFORM
            if(o.M < 0.) {
                ma = o.M + 2*M_PI;
            }
            else {
                ma = o.M;
            }
            double prefac = 0;
            if(ma > M_PI) {
                if(o.d < 0.919) prefac = 0.03869 * pow(o.d, -0.405);
            } else {
                prefac = 0.249124 * exp(-4.01221 * o.d);
            }
#else
#ifdef VARIATION
            if(o.M < 0.) {
                ma = o.M + 2*M_PI;
            }
            else {
                ma = o.M;
            }
            double prefac = (*alpha)*pow(od/(*r0), -m) * pow(1 + pow(od/(*r0), n), -k);
            prefac *= exp(-4.*ma/M_PI);
#else
            const double prefac = (*alpha)*pow(od/(*r0), -m) * pow(1 + pow(od/(*r0), n), -k);
#endif //VARIATION
#endif //NEWFORM

            const double dx = p.x - primary.x;
            const double dy = p.y - primary.y;
            const double dz = p.z - primary.z;

            const double dvx = p.vx - primary.vx;
            const double dvy = p.vy - primary.vy;
            const double dvz = p.vz - primary.vz;
            const double hx = (dy*dvz - dz*dvy); 					//angular momentum vector
            const double hy = (dz*dvx - dx*dvz);
            const double hz = (dx*dvy - dy*dvx);
            const double tx = (hy*dz - hz*dy)/o.h;
            const double ty = (hz*dx - hx*dz)/o.h;
            const double tz = (hx*dy - hy*dx)/o.h;
#ifdef NEWFORM
            // double ang = 0.5*omega_s*sim->dt;
            double ang = -M_PI/2.;
#endif

#ifdef DEBUG
            // fprintf(of, "%e %e %e %e %e %e %e %e ", sim->t, particles[i].ax, particles[i].ay, particles[i].az, prefac, dx, o.d, sim->G);
#endif

#ifdef NEWFORM
            // particles[i].ax += prefac*dx/o.d;
            // particles[i].ay += prefac*dy/o.d;
            // particles[i].az += prefac*dz/o.d;
            particles[i].ax += prefac*(cos(ang)*dx + sin(ang)*tx)/o.d;
            particles[i].ay += prefac*(cos(ang)*dy + sin(ang)*ty)/o.d;
            particles[i].az += prefac*(cos(ang)*dz + sin(ang)*tz)/o.d;
#else
            particles[i].ax += prefac*(*a1*dx + *a2*tx)/o.d;
            particles[i].ay += prefac*(*a1*dy + *a2*ty)/o.d;
            particles[i].az += prefac*(*a1*dz + *a2*tz)/o.d;
#endif

#ifdef DEBUG
                // fprintf(of, "%e %e %e\n", particles[i].ax, particles[i].ay, particles[i].az);
                // fprintf(of, "%e %e %e\n", prefac*dx/o.d, prefac*dy/o.d, prefac*dz/o.d);
                fprintf(of, "%e %e %e %e %e %e %e\n", sim->t, o.e, p.x, p.y, p.z, o.a, o.v);
                fclose(of);
#endif
        }
    }
}

// static double rebx_calculate_central_force_potential(struct reb_simulation* const sim, const double A, const double gamma, const int source_index){
//     const struct reb_particle* const particles = sim->particles;
// 	const int _N_real = sim->N - sim->N_var;
//     const struct reb_particle source = particles[source_index];
//     double H = 0.;
// 	for (int i=0;i<_N_real;i++){
// 		if(i == source_index){
//             continue;
//         }
//         const struct reb_particle p = particles[i];
//         const double dx = p.x - source.x;
//         const double dy = p.y - source.y;
//         const double dz = p.z - source.z;
//         const double r2 = dx*dx + dy*dy + dz*dz;

//         if (fabs(gamma+1.) < DBL_EPSILON){ // F propto 1/r
//             H -= p.m*A*log(sqrt(r2));
//         }
//         else{
//             H -= p.m*A*pow(r2, (gamma+1.)/2.)/(gamma+1.);
//         }
//     }		
//     return H;
// }

// double rebx_central_force_potential(struct rebx_extras* const rebx){
//     if (rebx->sim == NULL){
//         rebx_error(rebx, ""); // rebx_error gives meaningful err
//         return 0;
//     }
//     struct reb_simulation* sim = rebx->sim;
//     const int N_real = sim->N - sim->N_var;
//     struct reb_particle* const particles = sim->particles;
//     double Htot = 0.;
//     for (int i=0; i<N_real; i++){
//         const double* const Acentral = rebx_get_param(rebx, particles[i].ap, "Acentral");
//         if (Acentral != NULL){
//             const double* const gammacentral = rebx_get_param(rebx, particles[i].ap, "gammacentral");
//             if (gammacentral != NULL){
//                 Htot += rebx_calculate_central_force_potential(sim, *Acentral, *gammacentral, i);
//             }
//         }
//     }
//     return Htot;
// }

// double rebx_central_force_Acentral(const struct reb_particle p, const struct reb_particle primary, const double pomegadot, const double gamma){
//     struct reb_simulation* sim = p.sim;
//     const double G = sim->G;
//     const struct reb_orbit o = reb_tools_particle_to_orbit(G, p, primary);
//     if (fabs(gamma+2.) < DBL_EPSILON){  // precession goes to 0 at r^-2, so A diverges for gamma=-2
//         reb_error(sim, "Precession vanishes for force law varying as r^-2, so can't initialize Acentral from a precession rate for gamma=-2)\n");
//         return 0.;
//     }
//     return G*primary.m*pomegadot/(1.+gamma/2.)/pow(o.d, gamma+2.)/o.n;
// }
