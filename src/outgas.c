#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

static void rebx_calculate_outgas_forces(struct rebx_extras* const rebx, struct reb_simulation* const sim, double alpha, double n, double k, double m, double r0, double A01, double A02, double A03, const int source_index, struct reb_particle* const particles, const int N){
    const struct reb_particle source = particles[source_index];
    
    for (int i=0; i<N; i++){
        if(i == source_index) continue;
        
        struct reb_particle* p = &particles[i];
        
        double A1[3] = {source.x - p->x, source.y - p->y, source.z - p->z};
        double R = sqrt(A1[0]*A1[0] + A1[1]*A1[1] + A1[2]*A1[2]);
        
        // Normalize A1
        A1[0] /= R; A1[1] /= R; A1[2] /= R;
        
        double A2[3] = {p->vx, p->vy, p->vz};
        double norm_A2 = sqrt(A2[0]*A2[0] + A2[1]*A2[1] + A2[2]*A2[2]);
        
        // Normalize A2
        A2[0] /= norm_A2; A2[1] /= norm_A2; A2[2] /= norm_A2;
        
        double A3[3] = {
            A1[1]*A2[2] - A1[2]*A2[1],
            A1[2]*A2[0] - A1[0]*A2[2],
            A1[0]*A2[1] - A1[1]*A2[0]
        };
        A1[0] *= A01; A1[1] *= A01; A1[2] *= A01;
        A2[0] *= A02; A2[1] *= A02; A2[2] *= A02;
        A3[0] *= A03; A3[1] *= A03; A3[2] *= A03;
        
        // Scaling factor
        double scale = alpha * pow(R/r0, -m) * pow(1 + pow(R/r0, n), -k);
        
        // Apply scaling
        A1[0] *= scale; A1[1] *= scale; A1[2] *= scale;
        A2[0] *= scale; A2[1] *= scale; A2[2] *= scale;
        A3[0] *= scale; A3[1] *= scale; A3[2] *= scale;
        
        // Update particle acceleration
        p->ax += A1[0] + A2[0] + A3[0];
        p->ay += A1[1] + A2[1] + A3[1];
        p->az += A1[2] + A2[2] + A3[2];
    }
}

void rebx_outgas_forces(struct reb_simulation* const sim, struct rebx_force* const outgas_forces, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    
    double* alpha = rebx_get_param(rebx, outgas_forces->ap, "alpha");
    double* n = rebx_get_param(rebx, outgas_forces->ap, "n");
    double* k = rebx_get_param(rebx, outgas_forces->ap, "k");
    double* m = rebx_get_param(rebx, outgas_forces->ap, "m");
    double* r0 = rebx_get_param(rebx, outgas_forces->ap, "r0");
    
    if (!alpha || !n || !k || !m || !r0){
        reb_simulation_error(sim, "Missing parameters for outgas forces.\n");
        return;
    }
    
    int source_found = 0;
    for (int i = 0; i < N; i++){
        if (rebx_get_param(rebx, particles[i].ap, "outgas_source") != NULL){
            source_found = 1;
            rebx_calculate_outgas_forces(rebx, sim, *alpha, *n, *k, *m, *r0, i, particles, N);
        }
    }
    if (!source_found){
        rebx_calculate_outgas_forces(rebx, sim, *alpha, *n, *k, *m, *r0, 0, particles, N); // default to index 0
    }
}
