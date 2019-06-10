/* use as src/partition_fixed_costs.h */
#define HAVE_FIXED_COSTS 1
repartition_costs[1][0] =        277; /* sort/none */
repartition_costs[2][1] =       1913; /* self/density */
repartition_costs[2][3] =       1257; /* self/force */
repartition_costs[2][6] =         11; /* self/external_grav */
repartition_costs[2][23] =         23; /* self/bh_density */
repartition_costs[2][24] =         24; /* self/bh_swallow */
repartition_costs[2][25] =          3; /* self/do_swallow */
repartition_costs[2][26] =         18; /* self/bh_feedback */
repartition_costs[3][1] =        265; /* pair/density */
repartition_costs[3][3] =        184; /* pair/force */
repartition_costs[3][23] =         17; /* pair/bh_density */
repartition_costs[3][24] =         31; /* pair/bh_swallow */
repartition_costs[3][25] =          9; /* pair/do_swallow */
repartition_costs[3][26] =          8; /* pair/bh_feedback */
repartition_costs[9][0] =        169; /* ghost/none */
repartition_costs[12][0] =        209; /* drift_part/none */
repartition_costs[14][0] =         10; /* drift_bpart/none */
repartition_costs[15][0] =         28; /* drift_gpart/none */
repartition_costs[17][0] =         12; /* end_hydro_force/none */
repartition_costs[18][0] =         99; /* kick1/none */
repartition_costs[19][0] =        107; /* kick2/none */
repartition_costs[20][0] =        232; /* timestep/none */
repartition_costs[29][0] =         13; /* grav_end_force/none */
repartition_costs[43][0] =        315; /* bh_density_ghost/none */
