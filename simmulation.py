import numpy as np
import csv
from numerical_methods import shooting, continuation_parameter, augmented_frechet_matrix, analyse_point
from additional_functions import *


def bifurcation_analysis(time_steps, system, bc, bc_params, n_eq, params, cur_p, p_max,
                          step, max_step, file):
    eps = 10**(-4)
    y_approx = np.random.uniform(-1, 1, n_eq)
    Ys = []
    Ps = []
    fieldnames = make_fieldnames(n_eq)
    create_csv(file, fieldnames)
    
    while True:
        point_type = 'N'
        newton_steps, ys, det, sv = shooting(time_steps, y_approx, system, params, bc, bc_params)
        
        if newton_steps < 3:
            step = min(step * 2, max_step)
        elif newton_steps > 7:
            step = step / 2
        
        if np.abs(det) < eps:
            point_type = 'C'
            aug_F = augmented_frechet_matrix(time_steps, ys, system, params, cur_p, bc, bc_params)
            bif, dets = analyse_point(aug_F)
            if bif: point_type = 'B'

        info = make_info(fieldnames, params, det, sv, point_type)
        write_csv(file, info)
        
        y_approx = continuation_parameter(Ps, Ys, params[cur_p], y_approx, params[cur_p] + step)
        params[cur_p] += step
        if np.abs(params[cur_p]) > p_max: break
        
    
    