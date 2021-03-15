import matplotlib.pyplot as plt
from design_data import *
from math import sin, cos, asin, pi, radians, ceil
import numpy as np


e = m / 2 * (1 - x)
z2 = z1 + 1
dw = (m * (z1 + x - 2 * rc_star) - dx) / 2 + dx


def c_delta_beta(beta):
    denominator = 0
    for i in range(1, z2+1):
        deltai = i * (2 * pi) / z2 - beta
        gammai = asin(((1 - x) * sin(deltai)) / (
                      pow(1 + pow((1 - x), 2) - 2 * (1 - x) * cos(deltai), 1/2)))
        denominator += pow(e, 2) * z1 * pow(sin(deltai + gammai), 2)
    cdb = T1 / denominator
    return cdb


def ith_vertical_normal_force(i, beta):
    cdB = c_delta_beta(0)
    deltai = i * (2 * pi) / z2 - beta
    gammai = asin(((1 - x) * sin(deltai)) /
                  pow(1 + pow((1 - x), 2) - 2 * (1 - x) * cos(deltai), 1/2))
    Fnvi = cdB * e * z1 * pow(sin(deltai + gammai), 2)
    return Fnvi


def ith_horizontal_normal_force(i, beta):
    cdB = c_delta_beta(0)
    deltai = i * (2 * pi) / z2 - beta
    gammai = asin(((1 - x) * sin(deltai)) /
                  pow(1 + pow((1 - x), 2) - 2 * (1 - x) * cos(deltai), 1/2))
    Fnhi = cdB * e * z1 * sin(deltai + gammai) * cos(deltai + gammai)
    return Fnhi


def ith_normal_force(i, beta):
    Fn = pow(pow(ith_horizontal_normal_force(i, beta), 2) +
             pow(ith_vertical_normal_force(i, beta), 2), 1/2)
    return Fn


def phikj(j, beta):
    phikj = phik0 + (j - 1) * ((2 * pi) / zw) - beta / z1
    return phikj


def rkj(j, beta):
    rkj = pow((pow((dw/2), 2) - (dw/2) * (dp + 2 * e) *
               cos(phikj(j, beta) - beta) + pow(((dp + 2 * e) / 2), 2)), 1/2)
    return rkj


def psikj(j, beta):
    psikj = asin((dw/2) / rkj(j, beta) * sin(phikj(j, beta) - beta))
    return psikj


def ck_delta_beta(beta):
    denominator = 0
    for j in range(1, zw+1):
        denominator += pow(rkj(j, beta), 2) * pow(sin(psikj(j, beta)), 2)
    ckdb = (T1 * z1) / denominator
    return ckdb


def force_by_output_roller_j(j, beta):
    ckdb = ck_delta_beta(0)
    Fkj = ckdb * rkj(j, beta) * sin(psikj(j, beta))
    return Fkj


def vertical_load_on_bearing(beta):
    Fev = 0
    for i in range(1, ceil(z2/2+1)):
        Fev += ith_vertical_normal_force(i, beta)
    return Fev


def horizontal_load_on_bearing(beta):
    Feh = 0
    for i in range(1, ceil(z2/2+1)):
        Feh += ith_horizontal_normal_force(i, beta)
    for j in range(1, ceil(zw/2+1)):
        Feh += force_by_output_roller_j(j, beta)
    return Feh


def load_on_bearing(beta):
    Fe = pow(pow(horizontal_load_on_bearing(beta), 2) +
             pow(vertical_load_on_bearing(beta), 2), 1/2)
    return Fe


def plot(x, y):
    plt.axis(adjustable='box', anchor='C')
    plt.grid(linestyle='--')
    plt.xlabel('Hajtótengely forgási szöge (°)')
    plt.ylabel('Erő (N)')
    plt.plot(x, y)
    plt.show()


def coordinates_of_force_diagram_of_housing_roller_i(i):
    x = []
    y = []
    x_max = 0
    y_max = 0
    for deg in range(181):
        x.append(deg)
        y.append(round(ith_normal_force(i, radians(deg)), 1))
        if y[deg] > y_max:
            x_max = x[deg]
            y_max = y[deg]
    return x, y, x_max, y_max


def coordinates_of_force_diagram_of_output_roller_j(j):
    x = []
    y = []
    x_max = 0
    y_max = 0
    for num, deg in enumerate(range(175)):
        x.append(deg)
        y.append(round(force_by_output_roller_j(j, radians(deg)), 1))
        if y[num] > y_max:
            x_max = x[num]
            y_max = y[num]
    return x, y, x_max, y_max


def load_on_shaft_end_bearings():
    Fa = ((L1 + L2) * (load_on_bearing(0)/2) -
          L1 * (load_on_bearing(0)/2)) / (L1+L2+L3)
    return Fa


def minimum_output_pin_diameter():
    dpmin = pow((coordinates_of_force_diagram_of_output_roller_j(3)
                 [3]*(1.5*B+delta))/(0.1*1034), 1/3)
    return dpmin


def print_cycloid_profile_equation():
    print('\nCycloidal Disc Profile Equations:')
    print('X coordinates:\n', m/2, '*(', z1+1, '*sin(t)-', 1-x, '*sin(', z1+1,
          '*t)+(', 2*rc_star, '*(', 1-x, '*sin(', z1+1, '*t)-sin(t)))/sqrt(1-',
          round(2*(1-x), 5), '*cos(', z1, '*t)+', round((1-x)**2, 5), '))')
    print('Y coordinates:\n', m/2, '*(', z1+1, '*cos(t)-', 1-x, '*cos(', z1+1,
          '*t)+(', 2*rc_star, '*(', 1-x, '*cos(', z1+1, '*t)-cos(t)))/sqrt(1-',
          round(2*(1-x), 5), '*cos(', z1, '*t)+', round((1-x)**2, 5), '))')


def print_cycloid_disc_parameters():
    print('\nCycloidal Disc Parameters:')
    print('Pitch Diameter (d1):', m * z1, '[mm]')
    print('Output Pin Holes (dw):', round(dp+2*e, 1), '[mm]')
    print('Disc Thickness (B):', B, '[mm]')


def print_housing_parameters():
    print('\nHousing Parameters:')
    print('Pin Number (z2):', z1 + 1)
    print('Pin Diameter (dc):', 2 * rc_star * m, '[mm]')
    print('Pitch Diameter (d2):', m * (z1 + 1), '[mm]')
    print('Upper Limit of Dedendum Diameter (df2):',
          round(m*(z1 + 1 - 0.3 * rc_star), 5), '[mm]')
    print('Lower Limit of Dedendum Diameter (df2):',
          round(m*(3 - 2 * x + z1 - 2 * rc_star), 5), '[mm]')


def print_input_shaft_parameters():
    print('\nInput Shaft Parameters:')
    print('Eccentricity (e0):', round(m / 2 * (1 - x), 5), '[mm]')
    print('Eccentric Bearing Static Load (C0):', round(
        1.2 * load_on_bearing(0) * pow(60 * n * 10000 / 10**6, 3/10), 5), '[N]')
    print('Eccentric Bearing Working Lifespan (Lh1):', round(
        (10**6 / (60 * n)) * pow(eccentric_c0 / (1.2 * load_on_bearing(0)), 10/3)), '[h]')
    print('Shaft End Bearings Static Load (C0):', round(
        1.2 * load_on_shaft_end_bearings() * pow(60 * n * 10000 / 10**6, 3/10), 5), '[N]')
    print('Shaft End Bearing Working Lifespan (Lh1):', round(
        (10**6 / (60 * n)) * pow(shaft_end_c0 / (1.2 * load_on_shaft_end_bearings()), 10/3)), '[h]')


def print_output_shaft_parameters():
    print('\nOutput Shaft Parameters:')
    print('Output Pin Number (zw):', zw)
    print('Pitch Diameter (Dw):', round(dw), '[mm]')
    print('Output Pin Diameter (dp):', dp, '[mm]')
    print('Minimum Output Pin Diameter (dpmin):', round(
        minimum_output_pin_diameter(), 5), '[mm]')


def print_force_distribution():
    print('\nForce Distribution of the Cycloid Drive:')
    print('Vertical Force on Bearing at 0 degree (Fev):',
          round(vertical_load_on_bearing(0), 5), '[N]')
    print('Horizontal Force on Bearing at 0 degree (Feh):',
          round(horizontal_load_on_bearing(0), 5), '[N]')
    print('Force on Bearing at 0 degree (Fe):',
          round(load_on_bearing(0), 5), '[N]')
    print('Maximum Force on Single Housing Roller (Fn):', coordinates_of_force_diagram_of_housing_roller_i(15)[
        3], '[N]\tOn degree:', coordinates_of_force_diagram_of_housing_roller_i(15)[2], '[deg]')
    print('Maximum Force on Single Output Roller (Fk):', coordinates_of_force_diagram_of_output_roller_j(3)[
        3], '[N]\tOn degree:', coordinates_of_force_diagram_of_output_roller_j(3)[2], '[deg]')


if __name__ == "__main__":
    print_cycloid_profile_equation()
    print_cycloid_disc_parameters()
    print_housing_parameters()
    print_input_shaft_parameters()
    print_output_shaft_parameters()
    print_force_distribution()
    # plot(coordinates_of_force_diagram_of_housing_roller_i(15)
    #      [0], coordinates_of_force_diagram_of_housing_roller_i(15)[1])
    # plot(coordinates_of_force_diagram_of_output_roller_j(3)
    #      [0], coordinates_of_force_diagram_of_output_roller_j(3)[1])
