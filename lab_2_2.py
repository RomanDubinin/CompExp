from bokeh.plotting import figure, output_file, show
import bokeh
import itertools
import numpy as np
import math
import sys

N = 4
l = 2 + 0.1 * N

def F(x, y, d_y):
    return y + l*x*(1-x) + 2*l + 2

def G(x, y, d_y):
    return y

def real_func2(x):
    return (1/2)*e**(-x-1) * (2*e**(x+1) *(l*(x*x - x + 4) - 2) - 2*(l-1)*e**(2*x + 1) + e*(2-6*l) - e**(2*x) + e**(2*x + 2) - e**2 + 1)


def real_func(x):
    return 2.4*x*x - 2.4*x + 0.999999*e**x + e**(-x) - 2

def real_func_method(x_points):
    y_points = []
    for i in range(len(x_points)):
        y_points.append(real_func(x_points[i]))

    return y_points




def explicit_eiler(f, x_points, y_0, z_0, h):
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points)-1):
        new_y = y_points[i] + h * z_points[i]
        new_z = z_points[i] + h * f(x_points[i], y_points[i], z_points[i])
        y_points.append(new_y)
        z_points.append(new_z)

    return y_points

def recount_eiler(f, x_points, y_0, z_0, h):
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points) - 1):
        new_y = y_points[i] + (h/2) * (2*z_points[i] + h*f(x_points[i], y_points[i], z_points[i]))
        new_z = z_points[i] + (h/2) * (f(x_points[i], y_points[i], z_points[i]) + f(x_points[i+1], y_points[i] + h*z_points[i], z_points[i] + h*f(x_points[i], y_points[i], z_points[i])))
        y_points.append(new_y)
        z_points.append(new_z)

    return y_points


def recount_eiler_MY(f, x_points, y_0, z_0, h):
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points) - 1):
        new_y = y_points[i] + (h/2) * (z_points[i] + y_points[i] + h*z_points[i])
        new_z = z_points[i] + (h/2) * (f(x_points[i], y_points[i], z_points[i]) + x_points[i] + h*f(x_points[i], y_points[i], z_points[i]))
        y_points.append(new_y)
        z_points.append(new_z)

    return y_points

def runge_kutta(f, x_points, y_0, z_0, h):
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points) - 1):
        k1 = h*f(x_points[i],       y_points[i],        z_points[i])
        k2 = h*f(x_points[i] + h/2, y_points[i] + k1/2, z_points[i])
        k3 = h*f(x_points[i] + h/2, y_points[i] + k2/2, z_points[i])
        k4 = h*f(x_points[i] + h,   y_points[i] + k3,   z_points[i])

        new_y = y_points[i] + h*z_points[i]
        new_z = z_points[i] + (k1 + 2*k2 + 2*k3 + k4)/6

        y_points.append(new_y)
        z_points.append(new_z)

    return y_points


def get_next_M(x_points, F, method, y_0, dy_n, last_M, h):
    y_points = method(F, x_points, y_0, last_M, h)
    fi = ((y_points[-1] - y_points[-2]) / h) - dy_n

    y_points = method(G, x_points, 0, 1, h)
    d_fi = (y_points[-1] - y_points[-2]) / h

    return last_M - fi/d_fi

n = int(sys.argv[1])
e = 2.71828





x_points = np.linspace(0, 1, n+1, endpoint=True)
y_0 = 0
right_boundary = e - 1/e + l
h = 1/n

methods = [
    {"f": explicit_eiler, "name": "Эйлер"},
    {"f": recount_eiler,  "name": "Эйлер с пересчётом"},
    {"f": runge_kutta,    "name": "Рунге-Кутта"}
]
colors = ["blue", "red", "purple", "orange", "green"]
for i in range(len(methods)):
    p = figure(title="", plot_width=1100, plot_height=600)

    y_points = real_func_method(x_points)
    p.line(x=x_points, y=y_points, color="black", legend="real", line_width=2)

    last_M = 0
    y_points = [methods[i]["f"](F, x_points, y_0, last_M, h)]
    next_M = get_next_M(x_points, F, methods[i]["f"], y_0, right_boundary, last_M, 1/n)

    while abs(next_M - last_M) > 0.01:
        y_points.append(methods[i]["f"](F, x_points, y_0, next_M, h))
        last_M = next_M
        next_M = get_next_M(x_points, F, methods[i]["f"], y_0, right_boundary, last_M, 1/n)

    for j in range(len(y_points)):
        p.line(x=x_points, y=y_points[j], color="red", legend=methods[i]["name"])
    p.circle(x=x_points, y=y_points[-1], color="red", legend=methods[i]["name"])

    output_file("{0}_points.html".format(i))
    p.legend.orientation = "top_left"
    show(p)


