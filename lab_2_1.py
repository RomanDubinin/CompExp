from bokeh.plotting import figure, output_file, show
import bokeh
import itertools
import numpy as np
import math
import sys

y0 = 0.1

def function_to_integrate(x, y):
    return 30*y*(x - 0.2)*(x - 0.7)

def df_x(x, y):
    return y*(60*x - 27)

def df_xx(x, y):
    return y*60

def df_y(x, y):
    return 30*(x - 0.7)*(x - 0.2)

def df_yy(x, y):
    return 0

def df_xy(x, y):
    return 60*x - 27

def original_func(x):
    return 0.1 * math.exp(x*(10*x*x - 13.5*x + 4.2))

def explicit_eiler(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(len(x_points)):
        next_y = y_points[k] + f(x_points[k], y_points[k]) * h
        y_points.append(next_y)
    
    return y_points[:-1]

def implicit_eiler(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(len(x_points)-1):
        subsidiary_y = y_points[k] + f(x_points[k], y_points[k])*h
        next_y = y_points[k] + f(x_points[k+1], subsidiary_y) * h
        y_points.append(next_y)
    
    return y_points

def coshi(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(len(x_points)-1):
        subsidiary_y = y_points[k] + f(x_points[k], y_points[k])*h/2
        next_y = y_points[k] + f(x_points[k] + h/2, subsidiary_y) * h
        y_points.append(next_y)
    
    return y_points

def eiler_With_recount(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(len(x_points)-1):
        subsidiary_y = y_points[k] + f(x_points[k], y_points[k])*h
        next_y = y_points[k] + (f(x_points[k], y_points[k]) + f(x_points[k+1], subsidiary_y))*h/2
        y_points.append(next_y)
    
    return y_points

def runge_kutta(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(len(x_points)):
        k1 = h * f(x_points[k],       y_points[k])
        k2 = h * f(x_points[k] + h/2, y_points[k] + k1/2)
        k3 = h * f(x_points[k] + h/2, y_points[k] + k2/2)
        k4 = h * f(x_points[k] + h,   y_points[k] + k3)

        next_y = y_points[k] + (k1 + 2*k2 + 2*k3 + k4)/6
        y_points.append(next_y)
    
    return y_points[:-1]

def extrapolation_adams(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    y_points.append(y0 + h*f(x_points[0], y0))
    for k in range(1, len(x_points)):
        next_y = y_points[k] + (3/2 * f(x_points[k], y_points[k]) - 1/2 * f(x_points[k-1], y_points[k-1])) * h
        y_points.append(next_y)
    
    return y_points[:-1]

def teilor_3(f, df_x, df_y, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(0, len(x_points)):
        next_y = y_points[k] \
        + f(x_points[k], y_points[k])*h \
        + (df_x(x_points[k], y_points[k]) \
        + df_y(x_points[k], y_points[k]) * f(x_points[k], y_points[k]))*(h**2)/2
        y_points.append(next_y)
    
    return y_points[:-1]

def teilor_4(f, df_x, df_y, df_xx, df_yy, df_xy, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for k in range(len(x_points)):
        second_summand = f(x_points[k], y_points[k])*h
        third_summand = (df_x(x_points[k], y_points[k]) + df_y(x_points[k], y_points[k]) * f(x_points[k], y_points[k]))*(h**2)/2
        fourth_summand = (df_xx(x_points[k], y_points[k]) + \
                         2*f(x_points[k], y_points[k])*df_xy(x_points[k], y_points[k]) \
                         + df_yy(x_points[k], y_points[k])*(f(x_points[k], y_points[k])**2) \
                         + df_y(x_points[k], y_points[k])*(df_x(x_points[k], y_points[k]) \
                         + df_y(x_points[k], y_points[k])*f(x_points[k], y_points[k])))*(h**3)/3
        next_y = y_points[k] + second_summand + third_summand + fourth_summand
        y_points.append(next_y)
    
    return y_points[:-1]

n = int(sys.argv[1])

output_file("{0}_points.html".format(n))   

p = figure(title="30*y*(x - 0.2)*(x - 0.7), " + str(n) + "points", plot_width=1100, plot_height=600)

x_points = np.linspace(0, 1, n+1, endpoint=True)
#x_points = [x/n for x in range(n+1)]

y_points = explicit_eiler(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "red", legend="Эйлер")
p.circle(x=x_points, y=y_points, color = "red", legend="Эйлер")

y_points = coshi(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "blue", legend="Коши")

y_points = implicit_eiler(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "red", line_dash=[4, 4], legend="неявный Эйлер")

y_points = eiler_With_recount(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "red", legend="Эйлер с пересчётом")
p.square(x=x_points, y=y_points, color = "red", legend="Эйлер с пересчётом")

y_points = runge_kutta(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "purple", legend="Рунге-Кутта 4 го порядка")

y_points = extrapolation_adams(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "orange", legend="экстраполяционный метод Адамса(k=2)")

y_points = teilor_3(function_to_integrate, df_x, df_y, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "green", legend="Тейлор 3 го порядка")
p.triangle(x=x_points, y=y_points, color = "green", legend="Тейлор 3 го порядка")

y_points = teilor_4(function_to_integrate, df_x, df_y, df_xx, df_yy, df_xy, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "green", legend="Тейлор 4 го порядка")
p.square(x=x_points, y=y_points, color = "green", legend="Тейлор 4 го порядка")


y_points = [original_func(x) for x in x_points]
p.line(x=x_points, y=y_points, color = "black", legend="origin")
p.circle(x=x_points, y=y_points, color = "black", legend="origin")

p.legend.orientation = "top_left"
show(p)