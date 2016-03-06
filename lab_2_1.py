from bokeh.plotting import figure, output_file, show
import bokeh
import itertools
import numpy as np


y0 = 0.1
n = 100

def function_to_integrate(x, y):
    return 30*y*(x - 0.2)*(x - 0.7)

def explicit_eiler(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for i in range(len(x_points)):
        next_y = y_points[i] + f(x_points[i], y_points[i]) * h
        y_points.append(next_y)
    
    return y_points[:-1]

def coshi(f, x_points, y0, h):
    y_points = []
    y_points.append(y0)
    for i in range(len(x_points)):
        subsidiary_y = y_points[i] + f(x_points[i], y_points[i])*h/2
        next_y = y_points[i] + f(x_points[i] + h/2, subsidiary_y) * h
        y_points.append(next_y)
    
    return y_points[:-1]


output_file("lab_2_1.html")   
p = figure(title="30*y*(x - 0.2)*(x - 0.7), " + str(n) + "points", plot_width=500, plot_height=500)
    
x_points = np.linspace(0, 1, n, endpoint=True)

y_points = explicit_eiler(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "red")

y_points = coshi(function_to_integrate, x_points, y0, 1/n)
p.line(x=x_points, y=y_points, color = "blue")



show(p)