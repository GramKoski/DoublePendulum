import vpython

def animate(t, x1, y1, x2, y2):
    mass1 = vpython.sphere(color = vpython.color.red, radius = 0.3, make_trail = False)
    mass2 = vpython.sphere(color = vpython.color.green, radius = 0.3, make_trail = True)

    rod1 = vpython.cylinder(pos = vpython.vector(0,0,0), axis = vpython.vector(0,0,0), radius = 0.05)
    rod2 = vpython.cylinder(pos = vpython.vector(0,0,0), axis = vpython.vector(0,0,0), radius = 0.05)

    i = 0
    while t[i] <= 39.9:
        vpython.rate(25)
        mass1.pos= vpython.vector(x1[i], y1[i], 0)
        mass2.pos = vpython.vector(x2[i], y2[i], 0)
        rod1.axis = vpython.vector(x1[i], y1[i], 0)
        rod2.pos = vpython.vector(x1[i], y1[i], 0)
        rod2.axis = vpython.vector(x2[i]-x1[i], y2[i]-y1[i], 0)
        i += 1

