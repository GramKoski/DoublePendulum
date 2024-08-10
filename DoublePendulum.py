import turtle
import math

g = 9.8
M = 0.1
R = 0.1
pi = math.pi

THETA = pi / 4
THETADOT = 0

m = turtle.Turtle()
rod = turtle.Turtle()

s = turtle.Screen()

s.bgcolor('white')
s.setup(300,300)

m.shape('circle')
m.speed("fastest")
m.penup()
m.goto(R*1000*math.sin(THETA), -R*1000*math.cos(THETA))
m.pendown()

rod.penup()
rod.speed("fastest")
rod.shape('square')
rod.shapesize(0.1, R*100)
rod.setheading(math.degrees(THETA - pi/2))


def display(theta):
    m.goto(R*1000*math.sin(theta), -R*1000*math.cos(theta))
    rod.setheading(math.degrees(theta - pi/2))

def pendulumDynamics(theta, thetaDot):
    t = 0
    dt = 0.02
    while t < 4:
        thetaDot += -(g/R)* math.sin(theta)*dt
        theta += thetaDot*dt
        display(theta)
        print(theta)
        t += dt
    
pendulumDynamics(THETA, THETADOT)
turtle.done()


