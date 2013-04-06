# all of this code is completely useless.
# I'm just messing around right now.

class Point:
    long x = 0
    long y = 0

class ellipticCurve:
    long a = 1
    long b = 0

def addTwoPoints(ellipticCurve E, Point p1, Point p2):
    long slope = (p2.y-p1.y)/(p2.x-p1.x)
    long b = p2.y-slope*p2.x
    