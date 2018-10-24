import math


def atan3(p1, p2):
    '''
    Return the angle between the two points in the
    range between 0 and 2 * pi.

    The angle is first obtained from math.atan2. The if it is
    less 0, it is added to 2 * pi.

    @param p1: The first point
    @param p2: The second point
    @return: The angle between the two points.
    '''
    a = math.atan2(p1, p2)
    if a < 0:
        return 2 * math.pi + a
    else:
        return a


def clock_angle(a1, a2):
    '''
    The amount one needs to rotate from angle a1
    to angle a2 in a clockwise direction. I.e.
    increasing a1 until it reaches a2.

    @param a1: The first angle.
    @param a2: The second angle.
    '''
    if a2 < a1:
        a2 += 2. * math.pi
    return a2 - a1
