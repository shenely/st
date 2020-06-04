"""Common utilities"""

# built-in libraries
from math import *
from time import time
from collections import namedtuple
from struct import Struct

# external libraries
import numpy

# internal libraries
# ...

rect = namedtuple("rect", ("x", "y"))
basis = namedtuple("basis", ("o", "i", "j", "M"))


def appmag(x, xlim=(8, -4), ylim=(64, 192)):
    """Apparent magnitude to grayscale"""
    min_x, max_x = xlim
    min_y, max_y = ylim
    return min_y + (max_y - min_y) * (x - min_x) / (max_x - min_x)


def gcnav(lon0, lat0, lon1, lat1, max_gcd=180):
    """Great circle navigation"""
    if abs(lat1 - lat0) > max_gcd:
        raise

    lon0, lat0, lon1, lat1 = map(radians, (lon0, lat0, lon1, lat1))

    # precompute cos/sin
    cos_dlon = cos(lon1 - lon0)
    sin_dlon = sin(lon1 - lon0)
    cos_lat0 = cos(lat0)
    sin_lat0 = sin(lat0)
    cos_lat1 = cos(lat1)
    sin_lat1 = sin(lat1)

    # great circle distance
    gcd = degrees(acos(sin_lat0 * sin_lat1 + cos_lat0 * cos_lat1 * cos_dlon))
    if gcd > max_gcd:
        raise

    # azimuth/elevtion
    el = 90 - gcd
    az = degrees(atan2(cos_lat1 * sin_dlon,
                       cos_lat0 * sin_lat1 - sin_lat0 * cos_lat1 * cos_dlon))

    return az, el


def dualtri(xp, xn, p1, p2):
    """Complete dual triangles"""
    _, A = gcnav(*p1, *p2)
    _, B = gcnav(*xp, *xn)
    if round(A - B) != 0.0:  # sanity check
        raise

    # pre-compute sin/cos
    az1, el1 = map(radians, xp)
    az2, el2 = map(radians, xn)
    ra1, dec1 = map(radians, p1)
    ra2, dec2 = map(radians, p2)
    sin_A = sin(radians(A))
    cos_dec1 = cos(dec1)
    sin_dec1 = sin(dec1)
    cos_dec2 = cos(dec2)
    sin_dec2 = sin(dec2)
    sin_dra = sin(ra2 - ra1)
    cos_el1 = cos(el1)
    sin_el1 = sin(el1)
    cos_el2 = cos(el2)
    sin_el2 = sin(el2)
    sin_daz = sin(az1 - az2)

    # complete individual trianges
    th = atan2(cos_dec1 * cos_dec2 * sin_dra,
               sin_dec2 - sin_dec1 * sin_A)
    ph = atan2(cos_el1 * cos_el2 * sin_daz,
               sin_el2 - sin_el1 * sin_A)

    # calculate unknowns
    dec = pi / 2 - acos(sin_dec1 * sin_el1 +
                        cos_dec1 * cos_el1 * cos(th - ph))
    ra = ra1 + atan2(cos_el1 * cos_dec1 * sin(th - ph),
                     sin_el1 - sin_dec1 * sin(dec))
    r_bar = geo2vec(map(degrees, (ra, dec)))
    return r_bar


def topobase(p, left=False):
    """Topocentric coordinate system"""
    az, el = map(radians, p)
    
    r = cos(el)
    x = r * sin(az)
    y = r * cos(az)

    # XXX only one point is necessary since the other point is assumed
    # ... at the origin
    o = rect(x / 2, y / 2)  # origin
    i = rect(- x / r, - y / r)  # x-axis
    j = rect(*((- i.y, i.x)
               if not left else
               (i.y, - i.x)))  # y-axis
    M = 2 / r  # scale
    
    return basis(o, i, j, M)


def topo2rect(p, b, N=64):
    """Topocentric to cartesian coordinates"""
    az, el = map(radians, p)
    
    r = cos(el)
    x = (r * sin(az) - b.o.x) * b.M
    y = (r * cos(az) - b.o.y) * b.M
    
    return rect(int(N * (x * b.i.x + y * b.i.y)),
                int(N * (x * b.j.x + y * b.j.y)))


def geo2vec(p):
    """Geocentric coordinates to unit vector"""
    ra, dec = map(radians, p)
    return numpy.array([cos(dec) * cos(ra),
                        cos(dec) * sin(ra),
                        sin(dec)])

def vec2geo(r_hat):
    """Unit vector to geocentric coordinates"""
    ra = degrees(atan2(r_hat[1], r_hat[0]))
    dec = degrees(atan2(r_hat[2], sqrt(r_hat[0] ** 2 + r_hat[1] ** 2)))
    return ra, dec


# struct("!BiHBxxxB", ...)  # type, time, count, length, lrc
#          142124 11 = 16
# struct("!4s", ...)  # crc


def tabser(filename, body, data):
    """Serialize data to file"""
    # XXX checksums ignored
    head = Struct("!BiHBxxxB")
    body = Struct(body)
    # foot = Struct("!4s")

    buffer = bytearray([0] * (2 ** 16))
    head.pack_into(buffer, 0, 0, int(time()), len(data), body.size, 0),
    offset = head.size
    for row in data:
        body.pack_into(buffer, offset, *row, 0)
        offset += body.size
    else:
        print("write %d rows" % len(data))
        # offset = 2 ** 16 - foot.size
        # foot.pack_into(buffer, offset, bytes([0, 0, 0, 0]))
    with open(filename, "wb") as f:
        f.write(buffer)


def tabdes(filename, body):
    """Deserialize file to data"""
    # XXX checksums ignored
    head = Struct("!BiHBxxxB")
    body = Struct(body)
    # foot = Struct("!4s")

    data = []
    with open(filename, "rb") as f:
        buffer = f.read()
    _, _, count, length, _ = head.unpack_from(buffer, 0)
    offset = head.size
    for i in range(count):
        row = body.unpack_from(buffer, offset)
        data.append(row)
        offset += body.size
    else:
        print("read %d rows" % len(data))
        # offset = 2 ** 16 - foot.size
        # _, foot.unpack_from(buffer, offset))
        return data
