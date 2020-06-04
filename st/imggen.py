"""Generate star tracker image"""

# built-in libraries
from math import radians, cos, sin
from random import randrange
import argparse

# external libraries
import numpy.random
import matplotlib.pyplot

# internal libraries
from .common import appmag, gcnav, tabdes


def bg_noise(N, exp=16, std=4):
    """White background noise"""
    rng = numpy.random.default_rng()
    img = rng.normal(exp, std, (N, N))
    return img

def fg_noise(img, exp=128, std=64, lo=4, hi=8):
    """Point foreground noise"""
    N, N = img.shape
    rng = numpy.random.default_rng()
    for n in range(randrange(lo, hi)):
        x, y = rng.integers(0, N - 1, (2,))
        f = int(rng.normal(exp, std))
        img[x, y] = f


def cli():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-N", "--size", default=256, type=int,
                        help="image size")
    parser.add_argument("-r", "--ra", default=45, type=int,
                        help="right ascension")
    parser.add_argument("-d", "--dec", default=45, type=int,
                        help="declination")
    parser.add_argument("-f", "--field", default="star.nav", type=str,
                        help="nav stars")
    parser.add_argument("image", default="st.png", type=str,
                        help="output image")

    return parser.parse_args()


def main(fin, fout, N=256, ra=45, dec=45):
    img = bg_noise(N)
    fg_noise(img)
    
    data = tabdes(fin, "!HHhe2s3sxxx")
    star_map = [((lon, lat), amag) for (_, lon, lat, amag, _, _) in data]
        
    p0 = (ra, dec)
    for p1, amag in star_map:
        try:
            p = gcnav(*p0, *p1, 45)
        except:
            continue
        
        az, el = map(radians, p)
        r = cos(el)
        x = N // 2 + int(N * r * sin(az))
        y = N // 2 + int(N * r * cos(az))
        flag = (
            (0 <= x < N - 1) and
            (0 <= y < N - 1)
        )
        if not flag:
            continue

        f = appmag(amag)
        img[N - y - 1, N - x - 1] =  int(f)

    matplotlib.pyplot.imsave(fout, img, vmin=0, vmax=255, cmap="gray")


if __name__ == "__main__":
    ns = cli()
    main(ns.field, ns.image, ns.size, ns.ra, ns.dec)
