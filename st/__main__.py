"""A star tracker algorithm"""

# built-in libraries
from math import degrees, sqrt, acos, asin, atan2, log
import argparse

# external libraries
import numpy
import numpy.random
import scipy.linalg
import matplotlib.pyplot

# internal libraries
from .common import (rect, topobase, topo2rect,
                     geo2vec, vec2geo,
                     gcnav, dualtri, tabdes)


def filter_img(img, N=256, min_y=32, max_y=224):
    r, g, b, a = img[..., 0], img[..., 1], img[..., 2], img[..., 3]
    y = 255 * (0.2126 * r + 0.7152 * g + 0.0722 * b) * a
    point_list = []
    for y, x in zip(*numpy.where(numpy.logical_and(y > min_y, y < max_y))):
        r = sqrt((N // 2 - x - 1) ** 2 + (N // 2 - y - 1) ** 2)
        el = degrees(acos(r / N))
        az = degrees(atan2(N // 2 - x - 1, N // 2 - y - 1))
        point_list.append((az, el))
    else:
        return point_list


def find_matches(point_list, geo_hash):
    match_list = {}
    for xp in point_list:
        for xn in point_list:
            if xp is xn:
                continue
            
            try:
                p = gcnav(*xp, *xn)
            except:
                continue
            
            b = topobase(p, left=True)
            for p in point_list:
                if p is xp or p is xn:
                    continue
                
                try:
                    p = gcnav(*xp, *p)
                except:
                    continue

                p = topo2rect(p, b)
                pair = geo_hash.get(p)
                if pair is not None:
                    if (xp, xn) in match_list:
                        vote_hash = match_list[xp, xn]
                    elif (xn, xp) in match_list:
                        vote_hash = match_list[xn, xp]
                        pair = tuple(reversed(pair))
                    else:
                        vote_hash = match_list.setdefault((xp, xn), {})
                    vote_hash.setdefault(pair, 0)
                    vote_hash[pair] += 1
    return match_list


def tally_vote(match_list, star_map):
    bar_list = []
    for (xp, xn), vote_hash in match_list.items():
        for pair, votes in vote_hash.items():
            if pair is None:
                continue
            p1, star1 = star_map[pair[0]]
            p2, star2 = star_map[pair[1]]

            try:
                r_hat = dualtri(xp, xn, p1, p2)
            except:
                print("no votes for `%s %s` and `%s %s`" %
                      (*map(bytes.decode, star1),
                       *map(bytes.decode, star2)))
                continue

            print("voting %d for match on `%s %s` and `%s %s`" %
                  (votes,
                   *map(bytes.decode, star1),
                   *map(bytes.decode, star2)))
            bar_list.append((r_hat, votes))
    return bar_list


def calc_stats(bar_list, C=0.95, max_q=5):
    N = sum(n for (_, n) in bar_list)
    if N <= 1:
        return  # not enough measurements
    print("total votes %d" % N)
    
    r_bar = sum(n * _hat for (_hat, n) in bar_list) / N
    r = scipy.linalg.norm(r_bar)
    r_hat = r_bar / r
    
    d = 1 - sum(n * numpy.dot(_hat, r_hat) ** 2
                for (_hat, n) in bar_list) / N
    std_dev = sqrt(d / N) / r
    print("std dev %.2e [deg]" % degrees(std_dev))
    
    try:
        q = degrees(asin(sqrt(- log(1 - C)) * std_dev))
    except:
        return  # standard deviation too high
    if q > max_q:
        return  # confidence cone too wide

    ra, dec = vec2geo(r_hat)
    print("ra %.2f, dec %.2f (+/-%.2e) [deg]" % (ra, dec, q))

    return r_hat


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", default="star.nav", type=str,
                        help="nav stars")
    parser.add_argument("-m", "--metric", default="hash.geo", type=str,
                        help="geo hash")
    parser.add_argument("image", default="st.png", type=str,
                        help="input image")

    return parser.parse_args()


def main(navstar, geohash, fin):
    img = matplotlib.pyplot.imread(fin)
    point_list = filter_img(img)
    
    data = tabdes(geohash, "!bbHHxB")
    geo_hash = {rect(x, y): (i, j) for (x, y, i, j, _) in data}

    match_list = find_matches(point_list, geo_hash)
    
    data = tabdes(navstar, "!HHhe2s3sxxB")
    star_map = [((lon, lat), (x, abbr))
                for (_, lon, lat, _, x, abbr, _) in data]

    bar_list = tally_vote(match_list, star_map)

    calc_stats(bar_list)


if __name__ == "__main__":
    ns = cli()
    main(ns.field, ns.metric, ns.image)
