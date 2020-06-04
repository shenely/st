"""Build geometric hash table"""

# built-in libraries
from math import pi, radians
import argparse

# external libraries
# ...

# internal libraries
from .common import gcnav, topobase, topo2rect, tabser, tabdes

# struct("!bbHHxB", ...)  # x, y, 1st, 2nd, lrc
#          22222 3 111 = 8


def nn(star_map, max_gcd=45):
    """Nearest neighbor"""
    nn_map = {}
    for i, p0 in enumerate(star_map):
        for j, p1 in enumerate(star_map):
            if i == j:  # not neighbor of itself
                continue
            
            try:
                p = gcnav(*p0, *p1, max_gcd)
            except:
                continue
            
            nn_map.setdefault(i, {})[j] = p
    return nn_map


def geohash(nn_map):
    """Geometric hashing"""
    geo_hash = {}
    for i, star_map in nn_map.items():
        for j, p1 in star_map.items():
            b = topobase(p1)  # feature pair basis
            
            for k, p2 in star_map.items():
                if k == i or k == j:  # stars must be distinct
                    continue
                
                p = topo2rect(p2, b)
                # NOTE cut-off to store index in one bytes
                if abs(p.x) >= 128 or abs(p.y) >= 128:
                    continue
                
                if p in geo_hash:
                    print("dup (%d, %d)" % (*p,))
                print("hash (%d, %d) -> (%d, %d)" % (*p, i, j))
                geo_hash[p] = (i, j)
    return geo_hash


def cli():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", "--field", default="star.nav", type=str,
                        help="nav star")
    parser.add_argument("-m", "--metric", default="hash.geo", type=str,
                        help="geo hash")

    return parser.parse_args()

def main(fin, fout):
    data = tabdes(fin, "!HHhe2s3sxxB")
    star_map = [(lon, lat) for (_, lon, lat, _, _, _, _) in data]

    nn_map = nn(star_map)
    geo_hash = geohash(nn_map)
    
    data = [(p.x, p.y, i, j) for (p, (i, j)) in geo_hash.items()]
    tabser(fout, "!bbHHxB", data)


if __name__ == "__main__":
    ns = cli()
    main(ns.field, ns.metric)
