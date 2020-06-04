python3 -m st.navstar -f star.nav
python3 -m st.geohash -f star.nav -m hash.geo
python3 -m st.imggen -f star.nav --ra=45 --dec=45 test.png
python3 -m st -f star.nav -m hash.geo test.png
