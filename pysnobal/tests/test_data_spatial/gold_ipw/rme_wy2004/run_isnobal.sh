# RME 25 year run water year 2004
# inputs from /data/archive/RCEW/rme/25yr-run+Adam/wy04/spatial/data/input.04_v5.aw
# precip from /data/archive/RCEW/rme/25yr-run+Adam/wy04/spatial/data/precip/ppt/ppt.4b*
# init from /data/archive/RCEW/rme/25yr-run+Adam/wy04/spatial/data/init/init.ipw
# Ran in 16:55
time isnobal -v -t 60,15,1 -T 60,10,1 -n 8783 -i input/in -I init.ipw \
    -p ppt_desc.txt -O 24 -d 0.25 -e output/em -s output/snow \
    -P 1 -b 16 -M 0.01