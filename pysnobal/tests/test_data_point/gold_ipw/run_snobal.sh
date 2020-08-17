# normal output
snobal -z 2061 -m 0.01 -d 0.25 -t 60,60,15,1 -T 60,10,1 -s gold.snow.properties.input -h gold.inheight.input -p gold.snobal.ppt.input -i gold.snobal.data.input.short -O data -o gold.snobal.out -c

# all output
snobal -z 2061 -m 0.01 -d 0.25 -t 60,60,15,1 -T 60,10,1 -s gold.snow.properties.input -h gold.inheight.input -p gold.snobal.ppt.input -i gold.snobal.data.input.short -O data -o gold.snobal.out.all -c -O all

# run the full period of record
snobal -z 2061 -m 0.01 -d 0.25 -t 60,60,15,1 -T 60,10,1 -s gold.snow.properties.input -h gold.inheight.input -p gold.snobal.ppt.input -i gold.snobal.data.input -O data -o gold.snobal.out.25year -c