# step 2: separeate vertices by annotation region
rt=`pwd`			# project root directory
cd raw/VT			# vertex directory

mkdir -p 2			# create directory for step 2

# directories for 36 anatomy regions
for i in {0..36}; do mkdir -p 2/$i ; done

for i in {0..36}; do
    for f in 1/*.txt; do
	s=${f#*/}
	awk <1/$s "\$3==$i" >2/$i/$s
    done
done

cd 2
for i in {0..36}; do
    cd $i
    tar -zcf ../$i.tar.gz *
    cd ..
done
cd ..
