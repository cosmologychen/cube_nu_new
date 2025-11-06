# source module_load_t7920.sh

cd utilities/
make clean
make &&\
cd ../ &&\
./utilities/ic.x #&&\
# ./utilities/cicpower.x


#make clean
make
# ./utilities/ic.x && \
./main.x #&&\
./utilities/cicpower.x
