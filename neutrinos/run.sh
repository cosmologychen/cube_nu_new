export FOR_COARRAY_NUM_IMAGES=1
# mkdir tf_test
cd neutrinos/
# mkdir Pk_nu
make &&\
# chmod 777 ./neutrinos/test.x &&\
cd ../ &&\
# python ./neutrinos/Pk.py &&\
pwd &&\
./neutrinos/test.x > t2.txt
# ./test.x > t2.txt