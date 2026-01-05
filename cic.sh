# source module_load_t7920.sh

# cd utilities/
# # make clean
# make &&\
# cd ../ &&\
# # ./utilities/ic.x #&&\
# # ./utilities/cicpower.x


# # make clean
# make
# # ./utilities/ic.x && \
# ./main.x #&&\
# # ./utilities/cicpower.x


cd utilities/
# make clean
make &&\
cd ../ &&\
# make clean
make

rm run.log
# qsub pbsic.qsub
qsub pbscic.qsub
# qsub pbsmain.qsub
sleep 1
qstat -tan -u ChenBH
sleep 3
echo tail -f run.log
ls -lh run.log
tail -f run.log