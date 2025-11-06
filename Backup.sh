opath=$(pwd)
echo $opath 


cd
# cd cbh &&\
mkdir backup
rm -vr backup/cube_nu_new
mkdir backup/cube_nu_new
cp -v $opath/* ./backup/cube_nu_new
cp -vr $opath/Green ./backup/cube_nu_new
cp -vr $opath/tf_wmap9 ./backup/cube_nu_new
cp -vr $opath/velocity_conversion ./backup/cube_nu_new
cp -vr $opath/utilities ./backup/cube_nu_new

mkdir -p backup/cube_nu_new/neutrinos
cp -v $opath/neutrinos/*  ./backup/cube_nu_new/neutrinos

mkdir -p backup/cube_nu_new/reps-master
cp -v $opath/reps-master/*  ./backup/cube_nu_new/reps-master
cp -vr $opath/reps-master/EXAMPLES ./backup/cube_nu_new/reps-master
cp -vr $opath/reps-master/INCLUDE ./backup/cube_nu_new/reps-master
cp -vr $opath/reps-master/MODULES ./backup/cube_nu_new/reps-master
cp -vr $opath/reps-master/SOURCE ./backup/cube_nu_new/reps-master
cp -vr $opath/reps-master/tabulated_functions ./backup/cube_nu_new/reps-master

cd backup/cube_nu_new &&\
make clean &&\
cd neutrinos &&\
make clean &&\
cd ../utilities &&\
make clean &&\


echo
echo
echo
echo
echo ++++++++++
echo ++++++++++
cd 
pwd
# cd cbh &&\
cd backup
times=$(date +%-m-%-d-%-H-%-M)
tar -cvf cube_nu_$1_$times.tar cube_nu_new
ls -lh *tar
path=$(pwd)



echo
echo
echo
echo
echo ++++++++++
echo ++++++++++
echo time : $times
echo ++++++++++
echo ++++++++++
echo $path/cube_nu_$1_$times.tar
echo ++++++++++
echo ++++++++++

if [ "$#" -lt 3 ]; then
    if [ "$2" = "cp" ]; then
        cp $path/cube_nu_$1_$times.tar $opath
        echo  cp $path/cube_nu_$1_$times.tar $opath
    else
        echo "No command"
    fi
fi
