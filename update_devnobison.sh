rm -rf build_dev
git checkout dev
git pull
mkdir build_dev
cd build_dev
cmake ..
make
git checkout dev_nobison
cp src/*.c ../src/
