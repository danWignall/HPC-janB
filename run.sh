rm *.dat
g++ janB.cpp
./a.out $1
python draw.py $1
