g++ a0.cpp
./a.out > res.data0
python hist.py res.data0 res.hist0 0.

g++ a1.cpp
./a.out > res.data1
python hist.py res.data1 res.hist1 -0.1

g++ a2.cpp
./a.out > res.data2
python hist.py res.data2 res.hist2 0.1
