g++ a3.cpp

aval=1.
./a.out $aval > res.data
python hist.py res.data res.hist $aval

