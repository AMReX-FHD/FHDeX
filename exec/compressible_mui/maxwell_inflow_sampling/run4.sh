g++ a4.cpp

aval=1
./a.out $aval > res.data_a$aval
python hist4.py res.data_a$aval res.hist_a$aval $aval

