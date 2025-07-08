rm -rf _internal
rm -rf bin
rm MyESL_model_apply.exe
rm MyESL.exe
mkdir bin
cd src
..\..\mingw64\bin\g++ -std=c++17 preprocess_main.cpp preprocess.cpp -o preprocess -static -lstdc++fs
..\..\mingw64\bin\g++ -std=c++17 sg_lasso_main.cpp sg_lasso.cpp -static libopenblas.dll -static liblapack.dll -o sg_lasso -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 sg_lasso_leastr_main.cpp sg_lasso_leastr.cpp -static libopenblas.dll -static liblapack.dll -o sg_lasso_leastr -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 overlapping_sg_lasso_leastr_main.cpp overlapping_sg_lasso_leastr.cpp -static libopenblas.dll -static liblapack.dll -o overlapping_sg_lasso_leastr -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 overlapping_sg_lasso_logisticr_main.cpp overlapping_sg_lasso_logisticr.cpp -static libopenblas.dll -static liblapack.dll -o overlapping_sg_lasso_logisticr -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 slep_main.cpp gl_logisticr.cpp slep_utils.cpp -static libopenblas.dll -static liblapack.dll -o gl_logisticr -Iinclude -static -lpthread
mv preprocess.exe sg_lasso.exe sg_lasso_leastr.exe overlapping_sg_lasso_leastr.exe overlapping_sg_lasso_logisticr.exe gl_logisticr.exe ../bin
cp *.dll ..\bin
cp *.lib ..\bin
cd ..
pyinstaller MyESL.py
pyinstaller MyESL_model_apply.py
mv dist\MyESL\_internal .
mv dist\MyESL\MyESL.exe .
mv dist\MyESL_model_apply\MyESL_model_apply.exe .
rm -rf build
rm -rf dist
rm MyESL_model_apply.spec
rm MyESL.spec