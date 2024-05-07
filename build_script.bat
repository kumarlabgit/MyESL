rm -rf bin
mkdir bin
cd src
..\..\mingw64\bin\g++ -std=c++17 preprocess_main.cpp preprocess.cpp -o preprocess -static -lstdc++fs
..\..\mingw64\bin\g++ -std=c++17 sg_lasso_main.cpp sg_lasso.cpp -static libopenblas.dll -static liblapack.dll -o sg_lasso -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 sg_lasso_leastr_main.cpp sg_lasso_leastr.cpp -static libopenblas.dll -static liblapack.dll -o sg_lasso_leastr -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 overlapping_sg_lasso_leastr_main.cpp overlapping_sg_lasso_leastr.cpp -static libopenblas.dll -static liblapack.dll -o overlapping_sg_lasso_leastr -Iinclude -static -lpthread
..\..\mingw64\bin\g++ -std=c++17 overlapping_sg_lasso_logisticr_main.cpp overlapping_sg_lasso_logisticr.cpp -static libopenblas.dll -static liblapack.dll -o overlapping_sg_lasso_logisticr -Iinclude -static -lpthread
mv preprocess.exe sg_lasso.exe sg_lasso_leastr.exe overlapping_sg_lasso_leastr.exe overlapping_sg_lasso_logisticr.exe ../bin
cp *.dll ..\bin
cp *.lib ..\bin
cd ..
pyinstaller ESL_pipeline.py
pyinstaller ESL_model_apply.py
mv dist\ESL_pipeline\_internal .
mv dist\ESL_pipeline\ESL_pipeline.exe .
mv dist\ESL_model_apply\ESL_model_apply.exe .
rm -rf build
rm -rf dist