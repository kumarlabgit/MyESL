cd src
g++-8 -std=c++17 preprocess_main.cpp preprocess.cpp -o preprocess -static -lstdc++fs -lpthread
g++-8 -std=c++17 sg_lasso_main.cpp sg_lasso.cpp /usr/lib/x86_64-linux-gnu/libopenblas.a /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libpthread.a -o sg_lasso -Iinclude -static
g++-8 -std=c++17 sg_lasso_leastr_main.cpp sg_lasso_leastr.cpp /usr/lib/x86_64-linux-gnu/libopenblas.a /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libpthread.a -o sg_lasso_leastr -Iinclude -static
g++-8 -std=c++17 overlapping_sg_lasso_leastr_main.cpp overlapping_sg_lasso_leastr.cpp /usr/lib/x86_64-linux-gnu/libopenblas.a /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libpthread.a -o overlapping_sg_lasso_leastr -Iinclude -static
g++-8 -std=c++17 overlapping_sg_lasso_logisticr_main.cpp overlapping_sg_lasso_logisticr.cpp /usr/lib/x86_64-linux-gnu/libopenblas.a /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libpthread.a -o overlapping_sg_lasso_logisticr -Iinclude -static
g++-8 -std=c++17 slep_main.cpp gl_logisticr.cpp slep_utils.cpp /usr/lib/x86_64-linux-gnu/libopenblas.a /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libpthread.a -o gl_logisticr -Iinclude -static
mv preprocess sg_lasso sg_lasso_leastr overlapping_sg_lasso_leastr overlapping_sg_lasso_logisticr gl_logisticr ../bin
cd ..
