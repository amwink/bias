rm -rf z zlib*
wget http://sourceforge.net/projects/libpng/files/zlib/1.2.8/zlib-1.2.8.tar.gz
tar zxvf zlib-1.2.8.tar.gz
ln -sf zlib-1.2.8 z
cd z
wget --no-check-certificate https://raw.githubusercontent.com/Alexpux/MSYS2-packages/master/zlib/1.2.7-minizip-cygwin.patch
wget --no-check-certificate https://raw.githubusercontent.com/Alexpux/MSYS2-packages/master/zlib/1.2.7-zlib-symbols.patch
wget --no-check-certificate https://raw.githubusercontent.com/Alexpux/MSYS2-packages/master/zlib/zlib-1.2.8-msys2.patch
patch -p2 -i 1.2.7-minizip-cygwin.patch
patch -p2 -i 1.2.7-zlib-symbols.patch
patch -p1 -i zlib-1.2.8-msys2.patch 
perl -i -wpe s/zlib1/libz/ z/win32/Makefile*
perl -i -wpe s/zlib1/libz/ z/Makefile*
perl -i -wpe s/zlib1/libz/ z/configure
perl -i -wpe 's/-O2/-Ofast -march=native -fPIC/' `find . -name Makefile\*`
bash configure
make
cd ..
mkdir -p ../../../bin
mkdir -p ../../../lib
cp z/libz.so* ../../../bin/
cp z/libz.a* ../../../bin/
cp z/libz.dll* ../../../lib/
