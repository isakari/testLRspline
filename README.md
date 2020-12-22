## testLRspline
文献[1]の計算例を再現するための fortran コード. 

# 01example224
2.2.4 Example を再現する。

```usage
  cd src
  make
  cd ../bin
  ./a.out
  cd plot
  gnuplot plot.gp
```
で、lrmesh.png (refine前後のLRメッシュ) と surface.png (refine前後の適当な曲面) が見れる。

# 02diagonal
diagonal structured refinement を再現する。 main.f90 の
```
  do j=1,5 !level 5 まで
```
を変えるとlevel nまで diagonal refinement ができる。使い方は01example224と同じ。

# Reference
[1] Kjetil André, Johannessen, Trond Kvamsdal, and Tor Dokken, Isogeometric analysis using LR B-splines, Computer Methods in Applied Mechanics and Engineering, Volume 269, 1 February 2014, Pages 471-514
