[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)
[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)
[![unstable](http://badges.github.io/stability-badges/dist/unstable.svg)](http://github.com/badges/stability-badges)


# bubble
3-Dimensional Voronoi Diagram Divided by Cones.

## Definition of Voronoi Diagram


## Starting Point
Although bubbling is the elegant algorithm, the data structure and algorithm are complex and heavy.
However, contacted bubbles hinted starting point of an algorithm described later.

![bubbling](doc/fig/Voronoi_growth_euclidean.gif)



## Tessellation Algorithm
Our approach accepts a partial reconstruction.
Maybe a novel approach. 
Please contact me if you know a similar one.

### Discretization and Data Structure
Voronoi cells are represented as a set of cones in our discretization.
それぞれの円錐の底面の中心はセルの境界に一致する。
この記述ではギャップや円錐同士のオーバーラップが生じる。
分割数を増やすほど実際のセルに一致する性質は立方体のvoxel表現した場合と同様である。
情報の圧縮性は円錐を用いた場合の方が高い。
![bubbling](doc/fig/discretization.jpeg)
円錐の方向を決めるために、格子を必要とする。
格子点は球面にほぼ一様に分布する
また、Pkから伸びるそれぞれの円錐の底面までのベクトルの長さを決める計算の高速化のために、
これらのすべて格子点は繋留され、閉じた１つのネットワークが作成される。

## Contribution
1. Fork it ( https://github.com/toyaku-phys/bubble/fork )
2. Create your feature branch (git checkout -b my-new-feature)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin my-new-feature)
5. Create a new Pull Request to the bubble/master branch

Or write issue

## Versioning
We use [SemVer](http://semver.org/) for versioning. 
For the versions available, see the tags on this repository.

## Authors
* [**Hibiki Itoga**](https://github.com/misteltein) -Key programmer-
* [**yde**](https://github.com/master-yde) -Discussion partner-

## License
MIT-license
