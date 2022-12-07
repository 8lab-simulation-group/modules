ここでは、MBDynのAerodyn ver13に対応するモジュールの実装方法と、これを用いた解析方法について説明する。

1．モジュールの実装

1.1　モジュールを保存するフォルダを作成

Linux上のMBDynフォルダに移動し[~/mbdyn-1.7.3/modules/]に移動しここで[mkdir module-aerodyn13]を実行
実行後に [~/mbdyn-1.7.3/modules/]内に[module-aerodyn13]が作成されていることを確認する

作成したディレクトリの中に、さらに[mkdir source]で[source]というフォルダを作成
このフォルダはAerodyn本体のコンパイル等に使う


1.2　Aerodynをコンパイルする
共有フォルダ[\\as6104t-bed5\Public\FAST\FASTv7]内の

AD_v13.00.02a-bjj.exe
InflowWind_v1.02.00c-bjj.exe
NWTC_Lib_v1.07.00b-mlb.exe

の3つを自分のパソコンの適当な場所で展開する

展開した各フォルダの中にはそれぞれ[Sorce]というフォルダがあり、中には.f90ファイルがたくさん入っている
これらの.f90ファイルをすべて（3つ展開したものすべてから）コピーし、1.1で作った[Sorce]フォルダにペーストする

以下、1.1の[Sorce]ないで作業する
コピーしてきたプログラムをコンパイルするため、コマンド画面（黒い背景に白地）で以下のコマンドを実行

gfortran -ffixed-line-length-none -ffree-line-length-none -fPIC -O -c *.f90

これはフォルダ内のすべての.f90ファイルをコンパイルしろと命令している
しかし最初にこれを実行するとおそらくエラーが出る
これは上記コマンドを実行した場合アルファベットの早い順にコンパイルが行われていくわけだが、それだと先にコンパイルされ始めたプログラム内でまだコンパイルしてないプログラムを参照することができないからである
例えばa.f90とb.f90があって、a.f90の中でb.90のなかのサブルーチンを参照する場合、先にb.f90がコンパイルされていないとa.f90はコンパイルできないという具合である

従ってエラーが出たら、どのf90ファイルが「まだコンパイルされてないよ」とエラーを出されているのか確認し、適宜そのファイルからコンパイルしてやる必要がある
単一のプログラムをコンパイルする時は上記コマンドの「＊」の部分をプログラムファイル名に変えればよい


どうにかして地道にコンパイルを続け、すべてのf90ファイをコンパイルし終わるとたくさんの.oファイルと.modファイルが作成されているはずである
たくさんの.oファイルはリンカーとしてまとめる必要があるので、下記コマンドを実行

ar ru libAeroDyn.a *.o

実行後に[libAeroDyn.a]が作成されていることを確認する

確認出来たら、作成されたlibAeroDyn.aとさっきの.modファイルすべてをコピーし、ひとつ前のディレクトリに戻ってペーストする
（[~/mbdyn-1.7.3/modules/module-aerodyn13/]の中に大量の.modファイルとlibAeroDyn.aがあればOK)


1.2　Aerodyn13モジュールをコンパイル（ビルド）する
ここからはMBDynのrun-time moduleのビルド方法に則ってモジュールのビルドを行う

自作モジュールのコンパイルにおいて注意する点といえば、今回のようにメインソースとなるmodule-aerodyn.ccの他に参照するプログラム（mbdyn2ad.f90とか）がある場合
この場合はMakefile.incを作成し、ビルドの時に参照してもらう必要がある。
この中には

MODULE_DEPENDENCIES=mbdyn_ad.lo
MODULE_LINK=-L. -lAeroDyn

というコマンドが記載されており、mbdyn_ad.loのところを参照するプログラムの名前にしておいたり、AeroDynのところをメインソースファイルで定義しているクラス名にしたりしておく必要がある
またメインプログラムが参照するサブプログラムはメインより先にコンパイルされる必要があるため、ファイル名の文頭アルファベットをメインより早くしておく


コンパイルはいつも通り
./configure --enable-runtime-loading
./configure --with-module="aerodyn13"
LDFLAGS=-rdynamic
LIBS=/usr/lib/libltdl.a

make
sudo make install
で完了する


風車の解析においては線形ソルバーとしてUMFPACKが指定されることが多いので、こちらも実行可能にしておく

CPPFLAGS=-I/usr/include/suitesparse ./configure

make
sudo make install



2．解析方法

2.1　モデルを作る際の注意点
Aerodyn version13を使用するにあたって、モデル作成で注意すべきことは大きく3点

1点目はグローバル座標系の向き
風車の解析にあたって、必ず、風下方向にx軸、鉛直上向き方向にz軸をとる
つまりアップウィンド風車ならハブ位置はx軸負の値となるはずである


2点目は空力を計算する翼要素の座標の向き
翼要素が付随するnodeの要素座標系が翼端方向にz軸正、コードラインに沿って後縁方向にy軸正をとるようにする

特に、翼をbeam要素で定義する場合、beam要素は長手方向にxをとりがちなため、beam要素を定義するnodeをそのまま空力解析に使用すると大変なことになる
空力を与えるnodeはbeam定義用とは別に用意し、total jointでbeam定義のnodeと拘束するなどのモデル化が良いと思われる


3点目は空力計算に必要なnodeの作成
Aerodyn ver13はその仕様上、解析時にblade, hub, rotorfurl, nacelle, towerで示される5種類のnodeの位置、回転姿勢、並進速度、回転速度を必要とする
これらのnodeは設置されるべき位置と姿勢、運動の条件が決まっている
詳しくはAD_v13.00.02a-bjjフォルダの中にある「UserGuideAddendum_AeroDyn13Interface.pdf」を参照されたい

従って、「UserGuideAddendum_AeroDyn13Interface.pdf」に記載されている条件を満たしたnodeをモデル上に配置しておく必要がある
このnodeはあくまで解析に必要なデータを取得するためだけのもので力がかかったりすることはないのでdummy nodeなどを利用して設置するとよい


2.2　モジュールのインプットデータの作成
Aerodyn13モジュールを用いた解析に必要なインプットデータは以下の通りとなる

	Nacelle node label
	Hub node label
	RotorFurl node label
	Tower node label
	Blade length		(この定義もUserGuideAddendum_AeroDyn13Interface.pdf」に記載)
	Blade枚数
	Blade 1枚当たりに存在する空力を計算する翼要素の数
・	Blade node label	(,区切りで翼枚数分)
	Pitch角を考慮していない状態のBlade rootを参照するreference frameの定義
	Bladeを構成するnodeから見た、空力計算をする翼要素の中心位置と回転姿勢
	（force scale factor）
	AeroDynに関するアウトプットファイル名[.ad]
	AeroDynに関するインプットファイル名[.ipt]
	AeroDynに関するエレメントファイル名[.elm]

また.iptファイルの中のRELM（翼要素位置）はローターの回転軸からの距離での記載になっていることを確認する


以上をもとに作成したCART風車をモデルにしたサンプルがmodule-aerodyn13の中のdemoフォルダの中に格納されている

ただしこのサンプルはもともとあったAerodyn version12.58対応のAerodynモジュールのサンプルを参考にしている

何が言いたいかというと、R3年度荒木の修士論文で既存のAerodynamic Elementの空力計算が間違っていることが分かり、同様の結果を出すAerodynモジュール（ver.12.58）も適切な結果が得られないことが分かっている
しかしこのサンプルがこの間違った計算できれいな結果が出るように作られていたモデルであったようである
つまり、このdemoサンプルをそのまま回しても、途中でエラーが出て止まる
従って結果は参考にしないでいただきたい

あくまでdemoフォルダの中身は、「Aerodyn13モジュールを使う場合のインプットファイルの例」として参照すること




おまけ
ついでに荒木が使っていたMBDynの出力ファイルをエクセルに変換するプログラム「input.py」を残置しておく
mov, ine, jnt, frcなんでも対応可能
実行するとinput file name, data number, output file, time stepを聞かれるのでそれぞれ
入力ファイル名（.movとか拡張子まで込み）、欲しいデータ番号（node labelとか）、出力ファイル名（拡張子無し）、解析時間刻み幅を入れる
パッと散布図とかでグラフ作るときに便利なのでもしよければ




