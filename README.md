# AmedasPlot
* このプログラムは、アメダスの観測データ（気圧・気温・風）とレーダーGPVを重ね書きするプログラムです。
* 今後の目標としては、過去の事例を解析できるよう、下記の気象庁HPのJSON（10日分くらいしか残らない）以外からデータを取得したいと考えていますが、誰か代わりに改修してください。

## データ取得元
* アメダスデータ:気象庁HPのJSON（ https://www.jma.go.jp/bosai/amedas/data/map/{YYYY}{MM}{DD}{HH}{mm}00.json ）
* レーダーデータ:京都大学生存圏研究所（ http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original ）

## special thanks
* 黒良先生が開発・公開されているアメダスプロットコードをもとに改変しました。https://note.com/rkurora/n/n47e2099f74b0
* レーダーデータのプロットはこちらを参考にしました。https://hineken.com/blog/mete_arti/jmaradar_python/
* レーダーデータのカラーバーはこちらを参照しました。https://qiita.com/earth06/items/eb579d122bb67d964c40
