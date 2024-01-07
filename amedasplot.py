#!/usr/bin/env python
# coding: utf-8
# 
# AMeDAS 観測データプロット図 （地図領域は、緯度・経度のそれぞれ2度の範囲)
# 矢羽,気温(左上),湿球温度（右下）,露点温度(左下),海面更正気圧(右上)
# 黒良先生のプログラムにレーダーGPVを重ね合わせ
#   
# アメダス地点テーブルJSON  https://www.jma.go.jp/bosai/amedas/const/amedastable.json
# アメダス観測データJSON  YYYY/MM/DD HH:mm(JST)  https:https://www.jma.go.jp/bosai/amedas/data/map/{YYYY}{MM}{DD}{HH}{mm}00.json
# 生存圏研究所ダウンロード元サイト  http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import numpy as np
np.set_printoptions(threshold=np.inf)
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
import metpy.calc as mpcalc
import json
import urllib
import math
import matplotlib
import subprocess
import os, sys, glob
import pandas as pd
import pandas.tseries.offsets as offsets
from metpy.units import units
from datetime import datetime, timedelta
import pytz
from cartopy.mpl.ticker import LatitudeFormatter,LongitudeFormatter
import pybufrkit
import struct
from itertools import repeat
from scipy.interpolate import griddata,RectBivariateSpline
from scipy.ndimage import gaussian_filter

def parse_datetime(arg):
    try:
        # 引数が12桁の数字の場合、YYYYMMDDHHMM形式の文字列を解析
        if len(arg) == 12 and arg.isdigit():
            year = arg[:4]
            month = arg[4:6]
            day = arg[6:8]
            hour = arg[8:10]
            min = arg[10:12]
            dt = datetime(int(year), int(month), int(day), int(hour), int(min))
        elif len(arg) == 10 and arg.isdigit():
            year = arg[:4]
            month = arg[4:6]
            day = arg[6:8]
            hour = arg[8:10]
            min = 0
            dt = datetime(int(year), int(month), int(day), int(hour), int(min))        
        else: 
            raise ValueError()
        return dt        
    except ValueError:
        return None

def set_table(section5):
    max_level = struct.unpack_from('>H', section5, 15)[0]
    table = (
        -10, # define representative of level 0 (Missing Value)
        *struct.unpack_from('>'+str(max_level)+'H', section5, 18)
    )
    return np.array(table, dtype=np.int16)

def decode_runlength(code, hi_level):
    for raw in code:
        if raw <= hi_level:
            level = raw
            pwr = 0
            yield level
        else:
            length = (0xFF - hi_level)**pwr * (raw - (hi_level + 1))
            pwr += 1
            yield from repeat(level, length)

def load_jmara_grib2(file):
    with open(file, 'rb') as f:
        binary = f.read()

    len_ = {'sec0':16, 'sec1':21, 'sec3':72, 'sec4':82, 'sec6':6}

    end4 = len_['sec0'] + len_['sec1'] + len_['sec3'] + len_['sec4'] - 1
    len_['sec5'] = struct.unpack_from('>I', binary, end4+1)[0]
    section5 = binary[end4:(end4+len_['sec5']+1)]

    end6 = end4 + len_['sec5'] + len_['sec6']
    len_['sec7'] = struct.unpack_from('>I', binary, end6+1)[0]
    section7 = binary[end6:(end6+len_['sec7']+1)]

    highest_level = struct.unpack_from('>H', section5, 13)[0]
    level_table = set_table(section5)
    decoded = np.fromiter(
        decode_runlength(section7[6:], highest_level), dtype=np.int16
    ).reshape((3360, 2560))
    
    # 転置するとうまくいく
    transposed_flipped_data = np.flip(np.transpose(decoded), axis=1)

    # convert level to representative
    return level_table[transposed_flipped_data]

# 地点テーブル
# 読み込み設定
n_station_json='./amedastable.json'
# 描画指定：順に気圧(右上),気温(左上),湿球温度(右下),露点温度(左下))
npre_dispflag=False
temp_dispflag=False
wbt_dispflag=False
dp_dispflag=False

# マーカーサイズ
markersize_0 = 2
# 文字サイズ
char_size=16
# 矢羽の長さ
barb_length=8

# 地図の中心位置を指定
(lat_center, lon_center) = (36, 138)   # 関東付近
# 地図の描画範囲指定
#i_area = [lon_center - 2.0, lon_center + 2.0, lat_center - 2.0, lat_center + 2.0]
i_area = [lon_center - 3.5, lon_center + 3.5, lat_center - 2.5, lat_center + 2.5]
# 緯線・経線の指定
dlon,dlat=1,1   # 1度ごとに

# 緯度経度で指定したポイントの図上の座標などを取得する関数 transform_lonlat_to_figure() 
# 図法の座標 => pixel座標 => 図の座標　と3回の変換を行う
#  　pixel座標: plt.figureで指定した大きさxDPIに合わせ、左下を原点とするpixelで測った座標   
#  　図の座標: axesで指定した範囲を(0,1)x(0,1)とする座標
# 3つの座標を出力する
#    図の座標, Pixel座標, 図法の座標
def transform_lonlat_to_figure(lonlat, ax, proj):
    # lonlat:経度と緯度  (lon, lat) 
    # ax: Axes図の座標系    ex. fig.add_subplot()の戻り値
    # proj: axで指定した図法 
    #
    # 例 緯度経度をpointで与え、ステレオ図法る場合
    #    point = (140.0,35.0)
    #    proj= ccrs.Stereographic(central_latitude=60, central_longitude=140) 
    #    fig = plt.figure(figsize=(20,16))
    #    ax = fig.add_subplot(1, 1, 1, projection=proj)
    #    ax.set_extent([108, 156, 17, 55], ccrs.PlateCarree())
    #
    # 図法の変換
    # 参照  https://scitools.org.uk/cartopy/docs/v0.14/crs/index.html                    
    point_proj = proj.transform_point(*lonlat, ccrs.PlateCarree())
    #
    # pixel座標へ変換
    # 参照　https://matplotlib.org/stable/tutorials/advanced/transforms_tutorial.html
    point_pix = ax.transData.transform(point_proj)
    #
    # 図の座標へ変換                                                           
    point_fig = ax.transAxes.inverted().transform(point_pix)
    return point_fig, point_pix, point_proj

# 1地点のアメダスjsonデータから、elem要素で指定した値を返す(ただしFlagが 0以外は Noneとする)
# 要素:elem = 'temp','humidity','snow1h','snow6h','snow12h','snow24h','sun10m','sun1h',
#            'precipitation10m','precipitation1h','precipitation3h','precipitation24h'
#            'wind','windDirection'
def get_obs_value(amd_obs,elem):
    try:
        et = amd_obs[elem]
        if int(et[1]) != 0:
            return None
        return float(et[0])
    except Exception:
        return None

# 緯度・経度のtuple(度分形式)をtuple(度単位)に変換
def trans_tupleLL(lat_tuple, lon_tuple):
    lat = lat_tuple[0] + lat_tuple[1]/60.0
    lon = lon_tuple[0] + lon_tuple[1]/60.0
    return (lat, lon)

# 16方位風向と風速から、u,v作成   
def get_uv(ws, wdd16):
    if ws == None or wdd16 == None:
        u = None
        v = None
    else:
        wd = wdd16 / 8.0 * math.pi  # 1/8 = (360/16) / 360 * 2
        u = -1.0 * ws * math.sin(wd)
        v = -1.0 * ws * math.cos(wd)
    return u,v

# 観測データ日時
if len(sys.argv) == 2:
    arg = sys.argv[1]
    dt = parse_datetime(arg)
    if dt:
        year=dt.year
        month=dt.month
        day=dt.day
        hour=dt.hour
        min=0
        print("読み込み観測時刻 {:4d}/{:02d}/{:02d} {:02d}:{:02d}".format(year,month,day,hour,min))
    else:
        print('Usage: python script.py [YYYYMMDDHH(MM)]')
elif len(sys.argv) == 1:
    jst = pytz.timezone('Asia/Tokyo')
    dt = datetime.now(jst) - timedelta(minutes=30)
    year=dt.year
    month=dt.month
    day=dt.day
    hour=dt.hour
    min=0
    dt = datetime(int(year), int(month), int(day), int(hour), int(min))                
    print("読み込み観測時刻 {:4d}/{:02d}/{:02d} {:02d}:{:02d}".format(year,month,day,hour,min))        
else:
    print('Usage: python script.py [YYYYMMDDHH(MM)]')
    
# 観測データJSONの url作成
url_data_json= 'https://www.jma.go.jp/bosai/amedas/data/map/{:4d}{:02d}{:02d}{:02d}{:02d}00.json'
url_data_json=url_data_json.format(year,month,day,hour,min)
# 気象庁HPからアメダスデータを読み込む
response = urllib.request.urlopen(url_data_json)
content = response.read()
response.close()
data_json=content.decode()
dat_json = json.loads(data_json)

# アメダス地点Tableのurl
url_station_json="https://www.jma.go.jp/bosai/amedas/const/amedastable.json"
# アメダス地点Tableを読み込む
response = urllib.request.urlopen(url_station_json)
content = response.read()
response.close()
station_json=content.decode()
amd_json = json.loads(station_json)

# 気象庁全国合成レーダーGPVのダウンロード
def download_time(time):    
    # 生存圏研究所ダウンロード元サイト
    http  = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original"  
    # 保存先ディレクトリの指定
    Opath = "."
    day_dir = time.strftime("%Y/%m/%d")
    basename = "Z__C_RJTD_{}00_RDR_JMAGPV__grib2.tar".format(time.strftime("%Y%m%d%H%M"))
    fname    = "{}/{}/{}".format(Opath, day_dir, basename)
    # すでにファイルが存在しなければ、ダウンロードを行う
    if os.path.exists(fname):
        print("{} Already exists.".format(basename))
    else:
        os.makedirs("{}/{}".format(Opath, day_dir), exist_ok=True)
        url      = "{}/{}/{}".format(http,  day_dir, basename)
        # wgetコマンドでデータのダウンロード
        subprocess.run("wget {} -P {}/{}/".format(url, Opath, day_dir), shell=True)
        # tarコマンドでダウンロードした圧縮ファイルの解凍
        subprocess.run("tar -xvf {} -C {}/{}/".format(fname, Opath, day_dir), shell=True)   
    GgisFile = "{}/{}/Z__C_RJTD_{}00_RDR_JMAGPV_Ggis1km_Prr10lv_ANAL_grib2.bin".format(Opath, day_dir, time.strftime("%Y%m%d%H%M"))
    # wgrib2コマンドでgrib2データをバイナリファイルに変換
    #outfile  = "{}/{}/Z__C_RJTD_{}00_RDR_JMAGPV_Ggis1km_Prr10lv_ANAL.dat".format(Opath, day_dir, time.strftime("%Y%m%d%H%M"))
    #if not os.path.exists(outfile):
    #    subprocess.run("wgrib2 {} -d 1 -no_header -bin {}".format(GgisFile, outfile), shell=True)
    return GgisFile
    
# 描画する時間の指定(年,月,日,時,分)：データは10分ごと（前10分の雨量が記録されている）    
# アメダスデータと同じ時刻のレーダーGPVをダウンロード
time = pd.Timestamp(year,month,day,hour,min)
utc = time - offsets.Hour(9)
filepath = download_time(utc)
#rain = mkdata_time(utc)

# データを読む
rain = load_jmara_grib2(filepath) / 100 
        
# 図法指定                                                                             
proj = ccrs.PlateCarree()
latlon_proj = ccrs.PlateCarree()

# 図のSIZE指定inch                                                                        
fig = plt.figure(figsize=(20,15))
# 余白設定                                                                                
plt.subplots_adjust(left=0.04, right=1, bottom=0.06, top=0.98)                  
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent(i_area, latlon_proj)

# カラーバーの設定
#気象庁RGBカラー
jmacolors=np.array(
   [[0.95,0.95,0.95,1],#white
    [0.63,0.82,0.99,1],
    [0.13,0.55,0.99,1],
    [0.00,0.25,0.99,1],
    [0.98,0.96,0.00,1],
    [0.99,0.60,0.00,1],
    [0.99,0.16,0.00,1],
    [0.71,0.00,0.41,1]]
)
#等高線値
clevs = np.array([0,1,5,10,20,30,50,80]) 
   
## 入力ファイルの次元の指定（水平解像度約1km）
## 経度方向： 118°〜150° (1/80°間隔)
## 緯度方向： 20°〜48° (1/120°間隔)
mlon, mlat = 2560, 3360
slon, elon, rlon = 118, 150, 1/80
slat, elat, rlat = 20, 48, 1/120

lon = np.arange(slon, elon, rlon)
lat = np.arange(slat, elat, rlat)
LON, LAT = np.meshgrid(lon, lat)
LON, LAT = LON.T, LAT.T
cs = ax.contourf(LON, LAT, rain, colors=jmacolors, levels=clevs, extend="max")
cb = plt.colorbar(cs, orientation="vertical", ticks=clevs, shrink=0.6)    
cb.ax.tick_params(labelsize=15)

# グリッド線を引く                                                               
xticks=np.arange(-180,180,dlon)
yticks=np.arange(-90,90.1,dlat)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)

# 配列の宣言
lat_list_t = []
lon_list_t = []
temp_list = []
lat_list_p = []
lon_list_p = []
npre_list = []

# 地点プロット                                                                                                 
for stno,val in dat_json.items():
    # データ取得
    (wlat,wlon) = trans_tupleLL(amd_json[stno]['lat'], amd_json[stno]['lon'])
    walt = amd_json[stno]['alt']
    # 風
    wind_ok = True
    ws = get_obs_value(val,'wind')
    wd = get_obs_value(val,'windDirection')
    if ws is not None and wd is not None:
        (au,av) = get_uv(ws,wd)
    else:
        wind_ok = False
    # 気温
    temp = get_obs_value(val,'temp')
    if temp is None:
        temp = np.nan
    elif walt < 1500:
        # 配列に格納
        lat_list_t.append(wlat)
        lon_list_t.append(wlon)
        temp_list.append(temp)
    # 湿度
    hu = get_obs_value(val,'humidity')
    if hu is None:
        hu = -1.0
    # 露点温度
    dp_temp = -200.0
    if hu >= 0.0 and temp > -200.0:
        dp_temp = mpcalc.dewpoint_from_relative_humidity(temp * units.degC,hu/100.0).m
    # 更正気圧
    npre = get_obs_value(val,'normalPressure')
    if npre is None:
        npre = np.nan
    else:
        # 配列に格納
        lat_list_p.append(wlat)
        lon_list_p.append(wlon)
        npre_list.append(npre)
    
    # 気圧
    pre = get_obs_value(val,'pressure')
    if pre is None:
        pre = -1.0

    # 湿球温度
    wb_temp = -200.0
    if dp_temp > -200.0 and temp > -200.0 and pre > 0.0:
        wb_temp = mpcalc.wet_bulb_temperature(pre * units.hPa, temp * units.degC, dp_temp * units.degC).m
    
    ## プロット
    fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj) 
    if ( fig_z[0] > 0.01 and fig_z[0] < 0.99  and fig_z[1] > 0.01 and fig_z[1] < 0.99):
        ax.plot(wlon, wlat, marker='s' , markersize=markersize_0, color="brown", transform=latlon_proj)
        if wind_ok and au*au+av*av>4.0: # 矢羽プロット
            ax.barbs(wlon, wlat, 
                     (au * units('m/s')).to('kt').m, (av * units('m/s')).to('kt').m,
                     length=barb_length, transform=latlon_proj)
        if npre_dispflag and pre >= 0.0: # 気圧プロット
            ax.text(fig_z[0]+0.029, fig_z[1]+0.015,'{:6.1f}'.format(npre),size=char_size, color="black", transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")
        if temp_dispflag and temp > -200.0: # 気温プロット
            color_temp = "red"
            ax.text(fig_z[0]-0.025, fig_z[1]+0.015,'{:5.1f}'.format(temp),size=char_size, color=color_temp, transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")
        if wbt_dispflag and wb_temp > -200.0: # 湿球温度プロット
            if dp_temp < 0:
                color_temp = "purple"
            else:
                color_temp = "black"
            ax.text(fig_z[0]+0.025, fig_z[1]-0.003,'{:5.1f}'.format(wb_temp),size=char_size, color=color_temp, transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")
        if dp_dispflag and dp_temp > -200.0: # 露点温度プロット
            if dp_temp < 0:
                color_temp = "purple"
            else:
                color_temp = "black"
            ax.text(fig_z[0]-0.025, fig_z[1]-0.003,'{:5.1f}'.format(dp_temp),size=char_size, color=color_temp, transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")  

# 0.05度単位のグリッドを作成
grid_lon_t, grid_lat_t = np.meshgrid(np.arange(i_area[0], i_area[1], 0.05),
                                     np.arange(i_area[2], i_area[3], 0.05))

# 0.25度単位のグリッドを作成
grid_lon_p, grid_lat_p = np.meshgrid(np.arange(i_area[0], i_area[1], 0.25),
                                     np.arange(i_area[2], i_area[3], 0.25))

# 線形補間
grid_temp = griddata((lon_list_t, lat_list_t), temp_list, (grid_lon_t, grid_lat_t), method='linear')
grid_npre = griddata((lon_list_p, lat_list_p), npre_list, (grid_lon_p, grid_lat_p), method='linear')

# ガウシアンフィルタを適用
sigma = 1.0  # ガウス分布の標準偏差
#grid_temp = gaussian_filter(grid_temp, sigma=0.5)
#grid_npre = gaussian_filter(grid_npre, sigma=0.5)

# 等温線をプロット
levels = np.arange(-30, 60, 3)
cont = plt.contour(grid_lon_t, grid_lat_t, grid_temp, levels=levels, linewidths=2, linestyle='solid', colors='red')

# 等温線のラベルを付ける
plt.clabel(cont, fontsize=20)

# 等圧線をプロット
levels = np.arange(900, 1050, 1)
cont = plt.contour(grid_lon_p, grid_lat_p, grid_npre, levels=levels, linewidths=2, colors='black')

# 等圧線のラベルを付ける
plt.clabel(cont, fontsize=20)

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black') # 海岸線の解像度を上げる   
            
# 図の説明
plt.title('{}'.format("AMeDAS and RadarGPV"), loc='left',size=20)
plt.title('Valid Time: {}'.format(dt), loc='right',size=20);
#plt.savefig("{}.jpg".format(time.strftime("%Y%m%d%H%M")), format="jpg")
plt.savefig("latest.jpg", format="jpg")
