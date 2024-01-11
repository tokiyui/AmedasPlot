#!/usr/bin/env python
# coding: utf-8

# AMeDAS 観測データプロット図 （地図領域は、緯度・経度のそれぞれ2度の範囲)
# 矢羽,気温(左上),湿球温度（右下）,露点温度(左下),海面更正気圧(右上)
# 黒良先生のプログラムにレーダーGPVを重ね合わせ
# アメダス地点テーブルJSON  https://www.jma.go.jp/bosai/amedas/const/amedastable.json
# アメダス観測データJSON  YYYY/MM/DD HH:mm(JST)  https:https://www.jma.go.jp/bosai/amedas/data/map/{YYYY}{MM}{DD}{HH}{mm}00.json
# 生存圏研究所ダウンロード元サイト  http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original

import json, math, matplotlib, os, pygrib, pytz, struct, subprocess, sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import metpy.calc as mpcalc
import numpy as np
import pandas as pd
import pandas.tseries.offsets as offsets
from datetime import datetime, timedelta
from itertools import repeat
from metpy.units import units
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter, maximum_filter, minimum_filter
from urllib.request import urlopen
#from scipy import interpolate
from scipy.interpolate import interp2d, RectBivariateSpline, RegularGridInterpolator

## 極大/極小ピーク検出関数                                                             
def detect_peaks(image, filter_size, dist_cut, flag=0):
    # filter_size: この値xこの値 の範囲内の最大値のピークを検出                        
    # dist_cut: この距離内のピークは1つにまとめる                                      
    # flag:  0:maximum検出  0以外:minimum検出                                          
    if flag==0:
      local_max = maximum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant') 
      detected_peaks = np.ma.array(image, mask=~(image == local_max))
    else:
      local_min = minimum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant')
      detected_peaks = np.ma.array(image, mask=~(image == local_min))
    peaks_index = np.where((detected_peaks.mask != True))
    # peak間の距離行例を求める                                                         
    (x,y) = peaks_index
    size=y.size
    dist=np.full((y.size, y.size), -1.0)
    for i in range(size):
      for j in range(size):
        if i == j:
          dist[i][j]=0.0
        elif i>j:
          d = math.sqrt(((y[i] - y[j])*(y[i] - y[j])) + ((x[i] - x[j])*(x[i] - x[j])))
          dist[i][j]= d
          dist[j][i]= d
    # 距離がdist_cut内のpeaksの距離の和と、そのピーク番号を取得する 
    Kinrin=[]
    dSum=[]
    for i in range(size):
      tmpA=[]
      distSum=0.0
      for j in range(size):
        if dist[i][j] < dist_cut and dist[i][j] > 0.0:
          tmpA.append(j)
          distSum=distSum+dist[i][j]
      dSum.append(distSum)
      Kinrin.append(tmpA)
    # Peakから外すPeak番号を求める.  peak間の距離和が最も小さいものを残す              
    cutPoint=[]
    for i in range(size):
      val = dSum[i]
      val_i=image[x[i]][y[i]]
      for k in Kinrin[i]:
        val_k=image[x[k]][y[k]]
        if flag==0 and val_i < val_k:
            cutPoint.append(i)
            break
        if flag!=0 and val_i > val_k:
            cutPoint.append(i)
            break
        if val > dSum[k]:
            cutPoint.append(i)
            break
        if val == dSum[k] and i > k:
            cutPoint.append(i)
            break
    # 戻り値用に外すpeak番号を配列から削除                                             
    newx=[]
    newy=[]
    for i in range(size):
      if (i in cutPoint):
        continue
      newx.append(x[i])
      newy.append(y[i])
    peaks_index=(np.array(newx),np.array(newy))
    return peaks_index

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
    decoded = np.fromiter(decode_runlength(section7[6:], highest_level), dtype=np.int16).reshape((3360, 2560))

    transposed_flipped_data = np.flip(np.transpose(decoded), axis=1)
    return level_table[transposed_flipped_data] / 100
    
def read_msm(time):
    # 初期値から5時間後に配信される
    modeltime = time - offsets.Hour(5)
    # MSMは03シリーズ
    base_time = modeltime.replace(hour=modeltime.hour - (modeltime.hour % 3), minute=0, second=0)  
    # 対象時刻と初期値の時間差
    ft = (modeltime - base_time).total_seconds() // 3600
    # 生存圏研究所ダウンロード元サイト
    http  = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original"  
    # データファイル名の指定
    day_dir = base_time.strftime("%Y/%m/%d")
    basename = "Z__C_RJTD_{}00_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin".format(base_time.strftime("%Y%m%d%H%M"))
    fname    = "./{}".format(basename)
    # すでにファイルが存在しなければ、ダウンロードを行う
    if os.path.exists(fname):
        print("{} Already exists.".format(basename))
    else:
        url      = "{}/{}/{}".format(http,  day_dir, basename)
        # wgetコマンドでデータのダウンロード
        subprocess.run("wget {} -P ./".format(url), shell=True)
        
    grbs = pygrib.open(basename)
    # 海面気圧のデータを取得する
    prmsl_fc0 = grbs.select(parameterName='Pressure reduced to MSL', forecastTime=ft)[0]
    prmsl, lats, lons = prmsl_fc0.data()
    prmsl_flipped = np.flip(prmsl / 100, axis=0)
    # 気温のデータを取得する
    tmp_fc0 = grbs.select(parameterName='Temperature', forecastTime=ft)[0]
    tmp, lats, lons = tmp_fc0.data()
    tmp_flipped = np.flip(tmp - 273.15, axis=0)
    
    # ガウシアンフィルタを適用
    prmsl_flipped = gaussian_filter(prmsl_flipped, sigma=3.0) # sigmaはガウス分布の標準偏差
    tmp_flipped = gaussian_filter(tmp_flipped, sigma=1.0) # sigmaはガウス分布の標準偏差
    return prmsl_flipped, tmp_flipped
    
# 緯度経度で指定したポイントの図上の座標などを取得する関数 
# 図法の座標 => pixel座標 => 図の座標　と3回の変換を行う
# pixel座標: plt.figureで指定した大きさxDPIに合わせ、左下を原点とするpixelで測った座標   
# 図の座標: axesで指定した範囲を(0,1)x(0,1)とする座標
# 3つの座標（図の座標, Pixel座標, 図法の座標）を出力する 

def transform_lonlat_to_figure(lonlat, ax, proj):
    # lonlat:経度と緯度  (lon, lat) 
    # ax: Axes図の座標系    例：fig.add_subplot()の戻り値
    # proj: axで指定した図法 
    #
    # 例：緯度経度をpointで与え、ステレオ図法にする場合
    # point = (140.0,35.0)
    # proj= ccrs.Stereographic(central_latitude=60, central_longitude=140) 
    # fig = plt.figure(figsize=(20,16))
    # ax = fig.add_subplot(1, 1, 1, projection=proj)
    # ax.set_extent([108, 156, 17, 55], ccrs.PlateCarree())
    #
    # 図法の変換：参照  https://scitools.org.uk/cartopy/docs/v0.14/crs/index.html                    
    point_proj = proj.transform_point(*lonlat, ccrs.PlateCarree())
    # pixel座標へ変換：参照　https://matplotlib.org/stable/tutorials/advanced/transforms_tutorial.html
    point_pix = ax.transData.transform(point_proj)
    # 図の座標へ変換                                                           
    point_fig = ax.transAxes.inverted().transform(point_pix)
    return point_fig, point_pix, point_proj

# 1地点のアメダスjsonデータから、elem要素で指定した値を返す(ただしFlagが 0以外は Noneとする)
# 要素:elem = 'temp','humidity','snow1h','snow6h','snow12h','snow24h','sun10m','sun1h','precipitation10m','precipitation1h','precipitation3h','precipitation24h','wind','windDirection'
def get_obs_value(amd_obs,elem):
    try:
        et = amd_obs[elem]
        if int(et[1]) != 0:
            return None
        return float(et[0])
    except Exception:
        return None
    
# 気象庁全国合成レーダーGPVのダウンロード
def download_time(time):    
    # 生存圏研究所ダウンロード元サイト
    http  = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original"  
    # データファイル名の指定
    day_dir = time.strftime("%Y/%m/%d")
    basename = "Z__C_RJTD_{}00_RDR_JMAGPV__grib2.tar".format(time.strftime("%Y%m%d%H%M"))
    fname    = "./{}".format(basename)
    # すでにファイルが存在しなければ、ダウンロードを行う
    if os.path.exists(fname):
        print("{} Already exists.".format(basename))
    else:
        url      = "{}/{}/{}".format(http,  day_dir, basename)
        # wgetコマンドでデータのダウンロード
        subprocess.run("wget {} -P ./".format(url), shell=True)
        # tarコマンドでダウンロードした圧縮ファイルの解凍
        subprocess.run("tar -xvf {} -C ./".format(fname), shell=True)   
    GgisFile = "./Z__C_RJTD_{}00_RDR_JMAGPV_Ggis1km_Prr10lv_ANAL_grib2.bin".format(time.strftime("%Y%m%d%H%M"))
    return GgisFile

# 描画指定：順に気圧(右上),気温(左上),湿球温度(右下),露点温度(左下))
npre_dispflag=False
temp_dispflag=False
wbt_dispflag=False
dp_dispflag=False

markersize_0 = 2 # マーカーサイズ
char_size=16 # 文字サイズ
barb_length=8 # 矢羽の長さ
dlon,dlat=1,1   # 緯線・経線は1度ごと

# 描画地域と描画時刻の設定
if len(sys.argv) == 3 and (sys.argv[1] == '0' or sys.argv[1] == '1' or sys.argv[1] == '2' or sys.argv[1] == '3'): 
    area = sys.argv[1] #第一引数は描画地域（0:北海道、1:東日本、2:西日本、3:東北）
    arg = sys.argv[2] #第二引数は描画時刻
    dt = parse_datetime(arg)
elif len(sys.argv) == 2 and (sys.argv[1] == '0' or sys.argv[1] == '1' or sys.argv[1] == '2' or sys.argv[1] == '3'):
    area = sys.argv[1] #第二引数は省略可能で、その場合は30分以上前（データ取得可能時刻）の直近の正時
    jst = pytz.timezone('Asia/Tokyo')
    dt = datetime.now(jst) - timedelta(minutes=30)
elif len(sys.argv) == 2 and len(sys.argv[1]) >= 8:
    area = 1 #第一引数を省略する場合は、デフォルトで東日本を描画し、時刻は引数によって定める
    dt = parse_datetime(arg)
elif len(sys.argv) == 1:
    area = 1 #すべての引数を省略することも可能
    jst = pytz.timezone('Asia/Tokyo')
    dt = datetime.now(jst) - timedelta(minutes=30)
else:
    print('Usage: python script.py [areacode[0:North][1:East][2:West]] [YYYYMMDDHH(MM)]')
    exit()

print(area)

# 描画開始メッセージ    
if dt:
    year=dt.year
    month=dt.month
    day=dt.day
    hour=dt.hour
    min=0
    dt = datetime(int(year), int(month), int(day), int(hour), int(min))
    print("読み込み観測時刻 {:4d}/{:02d}/{:02d} {:02d}:{:02d}".format(year,month,day,hour,min))
else:
    print('Usage: python script.py [areacode[1:North][2:East][3:West]] [YYYYMMDDHH(MM)]')
    exit()
    
# 観測データJSONの url作成
url_data_json= 'https://www.jma.go.jp/bosai/amedas/data/map/{:4d}{:02d}{:02d}{:02d}{:02d}00.json'
url_data_json=url_data_json.format(year,month,day,hour,min)
# 気象庁HPからアメダスデータを読み込む
response = urlopen(url_data_json)
content = response.read()
response.close()
data_json=content.decode()
dat_json = json.loads(data_json)

# アメダス地点Tableのurl
url_station_json="https://www.jma.go.jp/bosai/amedas/const/amedastable.json"
# アメダス地点Tableを読み込む
response = urlopen(url_station_json)
content = response.read()
response.close()
station_json=content.decode()
amd_json = json.loads(station_json)
    
# 描画する時間の指定(年,月,日,時,分)：データは10分ごと（前10分の雨量が記録されている）    
# アメダスデータと同じ時刻のレーダーGPVをダウンロード
time = pd.Timestamp(year,month,day,hour,min)
utc = time - offsets.Hour(9)
filepath = download_time(utc)

# データを読む
rain = load_jmara_grib2(filepath) #レーダーGPV
prmsl, tmp = read_msm(utc) #MSM海面気圧

# 地形データ取得
data = np.fromfile("LANDSEA.MSM_5K", dtype=np.dtype('>f4'))  # 地形バイナリデータをBig Endianの単精度浮動小数点数として読み込む

# データサイズに応じてreshape（MSMのグリッドサイズに合わせる）
data = data.reshape(505, 481)

# メッシュグリッドの作成
grid_lon, grid_lat = np.meshgrid(np.arange(120, 150 + 0.0625, 0.0625), np.arange(22.4, 47.6, 0.05))
sealand = np.flip(data*10000, axis=0)

# ガウシアンフィルタを適用
sealand_filterd = gaussian_filter(sealand, sigma=10.0) # sigmaはガウス分布の標準偏差

# 図法指定                                                                             
proj = ccrs.PlateCarree()
latlon_proj = ccrs.PlateCarree()

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

# 地図の描画範囲指定
if (area == '0'):
    i_area = [139, 147, 40, 46]
    areaname="Hokkaido"
elif (area == '1'):
    i_area = [134, 142, 33, 39]
    areaname="East"
elif (area == '2'):
    i_area = [128, 136, 31, 37]
    areaname="West"
elif (area == '3'):
    i_area = [135, 143, 36, 42]
    areaname="Tohoku"
else:
    i_area = [134, 142, 33, 39]
    areaname="East"

# 図のSIZE指定inch                                                                        
fig = plt.figure(figsize=(16,12))
# 余白設定                                                                                
plt.subplots_adjust(left=0.04, right=1.1, bottom=0.0, top=1.0)                  
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent(i_area, latlon_proj)

# レーダーGPV描画
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
    # 緯度・経度のtuple(度分形式)をtuple(度単位)に変換
    wlat = amd_json[stno]['lat'][0] + amd_json[stno]['lat'][1]/60.0
    wlon = amd_json[stno]['lon'][0] + amd_json[stno]['lon'][1]/60.0
    walt = amd_json[stno]['alt']
    # 風
    wind_ok = True
    ws = get_obs_value(val,'wind')
    wd = get_obs_value(val,'windDirection')
    if ws is not None and wd is not None:
        # 16方位風向と風速から、u,v作成   
        if ws == None or wd == None:
            u = None
            v = None
        else:
            wd = wd / 8.0 * math.pi  # 1/8 = (360/16) / 360 * 2
            au = -1.0 * ws * math.sin(wd)
            av = -1.0 * ws * math.cos(wd)
    else:
        wind_ok = False
    # 気温
    temp = get_obs_value(val,'temp')
    if temp is None:
        temp = np.nan
    elif walt < 800: #標高の高い観測点は無視する
        # 配列に格納
        tempsl = temp
        lat_list_t.append(wlat)
        lon_list_t.append(wlon)
        temp_list.append(tempsl)
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

# 線形補間
grid_temp = griddata((lon_list_t, lat_list_t), temp_list, (grid_lon, grid_lat), method='linear')
grid_npre = griddata((lon_list_p, lat_list_p), npre_list, (grid_lon, grid_lat), method='linear')

# 海上のデータは観測がないためMSMで補正する
#grid_npre[sealand == 0] = (prmsl[sealand == 0] + grid_npre[sealand == 0]) / 2 #海上の格子はアメダスによる補外とMSM予報値の平均
grid_temp[sealand == 0] = (tmp[sealand == 0] + grid_temp[sealand == 0]) / 2
grid_npre[sealand_filterd <= 1] = prmsl[sealand_filterd <= 1] #陸地から十分離れた格子はMSM予報値をそのまま採用する
grid_temp[sealand_filterd <= 1] = tmp[sealand_filterd <= 1]

# データがない格子もMSM予報値をそのまま採用する
nan_indices_npre = np.isnan(grid_npre)
grid_npre[nan_indices_npre] = prmsl[nan_indices_npre]
nan_indices_temp = np.isnan(grid_temp)
grid_temp[nan_indices_temp] = tmp[nan_indices_temp]

# ガウシアンフィルタを適用
grid_temp = gaussian_filter(grid_temp, sigma=3.0) # sigmaはガウス分布の標準偏差
grid_npre = gaussian_filter(grid_npre, sigma=3.0) # sigmaはガウス分布の標準偏差

# 描画領域のデータを切り出す（等圧線のラベルを表示するためのおまじない）
lon_range = np.where((grid_lon[0, :] >= i_area[0] - 0.25) & (grid_lon[0, :] <= i_area[1] + 0.25))
lat_range = np.where((grid_lat[:, 0] >= i_area[2] - 0.25) & (grid_lat[:, 0] <= i_area[3] + 0.25))

# 切り出したい範囲のインデックスを取得
lon_indices = lon_range[0]
lat_indices = lat_range[0]

# 切り出し
grid_lon_sliced = grid_lon[lat_indices][:, lon_indices]
grid_lat_sliced = grid_lat[lat_indices][:, lon_indices]
psea_grid = grid_npre[lat_indices][:, lon_indices]
temp_grid = grid_temp[lat_indices][:, lon_indices]

# 等温線をプロット
levels = np.arange(-27, 57, 6)
cont = plt.contour(grid_lon_sliced, grid_lat_sliced, temp_grid, levels=levels, linewidths=2, linestyles='solid', colors='red')
levels2 = np.arange(-30, 60, 6)
cont2 = plt.contour(grid_lon_sliced, grid_lat_sliced, temp_grid, levels=levels2, linewidths=4, linestyles='solid', colors='red')

# 等温線のラベルを付ける
plt.clabel(cont, fontsize=20)
plt.clabel(cont2, fontsize=40)

# 等圧線をプロット
levels = np.arange(901, 1049, 2)
cont = plt.contour(grid_lon_sliced, grid_lat_sliced, psea_grid, levels=levels, linewidths=2, colors='black')
levels2 = np.arange(900, 1050, 2)
cont2 = plt.contour(grid_lon_sliced, grid_lat_sliced, psea_grid, levels=levels2, linewidths=4, colors='black')

# 等圧線のラベルを付ける
plt.clabel(cont, fontsize=20)
plt.clabel(cont2, fontsize=40)

## H stamp
maxid = detect_peaks(psea_grid, filter_size=40, dist_cut=10)
for i in range(len(maxid[0])):
    wlon = grid_lon_sliced[0][maxid[1][i]]
    wlat = grid_lat_sliced[maxid[0][i]][0]
    # 図の範囲内に座標があるか確認                                                                           
    fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj)
    if ( fig_z[0] > 0.05 and fig_z[0] < 0.95  and fig_z[1] > 0.05 and fig_z[1] < 0.95):
        ax.plot(wlon, wlat, marker='x' , markersize=32, color="blue", transform=latlon_proj)
        ax.text(wlon - 0.12, wlat + 0.12, 'H', size=60, color="blue", transform=latlon_proj)
        val = psea_grid[maxid[0][i]][maxid[1][i]]
        ival = int(val)
        ax.text(fig_z[0], fig_z[1] - 0.025, str(ival), size=48, color="blue", transform=ax.transAxes, verticalalignment="top", horizontalalignment="center")

## L stamp
minid = detect_peaks(psea_grid, filter_size=40, dist_cut=10, flag=1)
for i in range(len(minid[0])):
    wlon = grid_lon_sliced[0][minid[1][i]]
    wlat = grid_lat_sliced[minid[0][i]][0]
    # 図の範囲内に座標があるか確認                                                                           
    fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj)
    if ( fig_z[0] > 0.05 and fig_z[0] < 0.95  and fig_z[1] > 0.05 and fig_z[1] < 0.95):
        ax.plot(wlon, wlat, marker='x' , markersize=32, color="red", transform=latlon_proj)
        ax.text(wlon - 0.12, wlat + 0.12, 'L', size=60, color="red", transform=latlon_proj)
        val = psea_grid[minid[0][i]][minid[1][i]]
        ival = int(val)
        ax.text(fig_z[0], fig_z[1] - 0.025, str(ival), size=48, color="red", transform=ax.transAxes, verticalalignment="top", horizontalalignment="center")

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("AMeDAS and RadarGPV"), loc='left',size=20)
plt.title('Valid Time: {}'.format(dt), loc='right',size=20);
#plt.savefig("{}.jpg".format(time.strftime("%Y%m%d%H%M")), format="jpg")
plt.savefig("latest{}.jpg".format(areaname), format="jpg")
